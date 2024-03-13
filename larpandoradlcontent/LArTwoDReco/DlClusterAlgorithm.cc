/**
 *  @file   larpandoradlcontent/LArTwoDReco/DlClusterAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning clustering algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArTwoDReco/DlClusterAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include <Eigen/Dense>

#include <random>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

StatusCode DlClusterAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        if (pCaloHitList->empty())
            continue;
        std::cout << "Num hits: " << pCaloHitList->size() << std::endl;
        Eigen::MatrixXf hits(pCaloHitList->size(), 2);
        int i{0};
        for (const CaloHit *const pCaloHit : *pCaloHitList)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            hits(i, 0) = pos.GetX();
            hits(i, 1) = pos.GetZ();
            ++i;
        }
        // Edges can be double-counted, so use map of maps to avoid this
        std::map<const CaloHit *, std::map<const CaloHit *, bool>> edges;
        for (int r = 0; r < hits.rows(); ++r)
        {
            Eigen::RowVectorXf row(2);
            row << hits(r,0), hits(r,1);
            Eigen::MatrixXf norms((hits.rowwise() - row).array().pow(2).rowwise().sum());
            norms(r, 0) = std::numeric_limits<float>::max();
            Eigen::Index index1, index2;
            norms.col(0).minCoeff(&index1);
            auto iter0{pCaloHitList->begin()};
            std::advance(iter0, r);
            auto iter1{pCaloHitList->begin()};
            std::advance(iter1, index1);
            edges[*iter0][*iter1] = true;
            edges[*iter1][*iter0] = true;
            norms(index1, 0) = std::numeric_limits<float>::max();
            auto val2{norms.col(0).minCoeff(&index2)};
            auto val3{(hits.row(index1) - hits.row(index2)).array().pow(2).rowwise().sum()};
            if (val2 < val3(0))
            {
                auto iter2{pCaloHitList->begin()};
                std::advance(iter2, index2);
                edges[*iter0][*iter2] = true;
                edges[*iter2][*iter0] = true;
            }
        }

        for (const auto &[pCaloHit1, map] : edges)
        {
            const CartesianVector &pos1{pCaloHit1->GetPositionVector()};
            for (const auto &[pCaloHit2, dummy] : map)
            {
                const CartesianVector &pos2{pCaloHit2->GetPositionVector()};
                CartesianVector start(pos1.GetX(), 0, pos1.GetZ());
                CartesianVector end(pos2.GetX(), 0, pos2.GetZ());
                PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &start, &end, "e", BLUE, 2, 1));
            }
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlClusterAlgorithm::PrepareTrainingSample()
{
    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        if (pCaloHitList->empty())
            continue;
        CaloHitVector caloHitVector(pCaloHitList->begin(), pCaloHitList->end());
        std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

        std::map<const MCParticle *, CaloHitList> mcToHitsMap;
        for (const CaloHit *pCaloHit : caloHitVector)
        {
            try
            {
                const MCParticle *pMC{MCParticleHelper::GetMainMCParticle(pCaloHit)};
                mcToHitsMap[pMC].emplace_back(pCaloHit);
            }
            catch (const StatusCodeException &)
            {
            }
        }

        LArMvaHelper::MvaFeatureVector featureVector;
        featureVector.emplace_back(static_cast<double>(mcToHitsMap.size()));
        for (const auto & [pMC, caloHitList] : mcToHitsMap)
        {
            (void)pMC;
            featureVector.emplace_back(static_cast<double>(caloHitList.size()));
            for (const CaloHit *pCaloHit : *pCaloHitList)
            {
                const CartesianVector &position{pCaloHit->GetPositionVector()};
                const float x{position.GetX()}, z{position.GetZ()}, adc{pCaloHit->GetMipEquivalentEnergy()};
                featureVector.emplace_back(x);
                featureVector.emplace_back(z);
                featureVector.emplace_back(adc);
            }
        }
        const HitType view{pCaloHitList->front()->GetHitType()};
        std::string suffix{view == TPC_VIEW_U ? "_U.csv" : view == TPC_VIEW_V ? "_V.csv" : "_W.csv"};
        LArMvaHelper::ProduceTrainingExample(m_outputFilePrefix + suffix, true, featureVector);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlClusterAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFilePrefix", m_outputFilePrefix));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MaxEdgeLength", m_maxEdgeLengthSquared));
    m_maxEdgeLengthSquared *= m_maxEdgeLengthSquared;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MaxEdgeRatio", m_maxEdgeRatioSquared));
    m_maxEdgeRatioSquared *= m_maxEdgeRatioSquared;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Prune", m_prune));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
