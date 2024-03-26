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
#include "larpandoracontent/LArObjects/LArGraph.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include <Eigen/Dense>

#include <random>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlClusterAlgorithm::DlClusterAlgorithm() :
    m_fullyConnect{true},
    m_nSourceEdges{2},
    m_maxSecondaryCosine{0.996f},
    m_maxSecondaryDistance{3.f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlClusterAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
        if (pCaloHitList->empty())
            continue;
        LArGraph graph(m_fullyConnect, m_nSourceEdges, m_maxSecondaryCosine, m_maxSecondaryDistance);
        graph.MakeGraph(*pCaloHitList);
        const LArGraph::EdgeVector &edges{graph.GetEdges()};
        for (const LArGraph::Edge *const edge : edges)
        {
            const CartesianVector &start{edge->m_v0->GetPositionVector()};
            const CartesianVector &end{edge->m_v1->GetPositionVector()};
            PANDORA_MONITORING_API(AddLineToVisualization(this->GetPandora(), &start, &end, "e", BLUE, 2, 1));
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
        for (const auto &[pMC, caloHitList] : mcToHitsMap)
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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FullyConnect", m_fullyConnect));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NumSourceEdges", m_nSourceEdges));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSecondaryCosine", m_maxSecondaryCosine));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSecondaryDistance", m_maxSecondaryDistance));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
