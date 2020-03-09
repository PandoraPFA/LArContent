/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/MatchingBaseAlgorithm.cc
 *
 *  @brief  Implementation of the three dimension algorithm base class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/MatchingBaseAlgorithm.h"

using namespace pandora;

namespace lar_content
{

MatchingBaseAlgorithm::MatchingBaseAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MatchingBaseAlgorithm::~MatchingBaseAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MatchingBaseAlgorithm::CreateThreeDParticles(const ProtoParticleVector &protoParticleVector)
{
    bool particlesMade(false);
    const PfoList *pPfoList(nullptr); std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (const ProtoParticle &protoParticle : protoParticleVector)
    {
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        this->SetPfoParameters(protoParticle, pfoParameters);

        const ParticleFlowObject *pPfo(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
        particlesMade = true;
    }

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
    }

    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchingBaseAlgorithm::SetPfoParameters(const ProtoParticle &protoParticle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const
{
    this->SetPfoParticleId(pfoParameters);
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList.insert(pfoParameters.m_clusterList.end(), protoParticle.m_clusterListU.begin(), protoParticle.m_clusterListU.end());
    pfoParameters.m_clusterList.insert(pfoParameters.m_clusterList.end(), protoParticle.m_clusterListV.begin(), protoParticle.m_clusterListV.end());
    pfoParameters.m_clusterList.insert(pfoParameters.m_clusterList.end(), protoParticle.m_clusterListW.begin(), protoParticle.m_clusterListW.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchingBaseAlgorithm::SetPfoParticleId(PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const
{
    pfoParameters.m_particleId = E_MINUS; // Shower
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MatchingBaseAlgorithm::MakeClusterMerges(const ClusterMergeMap &clusterMergeMap)
{
    ClusterSet deletedClusters;

    ClusterList parentClusters;
    for (const auto &mapEntry : clusterMergeMap) parentClusters.push_back(mapEntry.first);
    parentClusters.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : parentClusters)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pParentCluster));
        const std::string &clusterListName(this->GetClusterListName(hitType));

        if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        ClusterList daughterClusters(clusterMergeMap.at(pParentCluster));
        daughterClusters.sort(LArClusterHelper::SortByNHits);

        for (const Cluster *const pDaughterCluster : daughterClusters)
        {
            if (deletedClusters.count(pParentCluster) || deletedClusters.count(pDaughterCluster))
                throw StatusCodeException(STATUS_CODE_FAILURE);

            this->UpdateUponDeletion(pDaughterCluster);
            this->UpdateUponDeletion(pParentCluster);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster, clusterListName, clusterListName));

            this->UpdateForNewCluster(pParentCluster);
            deletedClusters.insert(pDaughterCluster);
        }
    }

    return !(deletedClusters.empty());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MatchingBaseAlgorithm::MakeClusterSplits(const SplitPositionMap &splitPositionMap)
{
    bool changesMade(false);

    ClusterList splitClusters;
    for (const auto &mapEntry : splitPositionMap) splitClusters.push_back(mapEntry.first);
    splitClusters.sort(LArClusterHelper::SortByNHits);

    for (const Cluster *pCurrentCluster : splitClusters)
    {
        CartesianPointVector splitPositions(splitPositionMap.at(pCurrentCluster));
        std::sort(splitPositions.begin(), splitPositions.end(), MatchingBaseAlgorithm::SortSplitPositions);

        const HitType hitType(LArClusterHelper::GetClusterHitType(pCurrentCluster));
        const std::string &clusterListName(this->GetClusterListName(hitType));

        if (!((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType)))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

        for (const CartesianVector &splitPosition : splitPositions)
        {
            const Cluster *pLowXCluster(nullptr), *pHighXCluster(nullptr);
            this->UpdateUponDeletion(pCurrentCluster);

            if (this->MakeClusterSplit(splitPosition, pCurrentCluster, pLowXCluster, pHighXCluster))
            {
                changesMade = true;
                this->UpdateForNewCluster(pLowXCluster);
                this->UpdateForNewCluster(pHighXCluster);
                pCurrentCluster = pHighXCluster;
            }
            else
            {
                this->UpdateForNewCluster(pCurrentCluster);
            }
        }
    }

    return changesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MatchingBaseAlgorithm::MakeClusterSplit(const CartesianVector &splitPosition, const Cluster *&pCurrentCluster, const Cluster *&pLowXCluster,
    const Cluster *&pHighXCluster) const
{
    CartesianVector lowXEnd(0.f, 0.f, 0.f), highXEnd(0.f, 0.f, 0.f);

    try
    {
        LArPointingCluster pointingCluster(pCurrentCluster);
        const bool innerIsLowX(pointingCluster.GetInnerVertex().GetPosition().GetX() < pointingCluster.GetOuterVertex().GetPosition().GetX());
        lowXEnd = (innerIsLowX ? pointingCluster.GetInnerVertex().GetPosition() : pointingCluster.GetOuterVertex().GetPosition());
        highXEnd = (innerIsLowX ? pointingCluster.GetOuterVertex().GetPosition() : pointingCluster.GetInnerVertex().GetPosition());
    }
    catch (const StatusCodeException &) {return false;}

    const CartesianVector lowXUnitVector((lowXEnd - splitPosition).GetUnitVector());
    const CartesianVector highXUnitVector((highXEnd - splitPosition).GetUnitVector());

    CaloHitList caloHitList;
    pCurrentCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    std::string originalListName, fragmentListName;
    const ClusterList clusterList(1, pCurrentCluster);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, originalListName, fragmentListName));

    pLowXCluster = nullptr;
    pHighXCluster = nullptr;

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector unitVector((pCaloHit->GetPositionVector() - splitPosition).GetUnitVector());
        const float dotProductLowX(unitVector.GetDotProduct(lowXUnitVector));
        const float dotProductHighX(unitVector.GetDotProduct(highXUnitVector));

        const Cluster *&pClusterToModify((dotProductLowX > dotProductHighX) ? pLowXCluster : pHighXCluster);

        if (!pClusterToModify)
        {
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList.push_back(pCaloHit);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterToModify));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterToModify, pCaloHit));
        }
    }

    if (!pLowXCluster || !pHighXCluster)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, originalListName, fragmentListName));
        return false;
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MatchingBaseAlgorithm::SortSplitPositions(const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs)
{
    return (lhs.GetX() < rhs.GetX());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchingBaseAlgorithm::SelectAllInputClusters()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchingBaseAlgorithm::PreparationStep()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MatchingBaseAlgorithm::TidyUp()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MatchingBaseAlgorithm::Run()
{
    try
    {
        this->SelectAllInputClusters();
        this->PreparationStep();
        this->PerformMainLoop();
        this->ExamineOverlapContainer();
        this->TidyUp();
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();

        if (STATUS_CODE_SUCCESS != statusCodeException.GetStatusCode())
            throw statusCodeException;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MatchingBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
