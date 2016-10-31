/**
 *  @file   larpandoracontent/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArTwoDReco/LArSeedFinding/ClusterCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ClusterCharacterisationAlgorithm::ClusterCharacterisationAlgorithm() :
    m_slidingFitWindow(10),
    m_minHitsInCluster(20),
    m_maxLayerGapFraction(0.2f),
    m_maxWidthPerUnitLength(0.15f),
    m_maxShowerLength(1000.f),
    m_useDetectorGaps(true),
    m_overwriteExistingId(false),
    m_useUnavailableClusters(false),
    m_writeToTree(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ClusterCharacterisationAlgorithm::~ClusterCharacterisationAlgorithm()
{
    if (m_writeToTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::Run()
{
    m_clusterDirectionMap.clear();

    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "ClusterCharacterisationAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        for (const Cluster *const pCluster : *pClusterList)
        {
            if (!m_overwriteExistingId && (UNKNOWN_PARTICLE_TYPE != pCluster->GetParticleId()))
                continue;

            if (!m_useUnavailableClusters && !PandoraContentApi::IsAvailable(*this, pCluster))
                continue;

            PandoraContentApi::Cluster::Metadata metadata;

            if (this->IsClearTrack(pCluster, pClusterList))
            {
                metadata.m_particleId = MU_MINUS;
            }
            else
            {
                metadata.m_particleId = E_MINUS;
            }

            if (pCluster->GetParticleId() != metadata.m_particleId.Get())
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
        }
    }

    m_clusterDirectionMap.clear();
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster, const ClusterList *const pClusterList) const
{
    bool isTrueTrack(false);

    try
    {
        // ATTN Slightly curious definition of a clear track, but this is most-likely what is needed for shower-growing
        const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
        isTrueTrack = (PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId()));
    }
    catch (StatusCodeException &)
    {
    }

    // Tree variables here
    //--------------------------------------------------------------------------------------------------------------------------------------
    const int trueTrack(isTrueTrack ? 1 : 0);
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "trueTrack", trueTrack));

    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    const int hitTypeInt(static_cast<int>(hitType));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "hitType", hitTypeInt));

    const int nHits(pCluster->GetNCaloHits());
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHits", nHits));

    // vertexDistance
    float vertexDistance(-1.f);
    const VertexList *pVertexList = nullptr;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));
    const Vertex *const pSelectedVertex((pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);

    if (pSelectedVertex)
    {
        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSelectedVertex->GetPosition(), LArClusterHelper::GetClusterHitType(pCluster)));
        vertexDistance = LArClusterHelper::GetClosestDistance(vertexPosition2D, pCluster);
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vertexDistance", vertexDistance));

    // hit positions and energy
    FloatVector xPositions, zPositions;
    float mipEnergy(0.f);

    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        xPositions.push_back(pCaloHit->GetPositionVector().GetX());
        zPositions.push_back(pCaloHit->GetPositionVector().GetZ());
        mipEnergy += pCaloHit->GetMipEquivalentEnergy();
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "xPositions", &xPositions));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "zPositions", &zPositions));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mipEnergy", mipEnergy));

    // straight line length and integrated pathlength
    float straightLineLength(-1.f), integratedPathLength(-1.f);

    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, 5, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const CartesianVector globalMinLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
        const CartesianVector globalMaxLayerPosition(slidingFitResult.GetGlobalMaxLayerPosition());
        straightLineLength = (globalMaxLayerPosition - globalMinLayerPosition).GetMagnitude();

        integratedPathLength = 0.f;
        CartesianVector previousFitPosition(globalMinLayerPosition);
        const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());

        for (const auto &mapEntry : layerFitResultMap)
        {
            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);
            integratedPathLength += (thisFitPosition - previousFitPosition).GetMagnitude();
            previousFitPosition = thisFitPosition;
        }
    }
    catch (const StatusCodeException &)
    {
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "straightLineLength", straightLineLength));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "integratedPathLength", integratedPathLength));

    float straightLineLength10(-1.f), integratedPathLength10(-1.f);

    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, 10, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const CartesianVector globalMinLayerPosition(slidingFitResult.GetGlobalMinLayerPosition());
        const CartesianVector globalMaxLayerPosition(slidingFitResult.GetGlobalMaxLayerPosition());
        straightLineLength10 = (globalMaxLayerPosition - globalMinLayerPosition).GetMagnitude();

        integratedPathLength10 = 0.f;
        CartesianVector previousFitPosition(globalMinLayerPosition);
        const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());

        for (const auto &mapEntry : layerFitResultMap)
        {
            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);
            integratedPathLength10 += (thisFitPosition - previousFitPosition).GetMagnitude();
            previousFitPosition = thisFitPosition;
        }
    }
    catch (const StatusCodeException &)
    {
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "straightLineLength10", straightLineLength10));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "integratedPathLength10", integratedPathLength10));

    // shower fit width and gap length
    float showerFitWidth(-1.f), showerFitGapLength(-1.f);

    try
    {
        const TwoDSlidingShowerFitResult showerFitResult(pCluster, 10, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapP(showerFitResult.GetPositiveEdgeFitResult().GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapN(showerFitResult.GetNegativeEdgeFitResult().GetLayerFitResultMap());

        if (layerFitResultMapS.size() > 1)
        {
            CartesianVector globalMinLayerPositionOnAxis(0.f, 0.f, 0.f), globalMaxLayerPositionOnAxis(0.f, 0.f, 0.f);
            showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.begin()->second.GetL(), 0.f, globalMinLayerPositionOnAxis);
            showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.rbegin()->second.GetL(), 0.f, globalMaxLayerPositionOnAxis);
            const float straightLinePathLength((globalMaxLayerPositionOnAxis - globalMinLayerPositionOnAxis).GetMagnitude());

            if (straightLinePathLength > std::numeric_limits<float>::epsilon())
            {
                showerFitWidth = 0.f;
                showerFitGapLength = 0.f;
                CartesianVector previousLayerPosition(globalMinLayerPositionOnAxis);

                for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
                {
                    CartesianVector thisLayerPosition(0.f, 0.f, 0.f);
                    showerFitResult.GetShowerFitResult().GetGlobalPosition(iterS->second.GetL(), 0.f, thisLayerPosition);
                    const float thisGapLength((thisLayerPosition - previousLayerPosition).GetMagnitude());

                    if (thisGapLength > showerFitGapLength)
                    {
                        const float minZ(std::min(thisLayerPosition.GetZ(), previousLayerPosition.GetZ()));
                        const float maxZ(std::max(thisLayerPosition.GetZ(), previousLayerPosition.GetZ()));

                        if ((maxZ - minZ) < std::numeric_limits<float>::epsilon())
                            throw StatusCodeException(STATUS_CODE_FAILURE);

                        const float gapZ(LArGeometryHelper::CalculateGapDeltaZ(this->GetPandora(), minZ, maxZ, LArClusterHelper::GetClusterHitType(pCluster)));
                        const float correctedGapLength(thisGapLength * (1.f - gapZ / (maxZ - minZ)));

                        if (correctedGapLength > showerFitGapLength)
                            showerFitGapLength = correctedGapLength;
                    }

                    previousLayerPosition = thisLayerPosition;

                    LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(iterS->first);
                    LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(iterS->first);

                    if ((layerFitResultMapP.end() == iterP) || (layerFitResultMapN.end() == iterN))
                        continue;

                    showerFitWidth += std::fabs(iterP->second.GetFitT() - iterN->second.GetFitT());
                }
            }
        }
    }
    catch (StatusCodeException &)
    {
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerFitWidth", showerFitWidth));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "showerFitGapLength", showerFitGapLength));

    // nPointsOfContact, nHitsInPrimaryBranches
    int nPointsOfContact(0), nHitsInBranches(0);
    const CartesianVector vertexPosition2D(pSelectedVertex ? LArGeometryHelper::ProjectPosition(this->GetPandora(), pSelectedVertex->GetPosition(), hitType) : CartesianVector(0.f, 0.f, 0.f));

    ClusterVector candidateClusters;
    for (const Cluster *const pCandidateCluster : *pClusterList)
    {
        if ((pCandidateCluster == pCluster) || (pCandidateCluster->GetNCaloHits() < 5))
            continue;

        try
        {
            if (pSelectedVertex && this->IsVertexAssociated(LArPointingCluster(pCluster), vertexPosition2D))
                continue;
        }
        catch (const StatusCodeException &)
        {
        }

        candidateClusters.push_back(pCandidateCluster);
    }

    ClusterUsageMap forwardUsageMap, backwardUsageMap;
    this->FindAssociatedClusters(pCluster, candidateClusters, forwardUsageMap, backwardUsageMap);

    SeedAssociationList seedAssociationList;
    this->IdentifyClusterMerges(ClusterVector(1, pCluster), backwardUsageMap, seedAssociationList);

    for (const Cluster *const pBranchCluster : seedAssociationList.at(pCluster))
    {
        ++nPointsOfContact;
        nHitsInBranches += pBranchCluster->GetNCaloHits();
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nPointsOfContact", nPointsOfContact));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nHitsInBranches", nHitsInBranches));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

    //--------------------------------------------------------------------------------------------------------------------------------------
    // End tree variables calculations

    if (nHits < 6)
        return false;

    if (straightLineLength < std::numeric_limits<float>::epsilon())
        return false;

    if (straightLineLength > 80.f)
        return true;

    if (vertexDistance / straightLineLength > 0.4f)
        return false;

    if (nPointsOfContact > 4)
        return false;

    if (0 == nHits)
        return false;

    if (static_cast<float>(nHitsInBranches) / static_cast<float>(nHits) > 5.f)
        return false;

    if (integratedPathLength / straightLineLength > 1.3f)
        return false;

    if (showerFitWidth / straightLineLength > 1.7f)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsInCluster", m_minHitsInCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLayerGapFraction", m_maxLayerGapFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxWidthPerUnitLength", m_maxWidthPerUnitLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxShowerLength", m_maxShowerLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseDetectorGaps", m_useDetectorGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OverwriteExistingId", m_overwriteExistingId));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseUnavailableClusters", m_useUnavailableClusters));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    if (m_writeToTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputTree", m_treeName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputFile", m_fileName));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
