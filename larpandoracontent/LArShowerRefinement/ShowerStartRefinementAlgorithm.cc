/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the shower start refinement base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.h"

using namespace pandora;

namespace lar_content
{

ShowerStartRefinementAlgorithm::ShowerStartRefinementAlgorithm() : m_binSize(0.005)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ShowerStartRefinementAlgorithm::~ShowerStartRefinementAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ShowerDistribution", "ShowerDistribution.root", "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "THE LIMIT DOES NOT EXIST" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartRefinementAlgorithm::Run()
{
    PfoVector pfoVector;
    this->FillPfoVector(pfoVector);

    CartesianVector nuVertexPosition(0.f, 0.f, 0.f);
    if (this->GetNeutrinoVertex(nuVertexPosition) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    // run tools
    for (const ParticleFlowObject *const pPfo : pfoVector)
    {
        if (std::find(m_deletedPfos.begin(), m_deletedPfos.end(), pPfo) != m_deletedPfos.end())
            continue;

        for (ShowerStartRefinementBaseTool *const pShowerStartRefinementTool : m_algorithmToolVector)
            pShowerStartRefinementTool->Run(this, pPfo, nuVertexPosition);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementAlgorithm::FillPfoVector(PfoVector &pfoVector)
{
    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        if (!pPfoList || pPfoList->empty())
        {
            std::cout << "ShowerStartRefinementAlgorithm: unable to find pfo list " << pfoListName << std::endl;
            continue;
        }

        pfoVector.insert(pfoVector.begin(), pPfoList->begin(), pPfoList->end());
    }

    // This ordering is important.
    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartRefinementAlgorithm::GetNeutrinoVertex(CartesianVector &neutrinoVertex)
{
    const VertexList *pNuVertexList(nullptr);
    const StatusCode statusCode(PandoraContentApi::GetList(*this, m_neutrinoVertexListName, pNuVertexList));

    if (statusCode != STATUS_CODE_SUCCESS)
        return statusCode;

    if (!pNuVertexList || (pNuVertexList->size() != 1))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ShowerStartRefinementAlgorithm: unable to find vertex list " << m_neutrinoVertexListName << " if it does exist, it may have more than one nu vertex" << std::endl;

        return STATUS_CODE_NOT_INITIALIZED;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementAlgorithm::FillOwnershipMaps()
{
    m_hitToClusterMapU.clear(); m_hitToClusterMapV.clear(); m_hitToClusterMapW.clear();
    m_clusterToPfoMapU.clear(); m_clusterToPfoMapV.clear(); m_clusterToPfoMapW.clear();

    // First fill pfo maps
    PfoVector pfoVector;

    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        //PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, );
        PandoraContentApi::GetList(*this, pfoListName, pPfoList);

        if (!pPfoList || pPfoList->empty())
        {
            std::cout << "ShowerStartRefinementAlgorithm: unable to find pfo list " << pfoListName << std::endl;
            continue;
        }

        pfoVector.insert(pfoVector.begin(), pPfoList->begin(), pPfoList->end());
    }


    for (const ParticleFlowObject *const pPfo : pfoVector)
    {
        ClusterList twoDClusterList;
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_U, twoDClusterList);
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_V, twoDClusterList);
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, twoDClusterList);

        for (const Cluster *const pCluster : twoDClusterList)
        {
            CaloHitList caloHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

            CaloHitList isolated(pCluster->GetIsolatedCaloHitList());

            for (const CaloHit *const pIsolated : isolated)
            {
                if (std::find(caloHitList.begin(), caloHitList.end(), pIsolated) == caloHitList.end())
                    caloHitList.push_back(pIsolated);
            }

            const HitType hitType(caloHitList.front()->GetHitType());
            HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
            ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);

            for (const CaloHit *const pCaloHit : caloHitList)
                hitToClusterMap[pCaloHit] = pCluster;

            clusterToPfoMap[pCluster] = pPfo;
        }
    }

    // Now fill cluster maps
    StringVector clusterListNames;
    clusterListNames.push_back("ClustersU");
    clusterListNames.push_back("ClustersV");
    clusterListNames.push_back("ClustersW");

    ClusterList clusterList;
    for (const std::string &clusterListName : clusterListNames)
    {
        const ClusterList *pClusterList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            std::cout << "ShowerStartRefinementAlgorithm: No cluster list found, returning..." << std::endl;
            throw;
        }

        clusterList.insert(clusterList.end(), pClusterList->begin(), pClusterList->end());
    }

    for (const Cluster *const pCluster : clusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
        ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);

        if (clusterToPfoMap.find(pCluster) != clusterToPfoMap.end())
            continue;

        if (!pCluster->IsAvailable())
        {
            std::cout << "CLUSTER IS NOT AVAILABLE ISOBE" << std::endl;
            throw;
        }

        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        CaloHitList isolated(pCluster->GetIsolatedCaloHitList());

        for (const CaloHit *const pIsolated : isolated)
        {
            if (std::find(caloHitList.begin(), caloHitList.end(), pIsolated) == caloHitList.end())
                caloHitList.push_back(pIsolated);
        }

        for (const CaloHit *const pCaloHit : caloHitList)
            hitToClusterMap[pCaloHit] = pCluster;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::IsElectronPathway(const CaloHitList &hitsToAdd)
{
    int showerHitCount;

    for (const CaloHit *const pHitToAdd : hitsToAdd)
    {
        try
        {
            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pHitToAdd));
            const int pdg(std::abs(pMCParticle->GetParticleId()));

            if ((pdg == 11) || (pdg == 22))
                ++showerHitCount;
        }
        catch(...)
        {
            continue;
        }
    }

    const float showerProportion(static_cast<float>(showerHitCount) / static_cast<float>(hitsToAdd.size()));

    return (showerProportion > 0.5f); 
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerStartRefinementAlgorithm::AddElectronPathway(const ParticleFlowObject *const pShowerPfo, const CaloHitList &pathwayHitList)
{
    // This is so incredibly lazy isobel
    this->FillOwnershipMaps();

    const HitType hitType(pathwayHitList.front()->GetHitType());

    ClusterList showerClusters2D, showerClusters3D;
    LArPfoHelper::GetClusters(pShowerPfo, hitType, showerClusters2D);
    LArPfoHelper::GetClusters(pShowerPfo, TPC_3D, showerClusters3D);

    const ThreeDSlidingFitResult showerSlidingFit(showerClusters3D.front(), 1000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const CartesianVector &showerDirection3D(showerSlidingFit.GetGlobalMaxLayerPosition());

    std::map<const ParticleFlowObject*, int> showerHitCountMap;

    for (const CaloHit *const pPathwayHit : pathwayHitList)
    {
        const HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
        const ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);

        if (hitToClusterMap.find(pPathwayHit) == hitToClusterMap.end())
            continue;

        if (clusterToPfoMap.find(hitToClusterMap.at(pPathwayHit)) == clusterToPfoMap.end())
            continue;

        const ParticleFlowObject *const pPathwayShower(clusterToPfoMap.at(hitToClusterMap.at(pPathwayHit)));

        if (LArPfoHelper::IsTrack(pPathwayShower))
            continue;

        if (showerHitCountMap.find(pPathwayShower) == showerHitCountMap.end())
            showerHitCountMap[pPathwayShower] = 1;
        else
            showerHitCountMap[pPathwayShower] = showerHitCountMap[pPathwayShower] + 1;
    }

    PfoList significantShowersToMerge;

    for (const auto &entry : showerHitCountMap)
    {
        float contaminationRatio(static_cast<float>(entry.second) / static_cast<float>(pathwayHitList.size()));

        if (contaminationRatio < 0.3f)
            continue;

        ClusterList pathwayClusters3D;
        LArPfoHelper::GetClusters(entry.first, TPC_3D, pathwayClusters3D);

        const ThreeDSlidingFitResult pathwaySlidingFit(pathwayClusters3D.front(), 1000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const CartesianVector &pathwayDirection(pathwaySlidingFit.GetGlobalMaxLayerPosition());

        const float openingAngle(pathwayDirection.GetOpeningAngle(showerDirection3D) * 180 / M_PI);

        if (openingAngle > 5.f)
            significantShowersToMerge.push_back(entry.first);
    }

    // Add in hits first, then deal with merges
    for (const CaloHit *const pPathwayHit : pathwayHitList)
    {
        const HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
        const ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);
        std::string clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

        const Cluster *pParentCluster(nullptr);
        const ParticleFlowObject *pParentPfo(nullptr);

        if (hitToClusterMap.find(pPathwayHit) != hitToClusterMap.end())
        {
            pParentCluster = hitToClusterMap.at(pPathwayHit);

            if (clusterToPfoMap.find(pParentCluster) != clusterToPfoMap.end())
            {
                pParentPfo = clusterToPfoMap.at(pParentCluster);

                if (std::find(significantShowersToMerge.begin(), significantShowersToMerge.end(), pParentPfo) != significantShowersToMerge.end())
                    continue;
            }
        }

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

            CaloHitList clusterNormalHitList; const CaloHitList clusterIsolatedHitList(pParentCluster->GetIsolatedCaloHitList());
            pParentCluster->GetOrderedCaloHitList().FillCaloHitList(clusterNormalHitList);

            const bool isIsolated(std::find(clusterIsolatedHitList.begin(), clusterIsolatedHitList.end(), pPathwayHit) != clusterIsolatedHitList.end());

            if (!isIsolated && (clusterNormalHitList.size() == 1) && !(clusterIsolatedHitList.empty()))
            {
                const HitType isolatedHitType(LArClusterHelper::GetClusterHitType(pParentCluster));
                HitToClusterMap &isolatedHitToClusterMap(isolatedHitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);

                for (const CaloHit * const pIsolatedHit : clusterIsolatedHitList)
                {
                    isolatedHitToClusterMap.erase(pIsolatedHit);
                    const StatusCode isolatedStatusCode(PandoraContentApi::RemoveIsolatedFromCluster(*this, pParentCluster, pIsolatedHit));

                    if (isolatedStatusCode != STATUS_CODE_SUCCESS)
                    {
                        std::cout << "ISOBEL CANNOT REMOVE ISOLATED HIT?" << std::endl;
                        throw;
                    }
                }
            }

            const StatusCode statusCodeCluster(isIsolated ? PandoraContentApi::RemoveIsolatedFromCluster(*this, pParentCluster, pPathwayHit) : 
                PandoraContentApi::RemoveFromCluster(*this, pParentCluster, pPathwayHit));

            if (statusCodeCluster != STATUS_CODE_SUCCESS)
            {
                if (statusCodeCluster != STATUS_CODE_NOT_ALLOWED)
                {
                    std::cout << "ElectronStartRefinementTool: cluster jam" << std::endl;
                    throw StatusCodeException(statusCodeCluster);
                }

                if (pParentPfo)
                {
                    const StatusCode statusCodePfo(PandoraContentApi::RemoveFromPfo(*this, pParentPfo, pParentCluster));
                    const unsigned int nHits(LArPfoHelper::GetNumberOfTwoDHits(pParentPfo));

                    if (nHits == 0)
                        std::cout << "ElectronStartRefinementTool: ISOBEL - PFO HAS ZERO HITS" << std::endl;

                    if (statusCodePfo != STATUS_CODE_SUCCESS)
                    {
                        std::cout << "ElectronStartRefinementTool: pfo jam" << std::endl;
                        throw;
                    }
                }

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pParentCluster));
            }
        }

        if (!PandoraContentApi::IsAvailable(*this, pPathwayHit))
        {
            std::cout << "CALO HIT IS NOT AVAILABLE!!" << std::endl;
            throw;
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, showerClusters2D.front(), pPathwayHit));
    }

    // Now handle the shower merges
    for (const ParticleFlowObject *const pShowerToMerge : significantShowersToMerge)
    {
        m_deletedPfos.push_back(pShowerToMerge);
        this->MergeAndDeletePfos(pShowerPfo, pShowerToMerge);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerStartRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoVertexListName", m_neutrinoVertexListName));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "ShowerStartRefinementTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        ShowerStartRefinementBaseTool *const pShowerStartRefinementTool(dynamic_cast<ShowerStartRefinementBaseTool *>(*iter));

        if (!pShowerStartRefinementTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pShowerStartRefinementTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BinSize", m_binSize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
