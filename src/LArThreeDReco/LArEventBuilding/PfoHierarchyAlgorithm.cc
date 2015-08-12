/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/PfoHierarchyAlgorithm.cc
 * 
 *  @brief  Implementation of the pfo hierarchy algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArThreeDReco/LArEventBuilding/PfoHierarchyAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PfoHierarchyAlgorithm::PfoHierarchyAlgorithm() :
    m_halfWindowLayers(20),
    m_displayPfoInfoMap(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::SeparatePfos(const PfoInfoMap &pfoInfoMap, PfoList &assignedPfos, PfoList &unassignedPfos) const
{
    for (const PfoInfoMap::value_type mapIter : pfoInfoMap)
    {
        const PfoInfo *const pPfoInfo(mapIter.second);

        if (pPfoInfo->IsNeutrinoVertexAssociated() || pPfoInfo->GetParentPfo())
        {
            assignedPfos.insert(pPfoInfo->GetThisPfo());
        }
        else
        {
            unassignedPfos.insert(pPfoInfo->GetThisPfo());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoHierarchyAlgorithm::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_neutrinoListName,
        pPfoList));

    if (NULL == pPfoList)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "PfoHierarchyAlgorithm: pfo list unavailable." << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    for (const ParticleFlowObject *const pNeutrinoPfo : *pPfoList)
    {
        if (!LArPfoHelper::IsNeutrino(pNeutrinoPfo))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));

        PfoList daughterPfoList;
        LArPfoHelper::GetAllDownstreamPfos(pNeutrinoPfo, daughterPfoList);
        daughterPfoList.erase(pNeutrinoPfo);

        PfoInfoMap pfoInfoMap;

        try
        {
            this->GetInitialPfoInfoMap(daughterPfoList, pfoInfoMap);

            for (PfoRelationTool *const pPfoRelationTool : m_algorithmToolList)
                pPfoRelationTool->Run(this, pNeutrinoVertex, pfoInfoMap);

            this->ProcessPfoInfoMap(pNeutrinoPfo, pfoInfoMap);

            if (m_displayPfoInfoMap)
                this->DisplayPfoInfoMap(pNeutrinoPfo, pfoInfoMap);

            for (auto mapIter : pfoInfoMap) delete mapIter.second;
        }
        catch (StatusCodeException &statusCodeException)
        {
            std::cout << "PfoHierarchyAlgorithm: unable to process input pfo, " << statusCodeException.ToString() << std::endl;
            for (auto mapIter : pfoInfoMap) delete mapIter.second;
            throw statusCodeException;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::GetInitialPfoInfoMap(const PfoList &pfoList, PfoInfoMap &pfoInfoMap) const
{
    const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        PfoInfo *pPfoInfo(NULL);

        try
        {
            pPfoInfo = new PfoInfo(pPfo, m_halfWindowLayers, layerPitch);
            (void) pfoInfoMap.insert(PfoInfoMap::value_type(pPfo, pPfoInfo));
        }
        catch (StatusCodeException &)
        {
            delete pPfoInfo;
            std::cout << "Unable to calculate pfo information " << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::ProcessPfoInfoMap(const ParticleFlowObject *const pNeutrinoPfo, const PfoInfoMap &pfoInfoMap) const
{
    // Remove original neutrino->daughter links, where necessary
    const PfoList originalDaughters(pNeutrinoPfo->GetDaughterPfoList());

    for (const ParticleFlowObject *const pDaughterPfo : originalDaughters)
    {
        PfoInfoMap::const_iterator iter = pfoInfoMap.find(pDaughterPfo);

        if ((pfoInfoMap.end() == iter) || !(iter->second->IsNeutrinoVertexAssociated()))
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemovePfoParentDaughterRelationship(*this, pNeutrinoPfo, pDaughterPfo));
    }

    // Add new parent->daughter links
    for (const PfoInfoMap::value_type &iter : pfoInfoMap)
    {
        const PfoInfo *const pPfoInfo(iter.second);

        for (const ParticleFlowObject *const pDaughterPfo : pPfoInfo->GetDaughterPfoList())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPfoInfo->GetThisPfo(), pDaughterPfo));
    }

    // Deal with any remaining parent-less pfos
    for (const ParticleFlowObject *const pPfo : originalDaughters)
    {
        if (!pPfo->GetParentPfoList().empty())
            continue;

        // TODO Most appropriate decision - add as daughter of either i) nearest pfo, or ii) the neutrino (current approach)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNeutrinoPfo, pPfo));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::DisplayPfoInfoMap(const ParticleFlowObject *const pNeutrinoPfo, const PfoInfoMap &pfoInfoMap) const
{
    bool display(false);
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ));
    std::cout << "-Neutrino Pfo, nDaughters " << pNeutrinoPfo->GetDaughterPfoList().size() << ", nVertices " << pNeutrinoPfo->GetVertexList().size() << std::endl;

    for (const PfoInfoMap::value_type &iter : pfoInfoMap)
    {
        const PfoInfo *const pPfoInfo(iter.second);

        std::cout << "Pfo " << pPfoInfo->GetThisPfo() << ", vtxAssoc " << pPfoInfo->IsNeutrinoVertexAssociated()
                  << ", parent " << pPfoInfo->GetParentPfo() << ", nDaughters " << pPfoInfo->GetDaughterPfoList().size() << " (";

        for (const ParticleFlowObject *const pDaughterPfo : pPfoInfo->GetDaughterPfoList()) std::cout << pDaughterPfo << " ";
        std::cout << ") " << std::endl;

        if (pPfoInfo->IsNeutrinoVertexAssociated())
        {
            display = true;
            PfoList tempPfoList; tempPfoList.insert(pPfoInfo->GetThisPfo());
            PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &tempPfoList, "VertexPfo", RED, true, false));
        }
    }

    if (display)
    {
        PANDORA_MONITORING_API(VisualizeVertices(this->GetPandora(), &(pNeutrinoPfo->GetVertexList()), "NeutrinoVertex", ORANGE));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        display = false;
    }

    for (const PfoInfoMap::value_type &iter : pfoInfoMap)
    {
        const PfoInfo *const pPfoInfo(iter.second);

        if (!pPfoInfo->GetDaughterPfoList().empty())
        {
            display = true;
            PfoList tempPfoList; tempPfoList.insert(pPfoInfo->GetThisPfo());
            PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &tempPfoList, "ParentPfo", RED, true, false));
            PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &(pPfoInfo->GetDaughterPfoList()), "DaughterPfos", BLUE, true, false));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PfoHierarchyAlgorithm::PfoInfo::PfoInfo(const pandora::ParticleFlowObject *const pPfo, const unsigned int halfWindowLayers,
        const float layerPitch) :
    m_pThisPfo(pPfo),
    m_pCluster3D(NULL),
    m_pSlidingFitResult3D(NULL),
    m_isNeutrinoVertexAssociated(false),
    m_isInnerLayerAssociated(false),
    m_pParentPfo(NULL)
{
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusterList3D);

    if (1 != clusterList3D.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    m_pCluster3D = *(clusterList3D.begin());
    m_pSlidingFitResult3D = new ThreeDSlidingFitResult(m_pCluster3D, halfWindowLayers, layerPitch);
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoHierarchyAlgorithm::PfoInfo::PfoInfo(const PfoInfo &rhs) :
    m_pThisPfo(rhs.m_pThisPfo),
    m_pCluster3D(rhs.m_pCluster3D),
    m_pVertex3D(rhs.m_pVertex3D),
    m_pSlidingFitResult3D(NULL),
    m_isNeutrinoVertexAssociated(rhs.m_isNeutrinoVertexAssociated),
    m_isInnerLayerAssociated(rhs.m_isInnerLayerAssociated),
    m_pParentPfo(rhs.m_pParentPfo),
    m_daughterPfoList(rhs.m_daughterPfoList)
{
    if (!rhs.m_pSlidingFitResult3D)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    m_pSlidingFitResult3D = new ThreeDSlidingFitResult(m_pCluster3D, rhs.m_pSlidingFitResult3D->GetFirstFitResult().GetLayerFitHalfWindow(), 
        rhs.m_pSlidingFitResult3D->GetFirstFitResult().GetLayerPitch());
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoHierarchyAlgorithm::PfoInfo &PfoHierarchyAlgorithm::PfoInfo::operator=(const PfoInfo &rhs)
{
    if (this != &rhs)
    {
        if (!rhs.m_pSlidingFitResult3D)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        m_pThisPfo = rhs.m_pThisPfo;
        m_pCluster3D = rhs.m_pCluster3D;
        m_pVertex3D = rhs.m_pVertex3D;
        m_isNeutrinoVertexAssociated = rhs.m_isNeutrinoVertexAssociated;
        m_isInnerLayerAssociated = rhs.m_isInnerLayerAssociated;
        m_pParentPfo = rhs.m_pParentPfo;
        m_daughterPfoList = rhs.m_daughterPfoList;

        delete m_pSlidingFitResult3D;
        m_pSlidingFitResult3D = new ThreeDSlidingFitResult(m_pCluster3D, rhs.m_pSlidingFitResult3D->GetFirstFitResult().GetLayerFitHalfWindow(), 
            rhs.m_pSlidingFitResult3D->GetFirstFitResult().GetLayerPitch());
    }

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoHierarchyAlgorithm::PfoInfo::~PfoInfo()
{
    delete m_pSlidingFitResult3D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::PfoInfo::SetNeutrinoVertexAssociation(const bool isNeutrinoVertexAssociated)
{
    m_isNeutrinoVertexAssociated = isNeutrinoVertexAssociated;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::PfoInfo::SetInnerLayerAssociation(const bool isInnerLayerAssociated)
{
    m_isInnerLayerAssociated = isInnerLayerAssociated;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::PfoInfo::SetParentPfo(const pandora::ParticleFlowObject *const pParentPfo)
{
    if (m_pParentPfo)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

    m_pParentPfo = pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::PfoInfo::RemoveParentPfo()
{
    m_pParentPfo = NULL;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::PfoInfo::AddDaughterPfo(const pandora::ParticleFlowObject *const pDaughterPfo)
{
    if (!m_daughterPfoList.insert(pDaughterPfo).second)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyAlgorithm::PfoInfo::RemoveDaughterPfo(const pandora::ParticleFlowObject *const pDaughterPfo)
{
    if (!m_daughterPfoList.erase(pDaughterPfo))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoHierarchyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "PfoRelationTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        PfoRelationTool *const pPfoRelationTool(dynamic_cast<PfoRelationTool*>(*iter));

        if (NULL == pPfoRelationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pPfoRelationTool);
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "NeutrinoPfoListName", m_neutrinoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DisplayPfoInfoMap", m_displayPfoInfoMap));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
