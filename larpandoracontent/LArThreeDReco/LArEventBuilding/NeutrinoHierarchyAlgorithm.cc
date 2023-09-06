/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoHierarchyAlgorithm.cc
 *
 *  @brief  Implementation of the neutrino hierarchy algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoHierarchyAlgorithm.h"

using namespace pandora;

namespace lar_content
{

NeutrinoHierarchyAlgorithm::NeutrinoHierarchyAlgorithm() : m_halfWindowLayers(20), m_displayPfoInfoMap(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::SeparatePfos(const PfoInfoMap &pfoInfoMap, PfoVector &assignedPfos, PfoVector &unassignedPfos) const
{
    PfoVector sortedPfos;
    for (const auto &mapEntry : pfoInfoMap)
        sortedPfos.push_back(mapEntry.first);
    std::sort(sortedPfos.begin(), sortedPfos.end(), LArPfoHelper::SortByNHits);

    for (const Pfo *const pPfo : sortedPfos)
    {
        const PfoInfo *const pPfoInfo(pfoInfoMap.at(pPfo));

        if (pPfoInfo->IsNeutrinoVertexAssociated() || pPfoInfo->GetParentPfo())
        {
            assignedPfos.push_back(pPfoInfo->GetThisPfo());
        }
        else
        {
            unassignedPfos.push_back(pPfoInfo->GetThisPfo());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoHierarchyAlgorithm::Run()
{
    const ParticleFlowObject *pNeutrinoPfo(nullptr);
    PfoList candidateDaughterPfoList;

    try
    {
        this->GetNeutrinoPfo(pNeutrinoPfo);
        this->GetCandidateDaughterPfoList(candidateDaughterPfoList);
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoHierarchyAlgorithm: required input pfos unavailable." << std::endl;

        if (STATUS_CODE_NOT_FOUND != statusCodeException.GetStatusCode())
            throw statusCodeException;

        return STATUS_CODE_SUCCESS;
    }

    PfoInfoMap pfoInfoMap;

    try
    {
        if (!pNeutrinoPfo->GetVertexList().empty())
        {
            const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
            this->GetInitialPfoInfoMap(candidateDaughterPfoList, pfoInfoMap);

            for (PfoRelationTool *const pPfoRelationTool : m_algorithmToolVector)
                pPfoRelationTool->Run(this, pNeutrinoVertex, pfoInfoMap);
        }

        this->ProcessPfoInfoMap(pNeutrinoPfo, candidateDaughterPfoList, pfoInfoMap);
    }
    catch (StatusCodeException &statusCodeException)
    {
        std::cout << "NeutrinoHierarchyAlgorithm: unable to process input neutrino pfo, " << statusCodeException.ToString() << std::endl;
    }

    if (m_displayPfoInfoMap)
        this->DisplayPfoInfoMap(pNeutrinoPfo, pfoInfoMap);

    for (auto &mapIter : pfoInfoMap)
        delete mapIter.second;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::GetNeutrinoPfo(const ParticleFlowObject *&pNeutrinoPfo) const
{
    const PfoList *pPfoList = nullptr;
    PANDORA_THROW_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_neutrinoPfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "NeutrinoHierarchyAlgorithm: unable to find pfo list " << m_neutrinoPfoListName << std::endl;

        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    // ATTN Enforces that only one pfo, of neutrino-type, be in the specified input list
    pNeutrinoPfo = ((1 == pPfoList->size()) ? *(pPfoList->begin()) : nullptr);

    if (!pNeutrinoPfo || !LArPfoHelper::IsNeutrino(pNeutrinoPfo))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::GetCandidateDaughterPfoList(PfoList &candidateDaughterPfoList) const
{
    for (const std::string &daughterPfoListName : m_daughterPfoListNames)
    {
        const PfoList *pCandidatePfoList(nullptr);

        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, daughterPfoListName, pCandidatePfoList))
        {
            candidateDaughterPfoList.insert(candidateDaughterPfoList.end(), pCandidatePfoList->begin(), pCandidatePfoList->end());
        }
        else if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        {
            std::cout << "NeutrinoHierarchyAlgorithm: unable to find pfo list " << daughterPfoListName << std::endl;
        }
    }

    if (candidateDaughterPfoList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::GetInitialPfoInfoMap(const PfoList &pfoList, PfoInfoMap &pfoInfoMap) const
{
    const float pitchU{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U)};
    const float pitchV{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V)};
    const float pitchW{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W)};
    const float pitchMax{std::max({pitchU, pitchV, pitchW})};
    const float layerPitch(pitchMax);

    for (const ParticleFlowObject *const pPfo : pfoList)
    {
        PfoInfo *pPfoInfo(nullptr);

        try
        {
            pPfoInfo = new PfoInfo(pPfo, m_halfWindowLayers, layerPitch);
            (void)pfoInfoMap.insert(PfoInfoMap::value_type(pPfo, pPfoInfo));
        }
        catch (StatusCodeException &)
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "NeutrinoHierarchyAlgorithm: Unable to calculate pfo information " << std::endl;

            delete pPfoInfo;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::ProcessPfoInfoMap(const ParticleFlowObject *const pNeutrinoPfo, const PfoList &candidateDaughterPfoList,
    PfoInfoMap &pfoInfoMap, const unsigned int callDepth) const
{
    PfoVector candidateDaughterPfoVector(candidateDaughterPfoList.begin(), candidateDaughterPfoList.end());
    std::sort(candidateDaughterPfoVector.begin(), candidateDaughterPfoVector.end(), LArPfoHelper::SortByNHits);

    // Add neutrino->primary pfo links
    for (const ParticleFlowObject *const pDaughterPfo : candidateDaughterPfoVector)
    {
        PfoInfoMap::const_iterator iter = pfoInfoMap.find(pDaughterPfo);

        if ((pfoInfoMap.end() != iter) && (iter->second->IsNeutrinoVertexAssociated()))
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNeutrinoPfo, pDaughterPfo));
    }

    // Handle exceptional case where vertex selected far from any other particle passing creation thresholds
    if ((0 == callDepth) && (0 == pNeutrinoPfo->GetNDaughterPfos()))
    {
        this->AdjustVertexAndPfoInfo(pNeutrinoPfo, candidateDaughterPfoList, pfoInfoMap);
        return this->ProcessPfoInfoMap(pNeutrinoPfo, candidateDaughterPfoList, pfoInfoMap, callDepth + 1);
    }

    // Add primary pfo->daughter pfo links
    PfoVector sortedPfos;
    for (const auto &mapEntry : pfoInfoMap)
        sortedPfos.push_back(mapEntry.first);
    std::sort(sortedPfos.begin(), sortedPfos.end(), LArPfoHelper::SortByNHits);

    for (const Pfo *const pPfo : sortedPfos)
    {
        const PfoInfo *const pPfoInfo(pfoInfoMap.at(pPfo));

        PfoVector daughterPfos(pPfoInfo->GetDaughterPfoList().begin(), pPfoInfo->GetDaughterPfoList().end());
        std::sort(daughterPfos.begin(), daughterPfos.end(), LArPfoHelper::SortByNHits);

        for (const ParticleFlowObject *const pDaughterPfo : daughterPfos)
            PANDORA_THROW_RESULT_IF(
                STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pPfoInfo->GetThisPfo(), pDaughterPfo));
    }

    // Deal with any remaining parent-less pfos
    for (const ParticleFlowObject *const pRemainingPfo : candidateDaughterPfoVector)
    {
        if (!pRemainingPfo->GetParentPfoList().empty())
            continue;

        // TODO Most appropriate decision - add as daughter of either i) nearest pfo, or ii) the neutrino (current approach)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNeutrinoPfo, pRemainingPfo));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::AdjustVertexAndPfoInfo(
    const ParticleFlowObject *const pNeutrinoPfo, const PfoList &candidateDaughterPfoList, PfoInfoMap &pfoInfoMap) const
{
    PfoVector candidateDaughterPfoVector(candidateDaughterPfoList.begin(), candidateDaughterPfoList.end());
    std::sort(candidateDaughterPfoVector.begin(), candidateDaughterPfoVector.end(), LArPfoHelper::SortByNHits);

    // TODO Consider deleting the neutrino pfo if there are no daughter pfo candidates
    if (candidateDaughterPfoVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    ClusterList daughterClusterList3D;
    LArPfoHelper::GetThreeDClusterList(candidateDaughterPfoVector.front(), daughterClusterList3D);
    daughterClusterList3D.sort(LArClusterHelper::SortByNHits);

    // If no 3D hits, must leave vertex where it was
    if (daughterClusterList3D.empty())
        return;

    const Vertex *pOldNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
    const CartesianVector newVertexPosition(LArClusterHelper::GetClosestPosition(pOldNeutrinoVertex->GetPosition(), daughterClusterList3D.front()));

    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = newVertexPosition;
    parameters.m_vertexLabel = pOldNeutrinoVertex->GetVertexLabel();
    parameters.m_vertexType = pOldNeutrinoVertex->GetVertexType();
    parameters.m_x0 = pOldNeutrinoVertex->GetX0();

    std::string neutrinoVertexListName(m_neutrinoVertexListName);
    if (neutrinoVertexListName.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Vertex>(*this, neutrinoVertexListName));

    std::string temporaryVertexListName;
    const VertexList *pTemporaryVertexList(nullptr);
    const Vertex *pNewNeutrinoVertex(nullptr);
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTemporaryVertexList, temporaryVertexListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pNewNeutrinoVertex));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, neutrinoVertexListName));

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pNeutrinoPfo, pOldNeutrinoVertex));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pNeutrinoPfo, pNewNeutrinoVertex));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Vertex>(*this, pOldNeutrinoVertex, neutrinoVertexListName));

    for (auto &mapIter : pfoInfoMap)
        delete mapIter.second;

    pfoInfoMap.clear();

    for (PfoRelationTool *const pPfoRelationTool : m_algorithmToolVector)
        pPfoRelationTool->Run(this, pNewNeutrinoVertex, pfoInfoMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::DisplayPfoInfoMap(const ParticleFlowObject *const pNeutrinoPfo, const PfoInfoMap &pfoInfoMap) const
{
    bool display(false);
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    std::cout << "-Neutrino Pfo, nDaughters " << pNeutrinoPfo->GetDaughterPfoList().size() << ", nVertices "
              << pNeutrinoPfo->GetVertexList().size() << std::endl;

    PfoVector sortedPfos;
    for (const auto &mapEntry : pfoInfoMap)
        sortedPfos.push_back(mapEntry.first);
    std::sort(sortedPfos.begin(), sortedPfos.end(), LArPfoHelper::SortByNHits);

    for (const Pfo *const pPfo : sortedPfos)
    {
        const PfoInfo *const pPfoInfo(pfoInfoMap.at(pPfo));

        std::cout << "Pfo " << pPfoInfo->GetThisPfo() << ", vtxAssoc " << pPfoInfo->IsNeutrinoVertexAssociated() << ", parent "
                  << pPfoInfo->GetParentPfo() << ", nDaughters " << pPfoInfo->GetDaughterPfoList().size() << " (";

        for (const ParticleFlowObject *const pDaughterPfo : pPfoInfo->GetDaughterPfoList())
            std::cout << pDaughterPfo << " ";
        std::cout << ") " << std::endl;

        if (pPfoInfo->IsNeutrinoVertexAssociated())
        {
            display = true;
            const PfoList tempPfoList(1, pPfoInfo->GetThisPfo());
            PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &tempPfoList, "VertexPfo", RED, true, false));
        }
    }

    if (display)
    {
        PANDORA_MONITORING_API(VisualizeVertices(this->GetPandora(), &(pNeutrinoPfo->GetVertexList()), "NeutrinoVertex", ORANGE));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        display = false;
    }

    for (const Pfo *const pPfo : sortedPfos)
    {
        const PfoInfo *const pPfoInfo(pfoInfoMap.at(pPfo));

        if (!pPfoInfo->GetDaughterPfoList().empty())
        {
            display = true;
            const PfoList tempPfoList(1, pPfoInfo->GetThisPfo());
            PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &tempPfoList, "ParentPfo", RED, true, false));
            PANDORA_MONITORING_API(
                VisualizeParticleFlowObjects(this->GetPandora(), &(pPfoInfo->GetDaughterPfoList()), "DaughterPfos", BLUE, true, false));
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

NeutrinoHierarchyAlgorithm::PfoInfo::PfoInfo(const pandora::ParticleFlowObject *const pPfo, const unsigned int halfWindowLayers, const float layerPitch) :
    m_pThisPfo(pPfo),
    m_pCluster3D(nullptr),
    m_pVertex3D(nullptr),
    m_pSlidingFitResult3D(nullptr),
    m_isNeutrinoVertexAssociated(false),
    m_isInnerLayerAssociated(false),
    m_pParentPfo(nullptr)
{
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusterList3D);

    if (1 != clusterList3D.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    m_pCluster3D = *(clusterList3D.begin());
    m_pSlidingFitResult3D = new ThreeDSlidingFitResult(m_pCluster3D, halfWindowLayers, layerPitch);

    if (m_pSlidingFitResult3D->GetMinLayer() >= m_pSlidingFitResult3D->GetMaxLayer())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

NeutrinoHierarchyAlgorithm::PfoInfo::PfoInfo(const PfoInfo &rhs) :
    m_pThisPfo(rhs.m_pThisPfo),
    m_pCluster3D(rhs.m_pCluster3D),
    m_pVertex3D(rhs.m_pVertex3D),
    m_pSlidingFitResult3D(nullptr),
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

NeutrinoHierarchyAlgorithm::PfoInfo &NeutrinoHierarchyAlgorithm::PfoInfo::operator=(const PfoInfo &rhs)
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

NeutrinoHierarchyAlgorithm::PfoInfo::~PfoInfo()
{
    delete m_pSlidingFitResult3D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::PfoInfo::SetNeutrinoVertexAssociation(const bool isNeutrinoVertexAssociated)
{
    m_isNeutrinoVertexAssociated = isNeutrinoVertexAssociated;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::PfoInfo::SetInnerLayerAssociation(const bool isInnerLayerAssociated)
{
    m_isInnerLayerAssociated = isInnerLayerAssociated;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::PfoInfo::SetParentPfo(const pandora::ParticleFlowObject *const pParentPfo)
{
    if (m_pParentPfo)
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

    m_pParentPfo = pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::PfoInfo::RemoveParentPfo()
{
    m_pParentPfo = nullptr;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::PfoInfo::AddDaughterPfo(const pandora::ParticleFlowObject *const pDaughterPfo)
{
    if (m_daughterPfoList.end() != std::find(m_daughterPfoList.begin(), m_daughterPfoList.end(), pDaughterPfo))
        throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

    m_daughterPfoList.push_back(pDaughterPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoHierarchyAlgorithm::PfoInfo::RemoveDaughterPfo(const pandora::ParticleFlowObject *const pDaughterPfo)
{
    PfoList::iterator eraseIter = std::find(m_daughterPfoList.begin(), m_daughterPfoList.end(), pDaughterPfo);

    if (m_daughterPfoList.end() == eraseIter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    m_daughterPfoList.erase(eraseIter);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoHierarchyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "PfoRelationTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        PfoRelationTool *const pPfoRelationTool(dynamic_cast<PfoRelationTool *>(*iter));

        if (!pPfoRelationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pPfoRelationTool);
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DaughterPfoListNames", m_daughterPfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoVertexListName", m_neutrinoVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitHalfWindow", m_halfWindowLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DisplayPfoInfoMap", m_displayPfoInfoMap));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
