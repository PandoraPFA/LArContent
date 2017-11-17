/**
 *  @file   larpandoracontent/LArControlFlow/MasterAlgorithm.cc
 *
 *  @brief  Implementation of the master algorithm class.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArContent.h"

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

using namespace pandora;

namespace lar_content
{

MasterAlgorithm::MasterAlgorithm() :
    m_shouldRunAllHitsCosmicReco(true),
    m_shouldRunStitching(true),
    m_shouldRunCosmicHitRemoval(true),
    m_shouldRunSlicing(true),
    m_shouldRunNeutrinoRecoOption(true),
    m_shouldRunCosmicRecoOption(true),
    m_shouldIdentifyNeutrinoSlice(true),
    m_printOverallRecoStatus(false),
    m_pSlicingWorkerInstance(nullptr),
    m_pSliceNuWorkerInstance(nullptr),
    m_pSliceCRWorkerInstance(nullptr),
    m_pCosmicRayTaggingTool(nullptr),
    m_pNeutrinoIdTool(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Initialize()
{
    try
    {
        const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
        const DetectorGapList &gapList(this->GetPandora().GetGeometry()->GetDetectorGapList());

        for (const LArTPCMap::value_type &mapEntry : larTPCMap)
            m_crWorkerInstances.push_back(this->CreateWorkerInstance(*(mapEntry.second), gapList, m_crSettingsFile));

        if (m_shouldRunSlicing)
            m_pSlicingWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_slicingSettingsFile);

        if (m_shouldRunNeutrinoRecoOption)
            m_pSliceNuWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_nuSettingsFile);

        if (m_shouldRunNeutrinoRecoOption)
            m_pSliceCRWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_crSettingsFile);
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cout << "MasterAlgorithm: Exception during initialization of worker instances " << statusCodeException.ToString() << std::endl;
        return statusCodeException.GetStatusCode();
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Pandora *MasterAlgorithm::CreateWorkerInstance(const LArTPC &larTPC, const DetectorGapList &gapList, const std::string &settingsFile) const
{
    // The Pandora instance
    const Pandora *const pPandora(new Pandora());
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    MultiPandoraApi::AddDaughterPandoraInstance(&(this->GetPandora()), pPandora);

    // The LArTPC
    PandoraApi::Geometry::LArTPC::Parameters larTPCParameters;
    larTPCParameters.m_larTPCName = larTPC.GetLArTPCName();
    larTPCParameters.m_centerX = larTPC.GetCenterX();
    larTPCParameters.m_centerY = larTPC.GetCenterY();
    larTPCParameters.m_centerZ = larTPC.GetCenterZ();
    larTPCParameters.m_widthX = larTPC.GetWidthX();
    larTPCParameters.m_widthY = larTPC.GetWidthY();
    larTPCParameters.m_widthZ = larTPC.GetWidthZ();
    larTPCParameters.m_wirePitchU = larTPC.GetWirePitchU();
    larTPCParameters.m_wirePitchV = larTPC.GetWirePitchV();
    larTPCParameters.m_wirePitchW = larTPC.GetWirePitchW();
    larTPCParameters.m_wireAngleU = larTPC.GetWireAngleU();
    larTPCParameters.m_wireAngleV = larTPC.GetWireAngleV();
    larTPCParameters.m_sigmaUVW = larTPC.GetSigmaUVW();
    larTPCParameters.m_isDriftInPositiveX = larTPC.IsDriftInPositiveX();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, larTPCParameters));

    // The Gaps
    for (const DetectorGap *const pGap : gapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap*>(pGap));

        if (pLineGap && (pLineGap->GetLineEndX() > (larTPC.GetCenterX() - 0.5f * larTPC.GetWidthX())) && (pLineGap->GetLineStartX() > (larTPC.GetCenterX() + 0.5f * larTPC.GetWidthX())))
        {
            PandoraApi::Geometry::LineGap::Parameters lineGapParameters;
            lineGapParameters.m_lineGapType = pLineGap->GetLineGapType();
            lineGapParameters.m_lineStartX = pLineGap->GetLineStartX();
            lineGapParameters.m_lineEndX = pLineGap->GetLineEndX();
            lineGapParameters.m_lineStartZ = pLineGap->GetLineStartZ();
            lineGapParameters.m_lineEndZ = pLineGap->GetLineEndZ();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LineGap::Create(*pPandora, lineGapParameters));
        }
    }

    // Configuration
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, settingsFile));
    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Pandora *MasterAlgorithm::CreateWorkerInstance(const LArTPCMap &larTPCMap, const DetectorGapList &gapList, const std::string &settingsFile) const
{
    // The Pandora instance
    const Pandora *const pPandora(new Pandora());
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPandora, new lar_content::LArPseudoLayerPlugin));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetLArTransformationPlugin(*pPandora, new lar_content::LArRotationalTransformationPlugin));
    MultiPandoraApi::AddDaughterPandoraInstance(&(this->GetPandora()), pPandora);

    // The Parent LArTPC
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    const bool switchViews(pFirstLArTPC->IsDriftInPositiveX());

    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    PandoraApi::Geometry::LArTPC::Parameters larTPCParameters;
    larTPCParameters.m_larTPCName = "ParentLArTPC";
    larTPCParameters.m_centerX = 0.5f * (parentMaxX + parentMinX);
    larTPCParameters.m_centerY = 0.5f * (parentMaxY + parentMinY);
    larTPCParameters.m_centerZ = 0.5f * (parentMaxZ + parentMinZ);
    larTPCParameters.m_widthX = parentMaxX - parentMinX;
    larTPCParameters.m_widthY = parentMaxY - parentMinY;
    larTPCParameters.m_widthZ = parentMaxZ - parentMinZ;
    larTPCParameters.m_wirePitchU = std::max(pFirstLArTPC->GetWirePitchU(), pFirstLArTPC->GetWirePitchV());
    larTPCParameters.m_wirePitchV = std::max(pFirstLArTPC->GetWirePitchU(), pFirstLArTPC->GetWirePitchV());
    larTPCParameters.m_wirePitchW = pFirstLArTPC->GetWirePitchW();
    larTPCParameters.m_wireAngleU = switchViews ? -pFirstLArTPC->GetWireAngleV() : pFirstLArTPC->GetWireAngleU();
    larTPCParameters.m_wireAngleV = switchViews ? -pFirstLArTPC->GetWireAngleU() : pFirstLArTPC->GetWireAngleV();
    larTPCParameters.m_sigmaUVW = pFirstLArTPC->GetSigmaUVW();
    larTPCParameters.m_isDriftInPositiveX = switchViews ? !pFirstLArTPC->IsDriftInPositiveX() : pFirstLArTPC->IsDriftInPositiveX();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPandora, larTPCParameters));

    // The Gaps
    for (const DetectorGap *const pGap : gapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap*>(pGap));

        if (pLineGap)
        {
            PandoraApi::Geometry::LineGap::Parameters lineGapParameters;
            lineGapParameters.m_lineGapType = pLineGap->GetLineGapType();
            lineGapParameters.m_lineStartX = pLineGap->GetLineStartX();
            lineGapParameters.m_lineEndX = pLineGap->GetLineEndX();
            lineGapParameters.m_lineStartZ = pLineGap->GetLineStartZ();
            lineGapParameters.m_lineEndZ = pLineGap->GetLineEndZ();
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LineGap::Create(*pPandora, lineGapParameters));
        }
    }

    // Configuration
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPandora, settingsFile));
    return pPandora;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Run()
{
    //--------------------------------------------------------------------------------------------------------------------------------------
    // CR Reconstruction
    //--------------------------------------------------------------------------------------------------------------------------------------
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "Input", pCaloHitList)); // TODO
    const CaloHitList originalHitList(*pCaloHitList);

    for (const Pandora *const pCRWorker : m_crWorkerInstances)
    {
        const LArTPC &larTPC(pCRWorker->GetGeometry()->GetLArTPC());

        for (const CaloHit *const pCaloHit : originalHitList)
        {
            // TODO Whether to truncate hits at tpc boundaries - configurable parameter (replace map size check)
            if ((this->GetPandora().GetGeometry()->GetLArTPCMap().size() == 1) ||
                ((pCaloHit->GetPositionVector().GetX() > (larTPC.GetCenterX() - 0.5f * larTPC.GetWidthX())) &&
                (pCaloHit->GetPositionVector().GetX() < (larTPC.GetCenterX() + 0.5f * larTPC.GetWidthX()))))
            {
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(pCRWorker, pCaloHit));
            }
        }

        if (m_printOverallRecoStatus)
            std::cout << "Running cosmic-ray reconstruction worker instance" << std::endl;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pCRWorker));
    }

    //--------------------------------------------------------------------------------------------------------------------------------------
    // Recreate CR worker particles in master instance
    //--------------------------------------------------------------------------------------------------------------------------------------
    StitchingInfo stitchingInfo;

    for (const Pandora *const pCRWorker : m_crWorkerInstances)
    {
        const PfoList *pCRPfos(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*pCRWorker, pCRPfos));

        PfoList newPfoList;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Recreate(*pCRPfos, newPfoList));
        //const LArTPC &larTPC(pCRWorker->GetGeometry()->GetLArTPC());
    }

    //--------------------------------------------------------------------------------------------------------------------------------------
    // Stitching
    //--------------------------------------------------------------------------------------------------------------------------------------
    const PfoList *pRecreatedCRPfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(this->GetPandora(), pRecreatedCRPfos));
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), pRecreatedCRPfos, "pRecreatedCRPfos", GREEN);
PandoraMonitoringApi::ViewEvent(this->GetPandora());

    // Care with new hit creation within tools (logic for which should be encapsulated within tools)
    //    for (StitchingTool *const pStitchingTool : m_algorithmToolVector)
    //        pStitchingTool->Run(this, stitchingInfo);

    //--------------------------------------------------------------------------------------------------------------------------------------
    // CR tagging and hit removal
    //--------------------------------------------------------------------------------------------------------------------------------------
    PfoList parentCosmicRayPfos;
    for (const Pfo *const pPfo : *pRecreatedCRPfos)
    {
        if (pPfo->GetParentPfoList().empty())
            parentCosmicRayPfos.push_back(pPfo);
    }

    PfoList ambiguousPfos;
    m_pCosmicRayTaggingTool->FindAmbiguousPfos(parentCosmicRayPfos, ambiguousPfos);

    PfoList clearCosmicRayPfos;
    for (const Pfo *const pPfo : parentCosmicRayPfos)
    {
        if (ambiguousPfos.end() == std::find(ambiguousPfos.begin(), ambiguousPfos.end(), pPfo))
            clearCosmicRayPfos.push_back(pPfo);
    }

PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &clearCosmicRayPfos, "clearCosmicRayPfos", RED);
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &ambiguousPfos, "ambiguousPfos", BLUE);
PandoraMonitoringApi::ViewEvent(this->GetPandora());

    PfoList allPfosToDelete;
    LArPfoHelper::GetAllConnectedPfos(ambiguousPfos, allPfosToDelete);

    for (const Pfo *const pPfoToDelete : allPfosToDelete)
    {
        const ClusterList clusterList(pPfoToDelete->GetClusterList());
        const VertexList vertexList(pPfoToDelete->GetVertexList());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pPfoToDelete));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &clusterList));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, &vertexList));
    }

    //--------------------------------------------------------------------------------------------------------------------------------------
    // Slicing
    //--------------------------------------------------------------------------------------------------------------------------------------
    for (const CaloHit *const pCaloHit : originalHitList)
    {
        if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;

        const HitType hitType(pCaloHit->GetHitType());
        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            continue;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSlicingWorkerInstance, pCaloHit));
    }

    if (m_printOverallRecoStatus)
        std::cout << "Running slicing worker instance" << std::endl;

    const PfoList *pSlicePfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSlicingWorkerInstance));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSlicingWorkerInstance, pSlicePfos));

    if (m_printOverallRecoStatus)
        std::cout << "Identified " << pSlicePfos->size() << " slice(s)" << std::endl;
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), pSlicePfos, "slices", BLUE);
PandoraMonitoringApi::ViewEvent(this->GetPandora());

    //--------------------------------------------------------------------------------------------------------------------------------------
    // Slice hypotheses
    //--------------------------------------------------------------------------------------------------------------------------------------
    // TODO if no slicing, just give all non-CR hits to workers
    unsigned int sliceCounter(0);
    NeutrinoIdBaseTool::SliceHypotheses nuSliceHypotheses, crSliceHypotheses;

    for (const Pfo *const pSlicePfo : *pSlicePfos)
    {
        CaloHitList caloHitList;
        LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_U, caloHitList);
        LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_V, caloHitList);
        LArPfoHelper::GetCaloHits(pSlicePfo, TPC_VIEW_W, caloHitList);

        for (const CaloHit *const pSliceCaloHit : caloHitList)
        {
            const CaloHit *const pCaloHitInMaster(static_cast<const CaloHit*>(pSliceCaloHit->GetParentAddress()));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSliceNuWorkerInstance, pCaloHitInMaster));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Copy(m_pSliceCRWorkerInstance, pCaloHitInMaster));
        }

        if (m_printOverallRecoStatus)
            std::cout << "Running slice nu worker instance " << sliceCounter << " of " << pSlicePfos->size() << std::endl;

        const PfoList *pSliceNuPfos(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSliceNuWorkerInstance));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSliceNuWorkerInstance, pSliceNuPfos));
        nuSliceHypotheses.push_back(*pSliceNuPfos);

        if (m_printOverallRecoStatus)
            std::cout << "Running slice cr worker instance " << sliceCounter << " of " << pSlicePfos->size() << std::endl;

        const PfoList *pSliceCRPfos(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*m_pSliceCRWorkerInstance));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pSliceCRWorkerInstance, pSliceCRPfos));
        crSliceHypotheses.push_back(*pSliceCRPfos);

        ++sliceCounter;
    }

    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

PfoList allSliceNuOutcomes, allSliceCROutcomes;
for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
{
    allSliceNuOutcomes.insert(allSliceNuOutcomes.end(), nuSliceHypotheses.at(sliceIndex).begin(), nuSliceHypotheses.at(sliceIndex).end());
    allSliceCROutcomes.insert(allSliceCROutcomes.end(), crSliceHypotheses.at(sliceIndex).begin(), crSliceHypotheses.at(sliceIndex).end());
}
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &allSliceNuOutcomes, "allSliceNuOutcomes", GREEN);
PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &allSliceCROutcomes, "allSliceCROutcomes", RED);
PandoraMonitoringApi::ViewEvent(this->GetPandora());

    //--------------------------------------------------------------------------------------------------------------------------------------
    // Select best slice hypotheses
    //--------------------------------------------------------------------------------------------------------------------------------------
    // Recreate just the hypotheses that are required
    if (m_printOverallRecoStatus)
        std::cout << "Select best slice hypotheses" << std::endl;    

    PfoList selectedSlicePfos;
    m_pNeutrinoIdTool->SelectOutputPfos(nuSliceHypotheses, crSliceHypotheses, selectedSlicePfos);

    PfoList newSlicePfoList;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->Recreate(selectedSlicePfos, newSlicePfoList));

    //--------------------------------------------------------------------------------------------------------------------------------------
    // Tidy up
    //--------------------------------------------------------------------------------------------------------------------------------------
    for (const Pandora *const pCRWorker : m_crWorkerInstances)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pCRWorker));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSlicingWorkerInstance));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSliceNuWorkerInstance));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*m_pSliceCRWorkerInstance));

    if (m_printOverallRecoStatus)
        std::cout << "PatRec complete" << std::endl;  

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Copy(const Pandora *const pPandora, const CaloHit *const pCaloHit) const
{
    // TODO Useful protection to ensure that only calo hits owned by the master instance are copied
    PandoraApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = pCaloHit->GetPositionVector();
    parameters.m_expectedDirection = pCaloHit->GetExpectedDirection();
    parameters.m_cellNormalVector = pCaloHit->GetCellNormalVector();
    parameters.m_cellGeometry = pCaloHit->GetCellGeometry();
    parameters.m_cellSize0 = pCaloHit->GetCellSize0();
    parameters.m_cellSize1 = pCaloHit->GetCellSize1();
    parameters.m_cellThickness = pCaloHit->GetCellThickness();
    parameters.m_nCellRadiationLengths = pCaloHit->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pCaloHit->GetNCellInteractionLengths();
    parameters.m_time = pCaloHit->GetTime();
    parameters.m_inputEnergy = pCaloHit->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pCaloHit->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pCaloHit->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pCaloHit->GetHadronicEnergy();
    parameters.m_isDigital = pCaloHit->IsDigital();
    parameters.m_hitType = pCaloHit->GetHitType();
    parameters.m_hitRegion = pCaloHit->GetHitRegion();
    parameters.m_layer = pCaloHit->GetLayer();
    parameters.m_isInOuterSamplingLayer = pCaloHit->IsInOuterSamplingLayer();
    parameters.m_pParentAddress = static_cast<const void*>(pCaloHit);
    return PandoraApi::CaloHit::Create(*pPandora, parameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Recreate(const PfoList &inputPfoList, PfoList &newPfoList) const
{
    if (inputPfoList.empty())
        return STATUS_CODE_SUCCESS;

    // TODO if no pfos in input list is primary - raise exception

    std::string clusterListName;
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListName));

    std::string vertexListName;
    const VertexList *pVertexList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    std::string pfoListName;
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (const Pfo *const pPfo : inputPfoList)
    {
        if (pPfo->GetParentPfoList().empty())
            this->Recreate(pPfo, nullptr, newPfoList);
    }

    if (!pClusterList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, "RecreatedClusters"));// TODO
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, "RecreatedClusters"));
    }

    if (!pVertexList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, "RecreatedVertices"));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, "RecreatedVertices"));
    }

    if (!pPfoList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<ParticleFlowObject>(*this, "RecreatedPfos"));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*this, "RecreatedPfos"));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::Recreate(const ParticleFlowObject *const pInputPfo, const ParticleFlowObject *const pNewParentPfo, PfoList &newPfoList) const
{
    ClusterList inputClusterList2D, inputClusterList3D, newClusterList;
    LArPfoHelper::GetTwoDClusterList(pInputPfo, inputClusterList2D);
    LArPfoHelper::GetThreeDClusterList(pInputPfo, inputClusterList3D);

    for (const Cluster *const pInputCluster : inputClusterList2D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
            newCaloHitList.push_back(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
            newIsolatedCaloHitList.push_back(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));

        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    for (const Cluster *const pInputCluster : inputClusterList3D)
    {
        CaloHitList inputCaloHitList, newCaloHitList, newIsolatedCaloHitList;
        pInputCluster->GetOrderedCaloHitList().FillCaloHitList(inputCaloHitList);

        for (const CaloHit *const pInputCaloHit : inputCaloHitList)
        {
            const CaloHit *const pWorkerParentCaloHit(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));
            const CaloHit *const pMasterParentCaloHit(static_cast<const CaloHit*>(pWorkerParentCaloHit->GetParentAddress()));
            newCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit, pMasterParentCaloHit));
        }

        for (const CaloHit *const pInputCaloHit : pInputCluster->GetIsolatedCaloHitList())
        {
            const CaloHit *const pWorkerParentCaloHit(static_cast<const CaloHit*>(pInputCaloHit->GetParentAddress()));
            const CaloHit *const pMasterParentCaloHit(static_cast<const CaloHit*>(pWorkerParentCaloHit->GetParentAddress()));
            newIsolatedCaloHitList.push_back(this->CreateCaloHit(pInputCaloHit, pMasterParentCaloHit));
        }

        if (!newCaloHitList.empty())
            newClusterList.push_back(this->CreateCluster(pInputCluster, newCaloHitList, newIsolatedCaloHitList));
    }

    VertexList newVertexList;

    for (const Vertex *const pInputVertex : pInputPfo->GetVertexList())
        newVertexList.push_back(this->CreateVertex(pInputVertex));

    const ParticleFlowObject *const pNewPfo(this->CreatePfo(pInputPfo, newClusterList, newVertexList));
    newPfoList.push_back(pNewPfo);

    if (pNewParentPfo)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pNewParentPfo, pNewPfo))

    for (const ParticleFlowObject *const pInputDaughterPfo : pInputPfo->GetDaughterPfoList())
        this->Recreate(pInputDaughterPfo, pNewPfo, newPfoList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHit *MasterAlgorithm::CreateCaloHit(const CaloHit *const pInputCaloHit, const CaloHit *const pParentCaloHit) const
{
    PandoraContentApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = pInputCaloHit->GetPositionVector();
    parameters.m_expectedDirection = pInputCaloHit->GetExpectedDirection();
    parameters.m_cellNormalVector = pInputCaloHit->GetCellNormalVector();
    parameters.m_cellGeometry = pInputCaloHit->GetCellGeometry();
    parameters.m_cellSize0 = pInputCaloHit->GetCellSize0();
    parameters.m_cellSize1 = pInputCaloHit->GetCellSize1();
    parameters.m_cellThickness = pInputCaloHit->GetCellThickness();
    parameters.m_nCellRadiationLengths = pInputCaloHit->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pInputCaloHit->GetNCellInteractionLengths();
    parameters.m_time = pInputCaloHit->GetTime();
    parameters.m_inputEnergy = pInputCaloHit->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pInputCaloHit->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pInputCaloHit->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pInputCaloHit->GetHadronicEnergy();
    parameters.m_isDigital = pInputCaloHit->IsDigital();
    parameters.m_hitType = pInputCaloHit->GetHitType();
    parameters.m_hitRegion = pInputCaloHit->GetHitRegion();
    parameters.m_layer = pInputCaloHit->GetLayer();
    parameters.m_isInOuterSamplingLayer = pInputCaloHit->IsInOuterSamplingLayer();
    parameters.m_pParentAddress = static_cast<const void*>(pParentCaloHit);

    const CaloHit *pNewCaloHit(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pNewCaloHit));

    PandoraContentApi::CaloHit::Metadata metadata;
    metadata.m_isIsolated = pInputCaloHit->IsIsolated();
    metadata.m_isPossibleMip = pInputCaloHit->IsPossibleMip();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::AlterMetadata(*this, pNewCaloHit, metadata));

    return pNewCaloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *MasterAlgorithm::CreateCluster(const Cluster *const pInputCluster, const CaloHitList &newCaloHitList,
    const CaloHitList &newIsolatedCaloHitList) const
{
    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList = newCaloHitList;
    parameters.m_isolatedCaloHitList = newIsolatedCaloHitList;

    const Cluster *pNewCluster(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));

    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = pInputCluster->GetParticleId();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pNewCluster, metadata));

    return pNewCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Vertex *MasterAlgorithm::CreateVertex(const Vertex *const pInputVertex) const
{
    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = pInputVertex->GetPosition();
    parameters.m_vertexLabel = pInputVertex->GetVertexLabel();
    parameters.m_vertexType = pInputVertex->GetVertexType();

    const Vertex *pNewVertex(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pNewVertex));

    return pNewVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *MasterAlgorithm::CreatePfo(const ParticleFlowObject *const pInputPfo, const ClusterList &newClusterList,
    const VertexList &newVertexList) const
{
    PandoraContentApi::ParticleFlowObject::Parameters parameters;
    parameters.m_particleId = pInputPfo->GetParticleId();
    parameters.m_charge = pInputPfo->GetCharge();
    parameters.m_mass = pInputPfo->GetMass();
    parameters.m_energy = pInputPfo->GetEnergy();
    parameters.m_momentum = pInputPfo->GetMomentum();
    parameters.m_clusterList = newClusterList;
    parameters.m_trackList.clear();
    parameters.m_vertexList = newVertexList;

    const ParticleFlowObject *pNewPfo(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, parameters, pNewPfo));

    return pNewPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    ExternalSteeringParameters *pExternalParameters(nullptr);

    if (this->ExternalParametersPresent())
    {
        pExternalParameters = dynamic_cast<ExternalSteeringParameters*>(this->GetExternalParameters());

        if (!pExternalParameters)
            return STATUS_CODE_FAILURE;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunAllHitsCosmicReco, xmlHandle, "ShouldRunAllHitsCosmicReco", m_shouldRunAllHitsCosmicReco));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunStitching, xmlHandle, "ShouldRunStitching", m_shouldRunStitching));

    if (m_shouldRunStitching && !m_shouldRunAllHitsCosmicReco)
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldRunStitching requires ShouldRunAllHitsCosmicReco to be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldRunStitching)
    {
        // TODO
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunCosmicHitRemoval, xmlHandle, "ShouldRunCosmicHitRemoval", m_shouldRunCosmicHitRemoval));

    if (m_shouldRunCosmicHitRemoval && !m_shouldRunAllHitsCosmicReco)
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldRunCosmicHitRemoval requires ShouldRunAllHitsCosmicReco to be true" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldRunCosmicHitRemoval)
    {
        AlgorithmTool *pAlgorithmTool(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "CosmicRayTagging", pAlgorithmTool));
        m_pCosmicRayTaggingTool = dynamic_cast<CosmicRayTaggingBaseTool*>(pAlgorithmTool);

        if (!m_pCosmicRayTaggingTool)
            return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunSlicing, xmlHandle, "ShouldRunSlicing", m_shouldRunSlicing));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunNeutrinoRecoOption, xmlHandle, "ShouldRunNeutrinoRecoOption", m_shouldRunNeutrinoRecoOption));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunCosmicRecoOption, xmlHandle, "ShouldRunCosmicRecoOption", m_shouldRunCosmicRecoOption));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldIdentifyNeutrinoSlice, xmlHandle, "ShouldIdentifyNeutrinoSlice", m_shouldIdentifyNeutrinoSlice));

    if (m_shouldIdentifyNeutrinoSlice && (!m_shouldRunSlicing || !m_shouldRunNeutrinoRecoOption || !m_shouldRunCosmicRecoOption))
    {
        std::cout << "MasterAlgorithm::ReadSettings - ShouldIdentifyNeutrinoSlice requires ShouldRunSlicing and both neutrino and cosmic reconstruction options" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    if (m_shouldIdentifyNeutrinoSlice)
    {
        AlgorithmTool *pAlgorithmTool(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "NeutrinoId", pAlgorithmTool));
        m_pNeutrinoIdTool = dynamic_cast<NeutrinoIdBaseTool*>(pAlgorithmTool);

        if (!m_pNeutrinoIdTool)
            return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_printOverallRecoStatus, xmlHandle, "PrintOverallRecoStatus", m_printOverallRecoStatus));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CRSettingsFile", m_crSettingsFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NuSettingsFile", m_nuSettingsFile));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SlicingSettingsFile", m_slicingSettingsFile));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MasterAlgorithm::ReadExternalSettings(const ExternalSteeringParameters *const pExternalParameters, const InputBool inputBool,
    const TiXmlHandle xmlHandle, const std::string &xmlTag, bool &outputBool)
{
    if (pExternalParameters && inputBool.IsInitialized())
    {
        outputBool = inputBool.Get();
    }
    else
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, xmlTag, outputBool));
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
