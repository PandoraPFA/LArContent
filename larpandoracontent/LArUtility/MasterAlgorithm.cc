/**
 *  @file   larpandoracontent/LArUtility/MasterAlgorithm.cc
 *
 *  @brief  Implementation of the master algorithm class.
 *
 *  $Log: $
 */

#include "Api/PandoraApi.h"

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include "larpandoracontent/LArUtility/MasterAlgorithm.h"

#include "larpandoracontent/LArContent.h"

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
    m_pCosmicRayTaggingTool(nullptr),
    m_pEventSlicingTool(nullptr),
    m_pNeutrinoIdTool(nullptr),
    m_fastWorkerInstance(nullptr),
    m_sliceNuWorkerInstance(nullptr),
    m_sliceCrWorkerInstance(nullptr)
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
            m_fastWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_fastSettingsFile);

        if (m_shouldRunNeutrinoRecoOption)
            m_sliceNuWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_nuSettingsFile);

        if (m_shouldRunNeutrinoRecoOption)
            m_sliceCrWorkerInstance = this->CreateWorkerInstance(larTPCMap, gapList, m_crSettingsFile);
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
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "Input", pCaloHitList)); // TODO

    for (const Pandora *const pCRWorker : m_crWorkerInstances)
    {
        const LArTPC &larTPC(pCRWorker->GetGeometry()->GetLArTPC());

        for (const CaloHit *const pCaloHit : *pCaloHitList)
        {
            if ((pCaloHit->GetPositionVector().GetX() < (larTPC.GetCenterX() - 0.5f * larTPC.GetWidthX())) ||
                (pCaloHit->GetPositionVector().GetX() > (larTPC.GetCenterX() + 0.5f * larTPC.GetWidthX())))
            {
                continue;
            }

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
            parameters.m_pParentAddress = pCaloHit->GetParentAddress();
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*pCRWorker, parameters));
        }

        if (m_printOverallRecoStatus)
            std::cout << "Running cosmic-ray reconstruction worker instance" << std::endl;

        PandoraApi::ProcessEvent(*pCRWorker);
        PandoraApi::Reset(*pCRWorker);
    }

    return STATUS_CODE_SUCCESS;
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
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle,
            "CosmicRayTagging", pAlgorithmTool));

        m_pCosmicRayTaggingTool = dynamic_cast<CosmicRayTaggingBaseTool*>(pAlgorithmTool);

        if (!m_pCosmicRayTaggingTool)
            return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_shouldRunSlicing, xmlHandle, "ShouldRunSlicing", m_shouldRunSlicing));

    if (m_shouldRunSlicing)
    {
        AlgorithmTool *pAlgorithmTool(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle,
            "Slicing", pAlgorithmTool));

        m_pEventSlicingTool = dynamic_cast<EventSlicingBaseTool*>(pAlgorithmTool);

        if (!m_pEventSlicingTool)
            return STATUS_CODE_INVALID_PARAMETER;
    }

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
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle,
            "NeutrinoId", pAlgorithmTool));

        m_pNeutrinoIdTool = dynamic_cast<NeutrinoIdBaseTool*>(pAlgorithmTool);

        if (!m_pNeutrinoIdTool)
            return STATUS_CODE_INVALID_PARAMETER;
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReadExternalSettings(pExternalParameters, !pExternalParameters ? InputBool() :
        pExternalParameters->m_printOverallRecoStatus, xmlHandle, "PrintOverallRecoStatus", m_printOverallRecoStatus));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CRSettingsFile", m_crSettingsFile));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "NuSettingsFile", m_nuSettingsFile));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "FastSettingsFile", m_fastSettingsFile));

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
