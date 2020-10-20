/**
 *  @file   larpandoracontent/LArTrackShowerId/MvaPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the mva pfo characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArTrackShowerId/MvaPfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

template<typename T>
MvaPfoCharacterisationAlgorithm<T>::MvaPfoCharacterisationAlgorithm() :
    m_trainingSetMode(false),
    m_enableProbability(true),
    m_useThreeDInformation(true),
    m_minProbabilityCut(0.5f),
    m_minCaloHitsCut(5),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH"),
    m_writeToTree(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
template<typename T>
MvaPfoCharacterisationAlgorithm<T>::~MvaPfoCharacterisationAlgorithm()
{
    if (m_writeToTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "MvaPfoCharacterisationAlgorithm: Unable to write tree " << m_treeName << " to file " << m_fileName << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
bool MvaPfoCharacterisationAlgorithm<T>::IsClearTrack(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minCaloHitsCut)
        return false;

    const LArMvaHelper::MvaFeatureVector featureVector(LArMvaHelper::CalculateFeatures(m_featureToolVector, this, pCluster));

    if (m_trainingSetMode)
    {
        bool isTrueTrack(false);

        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
            isTrueTrack = ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId())));
        }
        catch (const StatusCodeException &) {}
        LArMvaHelper::ProduceTrainingExample(m_trainingOutputFile, isTrueTrack, featureVector);
        return isTrueTrack;
    }

    if (!m_enableProbability)
    {
        return LArMvaHelper::Classify(m_mva, featureVector);
    }
    else
    {
        return (LArMvaHelper::CalculateProbability(m_mva, featureVector) > m_minProbabilityCut);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
bool MvaPfoCharacterisationAlgorithm<T>::IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const
{
    if (!LArPfoHelper::IsThreeD(pPfo))
    {
        if (m_enableProbability)
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["TrackScore"] = -1.f;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
        }
        return (pPfo->GetParticleId() == MU_MINUS);
    }

    ClusterList wClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, wClusterList);

    //charge related features are only calculated using hits in W view
    // This won't work unless use 3D info is set to true - dev purposes only
    const PfoCharacterisationFeatureTool::FeatureToolVector &chosenFeatureToolVector(wClusterList.empty() ? m_featureToolVectorNoChargeInfo : m_featureToolVectorThreeD);
    // Purity, completeness
    // ATTN Assume your Pfos of interest are in a PfoList called myPfoList

    const LArMvaHelper::MvaFeatureVector featureVector(LArMvaHelper::CalculateFeatures(chosenFeatureToolVector, this, pPfo));

    if (m_trainingSetMode)
    {
        const PfoList myPfoList(1, pPfo);

        const MCParticleList *pMCParticleList = nullptr;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

        const CaloHitList *pCaloHitList = nullptr;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

        // Mapping target MCParticles -> truth associated Hits
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, m_primaryParameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);

        LArMCParticleHelper::MCContributionMapVector mcParticlesToGoodHitsMaps({targetMCParticleToHitsMap});

        LArMCParticleHelper::PfoContributionMap pfoToReconstructable2DHitsMap;
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(myPfoList, mcParticlesToGoodHitsMaps, pfoToReconstructable2DHitsMap, m_primaryParameters.m_foldBackHierarchy);

        LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;
        LArMCParticleHelper::MCParticleToPfoHitSharingMap mcParticleToPfoHitSharingMap;
        LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToReconstructable2DHitsMap, mcParticlesToGoodHitsMaps, pfoToMCParticleHitSharingMap, mcParticleToPfoHitSharingMap);

        const CaloHitList &allHitsInPfo(pfoToReconstructable2DHitsMap.at(pPfo));
        const int nHitsInPfoTotal(allHitsInPfo.size());
        int nHitsInBestMCParticleTotal(-1), bestMCParticlePdgCode(0);
        int nHitsSharedWithBestMCParticleTotal(-1);
        CartesianVector threeDVertexPosition(0.f, 0.f, 0.f);
        float hitsShower = 0;
        float hitsTrack = 0;
        const LArMCParticleHelper::MCParticleToSharedHitsVector &mcParticleToSharedHitsVector(pfoToMCParticleHitSharingMap.at(pPfo));

        for (const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleCaloHitListPair : mcParticleToSharedHitsVector)
        {
            const pandora::MCParticle *const pAssociatedMCParticle(mcParticleCaloHitListPair.first);
            const CaloHitList &allMCHits(targetMCParticleToHitsMap.at(pAssociatedMCParticle));
            const CaloHitList &associatedMCHits(mcParticleCaloHitListPair.second);

            if ((abs(pAssociatedMCParticle->GetParticleId()) == 11) || (pAssociatedMCParticle->GetParticleId()) == 22)
                hitsShower = hitsShower + associatedMCHits.size();

            else
                hitsTrack = hitsTrack + associatedMCHits.size();

            if (static_cast<int>(associatedMCHits.size()) > nHitsSharedWithBestMCParticleTotal)
            {
                nHitsSharedWithBestMCParticleTotal = associatedMCHits.size();
                nHitsInBestMCParticleTotal = allMCHits.size();
                bestMCParticlePdgCode = pAssociatedMCParticle->GetParticleId();
                threeDVertexPosition = pAssociatedMCParticle->GetVertex();
            }
        }

        float trackShowerHitsRatio;
        trackShowerHitsRatio = hitsTrack/(hitsTrack + hitsShower);
        int trueTrackInt = (trackShowerHitsRatio >= 0.5 ? 1 : 0);

        float xVertexPos = threeDVertexPosition.GetX();
        float yVertexPos = threeDVertexPosition.GetY();
        float zVertexPos = threeDVertexPosition.GetZ();

        const float completeness((nHitsInBestMCParticleTotal > 0) ? static_cast<float>(nHitsSharedWithBestMCParticleTotal) / static_cast<float>(nHitsInBestMCParticleTotal) : 0.f);
        const float purity((nHitsInPfoTotal > 0) ? static_cast<float>(nHitsSharedWithBestMCParticleTotal) / static_cast<float>(nHitsInPfoTotal) : 0.f);
        int pdgCode = bestMCParticlePdgCode;

        // End purity, completeness
        CaloHitList checkHitListW;
        CaloHitList checkHitListU;
        CaloHitList checkHitListV;
        CaloHitList checkHitListAll;

        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, checkHitListW);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, checkHitListU);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, checkHitListV);

        checkHitListAll.splice(checkHitListAll.end(), checkHitListW);
        checkHitListAll.splice(checkHitListAll.end(), checkHitListU);
        checkHitListAll.splice(checkHitListAll.end(), checkHitListV);

        LArMCParticleHelper::MCRelationMap mcPrimaryMap;
        LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;
        LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);

        LArMCParticleHelper::GetMCParticleToCaloHitMatches(&checkHitListAll, mcPrimaryMap, hitToMCMap, mcToTrueHitListMap);

        int showerCount(0);
        int mischaracterisedPfo(0);

        for (const CaloHit *pHit : checkHitListAll)
        {
            const MCParticle *pHitMCParticle(nullptr);

            try
            {
                pHitMCParticle = hitToMCMap.at(pHit);
            }
            catch (...) {continue;}

            if ((PHOTON == pHitMCParticle->GetParticleId()) || (E_MINUS == std::abs(pHitMCParticle->GetParticleId())))
            {
                ++showerCount;
            }
        }

        float showerProbability = (static_cast<float>(showerCount))/(static_cast<float>(hitToMCMap.size()));
        mischaracterisedPfo = ((((showerProbability < 0.5) && (trueTrackInt == 0)) || ((showerProbability > 0.5) && (trueTrackInt == 1))) ? 1 : 0);

        const bool isTrueTrack(1 == trueTrackInt);
        const bool isMainMCParticleSet(0 != pdgCode);

        if (isMainMCParticleSet)
        {
            if (completeness >= 0.0 && purity >= 0.0 && mischaracterisedPfo == 0 && (abs(xVertexPos) <= 340) && (abs(yVertexPos) <= 584) && (zVertexPos >= 200 && zVertexPos <= 1194))
            {
                std::string outputFile;
                outputFile.append(m_trainingOutputFile);
                const std::string end=((wClusterList.empty()) ? "noChargeInfo.txt" : ".txt");
                outputFile.append(end);
                LArMvaHelper::ProduceTrainingExample(outputFile, isTrueTrack, featureVector); // TODO Need this for sklearn training
            }
        }

        return isTrueTrack;
    }


    for (const LArMvaHelper::MvaFeature &featureValue : featureVector)
    {
        if (!featureValue.IsInitialized())
        {
            if (m_enableProbability)
            {
                object_creation::ParticleFlowObject::Metadata metadata;
                metadata.m_propertiesToAdd["TrackScore"] = -1.f;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
            }
            return (pPfo->GetParticleId() == MU_MINUS);
        }
    }

    //if no failures, proceed with MvaPfoCharacterisationAlgorithm classification
    if (!m_enableProbability)
    {
        return LArMvaHelper::Classify((wClusterList.empty() ? m_mvaNoChargeInfo : m_mva), featureVector);
    }
    else
    {
        const double score(LArMvaHelper::CalculateProbability((wClusterList.empty() ? m_mvaNoChargeInfo : m_mva), featureVector));
        object_creation::ParticleFlowObject::Metadata metadata;
        metadata.m_propertiesToAdd["TrackScore"] = score;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
        return (m_minProbabilityCut <= score);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode MvaPfoCharacterisationAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPrimaryGoodHits", m_primaryParameters.m_minPrimaryGoodHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsForGoodView", m_primaryParameters.m_minHitsForGoodView));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPrimaryGoodViews", m_primaryParameters.m_minPrimaryGoodViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectInputHits", m_primaryParameters.m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitSharingFraction", m_primaryParameters.m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPhotonPropagation", m_primaryParameters.m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FoldToPrimaries", m_primaryParameters.m_foldBackHierarchy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingSetMode", m_trainingSetMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsCut", m_minCaloHitsCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseThreeDInformation", m_useThreeDInformation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    // ATTN Support legacy XML configurations (note an order of precedence of XML keys exists)
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BdtFileName", m_mvaFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SvmFileName", m_mvaFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaFileName", m_mvaFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BdtName", m_mvaName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SvmName", m_mvaName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaName", m_mvaName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree)); // added by Mousam

    if (m_writeToTree)
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "OutputTree", m_treeName)); // added by Mousam
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
            "OutputFile", m_fileName)); //added by Mousam
    }

    if (m_useThreeDInformation)
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BdtFileNameNoChargeInfo", m_mvaFileNameNoChargeInfo));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SvmFileNameNoChargeInfo", m_mvaFileNameNoChargeInfo));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaFileNameNoChargeInfo", m_mvaFileNameNoChargeInfo));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BdtNameNoChargeInfo", m_mvaNameNoChargeInfo));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SvmNameNoChargeInfo", m_mvaNameNoChargeInfo));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaNameNoChargeInfo", m_mvaNameNoChargeInfo));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EnableProbability",  m_enableProbability));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinProbabilityCut", m_minProbabilityCut));

    if (m_trainingSetMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
    }
    else
    {
        if (m_mvaFileName.empty() || m_mvaName.empty())
        {
            std::cout << "MvaPfoCharacterisationAlgorithm: MvaFileName and MvaName must be set if in classification mode " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }

        const std::string fullMvaFileName(LArFileHelper::FindFileInPath(m_mvaFileName, m_filePathEnvironmentVariable));
        m_mva.Initialize(fullMvaFileName, m_mvaName);

        if (m_useThreeDInformation)
        {
            if (m_mvaFileNameNoChargeInfo.empty() || m_mvaNameNoChargeInfo.empty())
            {
                std::cout << "MvaPfoCharacterisationAlgorithm: MvaFileNameNoChargeInfo and MvaNameNoChargeInfo must be set if in classification mode for no charge info in 3D mode " << std::endl;
                return STATUS_CODE_INVALID_PARAMETER;
            }
            const std::string fullMvaFileNameNoChargeInfo(LArFileHelper::FindFileInPath(m_mvaFileNameNoChargeInfo, m_filePathEnvironmentVariable));
            m_mvaNoChargeInfo.Initialize(fullMvaFileNameNoChargeInfo, m_mvaNameNoChargeInfo);
        }
    }

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureTools", algorithmToolVector));

    if (m_useThreeDInformation)
    {
        AlgorithmToolVector algorithmToolVectorNoChargeInfo;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureToolsNoChargeInfo", algorithmToolVectorNoChargeInfo));
        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVectorThreeD));
        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVectorNoChargeInfo)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVectorNoChargeInfo));
    }
    else
    {
        for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVector));
    }

    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

template class MvaPfoCharacterisationAlgorithm<AdaBoostDecisionTree>;
template class MvaPfoCharacterisationAlgorithm<SupportVectorMachine>;

} // namespace lar_content
