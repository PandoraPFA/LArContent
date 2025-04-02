/**
 *  @file   larpandoracontent/LArCheating/CheatingMvaPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the mva pfo characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArCheating/CheatingMvaPfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

CheatingMvaPfoCharacterisationAlgorithm::CheatingMvaPfoCharacterisationAlgorithm() :
    m_persistFeatures(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingMvaPfoCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingMvaPfoCharacterisationAlgorithm::IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const
{
    if (!LArPfoHelper::IsThreeD(pPfo))
    {
        object_creation::ParticleFlowObject::Metadata metadata;
        metadata.m_propertiesToAdd["TrackScore"] = -1.f;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
        
        return (pPfo->GetParticleId() == MU_MINUS);
    }

    // Charge related features are only calculated using hits in W view
    ClusterList wClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, wClusterList);

    const PfoCharacterisationFeatureTool::FeatureToolMap &chosenFeatureToolMap(wClusterList.empty() ? m_featureToolMapNoChargeInfo : m_featureToolMapThreeD);
    const StringVector chosenFeatureToolOrder(wClusterList.empty() ? m_algorithmToolNamesNoChargeInfo : m_algorithmToolNames);
    StringVector featureOrder;
    const LArMvaHelper::MvaFeatureMap featureMap(
        LArMvaHelper::CalculateFeatures(chosenFeatureToolOrder, chosenFeatureToolMap, featureOrder, this, pPfo));

    for (auto const &[featureKey, featureValue] : featureMap)
    {
        (void)featureKey;

        if (!featureValue.IsInitialized())
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["TrackScore"] = -1.f;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
            return (pPfo->GetParticleId() == MU_MINUS);
        }
    }

    /* M. Sotgia: This flag enable a "cheating" mode for which the pfp features are computed, and their values 
      * are saved, but the track score is forced to be the cheated version of himself
      * 
      * TODO: make this a LArCheating indipendent module, with less overhead 
      * 
     */

    bool isTrueTrack(false);
    bool isMainMCParticleSet(false);

    std::cout << "** --> DEV : [CheatingMvaPfoCharacterisationAlgorithm::IsClearTrack] "
              << "Starting the cheating of the BDT" << std::endl;

    try
    {
        // This define the truth trackScore value for the Pfo: 1 tracks, 0 showers
        const MCParticle *const pMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
        isTrueTrack = ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId())));
        isMainMCParticleSet = (pMCParticle->GetParticleId() != 0);

        std::cout << "** --> DEV : [CheatingMvaPfoCharacterisationAlgorithm<T>::IsClearTrack] "
                  << "The MC particle found is with PID = " << pMCParticle->GetParticleId() << ", PHOTON == " << PHOTON
                  << ", E_MINUS == " << E_MINUS << ", assigned TrackScore == " << static_cast<int>(isTrueTrack) << " ("
                  << (isTrueTrack ? "track" : "shower") << ")" << std::endl;
    }
    catch (const StatusCodeException &)
    {
    }

    if (isMainMCParticleSet)
    {
        // I found a MC particle (the main MC particle is set) so I can return the correct value
        // SIDE NOTE: if this association do not happen, the track score is not present in the pfp metadata,
        // so CAFAna will assign a default value of -5

        object_creation::ParticleFlowObject::Metadata metadata;
        const double score = static_cast<int>(isTrueTrack);
        metadata.m_propertiesToAdd["TrackScore"] = score;

        std::cout << "** --> DEV : [CheatingMvaPfoCharacterisationAlgorithm<T>::IsClearTrack] "
                  << "score == " << score << std::endl;

        if (m_persistFeatures)
        {
            for (auto const &[name, value] : featureMap)
            {
                metadata.m_propertiesToAdd[name] = value.Get();
                std::cout << "** --> DEV : [CheatingMvaPfoCharacterisationAlgorithm<T>::IsClearTrack] "
                          << "addinng property : " << name << " with value " << value.Get() << std::endl;
            }
        }
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
        return isTrueTrack;
    }

    return isTrueTrack;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingMvaPfoCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PersistFeatures", m_persistFeatures));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));


    LArMvaHelper::AlgorithmToolMap algorithmToolMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        LArMvaHelper::ProcessAlgorithmToolListToMap(*this, xmlHandle, "FeatureTools", m_algorithmToolNames, algorithmToolMap));

    if (m_useThreeDInformation)
    {
        // and the map for NoChargeInfo
        LArMvaHelper::AlgorithmToolMap algorithmToolMapNoChargeInfo;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            LArMvaHelper::ProcessAlgorithmToolListToMap(
                *this, xmlHandle, "FeatureToolsNoChargeInfo", m_algorithmToolNamesNoChargeInfo, algorithmToolMapNoChargeInfo));

        for (auto const &[pAlgorithmToolName, pAlgorithmTool] : algorithmToolMap)
            PANDORA_RETURN_RESULT_IF(
                STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToMap(pAlgorithmTool, pAlgorithmToolName, m_featureToolMapThreeD));

        for (auto const &[pAlgorithmToolName, pAlgorithmTool] : algorithmToolMapNoChargeInfo)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                LArMvaHelper::AddFeatureToolToMap(pAlgorithmTool, pAlgorithmToolName, m_featureToolMapNoChargeInfo));
    }
    else
    {
        for (auto const &[pAlgorithmToolName, pAlgorithmTool] : algorithmToolMap)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToMap(pAlgorithmTool, pAlgorithmToolName, m_featureToolMap));
    }

    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
