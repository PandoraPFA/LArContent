/**
 *  @file   larpandoracontent/LArTrackShowerId/PfoCharacterisationBaseAlgorithm.cc
 *
 *  @brief  Implementation of the pfo characterisation base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArTrackShowerId/PfoCharacterisationBaseAlgorithm.h"

#include <set>

using namespace pandora;

namespace lar_content
{

PfoCharacterisationBaseAlgorithm::PfoCharacterisationBaseAlgorithm() :
    m_updateClusterIds(true), m_postBranchAddition(false), m_useThreeDInformation(true), m_minTrackLikeViews(2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoCharacterisationBaseAlgorithm::~PfoCharacterisationBaseAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCharacterisationBaseAlgorithm::Run()
{
    PfoList tracksToShowers, showersToTracks;

    for (const std::string &pfoListName : m_inputPfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

        if (!pPfoList || pPfoList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "PfoCharacterisationBaseAlgorithm: unable to find pfo list " << pfoListName << std::endl;

            continue;
        }

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            PandoraContentApi::ParticleFlowObject::Metadata pfoMetadata;
            const bool isTrackLike(m_useThreeDInformation ? this->IsClearTrack(pPfo) : this->IsClearTrack3x2D(pPfo));

            if (isTrackLike)
            {
                pfoMetadata.m_particleId = MU_MINUS;

                if (m_showerPfoListName == pfoListName)
                    showersToTracks.push_back(pPfo);
            }
            else
            {
                pfoMetadata.m_particleId = E_MINUS;

                if (m_trackPfoListName == pfoListName)
                    tracksToShowers.push_back(pPfo);
            }

            if (pPfo->GetParticleId() != pfoMetadata.m_particleId.Get())
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, pfoMetadata));

            if (!m_updateClusterIds)
                continue;

            ClusterList twoDClusterList;
            LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);

            for (const Cluster *const pCluster : twoDClusterList)
            {
                if (pCluster->GetParticleId() == pfoMetadata.m_particleId.Get())
                    continue;

                PandoraContentApi::Cluster::Metadata clusterMetadata;
                clusterMetadata.m_particleId = pfoMetadata.m_particleId.Get();
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, clusterMetadata));
            }
        }
    }

    if (!tracksToShowers.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_trackPfoListName, m_showerPfoListName, tracksToShowers));

    if (!showersToTracks.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_showerPfoListName, m_trackPfoListName, showersToTracks));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PfoCharacterisationBaseAlgorithm::IsClearTrack3x2D(const ParticleFlowObject *const pPfo) const
{
    ClusterList twoDClusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);

    typedef std::set<pandora::HitType> HitTypeSet;
    HitTypeSet hitTypeSet;

    unsigned int nTrackLikeViews(0);
    for (const Cluster *const pCluster : twoDClusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        if (!hitTypeSet.insert(hitType).second)
            continue;

        if (this->IsClearTrack(pCluster))
            ++nTrackLikeViews;

        if (nTrackLikeViews >= m_minTrackLikeViews)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoCharacterisationBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    m_inputPfoListNames.push_back(m_trackPfoListName);

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    m_inputPfoListNames.push_back(m_showerPfoListName);

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UpdateClusterIds", m_updateClusterIds));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinTrackLikeViews", m_minTrackLikeViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseThreeDInformation", m_useThreeDInformation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
