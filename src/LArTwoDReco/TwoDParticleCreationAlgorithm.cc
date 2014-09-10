/**
 *  @file   TwoDParticleCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the two dimensional particle creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/TwoDParticleCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode TwoDParticleCreationAlgorithm::Run()
{
    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    const ClusterList *pClusterListU = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListNameU, pClusterListU));

    if (NULL != pClusterListU)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreatePFOs(pClusterListU));
    }

    const ClusterList *pClusterListV = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListNameV, pClusterListV));

    if (NULL != pClusterListV)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreatePFOs(pClusterListV));
    }

    const ClusterList *pClusterListW = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListNameW, pClusterListW));

    if (NULL != pClusterListW)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CreatePFOs(pClusterListW));
    }

    if (!pPfoList->empty())
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDParticleCreationAlgorithm::CreatePFOs(const ClusterList *const pClusterList) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (pCluster->GetNCaloHits() < m_minHitsInCluster)
            continue;

        // TODO Finalize particle id and energy measurement here
        const float clusterEnergy(pCluster->GetElectromagneticEnergy());

        if (clusterEnergy < m_minClusterEnergy)
            continue;

        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;

        // TODO Finalize particle id here
        const ParticleType particleType(PandoraContentApi::GetPlugins(*this)->GetParticleId()->IsMuon(pCluster) ? MU_MINUS : PHOTON);

        const ClusterFitResult &fitToAllHitsResult(pCluster->GetFitToAllHitsResult());

        if (!fitToAllHitsResult.IsFitSuccessful())
            continue;

        // TODO Check remaining parameters
        pfoParameters.m_particleId = particleType;
        pfoParameters.m_charge = 0;
        pfoParameters.m_mass = 0.;
        pfoParameters.m_energy = clusterEnergy;
        pfoParameters.m_momentum = CartesianVector(fitToAllHitsResult.GetDirection() * clusterEnergy);
        pfoParameters.m_clusterList.insert(pCluster);

        ParticleFlowObject *pPfo(NULL);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDParticleCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameW", m_inputClusterListNameW));

    m_minHitsInCluster = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsInCluster", m_minHitsInCluster));

    m_minClusterEnergy = 0.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterEnergy", m_minClusterEnergy));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
