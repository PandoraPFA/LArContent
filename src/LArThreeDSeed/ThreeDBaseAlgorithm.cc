/**
 *  @file   LArContent/src/LArThreeDSeed/ThreeDBaseAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimension algorithm base class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArObjects/LArOverlapTensor.h"

#include "LArThreeDSeed/ThreeDBaseAlgorithm.h"

using namespace pandora;

namespace lar
{

template <typename T>
ThreeDBaseAlgorithm<T>::ThreeDBaseAlgorithm() :
    m_pInputClusterListU(NULL),
    m_pInputClusterListV(NULL),
    m_pInputClusterListW(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
ThreeDBaseAlgorithm<T>::~ThreeDBaseAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::SelectInputClusters()
{
    for (ClusterList::const_iterator iter = m_pInputClusterListU->begin(), iterEnd = m_pInputClusterListU->end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->IsAvailable())
            m_clusterVectorU.push_back(*iter);
    }

    for (ClusterList::const_iterator iter = m_pInputClusterListV->begin(), iterEnd = m_pInputClusterListV->end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->IsAvailable())
            m_clusterVectorV.push_back(*iter);
    }

    for (ClusterList::const_iterator iter = m_pInputClusterListW->begin(), iterEnd = m_pInputClusterListW->end(); iter != iterEnd; ++iter)
    {
        if ((*iter)->IsAvailable())
            m_clusterVectorW.push_back(*iter);
    }

    std::sort(m_clusterVectorU.begin(), m_clusterVectorU.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(m_clusterVectorV.begin(), m_clusterVectorV.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(m_clusterVectorW.begin(), m_clusterVectorW.end(), LArClusterHelper::SortByNOccupiedLayers);

//std::cout << "Clusters for 2D->3D matching " << std::endl;
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//ClusterList clusterListU; clusterListU.insert(m_clusterVectorU.begin(), m_clusterVectorU.end());
//PandoraMonitoringApi::VisualizeClusters(&clusterListU, "ClusterListU", RED);
//ClusterList clusterListV; clusterListV.insert(m_clusterVectorV.begin(), m_clusterVectorV.end());
//PandoraMonitoringApi::VisualizeClusters(&clusterListV, "ClusterListV", GREEN);
//ClusterList clusterListW; clusterListW.insert(m_clusterVectorW.begin(), m_clusterVectorW.end());
//PandoraMonitoringApi::VisualizeClusters(&clusterListW, "ClusterListW", BLUE);
//PandoraMonitoringApi::ViewEvent();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::PreparationStep()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::CreateThreeDParticles()
{
    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryPfoListAndSetCurrent(*this, pPfoList, pfoListName));

    for (typename ProtoParticleVector::const_iterator iter = m_protoParticleVector.begin(), iterEnd = m_protoParticleVector.end(); iter != iterEnd; ++iter)
    {
        // TODO - correct these placeholder parameters
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = 22;
        pfoParameters.m_charge = 0;
        pfoParameters.m_mass = 0.f;
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = CartesianVector(0., 0., 0.);
        pfoParameters.m_clusterList.insert(iter->m_clusterListU.begin(), iter->m_clusterListU.end());
        pfoParameters.m_clusterList.insert(iter->m_clusterListV.begin(), iter->m_clusterListV.end());
        pfoParameters.m_clusterList.insert(iter->m_clusterListW.begin(), iter->m_clusterListW.end());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters));
    }

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SavePfoList(*this, m_outputPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentPfoList(*this, m_outputPfoListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::UpdateTensor()
{
    ClusterList usedClusters;

    for (typename ProtoParticleVector::const_iterator iter = m_protoParticleVector.begin(), iterEnd = m_protoParticleVector.end(); iter != iterEnd; ++iter)
    {
        usedClusters.insert(iter->m_clusterListU.begin(), iter->m_clusterListU.end());
        usedClusters.insert(iter->m_clusterListV.begin(), iter->m_clusterListV.end());
        usedClusters.insert(iter->m_clusterListW.begin(), iter->m_clusterListW.end());
    }

    for (ClusterList::const_iterator iter = usedClusters.begin(), iterEnd = usedClusters.end(); iter != iterEnd; ++iter)
    {
        m_overlapTensor.RemoveCluster(*iter);
    }

    m_protoParticleVector.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm<T>::TidyUp()
{
    m_overlapTensor.Clear();
    m_protoParticleVector.clear();

    m_pInputClusterListU = NULL;
    m_pInputClusterListV = NULL;
    m_pInputClusterListW = NULL;

    m_clusterVectorU.clear();
    m_clusterVectorV.clear();
    m_clusterVectorW.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode ThreeDBaseAlgorithm<T>::Run()
{
    try
    {
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this,
            m_inputClusterListNameU, m_pInputClusterListU));
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this,
            m_inputClusterListNameV, m_pInputClusterListV));
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this,
            m_inputClusterListNameW, m_pInputClusterListW));

        if ((NULL == m_pInputClusterListU) || (NULL == m_pInputClusterListV) || (NULL == m_pInputClusterListW))
        {
            std::cout << "ThreeDBaseAlgorithm: one or more input cluster lists unavailable." << std::endl;
            throw StatusCodeException(STATUS_CODE_SUCCESS);
        }

        this->SelectInputClusters();
        this->PreparationStep();

        // Loop over selected modified input clusters and allow derived algorithm to populate tensor
        for (ClusterVector::const_iterator iterU = m_clusterVectorU.begin(), iterUEnd = m_clusterVectorU.end(); iterU != iterUEnd; ++iterU)
        {
            for (ClusterVector::const_iterator iterV = m_clusterVectorV.begin(), iterVEnd = m_clusterVectorV.end(); iterV != iterVEnd; ++iterV)
            {
                for (ClusterVector::const_iterator iterW = m_clusterVectorW.begin(), iterWEnd = m_clusterVectorW.end(); iterW != iterWEnd; ++iterW)
                    this->CalculateOverlapResult(*iterU, *iterV, *iterW);
            }
        }

        // Process results encoded in tensor
        while (this->ExamineTensor())
        {
            this->CreateThreeDParticles();
            this->UpdateTensor();
        }

        this->TidyUp();
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyUp();

        if (STATUS_CODE_SUCCESS != statusCodeException.GetStatusCode())
            throw statusCodeException;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode ThreeDBaseAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

template class ThreeDBaseAlgorithm<float>;
template class ThreeDBaseAlgorithm<TrackOverlapResult>;

} // namespace lar
