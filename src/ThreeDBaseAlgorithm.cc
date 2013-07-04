/**
 *  @file   ThreeDBaseAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimension algorithm base class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "ThreeDBaseAlgorithm.h"

using namespace pandora;

namespace lar
{

ThreeDBaseAlgorithm::ThreeDBaseAlgorithm() :
    m_pInputClusterListU(NULL),
    m_pInputClusterListV(NULL),
    m_pInputClusterListW(NULL)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDBaseAlgorithm::~ThreeDBaseAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDBaseAlgorithm::Run()
{
    try
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_inputClusterListNameU, m_pInputClusterListU));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_inputClusterListNameV, m_pInputClusterListV));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_inputClusterListNameW, m_pInputClusterListW));

        this->SelectInputClusters();
        this->ModifyInputClusters();

        // Derived algorithm creates the tensor instance (necessary if different algorithms store different objects in tensor)
        this->InitializeTensor();

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

        this->TidyMemberVariables();
    }
    catch (StatusCodeException &statusCodeException)
    {
        this->TidyMemberVariables();
        throw statusCodeException;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDBaseAlgorithm::CreateThreeDParticles()
{
    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryPfoListAndSetCurrent(*this, pPfoList, pfoListName));

    for (ProtoParticleVector::const_iterator iter = m_protoParticleVector.begin(), iterEnd = m_protoParticleVector.end(); iter != iterEnd; ++iter)
    {
        // TODO - correct these placeholder parameters
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = 22;
        pfoParameters.m_charge = 0;
        pfoParameters.m_mass = 0.f;
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = CartesianVector(0., 0., 0.);
        pfoParameters.m_clusterList.insert(iter->m_clusterVectorU.begin(), iter->m_clusterVectorU.end());
        pfoParameters.m_clusterList.insert(iter->m_clusterVectorV.begin(), iter->m_clusterVectorV.end());
        pfoParameters.m_clusterList.insert(iter->m_clusterVectorW.begin(), iter->m_clusterVectorW.end());
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters));
    }

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SavePfoList(*this, m_outputPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentPfoList(*this, m_outputPfoListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDBaseAlgorithm::TidyUp()
{
    m_protoParticleVector.clear();

    m_pInputClusterListU = NULL;
    m_pInputClusterListV = NULL;
    m_pInputClusterListW = NULL;

    m_clusterVectorU.clear();
    m_clusterVectorV.clear();
    m_clusterVectorW.clear();

    this->TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDBaseAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
