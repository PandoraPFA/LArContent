/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewTrackMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the three view track matching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTrackOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewTrackMatchingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

template<typename T>
ThreeViewTrackMatchingAlgorithm<T>::ThreeViewTrackMatchingAlgorithm() :
    m_slidingFitWindow(20),
    m_minClusterCaloHits(5),
    m_minClusterLengthSquared(3.f * 3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
ThreeViewTrackMatchingAlgorithm<T>::~ThreeViewTrackMatchingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
const TwoDSlidingFitResult &ThreeViewTrackMatchingAlgorithm<T>::GetCachedSlidingFitResult(const Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeViewTrackMatchingAlgorithm<T>::UpdateForNewCluster(const Cluster *const pNewCluster)
{
    try
    {
        this->AddToSlidingFitCache(pNewCluster);
    }
    catch (const StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;

        return;
    }

    ThreeViewMatchingAlgorithm<T>::UpdateForNewCluster(pNewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeViewTrackMatchingAlgorithm<T>::UpdateUponDeletion(const Cluster *const pDeletedCluster)
{
    this->RemoveFromSlidingFitCache(pDeletedCluster);
    ThreeViewMatchingAlgorithm<T>::UpdateUponDeletion(pDeletedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeViewTrackMatchingAlgorithm<T>::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    for (const Cluster *const pCluster : *pInputClusterList)
    {
        if (!pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
            continue;

        selectedClusterList.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeViewTrackMatchingAlgorithm<T>::SetPfoParameters(const ProtoParticle &protoParticle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const
{
    // TODO Correct these placeholder parameters
    pfoParameters.m_particleId = MU_MINUS; // Track
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList.insert(pfoParameters.m_clusterList.end(), protoParticle.m_clusterListU.begin(), protoParticle.m_clusterListU.end());
    pfoParameters.m_clusterList.insert(pfoParameters.m_clusterList.end(), protoParticle.m_clusterListV.begin(), protoParticle.m_clusterListV.end());
    pfoParameters.m_clusterList.insert(pfoParameters.m_clusterList.end(), protoParticle.m_clusterListW.begin(), protoParticle.m_clusterListW.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeViewTrackMatchingAlgorithm<T>::AddToSlidingFitCache(const Cluster *const pCluster)
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeViewTrackMatchingAlgorithm<T>::RemoveFromSlidingFitCache(const Cluster *const pCluster)
{
    TwoDSlidingFitResultMap::iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() != iter)
        m_slidingFitResultMap.erase(iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeViewTrackMatchingAlgorithm<T>::PreparationStep(ClusterList &clusterList)
{
    for (ClusterList::iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; )
    {
        try
        {
            this->AddToSlidingFitCache(*iter);
            ++iter;
        }
        catch (const StatusCodeException &statusCodeException)
        {
            clusterList.erase(iter++);

            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeViewTrackMatchingAlgorithm<T>::PreparationStep()
{
    this->PreparationStep(this->m_clusterListU);
    this->PreparationStep(this->m_clusterListV);
    this->PreparationStep(this->m_clusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
void ThreeViewTrackMatchingAlgorithm<T>::TidyUp()
{
    m_slidingFitResultMap.clear();
    return ThreeViewMatchingAlgorithm<T>::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
StatusCode ThreeViewTrackMatchingAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = std::sqrt(m_minClusterLengthSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    return ThreeViewMatchingAlgorithm<T>::ReadSettings(xmlHandle);
}

template class ThreeViewTrackMatchingAlgorithm<TransverseOverlapResult>;
template class ThreeViewTrackMatchingAlgorithm<LongitudinalOverlapResult>;
template class ThreeViewTrackMatchingAlgorithm<FragmentOverlapResult>;

} // namespace lar_content
