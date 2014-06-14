/**
 *  @file   LArContent/src/LArThreeDReco/LArShowerFragments/DeltaRayMatchingAlgorithm.cc
 * 
 *  @brief  Implementation of the delta ray shower matching algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArMCParticleHelper.h"
#include "LArHelpers/LArThreeDHelper.h"
#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDReco/LArShowerFragments/DeltaRayMatchingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode DeltaRayMatchingAlgorithm::Run()
{
    const PfoList *pPfoList = NULL;
    const StatusCode listStatusCode(PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    if (STATUS_CODE_SUCCESS != listStatusCode)
    {
        std::cout << "DeltaRayMatchingAlgorithm: Input pfo list unavailable " << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    m_slidingFitResultMap.clear();

    ClusterVector clustersU, clustersV, clustersW;
    this->GetInputClusters(m_inputClusterListNamesU, clustersU);
    this->GetInputClusters(m_inputClusterListNamesV, clustersV);
    this->GetInputClusters(m_inputClusterListNamesW, clustersW);

    this->ThreeViewMatching(clustersU, clustersV, clustersW);

    this->TwoViewMatching(clustersU, clustersV);
    this->TwoViewMatching(clustersV, clustersW);
    this->TwoViewMatching(clustersU, clustersW);

    this->OneViewMatching(clustersU);
    this->OneViewMatching(clustersV);
    this->OneViewMatching(clustersW);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::ThreeViewMatching(const ClusterVector &clustersU, const ClusterVector &clustersV, const ClusterVector &clustersW) const
{
    for (ClusterVector::const_iterator cIterW = clustersW.begin(), cIterWEnd = clustersW.end(); cIterW != cIterWEnd; ++cIterW)
    {
        Cluster *const pClusterW = *cIterW;
        float xminW(std::numeric_limits<float>::max()), xmaxW(-std::numeric_limits<float>::max());
        LArClusterHelper::GetClusterSpanX(pClusterW, xminW, xmaxW);

        for (ClusterVector::const_iterator cIterU = clustersU.begin(), cIterUEnd = clustersU.end(); cIterU != cIterUEnd; ++cIterU)
        {
            Cluster *const pClusterU = *cIterU;
            float xminU(std::numeric_limits<float>::max()), xmaxU(-std::numeric_limits<float>::max());
            LArClusterHelper::GetClusterSpanX(pClusterU, xminU, xmaxU);

            for (ClusterVector::const_iterator cIterV = clustersV.begin(), cIterVEnd = clustersV.end(); cIterV != cIterVEnd; ++cIterV)
            {
                Cluster *const pClusterV = *cIterV;

                float xminV(std::numeric_limits<float>::max()), xmaxV(-std::numeric_limits<float>::max());
                LArClusterHelper::GetClusterSpanX(pClusterV, xminV, xmaxV);

                if (!pClusterW->IsAvailable() || !pClusterU->IsAvailable() || !pClusterV->IsAvailable())
                    continue;

                const float xmin = std::max((std::max(xminU, xminV)), xminW);
                const float xmax = std::min((std::min(xmaxU, xmaxV)), xmaxW);

                if (xmax < xmin)
                    continue;

                const float pseudoChi2(this->CompareClusterTriplet(pClusterU, pClusterV, pClusterW));

                if (pseudoChi2 > m_chi2For3ViewMatching)
                    continue;

                ParticleFlowObject *pBestPFO = NULL;
                float distance(std::numeric_limits<float>::max());
                this->FindBestCosmicPFO(pClusterU, pClusterV, pClusterW, pBestPFO, distance);

                if ((distance > m_distanceFor3ViewMatching) || (NULL == pBestPFO))
                    continue;

                ClusterList clusterList;
                clusterList.insert(pClusterU);
                clusterList.insert(pClusterV);
                clusterList.insert(pClusterW);
                this->CreateDaughterPfo(clusterList, pBestPFO);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::TwoViewMatching(const ClusterVector &clusters1, const ClusterVector &clusters2) const
{
    for (ClusterVector::const_iterator cIter1 = clusters1.begin(), cIter1End = clusters1.end(); cIter1 != cIter1End; ++cIter1)
    {
        Cluster *const pCluster1 = *cIter1;
        float xmin1(std::numeric_limits<float>::max()), xmax1(-std::numeric_limits<float>::max());
        LArClusterHelper::GetClusterSpanX(pCluster1, xmin1, xmax1);

        for (ClusterVector::const_iterator cIter2 = clusters2.begin(), cIter2End = clusters2.end(); cIter2 != cIter2End; ++cIter2)
        {
            Cluster *const pCluster2 = *cIter2;
            float xmin2(std::numeric_limits<float>::max()), xmax2(-std::numeric_limits<float>::max());
            LArClusterHelper::GetClusterSpanX(pCluster2, xmin2, xmax2);

            if (std::min(xmax1, xmax2) < std::max(xmin1, xmin2))
                continue;

            ParticleFlowObject *pBestPFO = NULL;
            float distance(std::numeric_limits<float>::max());

            if (!pCluster1->IsAvailable() || !pCluster2->IsAvailable())
                continue;

            const HitType hitType1(LArThreeDHelper::GetClusterHitType(pCluster1));
            const HitType hitType2(LArThreeDHelper::GetClusterHitType(pCluster2));

            Cluster *pClusterU((TPC_VIEW_U == hitType1) ? pCluster1 : (TPC_VIEW_U == hitType2) ? pCluster2 : NULL);
            Cluster *pClusterV((TPC_VIEW_V == hitType1) ? pCluster1 : (TPC_VIEW_V == hitType2) ? pCluster2 : NULL);
            Cluster *pClusterW((TPC_VIEW_W == hitType1) ? pCluster1 : (TPC_VIEW_W == hitType2) ? pCluster2 : NULL);
            this->FindBestCosmicPFO(pClusterU, pClusterV, pClusterW, pBestPFO, distance);

            if ((distance > m_distanceFor2ViewMatching) || (NULL == pBestPFO))
                continue;

            ClusterList clusterList;
            clusterList.insert(pCluster1);
            clusterList.insert(pCluster2);
            this->CreateDaughterPfo(clusterList, pBestPFO);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::OneViewMatching(const ClusterVector &clusters) const
{
    for (ClusterVector::const_iterator cIter = clusters.begin(), cIterEnd = clusters.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *const pCluster = *cIter;
        ParticleFlowObject *pBestPFO = NULL;
        float distance(std::numeric_limits<float>::max());

        if (!pCluster->IsAvailable())
            continue;

        const HitType hitType(LArThreeDHelper::GetClusterHitType(pCluster));

        Cluster *pClusterU((TPC_VIEW_U == hitType) ? pCluster : NULL);
        Cluster *pClusterV((TPC_VIEW_V == hitType) ? pCluster : NULL);
        Cluster *pClusterW((TPC_VIEW_W == hitType) ? pCluster : NULL);
        this->FindBestCosmicPFO(pClusterU, pClusterV, pClusterW, pBestPFO, distance);

        if ((distance > m_distanceFor1ViewMatching) || (NULL == pBestPFO))
            continue;

        ClusterList clusterList;
        clusterList.insert(pCluster);
        this->CreateDaughterPfo(clusterList, pBestPFO);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::GetInputClusters(const pandora::StringVector &clusterListNames, pandora::ClusterVector &clusterVector)
{
    for (StringVector::const_iterator iter = clusterListNames.begin(), iterEnd = clusterListNames.end(); iter != iterEnd; ++iter)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, *iter, pClusterList));

        for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
        {
            Cluster *pCluster = *cIter;

            if (!pCluster->IsAvailable() || (pCluster->GetNCaloHits() < m_minCaloHitsPerCluster))
                continue;

            clusterVector.push_back(pCluster);
            this->AddToSlidingFitCache(pCluster);
        }
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DeltaRayMatchingAlgorithm::CompareClusterTriplet(Cluster *const pClusterU, Cluster *const pClusterV, Cluster *const pClusterW) const
{
    CartesianVector minimumCoordinatesU(0.f, 0.f, 0.f), maximumCoordinatesU(0.f, 0.f, 0.f);
    LArClusterHelper::GetClusterSpanXZ(pClusterU, minimumCoordinatesU, maximumCoordinatesU);
    const float xminU(minimumCoordinatesU.GetX()), xmaxU(maximumCoordinatesU.GetX());

    CartesianVector minimumCoordinatesV(0.f, 0.f, 0.f), maximumCoordinatesV(0.f, 0.f, 0.f);
    LArClusterHelper::GetClusterSpanXZ(pClusterV, minimumCoordinatesV, maximumCoordinatesV);
    const float xminV(minimumCoordinatesV.GetX()), xmaxV(maximumCoordinatesV.GetX());

    CartesianVector minimumCoordinatesW(0.f, 0.f, 0.f), maximumCoordinatesW(0.f, 0.f, 0.f);
    LArClusterHelper::GetClusterSpanXZ(pClusterW, minimumCoordinatesW, maximumCoordinatesW);
    const float xminW(minimumCoordinatesW.GetX()), xmaxW(maximumCoordinatesW.GetX());

    const float xmin = std::max((std::max(xminU, xminV)), xminW);
    const float xmax = std::min((std::min(xmaxU, xmaxV)), xmaxW);

    if (xmin >= xmax)
        return std::numeric_limits<float>::max();

    const float x((xmin + xmax) / 2.f);
    const float u(this->GetCoordinateAtX(pClusterU, x, xmin, xmax));
    const float v(this->GetCoordinateAtX(pClusterV, x, xmin, xmax));
    const float w(this->GetCoordinateAtX(pClusterW, x, xmin, xmax));
    const float uv2w(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_V, u, v));
    const float uw2v(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_W, u, w));
    const float vw2u(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_V, TPC_VIEW_W, v, w));

    const float pseudoChi2(((u - vw2u) * (u - vw2u) + (v - uw2v) * (v - uw2v) + (w - uv2w) * (w-uv2w)) / 3.f);
    return pseudoChi2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DeltaRayMatchingAlgorithm::GetCoordinateAtX(Cluster *const pCluster, const float x, const float xmin, const float xmax) const
{
    CartesianVector fitVector(0.f, 0.f, 0.f);

    try
    {
        const TwoDSlidingFitResult &slidingFitResult(this->GetCachedSlidingFitResult(pCluster));
        slidingFitResult.GetGlobalFitPositionAtX(x, fitVector);
    }
    catch (StatusCodeException &)
    {
        fitVector.SetValues(x, 0.f, LArClusterHelper::GetAverageZ(pCluster, xmin, xmax));
    }

    return fitVector.GetZ();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::FindBestCosmicPFO(Cluster *const pClusterU, Cluster *const pClusterV,
   Cluster *const pClusterW, ParticleFlowObject* &pBestPFO, float &distanceToBestPFO) const
{
    if ((NULL == pClusterU) && (NULL == pClusterV) && (NULL == pClusterW))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    distanceToBestPFO = std::numeric_limits<float>::max();

    const PfoList *pPfoList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
    {
        ParticleFlowObject *pPfo = *iter;

        float dU(std::numeric_limits<float>::max());
        float dV(std::numeric_limits<float>::max());
        float dW(std::numeric_limits<float>::max());

        bool isSubCluster(true);

        const ClusterList &pfoClusterList(pPfo->GetClusterList());
        for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const Cluster *const pPfoCluster = *cIter;

            const HitType pfoClusterHitType(LArThreeDHelper::GetClusterHitType(pPfoCluster));

            if ((TPC_VIEW_U != pfoClusterHitType) && (TPC_VIEW_V != pfoClusterHitType) && (TPC_VIEW_W != pfoClusterHitType))
            {
                if (TPC_3D == pfoClusterHitType)
                    continue;

                std::cout << "DeltaRayMatchingAlgorithm: Encountered unexpected hit type " << std::endl;
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
            }

            if ((pfoClusterHitType == TPC_VIEW_U) && (NULL != pClusterU))
            {
                dU = std::min(dU, ClusterHelper::GetDistanceToClosestHit(pPfoCluster, pClusterU));
                isSubCluster |= this->IsSubCluster(pPfoCluster, pClusterU);
            }

            if ((pfoClusterHitType == TPC_VIEW_V) && (NULL != pClusterV))
            {
                dV = std::min(dV, ClusterHelper::GetDistanceToClosestHit(pPfoCluster, pClusterV));
                isSubCluster |= this->IsSubCluster(pPfoCluster, pClusterV);
            }

            if ((pfoClusterHitType == TPC_VIEW_W) && (NULL != pClusterW))
            {
                dW = std::min(dW, ClusterHelper::GetDistanceToClosestHit(pPfoCluster, pClusterW));
                isSubCluster |= this->IsSubCluster(pPfoCluster, pClusterW);
            }
        }

        if (!isSubCluster)
            continue;


// --- BEGIN EVENT DISPLAY ---
// ClusterList tempList1, tempList2;
// for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
// {
// const Cluster *pPfoCluster = *cIter;
// tempList1.insert((Cluster*)pPfoCluster);
// }
// if (pClusterU) tempList2.insert((Cluster*)pClusterU);
// if (pClusterV) tempList2.insert((Cluster*)pClusterV);
// if (pClusterW) tempList2.insert((Cluster*)pClusterW);
// PandoraMonitoringApi::SetEveDisplayParameters(false, DETECTOR_VIEW_XZ);
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "PFO CLUSTER", RED);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "SUB CLUSTER", BLUE);
// PandoraMonitoringApi::ViewEvent();
// --- END EVENT DISPLAY ---

        float distanceSquaredSum(0.f);

        if (NULL != pClusterU)
            distanceSquaredSum += dU * dU;

        if (NULL != pClusterV)
            distanceSquaredSum += dV * dV;

        if (NULL != pClusterW)
            distanceSquaredSum += dW * dW;

        if (distanceSquaredSum < std::numeric_limits<float>::epsilon())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const float distanceToPFO = std::sqrt(distanceSquaredSum);

        if (distanceToPFO < distanceToBestPFO)
        {
            distanceToBestPFO = distanceToPFO;
            pBestPFO = pPfo;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DeltaRayMatchingAlgorithm::IsSubCluster(const Cluster *const pPFOCluster, const Cluster *const pSubCluster) const
{
    return (LArClusterHelper::GetLengthSquared(pPFOCluster) > 2.f * LArClusterHelper::GetLengthSquared(pSubCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::CreateDaughterPfo(const ClusterList &clusterList, ParticleFlowObject *const pParentPfo) const
{
    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    // TODO - correct these placeholder parameters
    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = 22;
    pfoParameters.m_charge = 0;
    pfoParameters.m_mass = 0.f;
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0., 0., 0.);
    pfoParameters.m_clusterList = clusterList;

    ParticleFlowObject *pDaughterPfo(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pDaughterPfo));

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pParentPfo, pDaughterPfo));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const TwoDSlidingFitResult &DeltaRayMatchingAlgorithm::GetCachedSlidingFitResult(Cluster *const pCluster) const
{
    TwoDSlidingFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayMatchingAlgorithm::AddToSlidingFitCache(Cluster *const pCluster)
{
    TwoDSlidingFitResult slidingFitResult;
    LArClusterHelper::LArTwoDSlidingFit(pCluster, m_slidingFitWindow, slidingFitResult);

    if (!m_slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, slidingFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "InputPfoListName", m_inputPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputPfoListName", m_outputPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesU", m_inputClusterListNamesU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesV", m_inputClusterListNamesV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesW", m_inputClusterListNamesW));

    m_chi2For3ViewMatching = 100.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Chi2For3ViewMatching", m_chi2For3ViewMatching));

    m_distanceFor3ViewMatching = 10.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceFor3ViewMatching", m_distanceFor3ViewMatching));

    m_distanceFor2ViewMatching = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceFor2ViewMatching", m_distanceFor2ViewMatching));

    m_distanceFor1ViewMatching = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceFor1ViewMatching", m_distanceFor1ViewMatching));

    m_minCaloHitsPerCluster = 2;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsPerCluster", m_minCaloHitsPerCluster));

    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
