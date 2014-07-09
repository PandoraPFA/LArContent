/**
 *  @file   LArContent/src/LArHelpers/LArPfoHelper.cc
 *
 *  @brief  Implementation of the pfo helper class.
 *
 *  $Log: $
 */

#include "Helpers/ClusterHelper.h"
#include "Helpers/XmlHelper.h"

#include "LArHelpers/LArPfoHelper.h"
#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar
{

void LArPfoHelper::GetClusters(const PfoVector &pfoVector, const HitType &hitType, ClusterVector &clusterVector)
{
    for (PfoVector::const_iterator pIter = pfoVector.begin(), pIterEnd = pfoVector.end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *pPfo = *pIter;
        LArPfoHelper::GetClusters(pPfo, hitType, clusterVector);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPfoHelper::GetClusters(const ParticleFlowObject *pPfo, const HitType &hitType, ClusterVector &clusterVector)
{
    const ClusterList &pfoClusterList = pPfo->GetClusterList();
    for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pPfoCluster = *cIter;

        if (hitType != LArClusterHelper::GetClusterHitType(pPfoCluster))
            continue;

        clusterVector.push_back(pPfoCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

float LArPfoHelper::GetTwoDLengthSquared(const ParticleFlowObject *const pPfo)
{
    float lengthSquared(0.f);

    const ClusterList &pfoClusterList = pPfo->GetClusterList();
    for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *pPfoCluster = *cIter;

        if (TPC_3D == LArClusterHelper::GetClusterHitType(pPfoCluster))
            continue;

        lengthSquared += LArClusterHelper::GetLengthSquared(pPfoCluster);
    }

    return lengthSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetThreeDLengthSquared(const ParticleFlowObject *const pPfo)
{
    float lengthSquared(0.f);

    const ClusterList &pfoClusterList = pPfo->GetClusterList();
    for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *pPfoCluster = *cIter;

        if (TPC_3D != LArClusterHelper::GetClusterHitType(pPfoCluster))
            continue;

        lengthSquared += LArClusterHelper::GetLengthSquared(pPfoCluster);
    }

    return lengthSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetClosestDistance(const ParticleFlowObject *const pPfo, const Cluster *const pCluster)
{
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

    ClusterVector clusterVector;
    LArPfoHelper::GetClusters(pPfo, hitType, clusterVector);

    if (clusterVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    float bestDistance(std::numeric_limits<float>::max());

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        const Cluster* pPfoCluster = *iter;
        const float thisDistance(LArClusterHelper::GetClosestDistance(pCluster, pPfoCluster));

        if (thisDistance < bestDistance)
            bestDistance = thisDistance;
    }

    return bestDistance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetTwoDSeparation(const ParticleFlowObject *const pPfo1, const ParticleFlowObject *const pPfo2)
{
    ClusterVector clusterVectorU1, clusterVectorV1, clusterVectorW1;
    ClusterVector clusterVectorU2, clusterVectorV2, clusterVectorW2;

    LArPfoHelper::GetClusters(pPfo1, TPC_VIEW_U, clusterVectorU1);
    LArPfoHelper::GetClusters(pPfo1, TPC_VIEW_V, clusterVectorV1);
    LArPfoHelper::GetClusters(pPfo1, TPC_VIEW_W, clusterVectorW1);

    LArPfoHelper::GetClusters(pPfo2, TPC_VIEW_U, clusterVectorU2);
    LArPfoHelper::GetClusters(pPfo2, TPC_VIEW_V, clusterVectorV2);
    LArPfoHelper::GetClusters(pPfo2, TPC_VIEW_W, clusterVectorW2);

    float numViews(0.f);
    float distanceSquared(0.f);

    if (!clusterVectorU1.empty() && !clusterVectorU2.empty())
    {
        distanceSquared += LArClusterHelper::GetClosestDistance(clusterVectorU1, clusterVectorU2);
        numViews += 1.f;
    }

    if (!clusterVectorV1.empty() && !clusterVectorV2.empty())
    {
        distanceSquared += LArClusterHelper::GetClosestDistance(clusterVectorV1, clusterVectorV2);
        numViews += 1.f;
    }

    if (!clusterVectorW1.empty() && !clusterVectorW2.empty())
    {
        distanceSquared += LArClusterHelper::GetClosestDistance(clusterVectorW1, clusterVectorW2);
        numViews += 1.f;
    }

    if (numViews < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return std::sqrt(distanceSquared / numViews);
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

float LArPfoHelper::GetThreeDSeparation(const ParticleFlowObject *const pPfo1, const ParticleFlowObject *const pPfo2)
{
    ClusterVector clusterVector1, clusterVector2;

    LArPfoHelper::GetClusters(pPfo1, TPC_3D, clusterVector1);
    LArPfoHelper::GetClusters(pPfo2, TPC_3D, clusterVector2);

    if (clusterVector1.empty() || clusterVector2.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return LArClusterHelper::GetClosestDistance(clusterVector1, clusterVector2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArPfoHelper::SortByNHits(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs)
{
    unsigned int nHitsLhs(0);
    for (ClusterList::const_iterator iter = pLhs->GetClusterList().begin(), iterEnd = pLhs->GetClusterList().end(); iter != iterEnd; ++iter)
    {
        const Cluster *pClusterLhs = *iter;

        if (TPC_3D != LArClusterHelper::GetClusterHitType(pClusterLhs))
            continue;

        nHitsLhs += pClusterLhs->GetNCaloHits();
    }

    unsigned int nHitsRhs(0);
    for (ClusterList::const_iterator iter = pRhs->GetClusterList().begin(), iterEnd = pRhs->GetClusterList().end(); iter != iterEnd; ++iter)
    {
        const Cluster *pClusterRhs = *iter;

        if (TPC_3D != LArClusterHelper::GetClusterHitType(pClusterRhs))
            continue;

        nHitsRhs += pClusterRhs->GetNCaloHits();
    }

    if (nHitsLhs != nHitsRhs)
        return (nHitsLhs > nHitsRhs);

    return (LArPfoHelper::GetTwoDLengthSquared(pLhs) > LArPfoHelper::GetTwoDLengthSquared(pRhs));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArPfoHelper::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
