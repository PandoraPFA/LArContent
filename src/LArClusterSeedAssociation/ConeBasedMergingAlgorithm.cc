/**
 *  @file   LArContent/src/LArClusterSeedAssociation/ConeBasedMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the cone based merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterSeedAssociation/ConeBasedMergingAlgorithm.h"

#include "LArHelpers/LArVertexHelper.h"

using namespace pandora;

namespace lar
{

StatusCode ConeBasedMergingAlgorithm::Run()
{
    // Input lists
    const ClusterList *pSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList));

    ClusterVector seedClusterVector(pSeedClusterList->begin(), pSeedClusterList->end());
    std::sort(seedClusterVector.begin(), seedClusterVector.end(), Cluster::SortByInnerLayer);

    const ClusterList *pNonSeedClusterList = NULL;
    const StatusCode statusCode(PandoraContentApi::GetClusterList(*this, m_nonSeedClusterListName, pNonSeedClusterList));

    if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_INITIALIZED != statusCode))
        return statusCode;

    if (STATUS_CODE_NOT_INITIALIZED == statusCode)
        return STATUS_CODE_SUCCESS;

    ClusterVector nonSeedClusterVector(pNonSeedClusterList->begin(), pNonSeedClusterList->end());
    std::sort(nonSeedClusterVector.begin(), nonSeedClusterVector.end(), Cluster::SortByInnerLayer);

    // Fit cones to parent clusters, using vertex for direction finding
    ConeParametersList coneParametersList;

    for (ClusterVector::const_iterator iter = seedClusterVector.begin(), iterEnd = seedClusterVector.end(); iter != iterEnd; ++iter)
    {
        try
        {
            if (!ParticleIdHelper::IsMuonFast(*iter))
                coneParametersList.push_back(ConeParameters(*iter));
        }
        catch (StatusCodeException &)
        {
        }
    }

    // Use cone fractions to match daughters to best parents
    ClusterMergeMap clusterMergeMap;

    for (ClusterVector::const_iterator iterI = nonSeedClusterVector.begin(), iterIEnd = nonSeedClusterVector.end(); iterI != iterIEnd; ++iterI)
    {
        Cluster *pDaughterCluster = *iterI;

        if (pDaughterCluster->GetNCaloHits() < 6)
            continue;

        MergeParameters bestMergeParameters;

        for (ConeParametersList::const_iterator iterJ = coneParametersList.begin(), iterJEnd = coneParametersList.end(); iterJ != iterJEnd; ++iterJ)
        {
            const ConeParameters &parentConeParameters(*iterJ);
            MergeParameters mergeParameters(pDaughterCluster, parentConeParameters);

            if (mergeParameters.GetCosThetaMax() < 0.6f)
                continue;

            if (mergeParameters.GetConeAxisProjection() < 0.f)
                continue;

            if (mergeParameters < bestMergeParameters)
                bestMergeParameters = mergeParameters;
        }

        if (bestMergeParameters.IsInitialized())
            clusterMergeMap[bestMergeParameters.GetParentCluster()].push_back(bestMergeParameters);
    }

    // Make the cluster merges
    for (ClusterMergeMap::const_iterator iter = clusterMergeMap.begin(), iterEnd = clusterMergeMap.end(); iter != iterEnd; ++iter)
    {
        MergeParametersList cosThetaMaxList(iter->second);

        std::sort(cosThetaMaxList.begin(), cosThetaMaxList.end(), ConeBasedMergingAlgorithm::SortByCosThetaMax);
        MergeParametersList coneAxisProjectionList;
//ClusterList angleClusterList;
        for (MergeParametersList::iterator mIter = cosThetaMaxList.begin(), mIterEnd = cosThetaMaxList.end(); mIter != mIterEnd; ++mIter)
        {
            MergeParameters &mergeParameters(*mIter);

            if (mergeParameters.GetCosThetaMax() > mergeParameters.GetParentCosConeHalfAngle())
            {
                coneAxisProjectionList.push_back(mergeParameters);
            }
//            else
//            {
//angleClusterList.insert(mergeParameters.GetDaughterCluster());
//            }
        }

        std::sort(coneAxisProjectionList.begin(), coneAxisProjectionList.end(), ConeBasedMergingAlgorithm::SortByConeAxisProjection);
        ClusterList finalDaughterList;
//ClusterList distanceClusterList;
        for (MergeParametersList::iterator mIter = coneAxisProjectionList.begin(), mIterEnd = coneAxisProjectionList.end(); mIter != mIterEnd; ++mIter)
        {
            MergeParameters &mergeParameters(*mIter);

            if (mergeParameters.GetConeAxisProjection() < 2.f * mergeParameters.GetParentConeLength())
            {
                finalDaughterList.insert(mergeParameters.GetDaughterCluster());
            }
//            else
//            {
//distanceClusterList.insert(mergeParameters.GetDaughterCluster());
//            }
        }
//ClusterList pTempList;
//pTempList.insert(const_cast<Cluster*>(iter->first));
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//PandoraMonitoringApi::VisualizeClusters(&pTempList, "CBMParent", RED);
//PandoraMonitoringApi::VisualizeClusters(&finalDaughterList, "CBMDaughters", BLUE);
//PandoraMonitoringApi::VisualizeClusters(&angleClusterList, "Angle", GREEN);
//PandoraMonitoringApi::VisualizeClusters(&distanceClusterList, "Distance", ORANGE);
//PandoraMonitoringApi::ViewEvent();
        for (ClusterList::const_iterator dauIter = finalDaughterList.begin(), dauIterEnd = finalDaughterList.end(); dauIter != dauIterEnd; ++dauIter)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, iter->first, *dauIter,
                m_seedClusterListName, m_nonSeedClusterListName));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ConeBasedMergingAlgorithm::ConeParameters::ConeParameters(Cluster *pCluster) :
    m_pCluster(pCluster),
    m_direction(0.f, 0.f, 0.f),
    m_apex(0.f, 0.f, 0.f),
    m_baseCentre(0.f, 0.f, 0.f),
    m_coneLength(0.f),
    m_coneCosHalfAngle(0.f),
    m_isForwardInZ(true)
{
    // Cluster fits and direction finding
    const ClusterHelper::ClusterFitResult &clusterFitResult(pCluster->GetFitToAllHitsResult());

    if (!clusterFitResult.IsFitSuccessful())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const bool isForwardInZ(LArVertexHelper::IsForwardInZ(pCluster));
    const bool isbackwardInZ(LArVertexHelper::IsBackwardInZ(pCluster));

    if (isForwardInZ == isbackwardInZ)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    // Cone axis, apex, base-centre
    const CartesianVector &intercept(clusterFitResult.GetIntercept());
    const CartesianVector &direction(clusterFitResult.GetDirection());
    const CartesianVector innerLayerCentre(intercept - direction * (pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()) - intercept).GetMagnitude());
    const CartesianVector outerLayerCentre(intercept + direction * (pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()) - intercept).GetMagnitude());

    m_isForwardInZ = isForwardInZ;
    m_direction = (isForwardInZ) ? direction : direction * -1.f;
    m_apex = (isForwardInZ) ? innerLayerCentre : outerLayerCentre;
    m_baseCentre = (isForwardInZ) ? outerLayerCentre : innerLayerCentre;
    m_coneLength = (m_baseCentre - m_apex).GetMagnitude();

    // Find representative parent cone cos half-angle
    FloatVector cosHalfAngleValues;
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
            cosHalfAngleValues.push_back(m_direction.GetCosOpeningAngle((*hitIter)->GetPositionVector() - m_apex));
    }

    std::sort(cosHalfAngleValues.begin(), cosHalfAngleValues.end());

    if (cosHalfAngleValues.empty())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    m_coneCosHalfAngle = cosHalfAngleValues.at(0.15f * cosHalfAngleValues.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ConeBasedMergingAlgorithm::ConeParameters::GetPositionCosHalfAngle(const CartesianVector &position) const
{
    return (m_direction.GetCosOpeningAngle(position - m_apex));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ConeBasedMergingAlgorithm::ConeParameters::GetConeAxisProjection(const CartesianVector &position) const
{
    return (m_direction.GetDotProduct(position - m_apex));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ConeBasedMergingAlgorithm::MergeParameters::MergeParameters(pandora::Cluster *const pDaughterCluster, const ConeParameters &parentConeParameters) :
    m_isInitialized(true),
    m_pDaughterCluster(pDaughterCluster),
    m_pParentCluster(parentConeParameters.GetCluster()),
    m_parentConeLength(parentConeParameters.GetConeLength()),
    m_parentCosConeHalfAngle(parentConeParameters.GetConeCosHalfAngle())
{
    const CartesianVector daughterInnerCentroid(pDaughterCluster->GetCentroid(pDaughterCluster->GetInnerPseudoLayer()));
    const CartesianVector daughterOuterCentroid(pDaughterCluster->GetCentroid(pDaughterCluster->GetOuterPseudoLayer()));
    const CartesianVector daughterMidpoint(daughterInnerCentroid + (daughterOuterCentroid - daughterInnerCentroid) * 0.5f);

    m_cosThetaInner = (parentConeParameters.GetPositionCosHalfAngle(daughterInnerCentroid));
    m_cosThetaOuter = (parentConeParameters.GetPositionCosHalfAngle(daughterOuterCentroid));
    m_cosThetaMidpoint = (parentConeParameters.GetPositionCosHalfAngle(daughterMidpoint));

    const CartesianVector &referencePosition(parentConeParameters.IsForwardInZ() ? daughterInnerCentroid : daughterOuterCentroid);
    m_coneAxisProjection = (parentConeParameters.GetConeAxisProjection(referencePosition));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConeBasedMergingAlgorithm::MergeParameters::operator< (const MergeParameters &rhs) const
{
    if (m_isInitialized != rhs.m_isInitialized)
        return m_isInitialized;

    if (m_cosThetaMidpoint != rhs.m_cosThetaMidpoint)
        return (m_cosThetaMidpoint > rhs.m_cosThetaMidpoint);

    return (this->GetCosThetaMax() > rhs.GetCosThetaMax());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

bool ConeBasedMergingAlgorithm::SortByCosThetaMax(const MergeParameters &lhs, const MergeParameters &rhs)
{
    return (lhs.GetCosThetaMax() > rhs.GetCosThetaMax());
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConeBasedMergingAlgorithm::SortByConeAxisProjection(const MergeParameters &lhs, const MergeParameters &rhs)
{
    return (lhs.GetConeAxisProjection() < rhs.GetConeAxisProjection());
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeBasedMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
