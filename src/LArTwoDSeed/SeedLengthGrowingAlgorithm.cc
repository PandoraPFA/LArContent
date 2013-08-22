/**
 *  @file   LArContent/src/LArTwoDSeed/SeedLengthGrowingAlgorithm.cc
 * 
 *  @brief  Implementation of the seed length growing algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArPointingClusterHelper.h"
#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArTwoDSeed/SeedLengthGrowingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode SeedLengthGrowingAlgorithm::Run()
{
    m_pointingClusterMap.clear();

    const ClusterList *pSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_seedClusterListName, pSeedClusterList));

    ClusterVector candidateClusters(pSeedClusterList->begin(),pSeedClusterList->end());

    const ClusterList *pNonSeedClusterList = NULL;
    const StatusCode statusCode(PandoraContentApi::GetClusterList(*this, m_nonSeedClusterListName, pNonSeedClusterList));

    if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_INITIALIZED != statusCode))
        return statusCode;

    if (STATUS_CODE_SUCCESS == statusCode)
        this->GetCandidateClusters(pNonSeedClusterList, candidateClusters);

    for (ClusterVector::const_iterator iter = candidateClusters.begin(), iterEnd = candidateClusters.end(); iter != iterEnd; ++iter)
    {
        if (!m_pointingClusterMap.insert(PointingClusterMap::value_type(*iter, LArPointingCluster(*iter))).second)
            return STATUS_CODE_FAILURE;
    }

    return SeedGrowingAlgorithm::Run();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedLengthGrowingAlgorithm::GetCandidateClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < 25.f)
            continue;

        if (pCluster->GetNCaloHits() < 25)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

SeedLengthGrowingAlgorithm::AssociationType SeedLengthGrowingAlgorithm::AreClustersAssociated(const Cluster *const pClusterSeed, const Cluster *const pCluster) const
{
    // Parent and daughter directions
    const ClusterDirection seedDirection(LArVertexHelper::GetDirectionInZ(pClusterSeed));
    const ClusterDirection clusterDirection(LArVertexHelper::GetDirectionInZ(pCluster));

    if ((DIRECTION_AMBIGUOUS == seedDirection) || (DIRECTION_UNKNOWN == seedDirection) || (DIRECTION_UNKNOWN == clusterDirection))
        return NONE;

    // Parent and daughter pointing information
    PointingClusterMap::const_iterator iterSeed = m_pointingClusterMap.find(pClusterSeed);
    PointingClusterMap::const_iterator iterCluster = m_pointingClusterMap.find(pCluster);

    if ((m_pointingClusterMap.end() == iterSeed) || (m_pointingClusterMap.end() == iterCluster))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const LArPointingCluster *const pPointingSeed(&(iterSeed->second));
    const LArPointingCluster *const pPointingCluster(&(iterCluster->second));

    // Investigate associations by considering pairs of vertices, note forbidden associations
    const bool isNodeII((BACKWARD != clusterDirection) && LArPointingClusterHelper::IsNode(pPointingSeed->GetInnerVertex().GetPosition(), pPointingCluster->GetInnerVertex()));
    const bool isNodeIO((FORWARD  != clusterDirection) && LArPointingClusterHelper::IsNode(pPointingSeed->GetInnerVertex().GetPosition(), pPointingCluster->GetOuterVertex()));
    const bool isNodeOI((BACKWARD != clusterDirection) && LArPointingClusterHelper::IsNode(pPointingSeed->GetOuterVertex().GetPosition(), pPointingCluster->GetInnerVertex()));
    const bool isNodeOO((FORWARD  != clusterDirection) && LArPointingClusterHelper::IsNode(pPointingSeed->GetOuterVertex().GetPosition(), pPointingCluster->GetOuterVertex()));

    const bool checkEmissionII((pPointingSeed->GetInnerVertex().GetPosition() - pPointingCluster->GetInnerVertex().GetPosition()).GetMagnitudeSquared() < 25.f);
    const bool checkEmissionIO((pPointingSeed->GetInnerVertex().GetPosition() - pPointingCluster->GetOuterVertex().GetPosition()).GetMagnitudeSquared() < 25.f);
    const bool checkEmissionOI((pPointingSeed->GetOuterVertex().GetPosition() - pPointingCluster->GetInnerVertex().GetPosition()).GetMagnitudeSquared() < 25.f);
    const bool checkEmissionOO((pPointingSeed->GetOuterVertex().GetPosition() - pPointingCluster->GetOuterVertex().GetPosition()).GetMagnitudeSquared() < 25.f);

    const bool isEmissionII(checkEmissionII && (BACKWARD == seedDirection) && (FORWARD  == clusterDirection) && LArPointingClusterHelper::IsEmission(pPointingSeed->GetInnerVertex().GetPosition(), pPointingCluster->GetInnerVertex()));
    const bool isEmissionIO(checkEmissionIO && (BACKWARD == seedDirection) && (BACKWARD == clusterDirection) && LArPointingClusterHelper::IsEmission(pPointingSeed->GetInnerVertex().GetPosition(), pPointingCluster->GetOuterVertex()));
    const bool isEmissionOI(checkEmissionOI && (FORWARD  == seedDirection) && (FORWARD  == clusterDirection) && LArPointingClusterHelper::IsEmission(pPointingSeed->GetOuterVertex().GetPosition(), pPointingCluster->GetInnerVertex()));
    const bool isEmissionOO(checkEmissionOO && (FORWARD  == seedDirection) && (BACKWARD == clusterDirection) && LArPointingClusterHelper::IsEmission(pPointingSeed->GetOuterVertex().GetPosition(), pPointingCluster->GetOuterVertex()));

    if (!isNodeII && !isNodeIO && !isNodeOI && !isNodeOO && !isEmissionII && !isEmissionIO && !isEmissionOI && !isEmissionOO)
        return NONE;

    // Find best association, with smallest turning angle
    float thetaII(std::acos(-1.f * pPointingSeed->GetInnerVertex().GetDirection().GetCosOpeningAngle(pPointingCluster->GetInnerVertex().GetDirection())));
    float thetaIO(std::acos(-1.f * pPointingSeed->GetInnerVertex().GetDirection().GetCosOpeningAngle(pPointingCluster->GetOuterVertex().GetDirection())));
    float thetaOI(std::acos(-1.f * pPointingSeed->GetOuterVertex().GetDirection().GetCosOpeningAngle(pPointingCluster->GetInnerVertex().GetDirection())));
    float thetaOO(std::acos(-1.f * pPointingSeed->GetOuterVertex().GetDirection().GetCosOpeningAngle(pPointingCluster->GetOuterVertex().GetDirection())));

    if (FORWARD == seedDirection)
        thetaII = M_PI - thetaII;

    if (FORWARD == seedDirection)
        thetaIO = M_PI - thetaIO;

    if (BACKWARD == seedDirection)
        thetaOI = M_PI - thetaOI;

    if (BACKWARD == seedDirection)
        thetaOO = M_PI - thetaOO;

    // Best node-related association type
    const float nodeThetaII(isNodeII ? thetaII : std::numeric_limits<float>::max());
    const float nodeThetaIO(isNodeIO ? thetaIO : std::numeric_limits<float>::max());
    const float nodeThetaOI(isNodeOI ? thetaOI : std::numeric_limits<float>::max());
    const float nodeThetaOO(isNodeOO ? thetaOO : std::numeric_limits<float>::max());
    const float bestNodeTheta(std::min(nodeThetaII, std::min(nodeThetaIO, std::min(nodeThetaOI, nodeThetaOO))));
    const AssociationType nodeAssociationType((bestNodeTheta < 0.174f) ? STRONG : (bestNodeTheta < 1.047f) ? SINGLE_ORDER : NONE);

    // Best emission-related association type
    const float emissionThetaII(isEmissionII ? thetaII : std::numeric_limits<float>::max());
    const float emissionThetaIO(isEmissionIO ? thetaIO : std::numeric_limits<float>::max());
    const float emissionThetaOI(isEmissionOI ? thetaOI : std::numeric_limits<float>::max());
    const float emissionThetaOO(isEmissionOO ? thetaOO : std::numeric_limits<float>::max());
    const float bestEmissionTheta(std::min(emissionThetaII, std::min(emissionThetaIO, std::min(emissionThetaOI, emissionThetaOO))));
    const AssociationType emissionAssociationType((bestEmissionTheta < 0.174f) ? STANDARD : NONE);

// if (std::max(nodeAssociationType, emissionAssociationType) > NONE)
// {
// std::cout << " nodeAssociationType " << nodeAssociationType << " bestNodeTheta " << bestNodeTheta << " nodeThetaII " << nodeThetaII << " nodeThetaIO " << nodeThetaIO << " nodeThetaOI " << nodeThetaOI << " nodeThetaOO " << nodeThetaOO << std::endl;
// std::cout << " emissionAssociationType " << emissionAssociationType << " bestEmissionTheta " << bestEmissionTheta << " emissionThetaII " << nodeThetaII << " emissionThetaIO " << emissionThetaIO << " emissionThetaOI " << emissionThetaOI << " emissionThetaOO " << emissionThetaOO << std::endl;
// ClusterList parent, daughter; parent.insert(const_cast<Cluster*>(pClusterSeed)); daughter.insert(const_cast<Cluster*>(pCluster));
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&parent, "parent", RED);
// PandoraMonitoringApi::VisualizeClusters(&daughter, "daughter", BLUE);
// PandoraMonitoringApi::ViewEvent();
// }

    return std::max(nodeAssociationType, emissionAssociationType);
}

} // namespace lar
