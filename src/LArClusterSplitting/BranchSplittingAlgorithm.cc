/**
 *  @file   LArContent/src/ClusterSplitting/BranchSplittingAlgorithm.cc
 * 
 *  @brief  Implementation of the branch splitting algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterSplitting/BranchSplittingAlgorithm.h"

#include "LArHelpers/LArPointingClusterHelper.h"

using namespace pandora;

namespace lar
{

StatusCode BranchSplittingAlgorithm::FindBestSplitPosition(const LArClusterHelper::TwoDSlidingFitResult &branchSlidingFit, const LArClusterHelper::TwoDSlidingFitResult &principalSlidingFit, CartesianVector &branchStartPosition, CartesianVector &branchStartDirection) const
{
   
    for (unsigned int principalForward = 0; principalForward < 2; ++principalForward)
    {
        const CartesianVector principalVertex(1==principalForward ? principalSlidingFit.GetGlobalMinLayerPosition() : principalSlidingFit.GetGlobalMaxLayerPosition());
        const CartesianVector principalDirection(1==principalForward ? principalSlidingFit.GetGlobalMinLayerDirection() : principalSlidingFit.GetGlobalMaxLayerDirection() * -1.f);
          
        float branchDistanceSquared(std::numeric_limits<float>::max());

        try{
	    branchStartPosition = LArPointingClusterHelper::GetProjectedPosition(principalVertex, principalDirection, branchSlidingFit.GetCluster());
            branchDistanceSquared = (branchStartPosition - principalVertex).GetMagnitudeSquared();
        }
        catch (StatusCodeException &)
        {
        }

        if (branchDistanceSquared > 200.f)
	    continue;


// ClusterList tempList1, tempList2;
// Cluster* pBranchCluster = (Cluster*)(branchSlidingFit.GetCluster());
// Cluster* pReplacementCluster = (Cluster*)(principalSlidingFit.GetCluster());
// tempList1.insert(pBranchCluster);
// tempList2.insert(pReplacementCluster);
// PANDORA_MONITORING_API(SetEveDisplayParameters(false, false, -1, 1));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "BranchCluster", GREEN));
// PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "PrincipalCluster", BLUE));
// PANDORA_MONITORING_API(AddMarkerToVisualization(&principalVertex, "PrincipalVertex", BLACK, 2.75));
// PANDORA_MONITORING_API(AddMarkerToVisualization(&branchStartPosition, "SplitPosition", RED, 2.75));
// PANDORA_MONITORING_API(ViewEvent());


    }
             
    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{ 
   

    return ClusterSplittingAndExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
