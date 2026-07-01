/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLThreeViewClearShowersTool.cc
 *
 *  @brief  Implementation of the three view clear showers tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLThreeViewClearShowersTool.h"


using namespace pandora;

namespace lar_dl_content
{

DLThreeViewClearShowersTool::DLThreeViewClearShowersTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLThreeViewClearShowersTool::Run(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // Get connected clusters
    DLMultiViewMatchingAlgorithm::ClusterGroupVector clusterGroupVector;
    pAlgorithm->GetConnectedGroups(clusterGroupVector);    

    // Loop over connected groups
    bool madeParticles(false);
    for (const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup : clusterGroupVector)
    {
        int nU(clusterGroup.m_clustersU.size()), nV(clusterGroup.m_clustersV.size()), nW(clusterGroup.m_clustersW.size());                
        
        // Look for 1:1:1 unambiguous match
        if ((nU * nV * nW) != 1)
            continue;
            
        this->CreateClearShowers(pAlgorithm, clusterGroup);
        madeParticles = true;
    }

    return madeParticles;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLThreeViewClearShowersTool::CreateClearShowers(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup)
{
    ClusterList pfoClusters({clusterGroup.m_clustersU.front()});
    pfoClusters.push_back(clusterGroup.m_clustersV.front());
    pfoClusters.push_back(clusterGroup.m_clustersW.front());
            
    pAlgorithm->CreatePfo(pfoClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLThreeViewClearShowersTool::ReadSettings([[maybe_unused]] const TiXmlHandle xmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
