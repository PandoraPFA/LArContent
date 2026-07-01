/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLTwoViewClearShowersTool.cc
 *
 *  @brief  Implementation of the two view clear showers tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLTwoViewClearShowersTool.h"

using namespace pandora;

namespace lar_dl_content
{

DLTwoViewClearShowersTool::DLTwoViewClearShowersTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLTwoViewClearShowersTool::Run(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix)
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
        int nView0(0), nView1(0);

        for (const ClusterList &clusterList : {clusterGroup.m_clustersU, clusterGroup.m_clustersV, clusterGroup.m_clustersW})
        {
            if (clusterList.size() == 0)
                ++nView0;
            else if (clusterList.size() == 1)
                ++nView1;
        }

        if ((nView0 == 1) && (nView1 == 2))
        {
            this->CreateClearShowers(pAlgorithm, clusterGroup);
            madeParticles = true;
        }
    }

    return madeParticles;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLTwoViewClearShowersTool::CreateClearShowers(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup)
{
    ClusterList pfoClusters;

    for (const ClusterList &clusterList : {clusterGroup.m_clustersU, clusterGroup.m_clustersV, clusterGroup.m_clustersW})
        if (clusterList.size() == 1)
            pfoClusters.push_back(clusterList.front());

    pAlgorithm->CreatePfo(pfoClusters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLTwoViewClearShowersTool::ReadSettings([[maybe_unused]] const TiXmlHandle xmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
