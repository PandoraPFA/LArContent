/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMergeTool.cc
 *
 *  @brief  Implementation of the delta ray merge tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewDeltaRayMergeTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

TwoViewDeltaRayMergeTool::TwoViewDeltaRayMergeTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewDeltaRayMergeTool::Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    bool mergesMade(false);

    this->MakeMerges(pAlgorithm, overlapMatrix, mergesMade);

    return mergesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

    void TwoViewDeltaRayMergeTool::MakeMerges(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix, bool &mergesMade) const
{
    bool mergeMade(true);

    while (mergeMade)
    {
        mergeMade = false;

        ClusterVector sortedKeyClusters;
        overlapMatrix.GetSortedKeyClusters(sortedKeyClusters);

        ClusterSet usedKeyClusters;
        for (const Cluster *const pKeyCluster : sortedKeyClusters)
        {
            if (usedKeyClusters.count(pKeyCluster))
                continue;

            ClusterSet checkedClusters;
            MatrixType::ElementList elementList;
            pAlgorithm->GetConnectedElements(pKeyCluster, true, elementList, checkedClusters);

            if (elementList.empty())
                continue;

            for (const MatrixType::Element &element : elementList)
	        {
                if (usedKeyClusters.count(element.GetCluster1()))
                    continue;

                usedKeyClusters.insert(element.GetCluster1());
            }

            if (this->PickOutGoodMatches(pAlgorithm, elementList))
            {
                mergeMade = true; mergesMade = true;
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewDeltaRayMergeTool::PickOutGoodMatches(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, const MatrixType::ElementList &elementList) const
{
    bool found(false);   
        
    float highestHitCount(-std::numeric_limits<float>::max()), bestChiSquared(0.f);
    MatrixType::Element bestElement(nullptr, nullptr, TrackTwoViewTopologyOverlapResult(TwoViewXOverlap(0.f, 0.f, 0.f, 0.f), PfoList(), nullptr, ClusterList(), 0.f));
    
    for (const MatrixType::Element &element : elementList)
    {            
      const Cluster *const pCluster1(element.GetCluster1()), *const pCluster2(element.GetCluster2());//, *const pCluster3(element.GetOverlapResult().GetBestMatchedCluster());

        //ATTN: Best matched cluster can be removed during pfo creation process and may not be replaceable
        //if (!pCluster3)
	//continue;

        const float chiSquared = element.GetOverlapResult().GetReducedChiSquared();
        const unsigned int hitSum(pCluster1->GetNCaloHits() + pCluster2->GetNCaloHits()); //+ (pCluster3 ? pCluster3->GetNCaloHits() : 0));

        if ((hitSum == highestHitCount) && (chiSquared < bestChiSquared))
        {
            found = true;
            bestChiSquared = chiSquared;
            highestHitCount = hitSum;
            bestElement = element;
                
            continue;
        }
            
        if (hitSum > highestHitCount)
        {
            found = true;
            bestChiSquared = chiSquared;
            highestHitCount = hitSum;
            bestElement = element;
        }
    }
    
    if (found)
    {
        pAlgorithm->CreatePfo(bestElement);
        return true;
    }
    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewDeltaRayMergeTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content
