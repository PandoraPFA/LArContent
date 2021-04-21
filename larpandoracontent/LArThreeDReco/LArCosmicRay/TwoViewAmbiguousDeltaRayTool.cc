/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewAmbiguousDeltaRayTool.cc
 *
 *  @brief  Implementation of the two view amgiuous delta ray tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/TwoViewAmbiguousDeltaRayTool.h"

using namespace pandora;

namespace lar_content
{

TwoViewAmbiguousDeltaRayTool::TwoViewAmbiguousDeltaRayTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewAmbiguousDeltaRayTool::Run(TwoViewDeltaRayMatchingAlgorithm *const pAlgorithm, MatrixType &overlapMatrix)
{
    m_pParentAlgorithm = pAlgorithm;
    
    if (PandoraContentApi::GetSettings(*m_pParentAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    this->ExamineConnectedElements(overlapMatrix);

    // ATTN: Prevent tensor tool loop running again
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewAmbiguousDeltaRayTool::ExamineConnectedElements(MatrixType &overlapMatrix) const
{
    bool particleCreated(true);

    while (particleCreated)
    {
        particleCreated = false;

        ClusterVector sortedKeyClusters;
        overlapMatrix.GetSortedKeyClusters(sortedKeyClusters);

        ClusterSet usedKeyClusters;

        for (const Cluster *const pKeyCluster : sortedKeyClusters)
        {
            if (usedKeyClusters.count(pKeyCluster))
                continue;

            MatrixType::ElementList elementList;
            overlapMatrix.GetConnectedElements(pKeyCluster, true, elementList);

            for (const MatrixType::Element &element : elementList)
                usedKeyClusters.insert(element.GetCluster1());

            if (this->PickOutGoodMatches(elementList))
            {
                particleCreated = true;
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewAmbiguousDeltaRayTool::PickOutGoodMatches(const MatrixType::ElementList &elementList) const
{
    unsigned int highestHitCount(0);
    float bestChiSquared(std::numeric_limits<float>::max());
    MatrixType::Element bestElement(nullptr, nullptr, TwoViewDeltaRayOverlapResult(TwoViewXOverlap(0.f, 0.f, 0.f, 0.f), PfoList(), nullptr, ClusterList(), 0.f));
    
    for (const MatrixType::Element &element : elementList)
    {            
        const float chiSquared = element.GetOverlapResult().GetReducedChiSquared();
        const Cluster *const pCluster1(element.GetCluster1()), *const pCluster2(element.GetCluster2());
        const unsigned int hitSum(pCluster1->GetNCaloHits() + pCluster2->GetNCaloHits());

        if ((hitSum > highestHitCount) || ((hitSum == highestHitCount) && (chiSquared < bestChiSquared)))
        {
            bestChiSquared = chiSquared;
            highestHitCount = hitSum;
            bestElement = element;
        }
    }

    if (bestElement.GetCluster1() && bestElement.GetCluster2())
    {
        m_pParentAlgorithm->CreatePfo(bestElement);
        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewAmbiguousDeltaRayTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
