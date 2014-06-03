/**
 *  @file   LArContent/src/LArThreeDReco/LArTransverseTrackMatching/TransverseTrackFragmentsTool.cc
 *
 *  @brief  Implementation of the transverse track fragments tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArTransverseTrackMatching/TransverseTrackFragmentsTool.h"

using namespace pandora;

namespace lar
{

bool TransverseTrackFragmentsTool::Run(ThreeDTransverseTrackFragmentsAlg */*pAlgorithm*/, TensorType &overlapTensor)
{
    if (PandoraSettings::ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << m_algorithmToolType << std::endl;

    bool particlesMade(false);
std::cout << " ********TransverseTrackFragmentsTool" << std::endl;
    for (TensorType::const_iterator iterU = overlapTensor.begin(), iterUEnd = overlapTensor.end(); iterU != iterUEnd; ++iterU)
    {
        if (!iterU->first->IsAvailable())
            continue;

        unsigned int nU(0), nV(0), nW(0);
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(iterU->first, true, elementList, nU, nV, nW);

std::cout << " nU " << nU << " nV " << nV << " nW " << nW << std::endl;
        for (TensorType::ElementList::const_iterator eIter = elementList.begin(); eIter != elementList.end(); ++eIter)
        {
std::cout << " msp " << eIter->GetOverlapResult().GetNMatchedSamplingPoints() << " nsp " << eIter->GetOverlapResult().GetNSamplingPoints() << " chi2Sum " << eIter->GetOverlapResult().GetChi2() << std::endl;
if (!(eIter->GetOverlapResult().GetCaloHitList().empty()))
    std::cout << " FragmentHitType " << (*(eIter->GetOverlapResult().GetCaloHitList().begin()))->GetHitType() << std::endl;
ClusterList tempList1, tempList2, tempList3;
tempList1.insert(eIter->GetClusterU());
tempList2.insert(eIter->GetClusterV());
tempList3.insert(eIter->GetClusterW());
PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "ClusterU", RED));
PANDORA_MONITORING_API(VisualizeClusters(&tempList2, "ClusterV", GREEN));
PANDORA_MONITORING_API(VisualizeClusters(&tempList3, "ClusterW", BLUE));
PANDORA_MONITORING_API(VisualizeClusters(&(eIter->GetOverlapResult().GetClusterList()), "matchedClusters", MAGENTA));
PANDORA_MONITORING_API(VisualizeCaloHits(&(eIter->GetOverlapResult().GetCaloHitList()), "matchedCaloHits", LIGHTGREEN));
PANDORA_MONITORING_API(ViewEvent());
        }
    }

    return particlesMade;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TransverseTrackFragmentsTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
