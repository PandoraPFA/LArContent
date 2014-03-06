/**
 *  @file   LArContent/src/LArThreeDReco/LArTrackMatching/ClearTracksTool.cc
 * 
 *  @brief  Implementation of the clear tracks tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArThreeDReco/LArTrackMatching/ClearTracksTool.h"

using namespace pandora;

namespace lar
{

StatusCode ClearTracksTool::Run(const SlidingFitResultMap &slidingFitResultMap, TensorType &overlapTensor, ProtoParticleVector &protoParticleVector)
{
    std::cout << "ClearTracksTool::Run() " << std::endl;

    TensorType::ElementList elementList;
    overlapTensor.GetUnambiguousElements(true, elementList);

    for (TensorType::ElementList::const_iterator iter = elementList.begin(), iterEnd = elementList.end(); iter != iterEnd; ++iter)
    {
        ProtoParticle protoParticle;
        protoParticle.m_clusterListU.insert(iter->GetClusterU());
        protoParticle.m_clusterListV.insert(iter->GetClusterV());
        protoParticle.m_clusterListW.insert(iter->GetClusterW());
        protoParticleVector.push_back(protoParticle);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ClearTracksTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    // Read settings from xml file here

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
