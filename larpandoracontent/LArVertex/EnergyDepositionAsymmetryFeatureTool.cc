/**
 *  @file   larpandoracontent/LArVertex/EnergyDepositionAsymmetryFeatureTool.cc
 *
 *  @brief  Implementation of the energy deposition asymmetry feature tool class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArVertex/EnergyDepositionAsymmetryFeatureTool.h"
#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content
{

EnergyDepositionAsymmetryFeatureTool::EnergyDepositionAsymmetryFeatureTool() :
    GlobalAsymmetryFeatureTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float EnergyDepositionAsymmetryFeatureTool::CalculateAsymmetry(const bool useEnergyMetrics, const CartesianVector &vertexPosition2D,
    const ClusterVector &clusterVector, const CartesianVector &localWeightedDirectionSum) const
{
    // Project every hit onto local event axis direction and record side of the projected vtx position on which it falls
    float beforeVtxEnergy(0.f), afterVtxEnergy(0.f);
    unsigned int beforeVtxHits(0), afterVtxHits(0);

    const CartesianVector localWeightedDirection(localWeightedDirectionSum.GetUnitVector());
    const float evtProjectedVtxPos(vertexPosition2D.GetDotProduct(localWeightedDirection));

    float minBeforeProjectedPos(std::numeric_limits<float>::max());
    float maxBeforeProjectedPos(-std::numeric_limits<float>::max());

    float minAfterProjectedPos(std::numeric_limits<float>::max());
    float maxAfterProjectedPos(-std::numeric_limits<float>::max());

    for (const Cluster *const pCluster : clusterVector)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
        std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPosition);

        for (const CaloHit *const pCaloHit : caloHitVector)
        {
            if (pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection) < evtProjectedVtxPos)
            {
                minBeforeProjectedPos = std::min(minBeforeProjectedPos, pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection));
                maxBeforeProjectedPos = std::max(maxBeforeProjectedPos, pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection));

                beforeVtxEnergy += pCaloHit->GetElectromagneticEnergy();
                ++beforeVtxHits;
            }

            else
            {
                minAfterProjectedPos = std::min(minAfterProjectedPos, pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection));
                maxAfterProjectedPos = std::max(maxAfterProjectedPos, pCaloHit->GetPositionVector().GetDotProduct(localWeightedDirection));

                afterVtxEnergy += pCaloHit->GetElectromagneticEnergy();
                ++afterVtxHits;
            }
        }
    }

    // Use energy metrics if possible, otherwise fall back on hit counting.
    const unsigned int totalHits(beforeVtxHits + afterVtxHits);

    const float beforeLength(std::fabs(maxBeforeProjectedPos - minBeforeProjectedPos));
    const float afterLength(std::fabs(maxAfterProjectedPos - minAfterProjectedPos));

    const float beforeVtxEnergyDeposition(beforeLength < std::numeric_limits<float>::epsilon() ? 0 : beforeVtxEnergy / beforeLength);

    const float afterVtxEnergyDeposition(afterLength < std::numeric_limits<float>::epsilon() ? 0 : afterVtxEnergy / afterLength);

    const float totalEnergyDeposition(beforeVtxEnergyDeposition + afterVtxEnergyDeposition);

    if (useEnergyMetrics && totalEnergyDeposition > std::numeric_limits<float>::epsilon())
        return std::fabs((afterVtxEnergyDeposition - beforeVtxEnergyDeposition)) / totalEnergyDeposition;

    if (0 == totalHits)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float beforeVtxHitDeposition(beforeLength < std::numeric_limits<float>::epsilon() ? 0 : beforeVtxHits / beforeLength);

    const float afterVtxHitDeposition(afterLength < std::numeric_limits<float>::epsilon() ? 0 : afterVtxHits / afterLength);

    const float totalHitDeposition(beforeVtxHitDeposition + afterVtxHitDeposition);

    if (totalHitDeposition > std::numeric_limits<float>::epsilon())
        return std::fabs((afterVtxHitDeposition - beforeVtxHitDeposition)) / totalHitDeposition;

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EnergyDepositionAsymmetryFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return GlobalAsymmetryFeatureTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
