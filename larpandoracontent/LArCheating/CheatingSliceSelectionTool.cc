/**
 *  @file   larpandoracontent/LArCheating/CheatingSliceSelectionTool.cc
 *
 *  @brief  Implementation of the cheating slice selection tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingSliceSelectionTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

using namespace pandora;

namespace lar_content
{

CheatingSliceSelectionTool::CheatingSliceSelectionTool() :
    m_maxSlices{1},
    m_threshold{-1.f},
    m_cutVariable{"completeness"}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void CheatingSliceSelectionTool::SelectSlices(const pandora::Algorithm *const /*pAlgorithm*/, const SliceVector &inputSliceVector, SliceVector &outputSliceVector)
{
    // ATTN Ensure this only runs if slicing enabled
    unsigned int sliceCounter{0};
    MetricSliceIndexMap reducedSliceVarIndexMap;

    FloatVector targetWeights(inputSliceVector.size());
    FloatVector otherWeights(inputSliceVector.size());

    // Calculate target and total weight for each slice
    for (const CaloHitList &sliceHits : inputSliceVector)
    {
        CaloHitList localHitList{};
        // ATTN Must ensure we copy the hit actually owned by master instance; access differs with/without slicing enabled
        for (const CaloHit *const pSliceCaloHit : sliceHits)
            localHitList.push_back(static_cast<const CaloHit *>(pSliceCaloHit->GetParentAddress()));

        for (const CaloHit *const pCaloHit : localHitList)
        {
            float thisTargetParticleWeight{0.f}, thisTotalWeight{0.f};
            const MCParticleWeightMap &hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());

            if (hitMCParticleWeightMap.empty())
                continue;

            MCParticleList mcParticleList;
            for (const auto &mapEntry : hitMCParticleWeightMap)
                mcParticleList.push_back(mapEntry.first);
            mcParticleList.sort(LArMCParticleHelper::SortByMomentum);

            for (const MCParticle *const pMCParticle : mcParticleList)
            {
                const float weight{hitMCParticleWeightMap.at(pMCParticle)};
                const MCParticle *const pParentMCParticle{LArMCParticleHelper::GetParentMCParticle(pMCParticle)};

                if (this->IsTarget(pParentMCParticle))
                    thisTargetParticleWeight += weight;

                thisTotalWeight += weight;
            }

            // ATTN normalise arbitrary input weights at this point
            if (thisTotalWeight > std::numeric_limits<float>::epsilon())
            {
                thisTargetParticleWeight *= 1.f / thisTotalWeight;
                thisTotalWeight = 1.f;
            }
            else
            {
                thisTargetParticleWeight = 0.f;
                thisTotalWeight = 0.f;
            }
            targetWeights[sliceCounter] += thisTargetParticleWeight;
            otherWeights[sliceCounter] += (1. - thisTargetParticleWeight);
        }

        ++sliceCounter;
    }

    float totalTargetWeight{0.f};
    for (const float weight : targetWeights)
        totalTargetWeight += weight;

    // Add slices passing cut threshold to map
    for (unsigned int i = 0; i < targetWeights.size(); ++i)
    {
        const float sliceWeight = targetWeights[i] + otherWeights[i];
        const float completeness{totalTargetWeight > std::numeric_limits<float>::epsilon() ? targetWeights[i] / totalTargetWeight : 0.f};
        const float purity{sliceWeight > std::numeric_limits<float>::epsilon() ? targetWeights[i] / sliceWeight : 0.f};
        // Already checked correctness of variable in ReadSettings
        if (m_cutVariable == "completeness")
        {
            if (completeness >= m_threshold)
                reducedSliceVarIndexMap.emplace(completeness, i);
        }
        else if (m_cutVariable == "purity")
        {
            if (purity >= m_threshold)
                reducedSliceVarIndexMap.emplace(purity, i);
        }
    }
    // Select the best m_maxSlices slices - prefix increment ensures all slices retained if m_maxSlices == 0
    std::vector<int> reducedSliceIndices{};
    int i = 0;
    for (const auto [cutVariable, index] : reducedSliceVarIndexMap)
    {                      // ATTN: Map is sorted on cut variable from max to min
        (void)cutVariable; // GCC 7 support, versions 8+ do not need this
        reducedSliceIndices.push_back(index);
        outputSliceVector.push_back(inputSliceVector[index]);
        if (++i == m_maxSlices)
            break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingSliceSelectionTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSlices", m_maxSlices));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Threshold", m_threshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CutVariable", m_cutVariable));
    std::transform(m_cutVariable.begin(), m_cutVariable.end(), m_cutVariable.begin(), [](unsigned char c) { return std::tolower(c); });
    if (m_cutVariable != "completeness" && m_cutVariable != "purity")
    {
        std::cout << "CheatingSliceSelectionTool::ReadSettings: Unknown cut variable \'" << m_cutVariable << "\'" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
