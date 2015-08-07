/**
 *  @file   LArContent/src/LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.cc
 * 
 *  @brief  Implementation of the branch associated pfos tool class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArEventBuilding/BranchAssociatedPfosTool.h"

using namespace pandora;

namespace lar_content
{

typedef PfoHierarchyAlgorithm::PfoInfo PfoInfo;
typedef PfoHierarchyAlgorithm::PfoInfoMap PfoInfoMap;

void BranchAssociatedPfosTool::Run(PfoHierarchyAlgorithm *const pAlgorithm, const Vertex *const /*pNeutrinoVertex*/, PfoInfoMap &pfoInfoMap)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this << ", " << this->GetType() << std::endl;

    bool associationsMade(true);

    while (associationsMade)
    {
        associationsMade = false;
        PfoList assignedPfos, unassignedPfos;
        pAlgorithm->SeparatePfos(pfoInfoMap, assignedPfos, unassignedPfos);

        if (unassignedPfos.empty())
            break;

        for (const ParticleFlowObject *const pParentPfo : assignedPfos)
        {
            PfoInfoMap::iterator parentMapIter(pfoInfoMap.find(pParentPfo));

            if (pfoInfoMap.end() == parentMapIter)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            PfoInfo *const pParentPfoInfo(parentMapIter->second);

            for (const ParticleFlowObject *const pPfo : unassignedPfos)
            {
                PfoInfoMap::iterator mapIter(pfoInfoMap.find(pPfo));

                if (pfoInfoMap.end() == mapIter)
                    throw StatusCodeException(STATUS_CODE_FAILURE);

                PfoInfo *const pPfoInfo(mapIter->second);

                if (false)
                {
                    associationsMade = true;
                    pParentPfoInfo->AddDaughterPfo(pPfoInfo->GetThisPfo());
                    pPfoInfo->SetParentPfo(pParentPfoInfo->GetThisPfo());
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BranchAssociatedPfosTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
