/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterCreation/ProvisionalClusteringAlgorithm.cc
 *
 *  @brief  Implementation of the cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include "larpandoracontent/LArTwoDReco/LArClusterCreation/ProvisionalClusteringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ProvisionalClusteringAlgorithm::ProvisionalClusteringAlgorithm() :
    m_maxGapSquared(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProvisionalClusteringAlgorithm::Run()
{
    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    this->PartitionHits(*pCaloHitList);
    this->ProcessPartition();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProvisionalClusteringAlgorithm::PartitionHits(const CaloHitList &caloHitList)
{
    m_apaHitMap.clear();

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pCaloHit));
        const ApaId apaId{pLArCaloHit ? pLArCaloHit->GetDaughterVolumeId() : 0};
        if (m_apaHitMap.find(apaId) == m_apaHitMap.end())
            m_apaHitMap[apaId] = CaloHitList();
        m_apaHitMap[apaId].emplace_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProvisionalClusteringAlgorithm::ProcessPartition()
{
    typedef KDTreeNodeInfoT<const pandora::CaloHit *, 2> KDNode;
    typedef std::vector<KDNode> KDNodeVector;
    typedef KDTreeLinkerAlgo<const pandora::CaloHit *, 2> KDTree;

    for (const auto &[apaId, caloHits] : m_apaHitMap)
    {
        if (caloHits.empty())
            continue;

        OrderedCaloHitList orderedCaloHits;
        orderedCaloHits.Add(caloHits);

        KDTree kdTree;
        KDNodeVector kdNodes;
        KDTreeBox kdRegion(fill_and_bound_2d_kd_tree(caloHits, kdNodes));
        kdTree.build(kdNodes, kdRegion);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProvisionalClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGapSquared", m_maxGapSquared));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
