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
#include "larpandoracontent/LArUtility/RecoTree.h"

#include "larpandoracontent/LArTwoDReco/LArClusterCreation/ProvisionalClusteringAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ProvisionalClusteringAlgorithm::ProvisionalClusteringAlgorithm() :
    m_maxGap(0.25f),
    m_maxGap2dSquared(2 * m_maxGap * m_maxGap)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProvisionalClusteringAlgorithm::Run()
{
    m_apaHitMap.clear();

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1, 1, 1));

    const CaloHitList *pCaloHitList{nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    this->PartitionHits(*pCaloHitList);
    this->ProcessPartition();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProvisionalClusteringAlgorithm::PartitionHits(const CaloHitList &caloHitList)
{
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
    for (const auto &[apaId, caloHits] : m_apaHitMap)
    {
        if (caloHits.empty())
            continue;

        //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHits, std::to_string(apaId), AUTOITER));
        //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

        KDTree kdTree;
        this->FillKDTree(caloHits, kdTree);

        std::cout << caloHits.size() << " hits in APA " << apaId << std::endl;
        OrderedCaloHitList provisionalOrderedCaloHits, orderedCaloHits;
        provisionalOrderedCaloHits.Add(caloHits);

        // Further sort each pseudo layer by position in x - avoids repeated sorts when bulding the reco tree
        for (OrderedCaloHitList::const_iterator iter = provisionalOrderedCaloHits.begin(); iter != provisionalOrderedCaloHits.end(); ++iter)
        {
            CaloHitVector caloHitVector(iter->second->begin(), iter->second->end());
            std::sort(caloHitVector.begin(), caloHitVector.end(), LArClusterHelper::SortHitsByPositionInX);
            if (caloHitVector.empty())
                continue;
            const CaloHitList sortedCaloHitList(caloHitVector.begin(), caloHitVector.end());
            orderedCaloHits.Add(sortedCaloHitList);
        }

        CaloHitSet ambiguousHits;
        this->TagAmbiguousHits(orderedCaloHits, ambiguousHits);

        CaloHitList unambiguousHitList, ambiguousHitList;
        for (OrderedCaloHitList::const_iterator iter = orderedCaloHits.begin(); iter != orderedCaloHits.end(); ++iter)
        {
            CaloHitVector currentLayerHits(iter->second->begin(), iter->second->end());
            for (const CaloHit *const pCaloHit : currentLayerHits)
            {
                if (ambiguousHits.count(pCaloHit) > 0)
                {
                    ambiguousHitList.push_back(pCaloHit);
                }
                else
                {
                    unambiguousHitList.push_back(pCaloHit);
                }
            }
        }
        if (!unambiguousHitList.empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &unambiguousHitList, std::to_string(apaId) + "_unambiguous", BLACK));
        }
        if (!ambiguousHitList.empty())
        {
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &ambiguousHitList, std::to_string(apaId) + "_ambiguous", RED));
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

        const float pitch{LArGeometryHelper::GetWirePitch(this->GetPandora(), caloHits.front()->GetHitType())};
        RecoTree recoTree(orderedCaloHits, ambiguousHits, pitch, this->GetPandora());
        recoTree.Populate();

        for (const auto &pNode : recoTree.GetRootNodes())
        {
            const CaloHitVector &nodeHits(pNode->GetHits());
            const CaloHitList hits(nodeHits.begin(), nodeHits.end());
            if (hits.empty())
                continue;

            // Visualize the hits in the reco tree
            PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits, std::to_string(apaId) + "_reco_tree", AUTOITER));
        }
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProvisionalClusteringAlgorithm::FillKDTree(const CaloHitList &caloHitList, KDTree &kdTree)
{
    typedef KDTreeNodeInfoT<const pandora::CaloHit *, 2> KDNode;
    typedef std::vector<KDNode> KDNodeVector;
    KDNodeVector kdNodes;
    KDTreeBox kdRegion(fill_and_bound_2d_kd_tree(caloHitList, kdNodes));
    kdTree.build(kdNodes, kdRegion);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProvisionalClusteringAlgorithm::TagAmbiguousHits(const OrderedCaloHitList &orderedCaloHits, CaloHitSet &ambiguousHits)
{
    for (OrderedCaloHitList::const_iterator iter = orderedCaloHits.begin(); iter != orderedCaloHits.end(); ++iter)
    {
        CaloHitVector currentLayerHits(iter->second->begin(), iter->second->end()), prevLayerHits, nextLayerHits;
        if (currentLayerHits.empty())
            continue;

        CaloHitList *pPrevLayerHits{nullptr}, *pNextLayerHits{nullptr};
        orderedCaloHits.GetCaloHitsInPseudoLayer(iter->first - 1, pPrevLayerHits);
        orderedCaloHits.GetCaloHitsInPseudoLayer(iter->first + 1, pNextLayerHits);
        if (pPrevLayerHits)
            prevLayerHits.insert(prevLayerHits.end(), pPrevLayerHits->begin(), pPrevLayerHits->end());
        if (pNextLayerHits)
            nextLayerHits.insert(nextLayerHits.end(), pNextLayerHits->begin(), pNextLayerHits->end());  

        for (const CaloHit *const pCurrentHit : currentLayerHits)
        {
            int neighbourhood[3][3]{
                {0, 0, 0},
                {0, 1, 0},
                {0, 0, 0}
            };
            // Consider current layer
            for (const CaloHit *const pOtherInLayerHit : currentLayerHits)
            {
                if (pCurrentHit == pOtherInLayerHit)
                    continue;

                if (pCurrentHit->GetPositionVector().GetX() >= pOtherInLayerHit->GetPositionVector().GetX())
                {
                    // Other hit is to the left
                    const float hit1Left{pCurrentHit->GetPositionVector().GetX() - 0.5f * pCurrentHit->GetCellSize1()};
                    const float hit2Right{pOtherInLayerHit->GetPositionVector().GetX() + 0.5f * pOtherInLayerHit->GetCellSize1()};
                    if (hit1Left - hit2Right < m_maxGap)
                    {
                        // Close enough to be a neighbour, check for blocked path, and if clear, add as neighbour
                        if (!LArClusterHelper::HasBlockedPath(currentLayerHits, pCurrentHit, pOtherInLayerHit))
                            ++neighbourhood[1][0];
                    }
                }
                else
                {
                    // Other hit is to the right
                    const float hit1Right{pCurrentHit->GetPositionVector().GetX() + 0.5f * pCurrentHit->GetCellSize1()};
                    const float hit2Left{pOtherInLayerHit->GetPositionVector().GetX() - 0.5f * pOtherInLayerHit->GetCellSize1()};
                    if (hit2Left - hit1Right < m_maxGap)
                    {
                        // Close enough to be a neighbour, check for blocked path, and if clear, add as neighbour
                        if (!LArClusterHelper::HasBlockedPath(currentLayerHits, pCurrentHit, pOtherInLayerHit))
                            ++neighbourhood[1][2];
                    }
                }
            }
            for (const CaloHit *const pOtherHit : prevLayerHits)
            {
                if (pCurrentHit == pOtherHit)
                    continue;

                if (pCurrentHit->GetPositionVector().GetX() >= pOtherHit->GetPositionVector().GetX())
                {
                    // Other hit is to the left
                    const float hit1Left{pCurrentHit->GetPositionVector().GetX() - 0.5f * pCurrentHit->GetCellSize1()};
                    const float hit2Right{pOtherHit->GetPositionVector().GetX() + 0.5f * pOtherHit->GetCellSize1()};
                    if (hit1Left - hit2Right < m_maxGap)
                    {
                        // Close enough to be a neighbour, check for blocked path, and if clear, add as neighbour
                        if (!(LArClusterHelper::HasBlockedPath(currentLayerHits, pCurrentHit, pOtherHit) ||
                              LArClusterHelper::HasBlockedPath(prevLayerHits, pCurrentHit, pOtherHit)))
                        {
                            const float overlap{hit2Right - hit1Left};
                            const float fractionOverlap{overlap / pOtherHit->GetCellSize1()};

                            if (fractionOverlap > 0.5f)
                            {
                                // Hit is below
                                ++neighbourhood[2][1];
                            }
                            else
                            {
                                // Hit is below and left
                                ++neighbourhood[2][0];
                            }
                        }
                    }
                }
                else
                {
                    // Other hit is to the right
                    const float hit1Right{pCurrentHit->GetPositionVector().GetX() + 0.5f * pCurrentHit->GetCellSize1()};
                    const float hit2Left{pOtherHit->GetPositionVector().GetX() - 0.5f * pOtherHit->GetCellSize1()};
                    if (hit2Left - hit1Right < m_maxGap)
                    {
                        // Close enough to be a neighbour, check for blocked path, and if clear, add as neighbour
                        if (!(LArClusterHelper::HasBlockedPath(currentLayerHits, pCurrentHit, pOtherHit) ||
                              LArClusterHelper::HasBlockedPath(prevLayerHits, pCurrentHit, pOtherHit)))
                        {
                            const float overlap{hit1Right - hit2Left};
                            const float fractionOverlap{overlap / pOtherHit->GetCellSize1()};

                            if (fractionOverlap > 0.5f)
                            {
                                // Hit is below
                                ++neighbourhood[2][1];
                            }
                            else
                            {
                                // Hit is below and right
                                ++neighbourhood[2][2];
                            }
                        }
                    }
                }
            }
            // Consider next layer
            for (const CaloHit *const pOtherHit : nextLayerHits)
            {
                if (pCurrentHit == pOtherHit)
                    continue;

                if (pCurrentHit->GetPositionVector().GetX() >= pOtherHit->GetPositionVector().GetX())
                {
                    // Other hit is to the left
                    const float hit1Left{pCurrentHit->GetPositionVector().GetX() - 0.5f * pCurrentHit->GetCellSize1()};
                    const float hit2Right{pOtherHit->GetPositionVector().GetX() + 0.5f * pOtherHit->GetCellSize1()};
                    if (hit1Left - hit2Right < m_maxGap)
                    {
                        // Close enough to be a neighbour, check for blocked path, and if clear, add as neighbour
                        if (!(LArClusterHelper::HasBlockedPath(currentLayerHits, pCurrentHit, pOtherHit) ||
                              LArClusterHelper::HasBlockedPath(nextLayerHits, pCurrentHit, pOtherHit)))
                        {
                            const float overlap{hit2Right - hit1Left};
                            const float fractionOverlap{overlap / pOtherHit->GetCellSize1()};

                            if (fractionOverlap > 0.5f)
                            {
                                // Hit is above
                                ++neighbourhood[0][1];
                            }
                            else
                            {
                                // Hit is above and left
                                ++neighbourhood[0][0];
                            }
                        }
                    }
                }
                else
                {
                    // Other hit is to the right
                    const float hit1Right{pCurrentHit->GetPositionVector().GetX() + 0.5f * pCurrentHit->GetCellSize1()};
                    const float hit2Left{pOtherHit->GetPositionVector().GetX() - 0.5f * pOtherHit->GetCellSize1()};
                    if (hit2Left - hit1Right < m_maxGap)
                    {
                        // Close enough to be a neighbour, check for blocked path, and if clear, add as neighbour
                        if (!(LArClusterHelper::HasBlockedPath(currentLayerHits, pCurrentHit, pOtherHit) ||
                              LArClusterHelper::HasBlockedPath(nextLayerHits, pCurrentHit, pOtherHit)))
                        {
                            const float overlap{hit1Right - hit2Left};
                            const float fractionOverlap{overlap / pOtherHit->GetCellSize1()};

                            if (fractionOverlap > 0.5f)
                            {
                                // Hit is above
                                ++neighbourhood[0][1];
                            }
                            else
                            {
                                // Hit is above and right
                                ++neighbourhood[0][2];
                            }
                        }
                    }
                }
            }

            float sum{0.f};
            FloatVector rowSums(3, 0.f), colSums(3, 0.f);
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j <3; ++j)
                {
                    sum += neighbourhood[i][j];
                    rowSums[i] += neighbourhood[i][j];
                    colSums[j] += neighbourhood[i][j];
                }
            }

            bool isAmbiguous{false};
            if (sum > 3)
            {
                isAmbiguous = true;
            }
            else if (sum == 3)
            {
                if (rowSums[1] == 1 && colSums[1] == 1)
                {
                    if (rowSums[0] == 0 || rowSums[2] == 0 || colSums[0] == 0 || colSums[2] == 0)
                    {
                        // Vertex-like, ambiguous
                        isAmbiguous = true;
                    }
                }
            }
            if (isAmbiguous)
            {
                ambiguousHits.insert(pCurrentHit);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProvisionalClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxGap", m_maxGap));
    m_maxGap2dSquared = 2 * m_maxGap * m_maxGap;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
