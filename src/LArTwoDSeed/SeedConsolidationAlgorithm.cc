/**
 *  @file   LArContent/src/LArTwoDSeed/SeedConsolidationAlgorithm.cc
 * 
 *  @brief  Implementation of the seed consolidation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArVertexHelper.h"

#include "LArTwoDSeed/SeedConsolidationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode SeedConsolidationAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterList(*this, pClusterList));

    // Create ordered list of shower-like clusters
    ClusterVector clusterVector;

    for (ClusterList::iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        if (!ParticleIdHelper::IsMuonFast(*iter))
            clusterVector.push_back(*iter);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

    // Examine associations between potential parent/daughter cluster pairs
    for (ClusterVector::iterator iterP = clusterVector.begin(), iterPEnd = clusterVector.end(); iterP != iterPEnd; ++iterP)
    {
        Cluster *pParentCluster = *iterP;

        if (NULL == pParentCluster)
            continue;

        ClusterIteratorList daughterIterators;

        for (ClusterVector::iterator iterD = clusterVector.begin(), iterDEnd = clusterVector.end(); iterD != iterDEnd; ++iterD)
        {
            Cluster *pDaughterCluster = *iterD;

            if ((NULL == pDaughterCluster) || (pParentCluster == pDaughterCluster))
                continue;

            if (pDaughterCluster->GetOrderedCaloHitList().size() > 2 * pParentCluster->GetOrderedCaloHitList().size())
                continue;

            if (this->PassesVertexAssociation(pParentCluster, pDaughterCluster))
                continue;

            try
            {
                 const ConeAssociation coneAssociation(pParentCluster, pDaughterCluster);

                 if (this->PassesConeAssociation(coneAssociation))
                 {
                     daughterIterators.push_back(iterD);
                     continue;
                 }

                 if (this->PassesRecoveryCuts(coneAssociation))
                 {
                     daughterIterators.push_back(iterD);
                     continue;
                }
            }
            catch (StatusCodeException &)
            {
            }

            try
            {
                const OverlapAssociation overlapAssociation(pParentCluster, pDaughterCluster);

                if (this->PassesOverlapAssociation(overlapAssociation))
                {
//std::cout << " SEEDCONSOLIDATIONALGORITHM WOULD HAVE MADE MERGE HERE!!!!!!!!! " << std::endl;
//std::cout << " overlapAssociation.GetNContactLayers() " << overlapAssociation.GetNContactLayers() << std::endl;
//std::cout << " overlapAssociation.GetLayerDistance50() " << overlapAssociation.GetLayerDistance50() << std::endl;
//std::cout << " overlapAssociation.GetNEnclosedLayers() " << overlapAssociation.GetNEnclosedLayers() << std::endl;
//std::cout << " overlapAssociation.GetNNonVtxContactGroups() " << overlapAssociation.GetNNonVtxContactGroups() << std::endl;
//std::cout << " overlapAssociation.GetNNonContactLayers() " << overlapAssociation.GetNNonContactLayers() << std::endl;
//std::cout << " overlapAssociation.GetNOverlapLayers() " << overlapAssociation.GetNOverlapLayers() << std::endl;
//ClusterList parent, daughter;
//parent.insert(pParentCluster);
//daughter.insert(pDaughterCluster);
//PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
//PandoraMonitoringApi::VisualizeClusters(&parent, "parent", RED);
//PandoraMonitoringApi::VisualizeClusters(&daughter, "daughter", GREEN);
//PandoraMonitoringApi::VisualizeClusters(pClusterList, "all", BLUE);
//PandoraMonitoringApi::ViewEvent();
                    daughterIterators.push_back(iterD);
                    continue;
                }
            }
            catch (StatusCodeException &)
            {
            }
        }

        for (ClusterIteratorList::const_iterator iterJ = daughterIterators.begin(), iterJEnd = daughterIterators.end(); iterJ != iterJEnd; ++iterJ)
        {
            Cluster *pDaughterCluster = *(*iterJ);
            *(*iterJ) = NULL;

// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// ClusterList tempList1, tempList2;
// tempList1.insert(pParentCluster);
// tempList2.insert(pDaughterCluster);
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "Parent", BLUE);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "Daughter", GREEN);
// PandoraMonitoringApi::ViewEvent();

            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pParentCluster, pDaughterCluster));
        }

        // Repeat with modified parent cluster
        if (!daughterIterators.empty())
            --iterP;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SeedConsolidationAlgorithm::PassesVertexAssociation( const Cluster* const pClusterParent, const Cluster* const pClusterDaughter )
{
    if ( LArVertexHelper::DoesCurrentVertexExist()==false ) return false;

    // Checks on proximity to vertex of parent and dughter clusters
    const float distanceToVertexSquaredParent = LArVertexHelper::GetDistanceSquaredToCurrentVertex(pClusterParent);
    const float distanceToVertexSquaredDaughter = LArVertexHelper::GetDistanceSquaredToCurrentVertex(pClusterDaughter);

    const float impactParameterToVertexParent = LArVertexHelper::GetImpactParameterToCurrentVertex(pClusterParent);
    const float impactParameterToVertexDaughter = LArVertexHelper::GetImpactParameterToCurrentVertex(pClusterDaughter);

    const CartesianVector parentInnerCentroid = pClusterParent->GetCentroid(pClusterParent->GetInnerPseudoLayer());
    const CartesianVector parentOuterCentroid = pClusterParent->GetCentroid(pClusterParent->GetOuterPseudoLayer());
    const CartesianVector parentDirection = (parentOuterCentroid-parentInnerCentroid).GetUnitVector();

    const CartesianVector daughterInnerCentroid = pClusterDaughter->GetCentroid(pClusterDaughter->GetInnerPseudoLayer());
    const CartesianVector daughterOuterCentroid = pClusterDaughter->GetCentroid(pClusterDaughter->GetOuterPseudoLayer());
    const CartesianVector daughterDirection = (daughterOuterCentroid-daughterInnerCentroid).GetUnitVector();

    const float cosineRelativeAngle = parentDirection.GetDotProduct(daughterDirection);
    const float parentLengthSquared = (parentOuterCentroid-parentInnerCentroid).GetMagnitudeSquared();

    if ( distanceToVertexSquaredParent < 25.0 && distanceToVertexSquaredDaughter < 25.0 ) 
        return true;

    if ( impactParameterToVertexParent < 2.5 && impactParameterToVertexDaughter < 2.5
      && distanceToVertexSquaredDaughter < 0.25 * parentLengthSquared && cosineRelativeAngle < 0.87 )
        return true;

    if ( distanceToVertexSquaredDaughter < distanceToVertexSquaredParent )
        return true; 

    // Checks on relative direction of parent and daughter clusters
    const bool isParentForwardInZ  = LArVertexHelper::IsForwardInZ(pClusterParent);
    const bool isParentBackwardInZ = LArVertexHelper::IsBackwardInZ(pClusterParent);

    const bool isDaughterForwardInZ = LArVertexHelper::IsForwardInZ(pClusterDaughter);
    const bool isDaughterBackwardInZ = LArVertexHelper::IsBackwardInZ(pClusterDaughter);

    if( true==isParentForwardInZ && true==isDaughterBackwardInZ )
        return true;

    if( true==isParentBackwardInZ && true==isDaughterForwardInZ )
        return true;

    if( true==isParentForwardInZ && pClusterParent->GetInnerPseudoLayer() + 10 > pClusterDaughter->GetInnerPseudoLayer() )
        return true;

    if( true==isParentBackwardInZ && pClusterDaughter->GetOuterPseudoLayer() + 10 > pClusterParent->GetOuterPseudoLayer() )
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
bool SeedConsolidationAlgorithm::PassesConeAssociation(const ConeAssociation &coneAssociation) const
{
    if (coneAssociation.GetEnclosedDaughterHitFraction() < 0.75f)
        return false;

    if (coneAssociation.GetEnclosedParentHitFraction() < 0.25f)
        return false;

    if (coneAssociation.GetCosRelativeAngle() < 0.87f)
        return false;

    if (coneAssociation.GetApexToDaughterDistance() < 20.f)
        return false;

    if (coneAssociation.GetApexToDaughterDistanceFraction() < 0.33f )
        return false;

    if (coneAssociation.GetParentToDaughterDistance() > 30.f)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SeedConsolidationAlgorithm::PassesRecoveryCuts(const ConeAssociation &coneAssociation) const
{
    if (coneAssociation.GetDaughterStartDistance() > 5.f)
        return false;

    if (coneAssociation.GetDaughterEndDistance() > 40.f )
        return false;

    if (coneAssociation.GetDaughterStartDistance() > coneAssociation.GetDaughterEndDistance())
        return false;

    if (coneAssociation.GetApexToDaughterDistance() < 40.f)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SeedConsolidationAlgorithm::PassesOverlapAssociation(const OverlapAssociation &overlapAssociation) const
{   
    if (overlapAssociation.GetNContactLayers() == 0)
        return false;

    if (overlapAssociation.GetLayerDistance50() > 4.f)
        return false;

    if (overlapAssociation.GetNEnclosedLayers() > 5)
        return true;

    if (overlapAssociation.GetNNonVtxContactGroups() > 1)
        return true;

    if ((overlapAssociation.GetNContactLayers() > 40) && (overlapAssociation.GetNNonContactLayers() < 0.5 * overlapAssociation.GetNOverlapLayers()))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

SeedConsolidationAlgorithm::ConeAssociation::ConeAssociation(const Cluster *const pParentCluster, const Cluster *const pDaughterCluster) :
    m_cosRelativeAngle(-1.f),
    m_daughterStartDistance(std::numeric_limits<float>::max()),
    m_daughterEndDistance(std::numeric_limits<float>::max()),
    m_apexToDaughterDistance(std::numeric_limits<float>::max()),
    m_parentToDaughterDistance(std::numeric_limits<float>::max()),
    m_coneApex(0.f, 0.f, 0.f),
    m_coneDirection(0.f, 0.f, 0.f),
    m_coneLength(0.f),
    m_isForwardsInZ(true),
    m_coneCosHalfAngleParent(1.f),
    m_coneCosHalfAngleDaughter(1.f),
    m_enclosedParentHitFraction(0.f),
    m_enclosedDaughterHitFraction(0.f)
{
    this->CalculateProperties(pParentCluster, pDaughterCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedConsolidationAlgorithm::ConeAssociation::CalculateProperties(const Cluster *const pParentCluster, const Cluster *const pDaughterCluster)
{
    const ClusterHelper::ClusterFitResult &parentFit(pParentCluster->GetFitToAllHitsResult());
    const ClusterHelper::ClusterFitResult &daughterFit(pDaughterCluster->GetFitToAllHitsResult());
    m_cosRelativeAngle = parentFit.GetDirection().GetCosOpeningAngle(daughterFit.GetDirection());

    // Define cone parameters and parent-daughter distances
    const float rInner(LArClusterHelper::GetClosestDistance(pDaughterCluster->GetCentroid(pDaughterCluster->GetInnerPseudoLayer()), pParentCluster));
    const float rOuter(LArClusterHelper::GetClosestDistance(pDaughterCluster->GetCentroid(pDaughterCluster->GetOuterPseudoLayer()), pParentCluster));
    const bool coneForwardsInZ(rOuter > rInner);
    m_daughterStartDistance = (coneForwardsInZ ? rInner : rOuter);
    m_daughterEndDistance = (coneForwardsInZ ? rOuter : rInner);
    m_coneDirection = coneForwardsInZ ? parentFit.GetDirection() : parentFit.GetDirection() * -1.f;

    const unsigned int parentStartLayer(coneForwardsInZ ? pParentCluster->GetInnerPseudoLayer() : pParentCluster->GetOuterPseudoLayer());
    const CartesianVector parentStartCentroid(pParentCluster->GetCentroid(parentStartLayer));
    m_coneApex = parentFit.GetIntercept() + m_coneDirection * (m_coneDirection.GetDotProduct(parentStartCentroid - parentFit.GetIntercept()));

    const unsigned int parentEndLayer(coneForwardsInZ ? pParentCluster->GetOuterPseudoLayer() : pParentCluster->GetInnerPseudoLayer());
    const CartesianVector parentEndCentroid(pParentCluster->GetCentroid(parentEndLayer));
    m_coneLength = m_coneDirection.GetDotProduct(parentEndCentroid - m_coneApex);

    const unsigned int daughterStartLayer(coneForwardsInZ ? pDaughterCluster->GetInnerPseudoLayer() : pDaughterCluster->GetOuterPseudoLayer());
    const CartesianVector daughterStartCentroid(pDaughterCluster->GetCentroid(daughterStartLayer));
    m_apexToDaughterDistance = m_coneDirection.GetDotProduct(daughterStartCentroid - m_coneApex);
    m_parentToDaughterDistance = m_coneDirection.GetDotProduct(daughterStartCentroid - parentEndCentroid);


    // Calculate 90% cone cos half angle of parent
    FloatVector coneCosHalfAnglesParent;
    const OrderedCaloHitList &parentCaloHitList(pParentCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = parentCaloHitList.begin(), iterEnd = parentCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
            coneCosHalfAnglesParent.push_back( m_coneDirection.GetDotProduct(((*hIter)->GetPositionVector() - m_coneApex).GetUnitVector()) );
    }

    std::sort(coneCosHalfAnglesParent.begin(), coneCosHalfAnglesParent.end());
    
    m_coneCosHalfAngleParent = coneCosHalfAnglesParent[static_cast<unsigned int>(0.10f * coneCosHalfAnglesParent.size())];


    // Calculate 10% cone cos half angle of daughter
    FloatVector coneCosHalfAnglesDaughter;
    const OrderedCaloHitList &daughterCaloHitList(pDaughterCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = daughterCaloHitList.begin(), iterEnd = daughterCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
            coneCosHalfAnglesDaughter.push_back( m_coneDirection.GetDotProduct(((*hIter)->GetPositionVector() - m_coneApex).GetUnitVector()) );
    } 

    std::sort(coneCosHalfAnglesDaughter.begin(), coneCosHalfAnglesDaughter.end());

    m_coneCosHalfAngleDaughter = coneCosHalfAnglesDaughter[static_cast<unsigned int>(0.90f * coneCosHalfAnglesDaughter.size())];


    // Calculate fraction of daughter hits bounded by the parent 90% cos half angle
    unsigned int nDaughterHitsBoundedByParent(0);
   
    for (OrderedCaloHitList::const_iterator iter = daughterCaloHitList.begin(), iterEnd = daughterCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
        {
            if (m_coneDirection.GetDotProduct(((*hIter)->GetPositionVector() - m_coneApex).GetUnitVector()) > m_coneCosHalfAngleParent)
                ++nDaughterHitsBoundedByParent;
        }
    }

    if (pDaughterCluster->GetNCaloHits() > 0)
        m_enclosedDaughterHitFraction = static_cast<float>(nDaughterHitsBoundedByParent) / static_cast<float>(pDaughterCluster->GetNCaloHits());
    

    // Calculate fraction of parent hits required to bound to daughter 10% cos half angle
    unsigned int ParentHitsBoundedByDaughter(0);
   
    for (OrderedCaloHitList::const_iterator iter = parentCaloHitList.begin(), iterEnd = parentCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
        {
            if (m_coneDirection.GetDotProduct(((*hIter)->GetPositionVector() - m_coneApex).GetUnitVector()) < m_coneCosHalfAngleDaughter)
                ++ParentHitsBoundedByDaughter;
        }
    }

    if (pParentCluster->GetNCaloHits() > 0)
        m_enclosedParentHitFraction = static_cast<float>(ParentHitsBoundedByDaughter) / static_cast<float>(pParentCluster->GetNCaloHits());


// CartesianVector p(0.,1.,0.);
// CartesianVector Marker1(m_coneApex);
// CartesianVector Marker2(m_coneApex + m_coneDirection * m_coneLength);
// CartesianVector Marker3(m_coneApex + m_coneDirection * m_coneLength + m_coneDirection.GetCrossProduct(p) * m_coneLength * sqrt ( 1.0 - m_coneCosHalfAngleParent * m_coneCosHalfAngleParent ) );
// CartesianVector Marker4(m_coneApex + m_coneDirection * m_coneLength - m_coneDirection.GetCrossProduct(p) * m_coneLength * sqrt ( 1.0 - m_coneCosHalfAngleParent * m_coneCosHalfAngleParent ) );

// ClusterList tempList1, tempList2;
// Cluster* tempC1 = (Cluster*)pParentCluster;
// Cluster* tempC2 = (Cluster*)pDaughterCluster;
// tempList1.insert(tempC1);
// tempList2.insert(tempC2);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "Parent", BLUE);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "Daughter", GREEN);
// PandoraMonitoringApi::AddMarkerToVisualization(&Marker1, "M1", AUTO, 1.);
// PandoraMonitoringApi::AddMarkerToVisualization(&Marker2, "M2", AUTO, 1.);
// PandoraMonitoringApi::AddMarkerToVisualization(&Marker3, "M3", AUTO, 1.);
// PandoraMonitoringApi::AddMarkerToVisualization(&Marker4, "M4", AUTO, 1.);
// PandoraMonitoringApi::ViewEvent();

}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

SeedConsolidationAlgorithm::OverlapAssociation::OverlapAssociation(const Cluster *const pParentCluster, const Cluster *const pDaughterCluster) :
    m_innerOverlapDistance(0.f),
    m_outerOverlapDistance(0.f),
    m_innerOverlapLayer(std::numeric_limits<unsigned int>::max()),
    m_outerOverlapLayer(0),
    m_nOverlapLayers(0),
    m_nContactLayers(0),
    m_nNonContactLayers(std::numeric_limits<unsigned int>::max()),
    m_nContactGroups(0),
    m_nNonVtxContactGroups(0),
    m_nEnclosedLayers(0),
    m_layerDistance50(std::numeric_limits<float>::max())
{
    this->CalculateProperties(pParentCluster, pDaughterCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedConsolidationAlgorithm::OverlapAssociation::CalculateProperties(const Cluster *const pParentCluster, const Cluster *const pDaughterCluster)
{
    if (pParentCluster->GetOuterPseudoLayer() < pDaughterCluster->GetInnerPseudoLayer())
        return;

    if (pDaughterCluster->GetOuterPseudoLayer() < pParentCluster->GetInnerPseudoLayer())
        return;

    // Overlap properties
    const Cluster *const pInnerCluster((pParentCluster->GetInnerPseudoLayer() < pDaughterCluster->GetInnerPseudoLayer()) ? pParentCluster : pDaughterCluster);
    const Cluster *const pInnerAltCluster((pInnerCluster == pParentCluster) ? pDaughterCluster : pParentCluster);
    const Cluster *const pOuterCluster = ((pParentCluster->GetOuterPseudoLayer() > pDaughterCluster->GetOuterPseudoLayer()) ? pParentCluster : pDaughterCluster);
    const Cluster *const pOuterAltCluster((pOuterCluster == pParentCluster) ? pDaughterCluster : pParentCluster);

    m_innerOverlapLayer = pInnerAltCluster->GetInnerPseudoLayer();
    m_outerOverlapLayer = pOuterAltCluster->GetOuterPseudoLayer();
    m_nOverlapLayers = m_outerOverlapLayer - m_innerOverlapLayer + 1;
    m_innerOverlapDistance = (pInnerCluster->GetCentroid(pInnerCluster->GetInnerPseudoLayer()) - pInnerAltCluster->GetCentroid(m_innerOverlapLayer)).GetMagnitude();
    m_outerOverlapDistance = (pOuterAltCluster->GetCentroid(m_outerOverlapLayer) - pOuterCluster->GetCentroid(pOuterCluster->GetOuterPseudoLayer())).GetMagnitude();

    // Layer proximities
    const OrderedCaloHitList &parentList(pParentCluster->GetOrderedCaloHitList());
    const OrderedCaloHitList &daughterList(pDaughterCluster->GetOrderedCaloHitList());

    FloatVector layerDistances;
    unsigned int layersSinceLastContact(std::numeric_limits<unsigned int>::max());

    for (unsigned int iLayer = m_innerOverlapLayer; iLayer <= m_outerOverlapLayer; ++iLayer)
    {
        const OrderedCaloHitList::const_iterator daughterListIter = daughterList.find(iLayer);

        if (daughterList.end() == daughterListIter)
        {
            ++layersSinceLastContact;
            continue;
        }

        bool enclosedHits(false);
        float closestHitDistance(std::numeric_limits<float>::max());
        this->EvaluateLayerContact(daughterListIter, parentList, closestHitDistance, enclosedHits);

        layerDistances.push_back(closestHitDistance);

        if (enclosedHits)
            ++m_nEnclosedLayers;

        if (closestHitDistance < 2.5f)
        {
            ++m_nContactLayers;

            if (layersSinceLastContact > 4)
            {
                ++m_nContactGroups;
                layersSinceLastContact = 0;
            }
        }
        else
        {
            ++layersSinceLastContact;
        }
    }

    m_nNonContactLayers = m_nOverlapLayers - m_nContactLayers;

    m_nNonVtxContactGroups = ((m_nNonVtxContactGroups > 0) && (std::min(m_innerOverlapDistance, m_outerOverlapDistance) < 5.f)) ?
        m_nContactGroups - 1 : m_nContactGroups;

    if(layerDistances.size()>0){
      std::sort(layerDistances.begin(), layerDistances.end());
      m_layerDistance50 = layerDistances[static_cast<unsigned int>(0.5 * layerDistances.size())];
    }
    else{
      m_layerDistance50 = 0.;     //Will fail anyway as there are no contact layers
    }

// ClusterList tempList1, tempList2;
// Cluster* tempC1 = (Cluster*)pParentCluster;
// Cluster* tempC2 = (Cluster*)pDaughterCluster;
// tempList1.insert(tempC1);
// tempList2.insert(tempC2);
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempList1, "Parent", BLUE);
// PandoraMonitoringApi::VisualizeClusters(&tempList2, "Daughter", GREEN);
// PandoraMonitoringApi::ViewEvent();

}

//------------------------------------------------------------------------------------------------------------------------------------------

void SeedConsolidationAlgorithm::OverlapAssociation::EvaluateLayerContact(const OrderedCaloHitList::const_iterator daughterListIter,
    const OrderedCaloHitList &parentList, float &closestHitDistance, bool &enclosedHits) const
{
    closestHitDistance = std::numeric_limits<float>::max();
    enclosedHits = false;

    const unsigned int daughterLayer(daughterListIter->first);
    const CaloHitList *const pDaughterHitList(daughterListIter->second);

    const unsigned int minLayer((daughterLayer > 2) ? daughterLayer - 2 : 0);
    const unsigned int maxLayer(daughterLayer + 2);

    bool distanceFound(false);
    float closestDistanceSquared(std::numeric_limits<float>::max());

    for (unsigned int iLayer = minLayer; iLayer <= maxLayer; ++iLayer)
    {
        const OrderedCaloHitList::const_iterator parentListIter = parentList.find(iLayer);

        if (parentList.end() == parentListIter)
            continue;

        const CaloHitList *const pParentHitList(parentListIter->second);

        for (CaloHitList::const_iterator dIter = pDaughterHitList->begin(), dIterEnd = pDaughterHitList->end(); dIter != dIterEnd; ++dIter)
        {
            CartesianVectorList cartesianVectorList;

            for (CaloHitList::const_iterator pIter = pParentHitList->begin(), pIterEnd = pParentHitList->end(); pIter != pIterEnd; ++pIter)
            {
                const CartesianVector displacement((*dIter)->GetPositionVector() - (*pIter)->GetPositionVector());
                const float distanceSquared(displacement.GetMagnitudeSquared());

                if (distanceSquared < closestDistanceSquared)
                {
                    closestDistanceSquared = distanceSquared;
                    distanceFound = true;
                }

                if (enclosedHits)
                    continue;

                for (CartesianVectorList::const_iterator vIter = cartesianVectorList.begin(), vIterEnd = cartesianVectorList.end(); vIter != vIterEnd; ++vIter)
                {
                    if (displacement.GetDotProduct(*vIter) < 0.f)
                    {
                        enclosedHits = true;
                        break;
                    }
                }

                cartesianVectorList.push_back(displacement);
            }
        }
    }

    if (distanceFound)
        closestHitDistance = std::sqrt(closestDistanceSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SeedConsolidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
