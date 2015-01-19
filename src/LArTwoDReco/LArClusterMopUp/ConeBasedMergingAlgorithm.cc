/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterMopUp/ConeBasedMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the cone based merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "LArTwoDReco/LArClusterMopUp/ConeBasedMergingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ConeBasedMergingAlgorithm::ConeBasedMergingAlgorithm() :
    m_slidingFitWindow(20),
    m_showerEdgeMultiplier(1.5f),
    m_coneAngleCentile(0.85f),
    m_maxConeLengthMultiplier(1.5f),
    m_minBoundedFraction(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeBasedMergingAlgorithm::ClusterMopUp(const ClusterList &pfoClusters, const ClusterList &remnantClusters,
    const ClusterToListNameMap &clusterToListNameMap) const
{
    ClusterAssociationMap clusterAssociationMap;
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    const VertexList *pVertexList(NULL);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));

    if ((pVertexList->size() != 1) || (VERTEX_3D != (*(pVertexList->begin()))->GetVertexType()))
        return;

    const Vertex *pVertex(*(pVertexList->begin()));

    for (ClusterList::const_iterator pIter = pfoClusters.begin(), pIterEnd = pfoClusters.end(); pIter != pIterEnd; ++pIter)
    {
        try
        {
            Cluster *pPfoCluster(*pIter);
            const TwoDSlidingShowerFitResult showerFitResult(pPfoCluster, m_slidingFitWindow, slidingFitPitch, m_showerEdgeMultiplier);

            const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());
            const LayerFitResultMap &layerFitResultMapP(showerFitResult.GetPositiveEdgeFitResult().GetLayerFitResultMap());
            const LayerFitResultMap &layerFitResultMapN(showerFitResult.GetNegativeEdgeFitResult().GetLayerFitResultMap());

            if (layerFitResultMapS.size() < 2)
                continue;

            // Cone direction
            CartesianVector minLayerPositionOnAxis(0.f, 0.f, 0.f), maxLayerPositionOnAxis(0.f, 0.f, 0.f);
            showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.begin()->second.GetL(), 0.f, minLayerPositionOnAxis);
            showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.rbegin()->second.GetL(), 0.f, maxLayerPositionOnAxis);

            const HitType hitType(LArClusterHelper::GetClusterHitType(pPfoCluster));
            const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType));
            const bool vertexAtMinL((vertexPosition2D - minLayerPositionOnAxis).GetMagnitudeSquared() < (vertexPosition2D - maxLayerPositionOnAxis).GetMagnitudeSquared());

            // Cone edges
            CoordinateList coordinateListP, coordinateListN;

            for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)
            {
                LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(iterS->first);
                LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(iterS->first);

                if (layerFitResultMapP.end() != iterP)
                    coordinateListP.push_back(Coordinate(iterP->second.GetL(), iterP->second.GetFitT()));

                if (layerFitResultMapN.end() != iterN)
                    coordinateListN.push_back(Coordinate(iterN->second.GetL(), iterN->second.GetFitT()));
            }

            if (coordinateListP.empty() || coordinateListN.empty())
                continue;

            std::sort(coordinateListP.begin(), coordinateListP.end(), ConeBasedMergingAlgorithm::SortCoordinates);
            std::sort(coordinateListN.begin(), coordinateListN.end(), ConeBasedMergingAlgorithm::SortCoordinates);

            const Coordinate maxP(coordinateListP.at(m_coneAngleCentile * coordinateListP.size()));
            const Coordinate minP(vertexAtMinL ? Coordinate(layerFitResultMapP.begin()->second.GetL(), layerFitResultMapP.begin()->second.GetFitT()) :
                Coordinate(layerFitResultMapP.rbegin()->second.GetL(), layerFitResultMapP.rbegin()->second.GetFitT()));

            const Coordinate maxN(coordinateListN.at((1.f - m_coneAngleCentile) * coordinateListN.size()));
            const Coordinate minN(vertexAtMinL ? Coordinate(layerFitResultMapN.begin()->second.GetL(), layerFitResultMapN.begin()->second.GetFitT()) :
                Coordinate(layerFitResultMapN.rbegin()->second.GetL(), layerFitResultMapN.rbegin()->second.GetFitT()));

            const float minL(layerFitResultMapS.begin()->second.GetL());
            const float maxL(minL + m_maxConeLengthMultiplier * std::fabs(layerFitResultMapS.rbegin()->second.GetL() - minL));

            if ((std::fabs(maxP.first - minP.first) < std::numeric_limits<float>::epsilon()) ||
                (std::fabs(maxN.first - minN.first) < std::numeric_limits<float>::epsilon()) ||
                (std::fabs(maxL - minL) < std::numeric_limits<float>::epsilon()))
            {
                continue;
            }

            // Bounded fraction calculation
            for (ClusterList::const_iterator rIter = remnantClusters.begin(), rIterEnd = remnantClusters.end(); rIter != rIterEnd; ++rIter)
            {
                Cluster *const pRemnantCluster(*rIter);
                const unsigned int nHits(pRemnantCluster->GetNCaloHits());

                unsigned int nMatchedHits(0);
                const OrderedCaloHitList &orderedCaloHitList(pRemnantCluster->GetOrderedCaloHitList());

                for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
                {
                    for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
                    {
                        float rL(0.f), rT(0.f);
                        showerFitResult.GetShowerFitResult().GetLocalPosition((*hIter)->GetPositionVector(), rL, rT);

                        if ((rL < minL) || (rL > maxL))
                            continue;

                        const float rTP(minP.second + (rL - minP.first) * ((maxP.second - minP.second) / (maxP.first - minP.first)));
                        const float rTN(minN.second + (rL - minN.first) * ((maxN.second - minN.second) / (maxN.first - minN.first)));

                        if ((rT > rTP) || (rT < rTN))
                            continue;

                        ++nMatchedHits;
                    }
                }

                const float boundedFraction((nHits > 0) ? static_cast<float>(nMatchedHits) / static_cast<float>(nHits) : 0.f);

                if (boundedFraction < m_minBoundedFraction)
                    continue;

                AssociationDetails &associationDetails(clusterAssociationMap[pRemnantCluster]);

                if (!associationDetails.insert(AssociationDetails::value_type(pPfoCluster, boundedFraction)).second)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
            }
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }

    this->MakeClusterMerges(clusterAssociationMap, clusterToListNameMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConeBasedMergingAlgorithm::SortCoordinates(const Coordinate &lhs, const Coordinate &rhs)
{
    return (lhs.second < rhs.second);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeBasedMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ShowerEdgeMultiplier", m_showerEdgeMultiplier));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "ConeAngleCentile", m_coneAngleCentile));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MaxConeLengthMultiplier", m_maxConeLengthMultiplier));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinBoundedFraction", m_minBoundedFraction));

    return ClusterMopUpAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
