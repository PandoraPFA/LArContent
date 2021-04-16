/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/DeltaRayExtensionAlgorithm.cc
 *
 *  @brief  Implementation of the delta ray extension algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArCosmicRay/DeltaRayExtensionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

DeltaRayExtensionAlgorithm::DeltaRayExtensionAlgorithm() :
    m_minClusterLength(1.f),
    m_maxClusterLength(10.f),
    m_maxLongitudinalDisplacement(2.5f),
    m_maxTransverseDisplacement(1.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayExtensionAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayExtensionAlgorithm::FillClusterAssociationMatrix(const ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    ClusterToCoordinateMap innerCoordinateMap, outerCoordinateMap;

    for (const Cluster *const pParentCluster : clusterVector)
    {
        for (const Cluster *const pDaughterCluster : clusterVector)
        {
            if (pParentCluster == pDaughterCluster)
                continue;

            this->FillClusterAssociationMatrix(pParentCluster, pDaughterCluster, innerCoordinateMap, outerCoordinateMap, clusterAssociationMatrix);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayExtensionAlgorithm::GetExtremalCoordinatesFromCache(const Cluster *const pCluster, ClusterToCoordinateMap &innerCoordinateMap,
    ClusterToCoordinateMap &outerCoordinateMap, CartesianVector &innerCoordinate, CartesianVector &outerCoordinate) const
{
    ClusterToCoordinateMap::const_iterator innerIter = innerCoordinateMap.find(pCluster);
    ClusterToCoordinateMap::const_iterator outerIter = outerCoordinateMap.find(pCluster);

    if ((innerCoordinateMap.end() == innerIter) || (outerCoordinateMap.end() == outerIter))
    {
        LArClusterHelper::GetExtremalCoordinates(pCluster, innerCoordinate, outerCoordinate);
        (void)innerCoordinateMap.insert(ClusterToCoordinateMap::value_type(pCluster, innerCoordinate));
        (void)outerCoordinateMap.insert(ClusterToCoordinateMap::value_type(pCluster, outerCoordinate));
    }
    else
    {
        innerCoordinate = innerIter->second;
        outerCoordinate = outerIter->second;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayExtensionAlgorithm::FillClusterAssociationMatrix(const Cluster *const pParentCluster, const Cluster *const pDaughterCluster,
    ClusterToCoordinateMap &innerCoordinateMap, ClusterToCoordinateMap &outerCoordinateMap, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    // Daughter cluster must be available for any association to proceed
    if (!PandoraContentApi::IsAvailable(*this, pDaughterCluster))
        return;

    // Weak association:   between parent cosmic-ray muon and daughter delta ray
    // Strong association: between parent and daughter fragments of delta ray
    // Figure of merit:    distance between parent and daughter clusters
    CartesianVector innerCoordinateP(0.f, 0.f, 0.f), outerCoordinateP(0.f, 0.f, 0.f);
    this->GetExtremalCoordinatesFromCache(pParentCluster, innerCoordinateMap, outerCoordinateMap, innerCoordinateP, outerCoordinateP);

    CartesianVector innerCoordinateD(0.f, 0.f, 0.f), outerCoordinateD(0.f, 0.f, 0.f);
    this->GetExtremalCoordinatesFromCache(pDaughterCluster, innerCoordinateMap, outerCoordinateMap, innerCoordinateD, outerCoordinateD);

    for (unsigned int useInnerD = 0; useInnerD < 2; ++useInnerD)
    {
        const CartesianVector daughterVertex(useInnerD == 1 ? innerCoordinateD : outerCoordinateD);
        const CartesianVector daughterEnd(useInnerD == 1 ? outerCoordinateD : innerCoordinateD);

        const float daughterLengthSquared((daughterEnd - daughterVertex).GetMagnitudeSquared());

        // Daughter cluster must be available and below a length cut for any association
        if (daughterLengthSquared > m_maxClusterLength * m_maxClusterLength)
            continue;

        const CartesianVector projectedVertex(LArClusterHelper::GetClosestPosition(daughterVertex, pParentCluster));

        const float daughterVertexDistanceSquared((projectedVertex - daughterVertex).GetMagnitudeSquared());
        const float daughterEndDistanceSquared((projectedVertex - daughterEnd).GetMagnitudeSquared());

        // Cut on proximity of daughter cluster to parent cluster
        if (daughterVertexDistanceSquared > std::min(m_maxLongitudinalDisplacement * m_maxLongitudinalDisplacement,
                                                std::min(daughterEndDistanceSquared, daughterLengthSquared)))
            continue;

        const ClusterAssociation::VertexType daughterVertexType(useInnerD == 1 ? ClusterAssociation::INNER : ClusterAssociation::OUTER);
        const float figureOfMerit(daughterVertexDistanceSquared);

        ClusterAssociation::AssociationType associationType(ClusterAssociation::WEAK);
        ClusterAssociation::VertexType parentVertexType(ClusterAssociation::UNDEFINED);

        for (unsigned int useInnerP = 0; useInnerP < 2; ++useInnerP)
        {
            const CartesianVector parentVertex(useInnerP == 1 ? innerCoordinateP : outerCoordinateP);
            const CartesianVector parentEnd(useInnerP == 1 ? outerCoordinateP : innerCoordinateP);

            const float parentVertexDistanceSquared((projectedVertex - parentVertex).GetMagnitudeSquared());
            const float parentEndDistanceSquared((projectedVertex - parentEnd).GetMagnitudeSquared());
            const float parentLengthSquared((parentEnd - parentVertex).GetMagnitudeSquared());

            if (parentVertexDistanceSquared < parentEndDistanceSquared)
                parentVertexType = (useInnerP == 1 ? ClusterAssociation::INNER : ClusterAssociation::OUTER);

            // Parent cluster must be available and below a length cut for a strong association
            if (!PandoraContentApi::IsAvailable(*this, pParentCluster) || parentLengthSquared > m_maxClusterLength * m_maxClusterLength)
                continue;

            // Require an end-to-end join between parent and daughter cluster
            if (parentVertexDistanceSquared > std::min(m_maxTransverseDisplacement * m_maxTransverseDisplacement,
                                                  std::min(parentEndDistanceSquared, daughterEndDistanceSquared)))
                continue;

            // Cut on pointing information
            const CartesianVector daughterDirection((daughterEnd - daughterVertex).GetUnitVector());
            const CartesianVector parentDirection((parentEnd - projectedVertex).GetUnitVector());

            const float forwardDistance(daughterDirection.GetDotProduct((daughterVertex - projectedVertex)));
            const float sidewaysDistance(daughterDirection.GetCrossProduct((daughterVertex - projectedVertex)).GetMagnitude());

            if (forwardDistance < 0.f || forwardDistance > m_maxLongitudinalDisplacement || sidewaysDistance > m_maxTransverseDisplacement)
                continue;

            if (-parentDirection.GetDotProduct(daughterDirection) < 0.25f)
                continue;

            associationType = ClusterAssociation::STRONG;
            break;
        }

        if (parentVertexType > ClusterAssociation::UNDEFINED)
        {
            (void)clusterAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(
                pDaughterCluster, ClusterAssociation(parentVertexType, daughterVertexType, associationType, figureOfMerit)));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayExtensionAlgorithm::FillClusterMergeMap(const ClusterAssociationMatrix &parentToDaughterMatrix, ClusterMergeMap &clusterMergeMap) const
{
    // Merge parent and daughter clusters if they are strongly associated
    // and the associations have the best figures of merit
    // (i.e. the P --> D association is the best P --> X association,
    //   and the P <-- D association is the best X <-- D association).
    ClusterAssociationMatrix daughterToParentMatrix;

    ClusterVector sortedParentClusters;
    for (const auto &mapEntry : parentToDaughterMatrix)
        sortedParentClusters.push_back(mapEntry.first);
    std::sort(sortedParentClusters.begin(), sortedParentClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : sortedParentClusters)
    {
        const ClusterAssociationMap &daughterToAssociationMap(parentToDaughterMatrix.at(pParentCluster));

        ClusterVector sortedLocalDaughterClusters;
        for (const auto &mapEntry : daughterToAssociationMap)
            sortedLocalDaughterClusters.push_back(mapEntry.first);
        std::sort(sortedLocalDaughterClusters.begin(), sortedLocalDaughterClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pDaughterCluster : sortedLocalDaughterClusters)
        {
            const ClusterAssociation &clusterAssociation(daughterToAssociationMap.at(pDaughterCluster));
            (void)daughterToParentMatrix[pDaughterCluster].insert(ClusterAssociationMap::value_type(pParentCluster, clusterAssociation));
        }
    }

    ClusterAssociationMatrix reducedParentToDaughterMatrix;

    ClusterVector sortedDaughterClusters;
    for (const auto &mapEntry : daughterToParentMatrix)
        sortedDaughterClusters.push_back(mapEntry.first);
    std::sort(sortedDaughterClusters.begin(), sortedDaughterClusters.end(), LArClusterHelper::SortByNHits);

    // Loop over parent clusters and select nearby daughter clusters that are closer than another parent cluster
    for (const Cluster *const pDaughterCluster : sortedDaughterClusters)
    {
        const ClusterAssociationMap &parentToAssociationMap(daughterToParentMatrix.at(pDaughterCluster));

        const Cluster *pBestInner(NULL);
        const Cluster *pBestOuter(NULL);

        float bestFomInner(std::numeric_limits<float>::max());
        float bestFomOuter(std::numeric_limits<float>::max());

        ClusterVector sortedLocalParentClusters;
        for (const auto &mapEntry : parentToAssociationMap)
            sortedLocalParentClusters.push_back(mapEntry.first);
        std::sort(sortedLocalParentClusters.begin(), sortedLocalParentClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pParentCluster : sortedLocalParentClusters)
        {
            const ClusterAssociation &clusterAssociation(parentToAssociationMap.at(pParentCluster));

            if (clusterAssociation.GetParent() == ClusterAssociation::INNER)
            {
                if (clusterAssociation.GetFigureOfMerit() < bestFomInner)
                {
                    bestFomInner = clusterAssociation.GetFigureOfMerit();

                    if (clusterAssociation.GetAssociation() == ClusterAssociation::STRONG)
                    {
                        pBestInner = pParentCluster;
                    }
                    else
                    {
                        pBestInner = NULL;
                    }
                }
            }

            if (clusterAssociation.GetParent() == ClusterAssociation::OUTER)
            {
                if (clusterAssociation.GetFigureOfMerit() < bestFomOuter)
                {
                    bestFomOuter = clusterAssociation.GetFigureOfMerit();

                    if (clusterAssociation.GetAssociation() == ClusterAssociation::STRONG)
                    {
                        pBestOuter = pParentCluster;
                    }
                    else
                    {
                        pBestOuter = NULL;
                    }
                }
            }
        }

        if (pBestInner)
        {
            ClusterAssociationMatrix::const_iterator iter3A = parentToDaughterMatrix.find(pBestInner);

            if (parentToDaughterMatrix.end() == iter3A)
                throw pandora::StatusCodeException(STATUS_CODE_FAILURE);

            const ClusterAssociationMap &parentToDaughterMap(iter3A->second);
            ClusterAssociationMap::const_iterator iter3B = parentToDaughterMap.find(pDaughterCluster);

            if (parentToDaughterMap.end() == iter3B)
                throw pandora::StatusCodeException(STATUS_CODE_FAILURE);

            const ClusterAssociation &bestAssociationInner(iter3B->second);
            (void)reducedParentToDaughterMatrix[pBestInner].insert(ClusterAssociationMap::value_type(pDaughterCluster, bestAssociationInner));
        }

        if (pBestOuter)
        {
            ClusterAssociationMatrix::const_iterator iter3A = parentToDaughterMatrix.find(pBestOuter);

            if (parentToDaughterMatrix.end() == iter3A)
                throw pandora::StatusCodeException(STATUS_CODE_FAILURE);

            const ClusterAssociationMap &parentToDaughterMap(iter3A->second);
            ClusterAssociationMap::const_iterator iter3B = parentToDaughterMap.find(pDaughterCluster);

            if (parentToDaughterMap.end() == iter3B)
                throw pandora::StatusCodeException(STATUS_CODE_FAILURE);

            const ClusterAssociation &bestAssociationOuter(iter3B->second);
            (void)reducedParentToDaughterMatrix[pBestOuter].insert(ClusterAssociationMap::value_type(pDaughterCluster, bestAssociationOuter));
        }
    }

    ClusterVector sortedReducedParentClusters;
    for (const auto &mapEntry : reducedParentToDaughterMatrix)
        sortedReducedParentClusters.push_back(mapEntry.first);
    std::sort(sortedReducedParentClusters.begin(), sortedReducedParentClusters.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pParentCluster : sortedReducedParentClusters)
    {
        const ClusterAssociationMap &daughterToAssociationMap(reducedParentToDaughterMatrix.at(pParentCluster));

        const Cluster *pBestInner(NULL);
        const Cluster *pBestOuter(NULL);

        float bestFomInner(std::numeric_limits<float>::max());
        float bestFomOuter(std::numeric_limits<float>::max());

        ClusterVector sortedLocalDaughterClusters;
        for (const auto &mapEntry : daughterToAssociationMap)
            sortedLocalDaughterClusters.push_back(mapEntry.first);
        std::sort(sortedLocalDaughterClusters.begin(), sortedLocalDaughterClusters.end(), LArClusterHelper::SortByNHits);

        for (const Cluster *const pDaughterCluster : sortedLocalDaughterClusters)
        {
            const ClusterAssociation &clusterAssociation(daughterToAssociationMap.at(pDaughterCluster));

            if (clusterAssociation.GetParent() == ClusterAssociation::INNER)
            {
                if (clusterAssociation.GetFigureOfMerit() < bestFomInner)
                {
                    bestFomInner = clusterAssociation.GetFigureOfMerit();

                    if (clusterAssociation.GetAssociation() == ClusterAssociation::STRONG)
                    {
                        pBestInner = pDaughterCluster;
                    }
                    else
                    {
                        pBestInner = NULL;
                    }
                }
            }

            if (clusterAssociation.GetParent() == ClusterAssociation::OUTER)
            {
                if (clusterAssociation.GetFigureOfMerit() < bestFomOuter)
                {
                    bestFomOuter = clusterAssociation.GetFigureOfMerit();

                    if (clusterAssociation.GetAssociation() == ClusterAssociation::STRONG)
                    {
                        pBestOuter = pDaughterCluster;
                    }
                    else
                    {
                        pBestOuter = NULL;
                    }
                }
            }
        }

        if (pBestInner)
        {
            ClusterList &parentList(clusterMergeMap[pParentCluster]);

            if (parentList.end() == std::find(parentList.begin(), parentList.end(), pBestInner))
                parentList.push_back(pBestInner);

            ClusterList &bestInnerList(clusterMergeMap[pBestInner]);

            if (bestInnerList.end() == std::find(bestInnerList.begin(), bestInnerList.end(), pParentCluster))
                bestInnerList.push_back(pParentCluster);
        }

        if (pBestOuter)
        {
            ClusterList &parentList(clusterMergeMap[pParentCluster]);

            if (parentList.end() == std::find(parentList.begin(), parentList.end(), pBestOuter))
                parentList.push_back(pBestOuter);

            ClusterList &bestOuterList(clusterMergeMap[pBestOuter]);

            if (bestOuterList.end() == std::find(bestOuterList.begin(), bestOuterList.end(), pParentCluster))
                bestOuterList.push_back(pParentCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxClusterLength", m_maxClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
