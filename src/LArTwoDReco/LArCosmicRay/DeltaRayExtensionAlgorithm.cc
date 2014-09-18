/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/DeltaRayExtensionAlgorithm.cc
 *
 *  @brief  Implementation of the delta ray extension algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArTwoDReco/LArCosmicRay/DeltaRayExtensionAlgorithm.h"

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
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayExtensionAlgorithm::FillClusterAssociationMatrix(const ClusterVector &clusterVector, ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterEndI = clusterVector.end(); iterI != iterEndI; ++iterI)
    {
        const Cluster *pParentCluster = *iterI;

        for (ClusterVector::const_iterator iterJ = clusterVector.begin(), iterEndJ = clusterVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            const Cluster *pDaughterCluster = *iterJ;

            if (pParentCluster == pDaughterCluster)
                continue;

            this->FillClusterAssociationMatrix(pParentCluster, pDaughterCluster, clusterAssociationMatrix);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayExtensionAlgorithm::FillClusterAssociationMatrix(const Cluster *const pParentCluster, const Cluster *const pDaughterCluster,
    ClusterAssociationMatrix &clusterAssociationMatrix) const
{
    // Daughter cluster must be available for any association to proceed
    if (!pDaughterCluster->IsAvailable())
        return;

    // Weak association:   between parent cosmic-ray muon and daughter delta ray
    // Strong association: between parent and daughter fragments of delta ray
    // Figure of merit:    distance between parent and daughter clusters
    CartesianVector innerCoordinateP(0.f,0.f,0.f), outerCoordinateP(0.f,0.f,0.f);
    CartesianVector innerCoordinateD(0.f,0.f,0.f), outerCoordinateD(0.f,0.f,0.f);
    LArClusterHelper::GetExtremalCoordinatesXZ(pParentCluster, innerCoordinateP, outerCoordinateP);
    LArClusterHelper::GetExtremalCoordinatesXZ(pDaughterCluster, innerCoordinateD, outerCoordinateD);

    for (unsigned int useInnerD = 0; useInnerD < 2; ++useInnerD)
    {
        const CartesianVector daughterVertex(useInnerD == 1 ? innerCoordinateD : outerCoordinateD);
        const CartesianVector daughterEnd(useInnerD == 1 ? outerCoordinateD : innerCoordinateD);
        const CartesianVector projectedVertex(LArClusterHelper::GetClosestPosition(daughterVertex, pParentCluster));

        const float daughterVertexDistanceSquared((projectedVertex - daughterVertex).GetMagnitudeSquared());
        const float daughterEndDistanceSquared((projectedVertex - daughterEnd).GetMagnitudeSquared());
        const float daughterLengthSquared((daughterEnd - daughterVertex).GetMagnitudeSquared());

        // Daughter cluster must be available and below a length cut for any association
        if (!pDaughterCluster->IsAvailable() || daughterLengthSquared > m_maxClusterLength * m_maxClusterLength)
            continue;

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
            if (!pParentCluster->IsAvailable() || parentLengthSquared > m_maxClusterLength * m_maxClusterLength)
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

            if (-parentDirection.GetDotProduct(daughterDirection) < 0.25)
                continue;

            associationType = ClusterAssociation::STRONG;
            break;
        }

        if (parentVertexType > ClusterAssociation::UNDEFINED)
        {
            (void) clusterAssociationMatrix[pParentCluster].insert(ClusterAssociationMap::value_type(pDaughterCluster,
                ClusterAssociation(parentVertexType, daughterVertexType, associationType, figureOfMerit)));
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

    for (ClusterAssociationMatrix::const_iterator iter1 = parentToDaughterMatrix.begin(), iterEnd1 = parentToDaughterMatrix.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster *pParentCluster(iter1->first);
        const ClusterAssociationMap &associationMap(iter1->second);

        for (ClusterAssociationMap::const_iterator iter2 = associationMap.begin(), iterEnd2 = associationMap.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster *pDaughterCluster(iter2->first);
            const ClusterAssociation &association(iter2->second);

            const ClusterAssociation::VertexType parentVertexType(association.GetParent());
            const ClusterAssociation::VertexType daughterVertexType(association.GetDaughter());
            const ClusterAssociation::AssociationType associationType(association.GetAssociation());
            const float figureOfMerit(association.GetFigureOfMerit());

            (void) daughterToParentMatrix[pDaughterCluster].insert(ClusterAssociationMap::value_type(pParentCluster,
                ClusterAssociation(daughterVertexType, parentVertexType, associationType, figureOfMerit)));
        }
    }


    ClusterAssociationMatrix reducedParentToDaughterMatrix;

    for (ClusterAssociationMatrix::const_iterator iter1 = daughterToParentMatrix.begin(), iterEnd1 = daughterToParentMatrix.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster *pDaughterCluster(iter1->first);
        const ClusterAssociationMap &associationMap(iter1->second);

        Cluster *pBestInner(NULL);
        Cluster *pBestOuter(NULL);

        float bestFomInner(std::numeric_limits<float>::max());
        float bestFomOuter(std::numeric_limits<float>::max());

        for (ClusterAssociationMap::const_iterator iter2 = associationMap.begin(), iterEnd2 = associationMap.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster *pParentCluster(iter2->first);
            const ClusterAssociation &association(iter2->second);

            if (association.GetParent() == ClusterAssociation::INNER)
            {
                if (association.GetFigureOfMerit() < bestFomInner)
                {
                    bestFomInner = association.GetFigureOfMerit();

                    if (association.GetAssociation() == ClusterAssociation::STRONG)
                    {
                        pBestInner = (Cluster*)pParentCluster;
                    }
                    else
                    {
                        pBestInner = NULL;
                    }
                }
            }

            if (association.GetParent() == ClusterAssociation::OUTER)
            {
                if (association.GetFigureOfMerit() < bestFomOuter)
                {
                    bestFomOuter = association.GetFigureOfMerit();

                    if (association.GetAssociation() == ClusterAssociation::STRONG)
                    {
                        pBestOuter = (Cluster*)pParentCluster;
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
            (void) reducedParentToDaughterMatrix[pBestInner].insert(ClusterAssociationMap::value_type(pDaughterCluster, bestAssociationInner));
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
            (void) reducedParentToDaughterMatrix[pBestOuter].insert(ClusterAssociationMap::value_type(pDaughterCluster, bestAssociationOuter));
        }
    }


    for (ClusterAssociationMatrix::const_iterator iter1 = reducedParentToDaughterMatrix.begin(), iterEnd1 = reducedParentToDaughterMatrix.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster *pParentCluster(iter1->first);
        const ClusterAssociationMap &associationMap(iter1->second);

        Cluster *pBestInner(NULL);
        Cluster *pBestOuter(NULL);

        float bestFomInner(std::numeric_limits<float>::max());
        float bestFomOuter(std::numeric_limits<float>::max());

        for (ClusterAssociationMap::const_iterator iter2 = associationMap.begin(), iterEnd2 = associationMap.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster *pDaughterCluster(iter2->first);
            const ClusterAssociation &association(iter2->second);

            if (association.GetParent() == ClusterAssociation::INNER)
            {
                if (association.GetFigureOfMerit() < bestFomInner)
                {
                    bestFomInner = association.GetFigureOfMerit();

                    if (association.GetAssociation() == ClusterAssociation::STRONG)
                    {
                        pBestInner = (Cluster*)pDaughterCluster;
                    }
                    else
                    {
                        pBestInner = NULL;
                    }
                }
            }

            if (association.GetParent() == ClusterAssociation::OUTER)
            {
                if (association.GetFigureOfMerit() < bestFomOuter)
                {
                    bestFomOuter = association.GetFigureOfMerit();

                    if (association.GetAssociation() == ClusterAssociation::STRONG)
                    {
                        pBestOuter = (Cluster*)pDaughterCluster;
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
            clusterMergeMap[pParentCluster].insert(pBestInner);
            clusterMergeMap[pBestInner].insert((Cluster*)pParentCluster);
        }

        if (pBestOuter)
        {
            clusterMergeMap[pParentCluster].insert(pBestOuter);
            clusterMergeMap[pBestOuter].insert((Cluster*)pParentCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxClusterLength", m_maxClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    return ClusterExtensionAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
