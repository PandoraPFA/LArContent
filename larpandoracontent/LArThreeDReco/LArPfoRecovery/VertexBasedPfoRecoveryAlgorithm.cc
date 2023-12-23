/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/VertexBasedPfoRecoveryAlgorithm.cc
 *
 *  @brief  Implementation of the vertex-based particle recovery algorithm
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArPfoRecovery/VertexBasedPfoRecoveryAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

using namespace pandora;

namespace lar_content
{

VertexBasedPfoRecoveryAlgorithm::VertexBasedPfoRecoveryAlgorithm() :
    m_slidingFitHalfWindow(10),
    m_maxLongitudinalDisplacement(5.f),
    m_maxTransverseDisplacement(2.f),
    m_twoViewChi2Cut(5.f),
    m_threeViewChi2Cut(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexBasedPfoRecoveryAlgorithm::Run()
{
    const VertexList *pVertexList = NULL;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetCurrentList(*this, pVertexList));

    const Vertex *const pSelectedVertex(
        (pVertexList && (pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : NULL);

    if (!pSelectedVertex)
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "VertexBasedPfoRecoveryAlgorithm: unable to find vertex in current list " << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    // Get the available clusters from each view
    ClusterVector availableClusters;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetAvailableClusters(m_inputClusterListNames, availableClusters));

    // Build a set of sliding fit results
    TwoDSlidingFitResultMap slidingFitResultMap;
    this->BuildSlidingFitResultMap(availableClusters, slidingFitResultMap);

    // Select seed clusters (adjacent to vertex)
    ClusterVector selectedClusters;
    this->SelectVertexClusters(pSelectedVertex, slidingFitResultMap, availableClusters, selectedClusters);
   
    //std::cout << "Selected Clusters !!!: " << selectedClusters.size() << std::endl;
    // Match the cluster end points
    ClusterSet vetoList;
    ParticleList particleList;
    //this->MatchThreeViews(pSelectedVertex, slidingFitResultMap, availableClusters, vetoList, particleList);
    this->MatchTwoViews(pSelectedVertex, slidingFitResultMap, availableClusters, vetoList, particleList);

    // Build new particles
    this->BuildParticles(particleList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexBasedPfoRecoveryAlgorithm::GetAvailableClusters(const StringVector inputClusterListNames, ClusterVector &clusterVector) const
{
    for (StringVector::const_iterator iter = inputClusterListNames.begin(), iterEnd = inputClusterListNames.end(); iter != iterEnd; ++iter)
    {
        const ClusterList *pClusterList(NULL);
        PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, *iter, pClusterList));

       // std::cout << "No. of Cluster: " << pClusterList->size() << std::endl;
	if (NULL == pClusterList)
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "VertexBasedPfoRecoveryAlgorithm: could not find cluster list " << *iter << std::endl;
            continue;
        }

        for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
        {
            const Cluster *const pCluster = *cIter;
           // std::cout << "Cluster Type: " << pCluster->GetNCaloHits() << std::endl;
	    if (!pCluster->IsAvailable())
                continue;

            if (pCluster->GetNCaloHits() <= 0)
                continue;

            if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
                continue;

	    clusterVector.push_back(pCluster);
        }
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

    return STATUS_CODE_SUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoRecoveryAlgorithm::BuildSlidingFitResultMap(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    //add in modified FitPitch value of 0.68 

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        if (slidingFitResultMap.end() == slidingFitResultMap.find(*iter))
        {
            try
            {
                const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(*iter)));
                const TwoDSlidingFitResult slidingFitResult(*iter, m_slidingFitHalfWindow, slidingFitPitch);
                const LArPointingCluster pointingCluster(slidingFitResult);

                if (pointingCluster.GetLengthSquared() < std::numeric_limits<float>::epsilon())
                    continue;

                if (!slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoRecoveryAlgorithm::SelectVertexClusters(const Vertex *const pVertex, const TwoDSlidingFitResultMap &slidingFitResultMap,
    const ClusterVector &inputClusters, ClusterVector &outputClusters) const
{
    const CartesianVector vertexU(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U));
    const CartesianVector vertexV(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V));
    const CartesianVector vertexW(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W));
    //std::cout << "Vertex Location: " << vertexU << vertexV << vertexW << std::endl;

    for (ClusterVector::const_iterator cIter = inputClusters.begin(), cIterEnd = inputClusters.end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster = *cIter;
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

      //  std::cout << "!!!HIT TYPE: " << hitType << std::endl;
        if (TPC_3D == hitType)
            continue;
      //  std::cout << "Number of Hits: " << pCluster->GetNCaloHits() << std::endl;
	 
        const CartesianVector vertexPosition((TPC_VIEW_U == hitType) ? vertexU : (TPC_VIEW_V == hitType) ? vertexV : vertexW);

        TwoDSlidingFitResultMap::const_iterator sIter = slidingFitResultMap.find(pCluster);
        if (slidingFitResultMap.end() == sIter)
            continue;

        const TwoDSlidingFitResult &slidingFitResult = sIter->second;
        const LArPointingCluster pointingCluster(slidingFitResult);

        for (unsigned int iVtx = 0; iVtx < 2; ++iVtx)
        {
            const LArPointingCluster::Vertex &pointingVertex((0 == iVtx) ? pointingCluster.GetInnerVertex() : pointingCluster.GetOuterVertex());

            float rL(0.f), rT(0.f);
            LArPointingClusterHelper::GetImpactParameters(pointingVertex, vertexPosition, rL, rT);
           // std::cout << "Impact Parameters: " << rL << " & " << rT << std::endl;
           // std::cout << "Cuts: " << m_maxLongitudinalDisplacement << " & " << m_maxTransverseDisplacement << std::endl;

            if (rL > -100.f && rL < m_maxLongitudinalDisplacement && rT < m_maxTransverseDisplacement)
            {
	//	std::cout << "Cluster Pass" << std::endl;
		outputClusters.push_back(pCluster);
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoRecoveryAlgorithm::MatchThreeViews(const Vertex *const pVertex, const TwoDSlidingFitResultMap &slidingFitResultMap,
    const ClusterVector &inputClusters, ClusterSet &vetoList, ParticleList &particleList) const
{
    while (true)
    {
        ClusterVector availableClusters, clustersU, clustersV, clustersW;
        this->SelectAvailableClusters(vetoList, inputClusters, availableClusters);
        this->SelectClusters(TPC_VIEW_U, availableClusters, clustersU);
        this->SelectClusters(TPC_VIEW_V, availableClusters, clustersV);
        this->SelectClusters(TPC_VIEW_W, availableClusters, clustersW);

        std::cout << clustersU.size() << clustersV.size() << clustersW.size() << std::endl;

	float chi2(m_threeViewChi2Cut);
        //std::cout << "three View Cut: " << m_threeViewChi2Cut << std::endl;
       	const Cluster *pCluster1(NULL);
        const Cluster *pCluster2(NULL);
        const Cluster *pCluster3(NULL);

        this->GetBestChi2(pVertex, slidingFitResultMap, clustersU, clustersV, clustersW, pCluster1, pCluster2, pCluster3, chi2);
        
        //std::cout << "U: " << clustersU.size() << " : " << pCluster1 << " V: " << clustersV.size() << " : " << pCluster2 << " W:" << clustersW.size() << " : " << pCluster3 << std::endl;
        if (NULL == pCluster1 || NULL == pCluster2 || NULL == pCluster3)
            return;

        const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
        const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
        const HitType hitType3(LArClusterHelper::GetClusterHitType(pCluster3));

        const Cluster *const pClusterU((TPC_VIEW_U == hitType1) ? pCluster1
                : (TPC_VIEW_U == hitType2)                      ? pCluster2
                : (TPC_VIEW_U == hitType3)                      ? pCluster3
                                                                : NULL);
        const Cluster *const pClusterV((TPC_VIEW_V == hitType1) ? pCluster1
                : (TPC_VIEW_V == hitType2)                      ? pCluster2
                : (TPC_VIEW_V == hitType3)                      ? pCluster3
                                                                : NULL);
        const Cluster *const pClusterW((TPC_VIEW_W == hitType1) ? pCluster1
                : (TPC_VIEW_W == hitType2)                      ? pCluster2
                : (TPC_VIEW_W == hitType3)                      ? pCluster3
                                                                : NULL);

        //std::cout << "UU: " << pClusterU << " VV: " << pClusterV << " WW:" << pClusterW << std::endl;
	particleList.push_back(Particle(pClusterU, pClusterV, pClusterW));

        vetoList.insert(pCluster1);
        vetoList.insert(pCluster2);
        vetoList.insert(pCluster3);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoRecoveryAlgorithm::MatchTwoViews(const Vertex *const pVertex, const TwoDSlidingFitResultMap &slidingFitResultMap,
    const ClusterVector &inputClusters, ClusterSet &vetoList, ParticleList &particleList) const
{
    while (true)
    {
        ClusterVector availableClusters, clustersU, clustersV, clustersW;
        this->SelectAvailableClusters(vetoList, inputClusters, availableClusters);
        this->SelectClusters(TPC_VIEW_U, availableClusters, clustersU);
        this->SelectClusters(TPC_VIEW_V, availableClusters, clustersV);
        this->SelectClusters(TPC_VIEW_W, availableClusters, clustersW);

	//std::cout << clustersU.size() << clustersV.size() << clustersW.size() << std::endl;

	float chi2(m_twoViewChi2Cut);
        const Cluster *pCluster1(NULL);
        const Cluster *pCluster2(NULL);

        this->GetBestChi2(pVertex, slidingFitResultMap, clustersU, clustersV, pCluster1, pCluster2, chi2);
        this->GetBestChi2(pVertex, slidingFitResultMap, clustersV, clustersW, pCluster1, pCluster2, chi2);
        this->GetBestChi2(pVertex, slidingFitResultMap, clustersW, clustersU, pCluster1, pCluster2, chi2);

	//if ( pCluster1 != NULL)
	//    std::cout << " *** " << pCluster1->GetNCaloHits() << std::endl;

	//if ( pCluster2 !=  NULL)
        //   std::cout << " *** " << pCluster2->GetNCaloHits() << std::endl;

       	if (NULL == pCluster1 || NULL == pCluster2)
            return;

        const HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1));
        const HitType hitType2(LArClusterHelper::GetClusterHitType(pCluster2));

        const Cluster *const pClusterU((TPC_VIEW_U == hitType1) ? pCluster1 : (TPC_VIEW_U == hitType2) ? pCluster2 : NULL);
        const Cluster *const pClusterV((TPC_VIEW_V == hitType1) ? pCluster1 : (TPC_VIEW_V == hitType2) ? pCluster2 : NULL);
        const Cluster *const pClusterW((TPC_VIEW_W == hitType1) ? pCluster1 : (TPC_VIEW_W == hitType2) ? pCluster2 : NULL);

	//std::cout << "We made one " << std::endl;
        particleList.push_back(Particle(pClusterU, pClusterV, pClusterW));
        
        vetoList.insert(pCluster1);
        vetoList.insert(pCluster2);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoRecoveryAlgorithm::GetBestChi2(const Vertex *const pVertex, const TwoDSlidingFitResultMap &slidingFitResultMap,
    const ClusterVector &clusters1, const ClusterVector &clusters2, const ClusterVector &clusters3, const Cluster *&pBestCluster1,
    const Cluster *&pBestCluster2, const Cluster *&pBestCluster3, float &bestChi2) const
{
    //std::cout << "1: " << clusters1.empty() << " 2: " << clusters2.empty() << " 3:" << clusters3.empty() << std::endl;	
    if (clusters1.empty() || clusters2.empty() || clusters3.empty())
        return;

    // First loop
    for (ClusterVector::const_iterator cIter1 = clusters1.begin(), cIterEnd1 = clusters1.end(); cIter1 != cIterEnd1; ++cIter1)
    {
        const Cluster *const pCluster1 = *cIter1;

        TwoDSlidingFitResultMap::const_iterator sIter1 = slidingFitResultMap.find(pCluster1);
        if (slidingFitResultMap.end() == sIter1)
            continue;

        const TwoDSlidingFitResult &slidingFitResult1 = sIter1->second;
        const LArPointingCluster pointingCluster1(slidingFitResult1);

        // Second loop
        for (ClusterVector::const_iterator cIter2 = clusters2.begin(), cIterEnd2 = clusters2.end(); cIter2 != cIterEnd2; ++cIter2)
        {
            const Cluster *const pCluster2 = *cIter2;

            TwoDSlidingFitResultMap::const_iterator sIter2 = slidingFitResultMap.find(pCluster2);
            if (slidingFitResultMap.end() == sIter2)
                continue;

            const TwoDSlidingFitResult &slidingFitResult2 = sIter2->second;
            const LArPointingCluster pointingCluster2(slidingFitResult2);

            // Third loop
            for (ClusterVector::const_iterator cIter3 = clusters3.begin(), cIterEnd3 = clusters3.end(); cIter3 != cIterEnd3; ++cIter3)
            {
                const Cluster *const pCluster3 = *cIter3;

                TwoDSlidingFitResultMap::const_iterator sIter3 = slidingFitResultMap.find(pCluster3);
                if (slidingFitResultMap.end() == sIter3)
                    continue;

                const TwoDSlidingFitResult &slidingFitResult3 = sIter3->second;
                const LArPointingCluster pointingCluster3(slidingFitResult3);

                // Calculate chi-squared
                const float thisChi2(this->GetChi2(pVertex, pointingCluster1, pointingCluster2, pointingCluster3));
                std::cout << "thisChi2: " << thisChi2 << "bestChi2: " << bestChi2 << std::endl;

                if (thisChi2 < bestChi2)
                {
                    bestChi2 = thisChi2;
                    pBestCluster1 = pCluster1;
                    pBestCluster2 = pCluster2;
                    pBestCluster3 = pCluster3;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoRecoveryAlgorithm::GetBestChi2(const Vertex *const pVertex, const TwoDSlidingFitResultMap &slidingFitResultMap,
    const ClusterVector &clusters1, const ClusterVector &clusters2, const Cluster *&pBestCluster1, const Cluster *&pBestCluster2, float &bestChi2) const
{
   // std::cout << "BestChi2 in two views: " << clusters1.size() << " clusters 2 size : " << clusters2.size() << std::endl;
    if (clusters1.empty() || clusters2.empty())
        return;

    // First loop
    for (ClusterVector::const_iterator cIter1 = clusters1.begin(), cIterEnd1 = clusters1.end(); cIter1 != cIterEnd1; ++cIter1)
    {
        const Cluster *const pCluster1 = *cIter1;

        TwoDSlidingFitResultMap::const_iterator sIter1 = slidingFitResultMap.find(pCluster1);
        if (slidingFitResultMap.end() == sIter1)
            continue;

        const TwoDSlidingFitResult &slidingFitResult1 = sIter1->second;
        const LArPointingCluster pointingCluster1(slidingFitResult1);
      //  std::cout << "Running first Bestchi2 Loop" << std::endl;
        // Second loop
        for (ClusterVector::const_iterator cIter2 = clusters2.begin(), cIterEnd2 = clusters2.end(); cIter2 != cIterEnd2; ++cIter2)
        {
            const Cluster *const pCluster2 = *cIter2;

            TwoDSlidingFitResultMap::const_iterator sIter2 = slidingFitResultMap.find(pCluster2);
        //    std::cout << "Running second Bestchi2 Loop" << std::endl;
	    if (slidingFitResultMap.end() == sIter2)
                continue;

            const TwoDSlidingFitResult &slidingFitResult2 = sIter2->second;
            const LArPointingCluster pointingCluster2(slidingFitResult2);

            // Calculate chi-squared
            const float thisChi2(this->GetChi2(pVertex, pointingCluster1, pointingCluster2));
          //  std::cout << "thisChi2: " << thisChi2 << "bestChi2: " << bestChi2 << std::endl;

            if (thisChi2 < bestChi2)
            {
                bestChi2 = thisChi2;
                pBestCluster1 = pCluster1;
                pBestCluster2 = pCluster2;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexBasedPfoRecoveryAlgorithm::GetChi2(
    const Vertex *const pVertex, const LArPointingCluster &pointingCluster1, const LArPointingCluster &pointingCluster2) const
{
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pointingCluster1.GetCluster()));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pointingCluster2.GetCluster()));

    if (hitType1 == hitType2)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const CartesianVector vertex1(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType1));
    const CartesianVector vertex2(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType2));

    const LArPointingCluster::Vertex &pointingVertex1(this->GetOuterVertex(vertex1, pointingCluster1));
    const LArPointingCluster::Vertex &pointingVertex2(this->GetOuterVertex(vertex2, pointingCluster2));

    float chi2(0.f);
    CartesianVector mergedPosition(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions3D(
        this->GetPandora(), hitType1, hitType2, pointingVertex1.GetPosition(), pointingVertex2.GetPosition(), mergedPosition, chi2);

    return chi2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float VertexBasedPfoRecoveryAlgorithm::GetChi2(const Vertex *const pVertex, const LArPointingCluster &pointingCluster1,
    const LArPointingCluster &pointingCluster2, const LArPointingCluster &pointingCluster3) const
{
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pointingCluster1.GetCluster()));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pointingCluster2.GetCluster()));
    const HitType hitType3(LArClusterHelper::GetClusterHitType(pointingCluster3.GetCluster()));

    if ((hitType1 == hitType2) || (hitType2 == hitType3) || (hitType3 == hitType1))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const CartesianVector vertex1(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType1));
    const CartesianVector vertex2(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType2));
    const CartesianVector vertex3(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), hitType3));

    const LArPointingCluster::Vertex &pointingVertex1(this->GetOuterVertex(vertex1, pointingCluster1));
    const LArPointingCluster::Vertex &pointingVertex2(this->GetOuterVertex(vertex2, pointingCluster2));
    const LArPointingCluster::Vertex &pointingVertex3(this->GetOuterVertex(vertex3, pointingCluster3));

    float chi2(0.f);
    CartesianVector mergedPosition(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), hitType1, hitType2, hitType3, pointingVertex1.GetPosition(),
        pointingVertex2.GetPosition(), pointingVertex3.GetPosition(), mergedPosition, chi2);

    return chi2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoRecoveryAlgorithm::SelectAvailableClusters(const ClusterSet &vetoList, const ClusterVector &inputVector, ClusterVector &outputVector) const
{
    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        if (0 == vetoList.count(*iter))
	    outputVector.push_back(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoRecoveryAlgorithm::SelectClusters(const HitType hitType, const ClusterVector &inputVector, ClusterVector &outputVector) const
{
    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        if (hitType == LArClusterHelper::GetClusterHitType(*iter))
            outputVector.push_back(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArPointingCluster::Vertex &VertexBasedPfoRecoveryAlgorithm::GetInnerVertex(const CartesianVector &vertex, const LArPointingCluster &cluster) const
{
    const float innerDistance((vertex - cluster.GetInnerVertex().GetPosition()).GetMagnitudeSquared());
    const float outerDistance((vertex - cluster.GetOuterVertex().GetPosition()).GetMagnitudeSquared());

    if (innerDistance < outerDistance)
        return cluster.GetInnerVertex();
    else
        return cluster.GetOuterVertex();
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArPointingCluster::Vertex &VertexBasedPfoRecoveryAlgorithm::GetOuterVertex(const CartesianVector &vertex, const LArPointingCluster &cluster) const
{
    const LArPointingCluster::Vertex &innerVertex = this->GetInnerVertex(vertex, cluster);

    if (innerVertex.IsInnerVertex())
        return cluster.GetOuterVertex();
    else
        return cluster.GetInnerVertex();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void VertexBasedPfoRecoveryAlgorithm::BuildParticles(const ParticleList &particleList)
{
     
    std::cout << "*** : " << particleList.size() << std::endl;	
    if (particleList.empty())
        return;
   
    const PfoList *pPfoList = NULL;
    std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (ParticleList::const_iterator iter = particleList.begin(), iterEnd = particleList.end(); iter != iterEnd; ++iter)
    {
        const Particle &particle = *iter;

        ClusterList clusterList;
        const Cluster *const pClusterU = particle.m_pClusterU;
        const Cluster *const pClusterV = particle.m_pClusterV;
        const Cluster *const pClusterW = particle.m_pClusterW;

        const bool isAvailableU((NULL != pClusterU) ? pClusterU->IsAvailable() : true);
        const bool isAvailableV((NULL != pClusterV) ? pClusterV->IsAvailable() : true);
        const bool isAvailableW((NULL != pClusterW) ? pClusterW->IsAvailable() : true);

        if (!(isAvailableU && isAvailableV && isAvailableW))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        if (pClusterU)
            clusterList.push_back(pClusterU);
        if (pClusterV)
            clusterList.push_back(pClusterV);
        if (pClusterW)
            clusterList.push_back(pClusterW);

        // TODO Correct these placeholder parameters
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = MU_MINUS; // TRACK
        pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
        pfoParameters.m_clusterList = clusterList;

        const ParticleFlowObject *pPfo(NULL);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
    }

    if (!pPfoList->empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

VertexBasedPfoRecoveryAlgorithm::Particle::Particle(const Cluster *const pClusterU, const Cluster *const pClusterV, const Cluster *const pClusterW) :
    m_pClusterU(pClusterU),
    m_pClusterV(pClusterV),
    m_pClusterW(pClusterW)
{
    if (NULL == m_pClusterU && NULL == m_pClusterV && NULL == m_pClusterW)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const HitType hitTypeU(NULL == m_pClusterU ? TPC_VIEW_U : LArClusterHelper::GetClusterHitType(m_pClusterU));
    const HitType hitTypeV(NULL == m_pClusterV ? TPC_VIEW_V : LArClusterHelper::GetClusterHitType(m_pClusterV));
    const HitType hitTypeW(NULL == m_pClusterW ? TPC_VIEW_W : LArClusterHelper::GetClusterHitType(m_pClusterW));

    if (!(TPC_VIEW_U == hitTypeU && TPC_VIEW_V == hitTypeV && TPC_VIEW_W == hitTypeW))
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexBasedPfoRecoveryAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitHalfWindow", m_slidingFitHalfWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TwoViewChi2Cut", m_twoViewChi2Cut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ThreeViewChi2Cut", m_threeViewChi2Cut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
