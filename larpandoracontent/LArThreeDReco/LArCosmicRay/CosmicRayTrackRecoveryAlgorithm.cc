/**
 *  @file   larpandoracontent/LArTwoDReco/LArCosmicRay/CosmicRayTrackRecoveryAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray longitudinal track recovery algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/LArCosmicRay/CosmicRayTrackRecoveryAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

using namespace pandora;

namespace lar_content
{

CosmicRayTrackRecoveryAlgorithm::CosmicRayTrackRecoveryAlgorithm() :
    m_clusterMinLength(10.f),
    m_clusterMinSpanZ(2.f),
    m_clusterMinOverlapX(6.f),
    m_clusterMaxDeltaX(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackRecoveryAlgorithm::Run()
{
    // Get the available clusters for each view
    ClusterVector availableClustersU, availableClustersV, availableClustersW;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetAvailableClusters(m_inputClusterListNameU, availableClustersU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetAvailableClusters(m_inputClusterListNameV, availableClustersV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetAvailableClusters(m_inputClusterListNameW, availableClustersW));

    // Select clean clusters in each view
    ClusterVector cleanClustersU, cleanClustersV, cleanClustersW;
    this->SelectCleanClusters(availableClustersU, cleanClustersU);
    this->SelectCleanClusters(availableClustersV, cleanClustersV);
    this->SelectCleanClusters(availableClustersW, cleanClustersW);

    // Calculate sliding fit results for clean clusters
    TwoDSlidingFitResultMap slidingFitResultMap;
    this->BuildSlidingFitResultMap(cleanClustersU, slidingFitResultMap);
    this->BuildSlidingFitResultMap(cleanClustersV, slidingFitResultMap);
    this->BuildSlidingFitResultMap(cleanClustersW, slidingFitResultMap);

    // Match clusters between pairs of views (using start/end information)
    ClusterAssociationMap matchedClustersUV, matchedClustersVW, matchedClustersWU;
    this->MatchViews(cleanClustersU, cleanClustersV, slidingFitResultMap, matchedClustersUV);
    this->MatchViews(cleanClustersV, cleanClustersW, slidingFitResultMap, matchedClustersVW);
    this->MatchViews(cleanClustersW, cleanClustersU, slidingFitResultMap, matchedClustersWU);

    // Create candidate particles using one, two and three primary views
    ParticleList candidateParticles;

    this->MatchThreeViews(cleanClustersU, cleanClustersV, cleanClustersW,
        matchedClustersUV, matchedClustersVW, matchedClustersWU, candidateParticles);
    
    this->MatchTwoViews(cleanClustersU, cleanClustersV, cleanClustersW,
        matchedClustersUV, matchedClustersVW, matchedClustersWU, candidateParticles);
   
    this->MatchOneView(cleanClustersU, cleanClustersV, cleanClustersW,
        matchedClustersUV, matchedClustersVW, matchedClustersWU, candidateParticles);
 
    // Build particle flow objects from candidate particles
    this->BuildParticles(candidateParticles);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackRecoveryAlgorithm::GetAvailableClusters(const std::string &inputClusterListName, ClusterVector &clusterVector) const
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,
        inputClusterListName, pClusterList))

    if (!pClusterList || pClusterList->empty())
    { 
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CosmicRayBaseMatchingAlgorithm: unable to find cluster list " << inputClusterListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
    {
        const Cluster *const pCluster = *cIter;

        if (!pCluster->IsAvailable())
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRecoveryAlgorithm::SelectCleanClusters(const ClusterVector &inputVector, ClusterVector &outputVector) const
{
    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        // Remove clusters below a minimum length
        if (LArClusterHelper::GetLengthSquared(pCluster) < m_clusterMinLength * m_clusterMinLength)
            continue;

        // Remove clusters nearly parallel to Z or X
        CartesianVector minCoordinate(0.f, 0.f, 0.f);
        CartesianVector maxCoordinate(0.f, 0.f, 0.f);
        LArClusterHelper::GetClusterBoundingBox(pCluster, minCoordinate, maxCoordinate);

        const CartesianVector deltaCoordinate(maxCoordinate - minCoordinate);
        if (std::fabs(deltaCoordinate.GetZ()) < m_clusterMinSpanZ || std::fabs(deltaCoordinate.GetX()) < m_clusterMinOverlapX)
            continue;

        outputVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRecoveryAlgorithm::BuildSlidingFitResultMap(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap) const
{
    const unsigned int m_halfWindowLayers(25);
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        if (slidingFitResultMap.end() == slidingFitResultMap.find(*iter))
        {
            try
            {
                const TwoDSlidingFitResult slidingFitResult(*iter, m_halfWindowLayers, slidingFitPitch);

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

void CosmicRayTrackRecoveryAlgorithm::MatchViews(const ClusterVector &clusterVector1, const ClusterVector &clusterVector2, 
    const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterAssociationMap &clusterAssociationMap) const
{
    for (ClusterVector::const_iterator iter1 = clusterVector1.begin(), iterEnd1 = clusterVector1.end(); iter1 != iterEnd1; ++iter1)
        this->MatchClusters(*iter1, clusterVector2, slidingFitResultMap, clusterAssociationMap);

    for (ClusterVector::const_iterator iter2 = clusterVector2.begin(), iterEnd2 = clusterVector2.end(); iter2 != iterEnd2; ++iter2)
        this->MatchClusters(*iter2, clusterVector1, slidingFitResultMap, clusterAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRecoveryAlgorithm::MatchClusters(const Cluster* const pSeedCluster, const ClusterVector &targetClusters, 
    const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterAssociationMap &clusterAssociationMap) const
{
    // Match seed cluster to target clusters according to alignment in X position of start/end positions
    // Two possible matches: (a) one-to-one associations where both the track start and end positions match up
    //                       (b) one-to-two associations where the track is split into two clusters in one view
    // Require overlap in X (according to clusterMinOverlapX) and alignment in X (according to clusterMaxDeltaX)

    TwoDSlidingFitResultMap::const_iterator fsIter = slidingFitResultMap.find(pSeedCluster);

    if (slidingFitResultMap.end() == fsIter)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const TwoDSlidingFitResult &slidingFitResult1(fsIter->second);
    const CartesianVector &innerVertex1(slidingFitResult1.GetGlobalMinLayerPosition());
    const CartesianVector &outerVertex1(slidingFitResult1.GetGlobalMaxLayerPosition());
    const float xSpan1(std::fabs(outerVertex1.GetX() - innerVertex1.GetX()));

    const Cluster *pBestClusterInner(NULL);
    const Cluster *pBestClusterOuter(NULL);
    const Cluster *pBestCluster(NULL);

    float bestDisplacementInner(m_clusterMaxDeltaX);
    float bestDisplacementOuter(m_clusterMaxDeltaX);
    float bestDisplacement(2.f * m_clusterMaxDeltaX);

    for (ClusterVector::const_iterator tIter = targetClusters.begin(), tIterEnd = targetClusters.end(); tIter != tIterEnd; ++tIter)
    {
        const Cluster *const pTargetCluster = *tIter;

        if (LArClusterHelper::GetClusterHitType(pSeedCluster) == LArClusterHelper::GetClusterHitType(pTargetCluster))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        TwoDSlidingFitResultMap::const_iterator ftIter = slidingFitResultMap.find(*tIter);

        if (slidingFitResultMap.end() == ftIter)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const TwoDSlidingFitResult &slidingFitResult2(ftIter->second);
        const CartesianVector &innerVertex2(slidingFitResult2.GetGlobalMinLayerPosition());
        const CartesianVector &outerVertex2(slidingFitResult2.GetGlobalMaxLayerPosition());
        const float xSpan2(std::fabs(outerVertex2.GetX() - innerVertex2.GetX()));

        if (xSpan2 > 1.5f * xSpan1)
            continue;

        const float xMin1(std::min(innerVertex1.GetX(), outerVertex1.GetX()));
        const float xMax1(std::max(innerVertex1.GetX(), outerVertex1.GetX()));
        const float xMin2(std::min(innerVertex2.GetX(), outerVertex2.GetX()));
        const float xMax2(std::max(innerVertex2.GetX(), outerVertex2.GetX()));
        const float xOverlap(std::min(xMax1,xMax2) - std::max(xMin1,xMin2));

        if (xOverlap < m_clusterMinOverlapX)
            continue;
            
        const float dxMin(std::fabs(xMin2 - xMin1));
        const float dxMax(std::fabs(xMax2 - xMax1));
           
        if (dxMin < bestDisplacementInner)
        {
            pBestClusterInner = pTargetCluster;
            bestDisplacementInner = dxMin;
        }

        if (dxMax < bestDisplacementOuter)
        {
            pBestClusterOuter = pTargetCluster;
            bestDisplacementOuter = dxMax;
        }

        if (dxMin + dxMax < bestDisplacement)
        {
            pBestCluster = pTargetCluster;
            bestDisplacement = dxMin + dxMax;
        }
    }

    if (pBestCluster)
    {
        clusterAssociationMap[pSeedCluster].insert(pBestCluster);
        clusterAssociationMap[pBestCluster].insert(pSeedCluster);
    }
    else if (pBestClusterInner && pBestClusterOuter)
    {
        TwoDSlidingFitResultMap::const_iterator iterInner = slidingFitResultMap.find(pBestClusterInner);  
        TwoDSlidingFitResultMap::const_iterator iterOuter = slidingFitResultMap.find(pBestClusterOuter);

        if (slidingFitResultMap.end() == iterInner || slidingFitResultMap.end() == iterOuter)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const LArPointingCluster pointingClusterInner(iterInner->second);
        const LArPointingCluster pointingClusterOuter(iterOuter->second);

        LArPointingCluster::Vertex pointingVertexInner, pointingVertexOuter;
        LArPointingClusterHelper::GetClosestVertices(pointingClusterInner, pointingClusterOuter, pointingVertexInner, pointingVertexOuter);

        const LArPointingCluster::Vertex pointingEndInner(pointingVertexInner.IsInnerVertex() ? pointingClusterInner.GetOuterVertex() : pointingClusterInner.GetInnerVertex());
        const LArPointingCluster::Vertex pointingEndOuter(pointingVertexOuter.IsInnerVertex() ? pointingClusterOuter.GetOuterVertex() : pointingClusterOuter.GetInnerVertex());

        const float rSpan((pointingEndInner.GetPosition() - pointingEndOuter.GetPosition()).GetMagnitude());

        if (LArPointingClusterHelper::IsEmission(pointingVertexInner.GetPosition(), pointingVertexOuter, -1.f, 0.75f * rSpan, 5.f, 10.f) &&
            LArPointingClusterHelper::IsEmission(pointingVertexOuter.GetPosition(), pointingVertexInner, -1.f, 0.75f * rSpan, 5.f, 10.f))
        {
            clusterAssociationMap[pBestClusterInner].insert(pSeedCluster);
            clusterAssociationMap[pBestClusterOuter].insert(pSeedCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRecoveryAlgorithm::MatchThreeViews(const ClusterVector &clusterVectorU,  const ClusterVector &clusterVectorV,
    const ClusterVector &clusterVectorW, const ClusterAssociationMap &matchedClustersUV, const ClusterAssociationMap &matchedClustersVW, 
    const ClusterAssociationMap &matchedClustersWU, ParticleList &particleList) const
{  
    ClusterList vetoList;
    this->BuildVetoList(particleList, vetoList);

    ParticleList newParticleList;

    const ClusterVector &clusterVector1(clusterVectorU);
    const ClusterVector &clusterVector2(clusterVectorV);
    const ClusterVector &clusterVector3(clusterVectorW);
 
    const ClusterAssociationMap &matchedClusters12(matchedClustersUV);
    const ClusterAssociationMap &matchedClusters23(matchedClustersVW);
    const ClusterAssociationMap &matchedClusters31(matchedClustersWU);

    for (ClusterVector::const_iterator iter1 = clusterVector1.begin(), iterEnd1 = clusterVector1.end(); iter1 != iterEnd1; ++iter1)
    {
        const Cluster *const pCluster1 = *iter1;

        if (vetoList.count(pCluster1))
            continue;

        const ClusterAssociationMap::const_iterator iter311 = matchedClusters31.find(pCluster1);
        const ClusterList matchedClusters31_pCluster1(iter311 != matchedClusters31.end() ? iter311->second : ClusterList());

        const ClusterAssociationMap::const_iterator iter121 = matchedClusters12.find(pCluster1);
        const ClusterList matchedClusters12_pCluster1(iter121 != matchedClusters12.end() ? iter121->second : ClusterList());

        for (ClusterVector::const_iterator iter2 = clusterVector2.begin(), iterEndV = clusterVector2.end(); iter2 != iterEndV; ++iter2)
        { 
            const Cluster *const pCluster2 = *iter2;

            if (vetoList.count(pCluster2))
                continue;

            const ClusterAssociationMap::const_iterator iter122 = matchedClusters12.find(pCluster2);
            const ClusterList matchedClusters12_pCluster2(iter122 != matchedClusters12.end() ? iter122->second : ClusterList());

            const ClusterAssociationMap::const_iterator iter232 = matchedClusters23.find(pCluster2);
            const ClusterList matchedClusters23_pCluster2(iter232 != matchedClusters23.end() ? iter232->second : ClusterList());

            for (ClusterVector::const_iterator iter3 = clusterVector3.begin(), iterEnd3 = clusterVector3.end(); iter3 != iterEnd3; ++iter3)
            { 
                const Cluster *const pCluster3 = *iter3;

                if (vetoList.count(pCluster3))
                    continue;

                const ClusterAssociationMap::const_iterator iter233 = matchedClusters23.find(pCluster3);
                const ClusterList matchedClusters23_pCluster3(iter233 != matchedClusters23.end() ? iter233->second : ClusterList());

                const ClusterAssociationMap::const_iterator iter313 = matchedClusters31.find(pCluster3);
                const ClusterList matchedClusters31_pCluster3(iter313 != matchedClusters31.end() ? iter313->second : ClusterList());

                const bool match12((matchedClusters12_pCluster1.size() +  matchedClusters12_pCluster2.size() > 0) &&
                    ((matchedClusters12_pCluster1.size() == 1 && matchedClusters12_pCluster1.find(pCluster2) != matchedClusters12_pCluster1.end()) || 
                     (matchedClusters12_pCluster1.size() == 0)) &&
                    ((matchedClusters12_pCluster2.size() == 1 && matchedClusters12_pCluster2.find(pCluster1) != matchedClusters12_pCluster2.end()) || 
                     (matchedClusters12_pCluster2.size() == 0)));

                const bool match23((matchedClusters23_pCluster2.size() +  matchedClusters23_pCluster3.size() > 0) &&
                    ((matchedClusters23_pCluster2.size() == 1 && matchedClusters23_pCluster2.find(pCluster3) != matchedClusters23_pCluster2.end()) || 
                     (matchedClusters23_pCluster2.size() == 0)) &&
                    ((matchedClusters23_pCluster3.size() == 1 && matchedClusters23_pCluster3.find(pCluster2) != matchedClusters23_pCluster3.end()) || 
                    (matchedClusters23_pCluster3.size() == 0)));

                const bool match31((matchedClusters31_pCluster3.size() +  matchedClusters31_pCluster1.size() > 0) &&
                    ((matchedClusters31_pCluster3.size() == 1 && matchedClusters31_pCluster3.find(pCluster1) != matchedClusters31_pCluster3.end()) || 
                    (matchedClusters31_pCluster3.size() == 0)) &&
                    ((matchedClusters31_pCluster1.size() == 1 && matchedClusters31_pCluster1.find(pCluster3) != matchedClusters31_pCluster1.end()) || 
                    (matchedClusters31_pCluster1.size() == 0)));

                if (match12 && match23 && match31)
                {
                    Particle newParticle;
                    newParticle.m_clusterList.insert(pCluster1);
                    newParticle.m_clusterList.insert(pCluster2);
                    newParticle.m_clusterList.insert(pCluster3);
                    newParticleList.push_back(newParticle);
                }
            }
        }
    }

    return this->RemoveAmbiguities(newParticleList, particleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRecoveryAlgorithm::MatchTwoViews(const ClusterVector &clusterVectorU,  const ClusterVector &clusterVectorV,
    const ClusterVector &clusterVectorW, const ClusterAssociationMap &matchedClustersUV, const ClusterAssociationMap &matchedClustersVW, 
    const ClusterAssociationMap &matchedClustersWU, ParticleList &particleList) const
{
    ClusterList vetoList;
    this->BuildVetoList(particleList, vetoList);

    ParticleList newParticleList;

    for (unsigned int iView = 0; iView < 3; ++iView)
    {
        const ClusterVector &clusterVector1((0 == iView) ? clusterVectorU : (1 == iView) ? clusterVectorV : clusterVectorW);
        const ClusterVector &clusterVector2((0 == iView) ? clusterVectorV : (1 == iView) ? clusterVectorW : clusterVectorU);    
        const ClusterVector &clusterVector3((0 == iView) ? clusterVectorW : (1 == iView) ? clusterVectorU : clusterVectorV);
 
        const ClusterAssociationMap &matchedClusters12(((0 == iView) ? matchedClustersUV : (1 == iView) ? matchedClustersVW : matchedClustersWU));
        const ClusterAssociationMap &matchedClusters23(((0 == iView) ? matchedClustersVW : (1 == iView) ? matchedClustersWU : matchedClustersUV));
        const ClusterAssociationMap &matchedClusters31(((0 == iView) ? matchedClustersWU : (1 == iView) ? matchedClustersUV : matchedClustersVW));

        for (ClusterVector::const_iterator iter1 = clusterVector1.begin(), iterEnd1 = clusterVector1.end(); iter1 != iterEnd1; ++iter1)
        {
            const Cluster *const pCluster1 = *iter1;

            if (vetoList.count(pCluster1))
                continue;

            const ClusterAssociationMap::const_iterator iter311 = matchedClusters31.find(pCluster1);
            const ClusterList matchedClusters31_pCluster1(iter311 != matchedClusters31.end() ? iter311->second : ClusterList());

            const ClusterAssociationMap::const_iterator iter121 = matchedClusters12.find(pCluster1);
            const ClusterList matchedClusters12_pCluster1(iter121 != matchedClusters12.end() ? iter121->second : ClusterList());

            for (ClusterVector::const_iterator iter2 = clusterVector2.begin(), iterEnd2 = clusterVector2.end(); iter2 != iterEnd2; ++iter2)
            { 
                const Cluster *const pCluster2 = *iter2;

                if (vetoList.count(pCluster2))
                    continue;

                const ClusterAssociationMap::const_iterator iter122 = matchedClusters12.find(pCluster2);
                const ClusterList matchedClusters12_pCluster2(iter122 != matchedClusters12.end() ? iter122->second : ClusterList());

                const ClusterAssociationMap::const_iterator iter232 = matchedClusters23.find(pCluster2);
                const ClusterList matchedClusters23_pCluster2(iter232 != matchedClusters23.end() ? iter232->second : ClusterList());

                const bool match12((matchedClusters12_pCluster1.size() == 1 && matchedClusters12_pCluster1.find(pCluster2) != matchedClusters12_pCluster1.end()) &&
                    (matchedClusters12_pCluster2.size() == 1 && matchedClusters12_pCluster2.find(pCluster1) != matchedClusters12_pCluster2.end()) &&
                    (matchedClusters23_pCluster2.size() == 0 && matchedClusters31_pCluster1.size() == 0));

                if (!match12)
                continue;

                Particle newParticle;
                newParticle.m_clusterList.insert(pCluster1);
                newParticle.m_clusterList.insert(pCluster2);
                    
                for (ClusterVector::const_iterator iter3 = clusterVector3.begin(), iterEnd3 = clusterVector3.end(); iter3 != iterEnd3; ++iter3)
                { 
                    const Cluster *const pCluster3 = *iter3;
  
                    if (vetoList.count(pCluster3))
                        continue;
 
                    const ClusterAssociationMap::const_iterator iter233 = matchedClusters23.find(pCluster3);
                    const ClusterList matchedClusters23_pCluster3(iter233 != matchedClusters23.end() ? iter233->second : ClusterList());

                    const ClusterAssociationMap::const_iterator iter313 = matchedClusters31.find(pCluster3);
                    const ClusterList matchedClusters31_pCluster3(iter313 != matchedClusters31.end() ? iter313->second : ClusterList());

                    const bool match3((matchedClusters31_pCluster3.size() +  matchedClusters23_pCluster3.size() > 0) &&
                        ((matchedClusters31_pCluster3.size() == 1 && matchedClusters31_pCluster3.find(pCluster1) != matchedClusters31_pCluster3.end()) || 
                        (matchedClusters31_pCluster3.size() == 0)) &&
                        ((matchedClusters23_pCluster3.size() == 1 && matchedClusters23_pCluster3.find(pCluster2) != matchedClusters23_pCluster3.end()) || 
                        (matchedClusters23_pCluster3.size() == 0)));

                    if (match3)
                        newParticle.m_clusterList.insert(pCluster3);
                }

                newParticleList.push_back(newParticle);
            }
        }
    }   

    return this->RemoveAmbiguities(newParticleList, particleList);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRecoveryAlgorithm::MatchOneView(const ClusterVector &clusterVectorU,  const ClusterVector &clusterVectorV,
    const ClusterVector &clusterVectorW, const ClusterAssociationMap &matchedClustersUV, const ClusterAssociationMap &matchedClustersVW, 
    const ClusterAssociationMap &matchedClustersWU, ParticleList &particleList) const
{
    ClusterList vetoList;
    this->BuildVetoList(particleList, vetoList);

    ParticleList newParticleList;

    for (unsigned int iView = 0; iView < 3; ++iView)
    {
        const ClusterVector &clusterVector1((0 == iView) ? clusterVectorU : (1 == iView) ? clusterVectorV : clusterVectorW);
        const ClusterVector &clusterVector2((0 == iView) ? clusterVectorV : (1 == iView) ? clusterVectorW : clusterVectorU);
        const ClusterVector &clusterVector3((0 == iView) ? clusterVectorW : (1 == iView) ? clusterVectorU : clusterVectorV);
 
        const ClusterAssociationMap &matchedClusters12(((0 == iView) ? matchedClustersUV : (1 == iView) ? matchedClustersVW : matchedClustersWU));
        const ClusterAssociationMap &matchedClusters23(((0 == iView) ? matchedClustersVW : (1 == iView) ? matchedClustersWU : matchedClustersUV));
        const ClusterAssociationMap &matchedClusters31(((0 == iView) ? matchedClustersWU : (1 == iView) ? matchedClustersUV : matchedClustersVW));

        for (ClusterVector::const_iterator iter1 = clusterVector1.begin(), iterEnd1 = clusterVector1.end(); iter1 != iterEnd1; ++iter1)
        {
            const Cluster *const pCluster1 = *iter1;

            if (vetoList.count(pCluster1))
                continue;

            const ClusterAssociationMap::const_iterator iter311 = matchedClusters31.find(pCluster1);
            const ClusterList matchedClusters31_pCluster1(iter311 != matchedClusters31.end() ? iter311->second : ClusterList());

            const ClusterAssociationMap::const_iterator iter121 = matchedClusters12.find(pCluster1);
            const ClusterList matchedClusters12_pCluster1(iter121 != matchedClusters12.end() ? iter121->second : ClusterList());

            if (matchedClusters12_pCluster1.size() + matchedClusters31_pCluster1.size() > 0)
                continue;

            Particle newParticle;
            newParticle.m_clusterList.insert(pCluster1);

            for (ClusterVector::const_iterator iter2 = clusterVector2.begin(), iterEnd2 = clusterVector2.end(); iter2 != iterEnd2; ++iter2)
            {
                const Cluster *const pCluster2 = *iter2;

                if (vetoList.count(pCluster2))
                    continue;

                const ClusterAssociationMap::const_iterator iter122 = matchedClusters12.find(pCluster2);
                const ClusterList matchedClusters12_pCluster2(iter122 != matchedClusters12.end() ? iter122->second : ClusterList());

                const ClusterAssociationMap::const_iterator iter232 = matchedClusters23.find(pCluster2);
                const ClusterList matchedClusters23_pCluster2(iter232 != matchedClusters23.end() ? iter232->second : ClusterList());

                if (matchedClusters12_pCluster2.size() == 1 && matchedClusters12_pCluster2.find(pCluster1) != matchedClusters12_pCluster2.end() &&
                    matchedClusters23_pCluster2.size() == 0)
                    newParticle.m_clusterList.insert(pCluster2);
            }

            for (ClusterVector::const_iterator iter3 = clusterVector3.begin(), iterEnd3 = clusterVector3.end(); iter3 != iterEnd3; ++iter3)
            {
                const Cluster *const pCluster3 = *iter3;
 
                if (vetoList.count(pCluster3))
                    continue;

                const ClusterAssociationMap::const_iterator iter233 = matchedClusters23.find(pCluster3);
                const ClusterList matchedClusters23_pCluster3(iter233 != matchedClusters23.end() ? iter233->second : ClusterList());

                const ClusterAssociationMap::const_iterator iter313 = matchedClusters31.find(pCluster3);
                const ClusterList matchedClusters31_pCluster3(iter313 != matchedClusters31.end() ? iter313->second : ClusterList());

                if (matchedClusters31_pCluster3.size() == 1 && matchedClusters31_pCluster3.find(pCluster1) != matchedClusters31_pCluster3.end() &&
                    matchedClusters23_pCluster3.size() == 0)
                    newParticle.m_clusterList.insert(pCluster3);
            }

            if (newParticle.m_clusterList.size() > 1)
                newParticleList.push_back(newParticle);
        }
    }

    return this->RemoveAmbiguities(newParticleList, particleList);
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRecoveryAlgorithm::BuildVetoList(const ParticleList &particleList, ClusterList &vetoList) const
{
    for (ParticleList::const_iterator pIter = particleList.begin(), pIterEnd = particleList.end(); pIter != pIterEnd; ++pIter)
    {
        const Particle &particle = *pIter;
        vetoList.insert(particle.m_clusterList.begin(), particle.m_clusterList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void CosmicRayTrackRecoveryAlgorithm::RemoveAmbiguities(const ParticleList &inputParticleList, ParticleList &outputParticleList) const
{
    for (ParticleList::const_iterator pIter1 = inputParticleList.begin(), pIterEnd1 = inputParticleList.end(); pIter1 != pIterEnd1; ++pIter1)
    {
        const Particle &particle1 = *pIter1;
        const ClusterList &clusterList1 = particle1.m_clusterList;

        try
        {
            bool isUnique(true);

            for (ParticleList::const_iterator pIter2 = pIter1, pIterEnd2 = inputParticleList.end(); pIter2 != pIterEnd2; ++pIter2)
            {
                const Particle &particle2 = *pIter2;
                const ClusterList &clusterList2 = particle2.m_clusterList;

                ClusterList duplicateList;

                for (ClusterList::const_iterator cIter1 = clusterList1.begin(), cIterEnd1 = clusterList1.end(); cIter1 != cIterEnd1; ++cIter1)
                {
                    const Cluster *pCluster = *cIter1;

                    if (clusterList2.count(pCluster))
                        duplicateList.insert(pCluster);
                }

                if (duplicateList.size() > 0 && clusterList1.size() != clusterList2.size())
                {
                    isUnique = false;
                    break;
                }
            }

            if (!isUnique)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            outputParticleList.push_back(particle1);
        }
        catch (StatusCodeException &)
        {
            std::cout << " Warning in CosmicRayTrackRecoveryAlgorithm: found duplicate particles in candidate list " << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackRecoveryAlgorithm::BuildParticles(const ParticleList &particleList)
{
    if (particleList.empty())
        return;

    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (ParticleList::const_iterator iter = particleList.begin(), iterEnd = particleList.end(); iter != iterEnd; ++iter)
    {
        const Particle &particle = *iter;

        ClusterList clusterList;
        this->MergeClusters(particle.m_clusterList, clusterList);

        if (clusterList.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

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

void CosmicRayTrackRecoveryAlgorithm::MergeClusters(const ClusterList &inputClusterList, ClusterList &outputClusterList) const
{
    ClusterList clusterListU, clusterListV, clusterListW;
    LArClusterHelper::GetClustersByHitType(inputClusterList, TPC_VIEW_U, clusterListU);
    LArClusterHelper::GetClustersByHitType(inputClusterList, TPC_VIEW_V, clusterListV);
    LArClusterHelper::GetClustersByHitType(inputClusterList, TPC_VIEW_W, clusterListW);

    for (unsigned int iView = 0; iView < 3; ++iView)
    {
        const ClusterList clusterList((0 == iView) ? clusterListU : (1 == iView) ?  clusterListV : clusterListW);
        const std::string inputClusterListName((0 == iView) ? m_inputClusterListNameU : (1 == iView) ? m_inputClusterListNameV : m_inputClusterListNameW); 

        if (clusterList.empty())
            continue;

        const Cluster *const pSeedCluster = *(clusterList.begin());

        if (!pSeedCluster->IsAvailable())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
        {
            const Cluster *const pAssociatedCluster = *iter;

            if (pAssociatedCluster == pSeedCluster)
                continue;

            if (!pAssociatedCluster->IsAvailable())
                throw StatusCodeException(STATUS_CODE_FAILURE);

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pSeedCluster, pAssociatedCluster,
                inputClusterListName, inputClusterListName));
        }

        outputClusterList.insert(pSeedCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackRecoveryAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinLength", m_clusterMinLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinSpanZ", m_clusterMinSpanZ));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinOverlapX", m_clusterMinOverlapX));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMaxDeltaX", m_clusterMaxDeltaX));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
