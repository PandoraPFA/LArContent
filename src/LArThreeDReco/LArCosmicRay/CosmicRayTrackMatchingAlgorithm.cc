/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the cosmic ray splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArThreeDReco/LArCosmicRay/CosmicRayTrackMatchingAlgorithm.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

using namespace pandora;

namespace lar
{

StatusCode CosmicRayTrackMatchingAlgorithm::Run()
{
    std::cout << " --- CosmicRayTrackMatchingAlgorithm::Run() --- " << std::endl;

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

    // Build a map of sliding linear fit results
    TwoDSlidingFitResultMap slidingFitCache;
    this->AddToSlidingFitCache(cleanClustersU, slidingFitCache);
    this->AddToSlidingFitCache(cleanClustersV, slidingFitCache);
    this->AddToSlidingFitCache(cleanClustersW, slidingFitCache);

    // Build associations between pairs of views
    ClusterAssociationMap matchedClusterUV, matchedClusterVW, matchedClusterWU;
    this->MatchTracks(slidingFitCache, cleanClustersU, cleanClustersV, matchedClusterUV);
    this->MatchTracks(slidingFitCache, cleanClustersV, cleanClustersW, matchedClusterVW);
    this->MatchTracks(slidingFitCache, cleanClustersW, cleanClustersU, matchedClusterWU);


    ParticleList particleList;
    this->MatchThreeViews(matchedClusterUV, matchedClusterVW, matchedClusterWU, particleList);
    this->MatchTwoViews(matchedClusterUV, matchedClusterVW, matchedClusterWU, particleList);
    this->BuildParticles(particleList);


 
// --- END EVENT DISPLAY ---
ClusterList tempListU, tempListV, tempListW;
for (ParticleList::const_iterator iter = particleList.begin(), iterEnd = particleList.end(); iter != iterEnd; ++iter)
{
const Particle &particle = *iter;
Cluster* pClusterU = const_cast<Cluster*>(particle.m_pClusterU);
Cluster* pClusterV = const_cast<Cluster*>(particle.m_pClusterV);
Cluster* pClusterW = const_cast<Cluster*>(particle.m_pClusterW);
std::cout << " pU=" << pClusterU << " pV=" << pClusterV << " pW=" << pClusterW << std::endl;
if(pClusterU) tempListU.insert(pClusterU);
if(pClusterV) tempListV.insert(pClusterV);
if(pClusterW) tempListW.insert(pClusterW);
}
PandoraMonitoringApi::SetEveDisplayParameters(false, DETECTOR_VIEW_XZ);
PandoraMonitoringApi::VisualizeClusters(&tempListU, "MatchedClusterU", RED);
PandoraMonitoringApi::VisualizeClusters(&tempListV, "MatchedClusterV", BLUE);
PandoraMonitoringApi::VisualizeClusters(&tempListW, "MatchedClusterW", GREEN);
PandoraMonitoringApi::ViewEvent();
// --- END EVENT DISPLAY ---



    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackMatchingAlgorithm::GetAvailableClusters(const std::string inputClusterListNames, ClusterVector &clusterVector) const
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, 
        inputClusterListNames, pClusterList))

    if (NULL == pClusterList)
    {
        std::cout << "CosmicRayTrackMatchingAlgorithm: could not find cluster list " << inputClusterListNames << std::endl;  
        return STATUS_CODE_SUCCESS;
    }

    for (ClusterList::const_iterator cIter = pClusterList->begin(), cIterEnd = pClusterList->end(); cIter != cIterEnd; ++cIter)
    {
        Cluster *pCluster = *cIter;
        if (!pCluster->IsAvailable())
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::SelectCleanClusters(const ClusterVector &inputVector, ClusterVector &outputVector) const
{
    for (ClusterVector::const_iterator iter = inputVector.begin(), iterEnd = inputVector.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_clusterMinLength * m_clusterMinLength)
            continue;

        outputVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::AddToSlidingFitCache(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitCache) const
{
    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        if (slidingFitCache.end() == slidingFitCache.find(*iter))
        {
            TwoDSlidingFitResult slidingFitResult;
            LArClusterHelper::LArTwoDSlidingFit(*iter, m_halfWindowLayers, slidingFitResult);

            if (!slidingFitCache.insert(TwoDSlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
                throw StatusCodeException(STATUS_CODE_FAILURE);
        }
    }
}
 
//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::MatchTracks(const TwoDSlidingFitResultMap &slidingFitCache, const ClusterVector &clusterVector1,
    const ClusterVector &clusterVector2, ClusterAssociationMap &matchedClusters12) const
{
    // Check that there are input clusters from both views
    if (clusterVector1.empty() || clusterVector2.empty())
        return;

    const HitType hitType1(LArThreeDHelper::GetClusterHitType(*clusterVector1.begin()));
    const HitType hitType2(LArThreeDHelper::GetClusterHitType(*clusterVector2.begin()));

    if (hitType1 == hitType2)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    for (ClusterVector::const_iterator iter1 = clusterVector1.begin(), iterEnd1 = clusterVector1.end(); iter1 != iterEnd1; ++iter1)
    {
        Cluster* pCluster1 = *iter1;

        TwoDSlidingFitResultMap::const_iterator sIter1 = slidingFitCache.find(pCluster1);
        if (slidingFitCache.end() == sIter1)
            continue;

        const TwoDSlidingFitResult &slidingFitResult1(sIter1->second);

        for (ClusterVector::const_iterator iter2 = clusterVector2.begin(), iterEnd2 = clusterVector2.end(); iter2 != iterEnd2; ++iter2)
        {
            Cluster* pCluster2 = *iter2;

            TwoDSlidingFitResultMap::const_iterator sIter2 = slidingFitCache.find(pCluster2);
            if (slidingFitCache.end() == sIter2)
                continue;

            const TwoDSlidingFitResult &slidingFitResult2(sIter2->second);

            if (this->MatchTracks(slidingFitResult1, slidingFitResult2))
                matchedClusters12[pCluster1].insert(pCluster2);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CosmicRayTrackMatchingAlgorithm::MatchTracks(const TwoDSlidingFitResult &slidingFit1, const TwoDSlidingFitResult &slidingFit2) const
{
    // Require a good X overlap between clusters in the first and second views
    const Cluster* pCluster1(slidingFit1.GetCluster());
    const Cluster* pCluster2(slidingFit2.GetCluster());

    float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
    LArClusterHelper::GetClusterSpanX(pCluster1, xMin1, xMax1);
    LArClusterHelper::GetClusterSpanX(pCluster2, xMin2, xMax2);

    const float xOverlap(std::min(xMax1,xMax2) - std::max(xMin1,xMin2));
    const float xSpan(std::max(xMax1,xMax2) - std::min(xMin1,xMin2));

    if (xOverlap < m_minXOverlap || xOverlap/xSpan < m_minXOverlapFraction)
        return false;

// --- BEGIN EVENT DISPLAY ---
ClusterList tempList1, tempList2;
tempList1.insert((Cluster*)pCluster1);
tempList2.insert((Cluster*)pCluster2);
PandoraMonitoringApi::SetEveDisplayParameters(false, DETECTOR_VIEW_XZ);
PandoraMonitoringApi::VisualizeClusters(&tempList1, "Cluster1", RED);
PandoraMonitoringApi::VisualizeClusters(&tempList2, "Cluster2", BLUE);
PandoraMonitoringApi::ViewEvent();
// --- END EVENT DISPLAY ---

    return true;
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::MatchThreeViews(const ClusterAssociationMap &matchedClusters12, 
    const ClusterAssociationMap &matchedClusters23,  const ClusterAssociationMap &matchedClusters31, ParticleList &matchedParticles) const
{
    if (matchedClusters12.empty() || matchedClusters23.empty() || matchedClusters31.empty())
        return;

    ParticleList candidateParticles;

    for (ClusterAssociationMap::const_iterator iter12 = matchedClusters12.begin(), iterEnd12 = matchedClusters12.end(); iter12 != iterEnd12; 
        ++iter12)
    {
        const Cluster* pCluster1 = iter12->first;
        const ClusterList &clusterList2 = iter12->second;

        for(ClusterList::const_iterator iter2 = clusterList2.begin(), iterEnd2 = clusterList2.end(); iter2 != iterEnd2; ++iter2)
        {
            const Cluster* pCluster2 = *iter2;  

            ClusterAssociationMap::const_iterator iter23 = matchedClusters23.find(pCluster2);
            if (matchedClusters23.end() == iter23)
                continue; 
 
            const ClusterList &clusterList3 = iter23->second;

            for(ClusterList::const_iterator iter3 = clusterList3.begin(), iterEnd3 = clusterList3.end(); iter3 != iterEnd3; ++iter3)
            {
                const Cluster* pCluster3 = *iter3;  
 
                ClusterAssociationMap::const_iterator iter31 = matchedClusters31.find(pCluster3);
                if (matchedClusters31.end() == iter31)
                    continue; 

                const ClusterList &clusterList1 = iter31->second; 
                Cluster* pCluster1check = const_cast<Cluster*>(pCluster1);
                ClusterList::const_iterator iter1 = clusterList1.find(pCluster1check);

                if (clusterList1.end() == iter1)
                    continue;

                const HitType hitType1(LArThreeDHelper::GetClusterHitType(pCluster1));
                const HitType hitType2(LArThreeDHelper::GetClusterHitType(pCluster2));
                const HitType hitType3(LArThreeDHelper::GetClusterHitType(pCluster3));

                const Cluster* pClusterU((TPC_VIEW_U == hitType1) ? pCluster1 : (TPC_VIEW_U == hitType2) ? pCluster2 : 
                                         (TPC_VIEW_U == hitType3) ? pCluster3 : NULL);
                const Cluster* pClusterV((TPC_VIEW_V == hitType1) ? pCluster1 : (TPC_VIEW_V == hitType2) ? pCluster2 : 
                                         (TPC_VIEW_V == hitType3) ? pCluster3 : NULL);
                const Cluster* pClusterW((TPC_VIEW_W == hitType1) ? pCluster1 : (TPC_VIEW_W == hitType2) ? pCluster2 : 
                                         (TPC_VIEW_W == hitType3) ? pCluster3 : NULL);

                candidateParticles.push_back(Particle(pClusterU, pClusterV, pClusterW));
            }
        }
    }

    return this->ResolveAmbiguities(candidateParticles, matchedParticles);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::MatchTwoViews(const ClusterAssociationMap &matchedClusters12, 
    const ClusterAssociationMap &matchedClusters23, const ClusterAssociationMap &matchedClusters31, ParticleList &matchedParticles) const
{
    ParticleList candidateParticles;
    this->MatchTwoViews(matchedClusters12, candidateParticles);
    this->MatchTwoViews(matchedClusters23, candidateParticles); 
    this->MatchTwoViews(matchedClusters31, candidateParticles);

    return this->ResolveAmbiguities(candidateParticles, matchedParticles);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::MatchTwoViews(const ClusterAssociationMap &matchedClusters12, ParticleList &matchedParticles) const
{ 
    if (matchedClusters12.empty())
        return;
    
    for (ClusterAssociationMap::const_iterator iter12 = matchedClusters12.begin(), iterEnd12 = matchedClusters12.end(); iter12 != iterEnd12; 
        ++iter12)
    {
        const Cluster* pCluster1 = iter12->first;
        const ClusterList &clusterList2 = iter12->second;      

        for (ClusterList::const_iterator iter2 = clusterList2.begin(), iterEnd2 = clusterList2.end() ; iter2 != iterEnd2; ++iter2)
        {
            const Cluster* pCluster2 = *iter2;
            
            const HitType hitType1(LArThreeDHelper::GetClusterHitType(pCluster1));
            const HitType hitType2(LArThreeDHelper::GetClusterHitType(pCluster2));
                
            const Cluster* pClusterU((TPC_VIEW_U == hitType1) ? pCluster1 : (TPC_VIEW_U == hitType2) ? pCluster2 : NULL);
            const Cluster* pClusterV((TPC_VIEW_V == hitType1) ? pCluster1 : (TPC_VIEW_V == hitType2) ? pCluster2 : NULL);
            const Cluster* pClusterW((TPC_VIEW_W == hitType1) ? pCluster1 : (TPC_VIEW_W == hitType2) ? pCluster2 : NULL);

            matchedParticles.push_back(Particle(pClusterU, pClusterV, pClusterW));
        }
    }
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::ResolveAmbiguities(const ParticleList &candidateParticles, ParticleList &matchedParticles) const
{
    for (ParticleList::const_iterator iter1 = candidateParticles.begin(), iterEnd1 = candidateParticles.end(); iter1 != iterEnd1; ++iter1)
    {
        const Particle &particle1 = *iter1;

        bool isGoodMatch(true);

        for (ParticleList::const_iterator iter2 = candidateParticles.begin(), iterEnd2 = candidateParticles.end(); iter2 != iterEnd2; ++iter2)
        {
            const Particle &particle2 = *iter2;
 
            const bool commonU(particle1.m_pClusterU == particle2.m_pClusterU);
            const bool commonV(particle1.m_pClusterV == particle2.m_pClusterV);
            const bool commonW(particle1.m_pClusterW == particle2.m_pClusterW);

            const bool ambiguousU(commonU && NULL != particle1.m_pClusterU);
            const bool ambiguousV(commonV && NULL != particle1.m_pClusterV);
            const bool ambiguousW(commonW && NULL != particle1.m_pClusterW);
 
            if (commonU && commonV && commonW)
                continue;

            if (ambiguousU || ambiguousV || ambiguousW)
            {
                isGoodMatch = false;
                break;   
            }
        }

        if (isGoodMatch)
            matchedParticles.push_back(particle1);
    }
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicRayTrackMatchingAlgorithm::BuildParticles(const ParticleList &particleList)
{
    if (particleList.empty())
        return;

    const PfoList *pPfoList = NULL; std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    for (ParticleList::const_iterator iter = particleList.begin(), iterEnd = particleList.end(); iter != iterEnd; ++iter)
    {
        const Particle &particle = *iter;
        
        ClusterList clusterList;
        Cluster *pClusterU = const_cast<Cluster*>(particle.m_pClusterU);
        Cluster *pClusterV = const_cast<Cluster*>(particle.m_pClusterV);
        Cluster *pClusterW = const_cast<Cluster*>(particle.m_pClusterW);

        const bool isAvailableU((NULL != pClusterU) ? pClusterU->IsAvailable() : true);
        const bool isAvailableV((NULL != pClusterV) ? pClusterV->IsAvailable() : true);
        const bool isAvailableW((NULL != pClusterW) ? pClusterW->IsAvailable() : true);

        if(!(isAvailableU && isAvailableV && isAvailableW))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        if (pClusterU) clusterList.insert(pClusterU);
        if (pClusterV) clusterList.insert(pClusterV);
        if (pClusterW) clusterList.insert(pClusterW);
       
        // TODO - correct these placeholder parameters
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = 22;
        pfoParameters.m_charge = 0;
        pfoParameters.m_mass = 0.f;
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = CartesianVector(0., 0., 0.);
        pfoParameters.m_clusterList = clusterList;

        ParticleFlowObject *pPfo(NULL);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
    }

    if (!pPfoList->empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------
  
CosmicRayTrackMatchingAlgorithm::Particle::Particle(const Cluster *pClusterU, const Cluster *pClusterV, const Cluster *pClusterW) :
    m_pClusterU(pClusterU),
    m_pClusterV(pClusterV),
    m_pClusterW(pClusterW)
{
    if (NULL == m_pClusterU && NULL == m_pClusterV && NULL == m_pClusterW)
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const HitType hitTypeU(NULL == m_pClusterU ? TPC_VIEW_U : LArThreeDHelper::GetClusterHitType(m_pClusterU));   
    const HitType hitTypeV(NULL == m_pClusterV ? TPC_VIEW_V : LArThreeDHelper::GetClusterHitType(m_pClusterV));
    const HitType hitTypeW(NULL == m_pClusterW ? TPC_VIEW_W : LArThreeDHelper::GetClusterHitType(m_pClusterW));      

    if (!(TPC_VIEW_U == hitTypeU && TPC_VIEW_V == hitTypeV && TPC_VIEW_W == hitTypeW))
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayTrackMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));

 

    m_clusterMinLength = 10.f; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinLength", m_clusterMinLength));

    m_halfWindowLayers = 15;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitHalfWindow", m_halfWindowLayers));

    m_minXOverlap = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlap", m_minXOverlap));

    m_minXOverlapFraction = 0.8f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

  

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
