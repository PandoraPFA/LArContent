/**
 *  @file   LArContent/src/LArTwoDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.cc
 * 
 *  @brief  Implementation of the cosmic ray shower matching algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "LArHelpers/LArMCParticleHelper.h"
#include "LArHelpers/LArThreeDHelper.h"
#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArTwoDReco/LArCosmicRay/CosmicRayShowerMatchingAlgorithm.h"

using namespace pandora;

namespace lar
{
    
    StatusCode CosmicRayShowerMatchingAlgorithm::Run()
    {

        ctopmap_t clusterToPfoMap;
        ptocmultimap_t pfoAssociatedClusterMap;

        // check the input PFO lists exists
        const PfoList *pPfoList = NULL;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

        const PfoList *pCosmicPfoList = NULL;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCosmicRayPfoListName, pCosmicPfoList));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CosmicRay3DShowerMatching(clusterToPfoMap,pfoAssociatedClusterMap));
        if(m_visualiseAllMatches||m_visualiseBadMatches)this->VisualiseMatches(clusterToPfoMap,pfoAssociatedClusterMap);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CosmicRayShowerMatching(clusterToPfoMap,pfoAssociatedClusterMap));

        return STATUS_CODE_SUCCESS;

    }

    //-------------------------------------------------------------------------------------------------------------------------------------------

    StatusCode CosmicRayShowerMatchingAlgorithm::CosmicRayShowerMatching(ctopmap_t &clusterToPfoMap, ptocmultimap_t &pfoAssociatedClusterMap) const
    {
    
        
        // associated clusters to matched cosmics
        
        const PfoList *pCosmicPfoList = NULL;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCosmicRayPfoListName, pCosmicPfoList));
        for (PfoList::const_iterator iterPfo = pCosmicPfoList->begin(), iterPfoEnd = pCosmicPfoList->end(); iterPfo != iterPfoEnd; ++iterPfo)
        {
            ParticleFlowObject *pPfo = *iterPfo;
            
            
            bool neutrinoPfo(false);
            //bool foundMCParticle(false); // --- UNUSED ---
            ClusterList pfoClusterList(pPfo->GetClusterList());
            for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
            {
                const Cluster *const pPfoCluster = *cIter;
                try {
                    const MCParticle *pPfoMCParticle(MCParticleHelper::GetMainMCParticle(pPfoCluster));
                    const MCParticle *pPfoParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pPfoMCParticle));
                    if(abs(pPfoParentMCParticle->GetParticleId())==12 || abs(pPfoParentMCParticle->GetParticleId())==14 )neutrinoPfo = true;
                    //foundMCParticle= true; // --- UNUSED ---
                } catch (StatusCodeException &) {                }
            }
            
            
            
            if(pfoAssociatedClusterMap.count(pPfo)!=0){
                std::pair<ptocIter,ptocIter> pos;
                pos = pfoAssociatedClusterMap.equal_range(pPfo);
                for(ptocIter it =pos.first; it !=pos.second; it++){
                    Cluster* pCluster = const_cast<Cluster*>(it->second);
                    bool neutrinoCluster(false);
                    try {
                        const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
                        const MCParticle *pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
                        if(abs(pParentMCParticle->GetParticleId())==12 || abs(pParentMCParticle->GetParticleId())==14 )neutrinoCluster = true;
                    } catch (StatusCodeException &) {     }
                    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo, pCluster));
                    if(m_debugPrintingOn){
                        if(neutrinoCluster && !neutrinoPfo)std::cout << " Error   - part of neutrino matched to COSMIC? " << std::endl;
                        if(neutrinoCluster && neutrinoPfo)std::cout  << " success - part of neutrino matched to neutrino " << std::endl;
                    }
                }
            }
        }

        return STATUS_CODE_SUCCESS;
    }
    
    
    //------------------------------------------------------------------------------------------------------------------------------------------
    
    StatusCode CosmicRayShowerMatchingAlgorithm::VisualiseMatches(ctopmap_t &clusterToPfoMap, ptocmultimap_t &pfoAssociatedClusterMap) const
    {
        
        
        // visulaise all matches - mote just the cosmics
        
        PfoList tempPfoList;
        this->GetAllPfos(tempPfoList);

       
        for (PfoList::const_iterator iterPfo = tempPfoList.begin(), iterPfoEnd = tempPfoList.end(); iterPfo != iterPfoEnd; ++iterPfo)
        {
            const ParticleFlowObject *pPfo = *iterPfo;
            if(pfoAssociatedClusterMap.count(pPfo)!=0){
                ClusterList tempList;
                ClusterList tempList1;
                bool foundNeutrinoCluster(false);
                std::pair<ptocIter,ptocIter> pos;
                pos = pfoAssociatedClusterMap.equal_range(pPfo);
                for(ptocIter it =pos.first; it !=pos.second; it++){
                    Cluster* pCluster = const_cast<Cluster*>(it->second);
                    bool neutrinoCluster(false);
                    try {
                        const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
                        const MCParticle *pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
                        if(abs(pParentMCParticle->GetParticleId())==12 || abs(pParentMCParticle->GetParticleId())==14 )neutrinoCluster = true;

                    } catch (StatusCodeException &) {  }
                    if(neutrinoCluster)foundNeutrinoCluster=true;
                    neutrinoCluster ? tempList1.insert(pCluster) : tempList.insert(pCluster);
                }
  
                
                
                
                bool neutrinoPfo(false);
                bool foundMCParticle(false);
                ClusterList pfoClusterList(pPfo->GetClusterList());
                for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
                {
                    const Cluster *const pPfoCluster = *cIter;
                    try {
                        const MCParticle *pPfoMCParticle(MCParticleHelper::GetMainMCParticle(pPfoCluster));
                        const MCParticle *pPfoParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pPfoMCParticle));
                        if(abs(pPfoParentMCParticle->GetParticleId())==12 || abs(pPfoParentMCParticle->GetParticleId())==14 )neutrinoPfo = true;
                        foundMCParticle= true;
                    } catch (StatusCodeException &) {  }
                }
                bool badMatch(false);
                if(neutrinoPfo && !foundNeutrinoCluster)badMatch=true;
                if(!neutrinoPfo && foundNeutrinoCluster)badMatch=true;
                if(m_visualiseAllMatches||(m_visualiseBadMatches&& badMatch)){
                    if(foundMCParticle){
                        if(neutrinoPfo)PANDORA_MONITORING_API(VisualizeClusters(&pfoClusterList, "PfoCluster", BLACK));
                        if(!neutrinoPfo)PANDORA_MONITORING_API(VisualizeClusters(&pfoClusterList, "PfoCluster", GREEN));
                    }else{
                        PANDORA_MONITORING_API(VisualizeClusters(&pfoClusterList, "PfoCluster", MAGENTA));
                    }
                    PANDORA_MONITORING_API(VisualizeClusters(&tempList, "Cluster", BLUE));
                    PANDORA_MONITORING_API(VisualizeClusters(&tempList1, "Neutrino Cluster", RED));
                    //  PANDORA_MONITORING_API(AddMarkerToVisualization(&splitPosition, "SplitPosition", RED, 2.75));
                    PANDORA_MONITORING_API(ViewEvent());
                }
            }
            
        }
        
        return STATUS_CODE_SUCCESS;
        
    }
    
    //------------------------------------------------------------------------------------------------------------------------------------------
    
    StatusCode CosmicRayShowerMatchingAlgorithm::CosmicRay3DShowerMatching(ctopmap_t &clusterToPfoMap, ptocmultimap_t &pfoAssociatedClusterMap) const
    {
    
        std::vector<const Cluster*> UClusters, VClusters, WClusters;
        
        // Loop over W clusters
        for (StringVector::const_iterator sIterW = m_inputClusterListNamesW.begin(), sIterEndW = m_inputClusterListNamesW.end(); sIterW != sIterEndW; ++sIterW){
            const ClusterList *pClusterListW = NULL;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, *sIterW, pClusterListW))
            for (ClusterList::const_iterator cIterW = pClusterListW->begin(), cIterWEnd = pClusterListW->end(); cIterW != cIterWEnd; ++cIterW)WClusters.push_back(*cIterW);
        }
        
        // Loop over U clusters
        for (StringVector::const_iterator sIterU = m_inputClusterListNamesU.begin(), sIterEndU = m_inputClusterListNamesU.end(); sIterU != sIterEndU; ++sIterU){
            const ClusterList *pClusterListU = NULL;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, *sIterU, pClusterListU));
            for(ClusterList::const_iterator cIterU = pClusterListU->begin(), cIterUEnd = pClusterListU->end(); cIterU != cIterUEnd; ++cIterU)UClusters.push_back(*cIterU);
        }
    
        // Loop over V clusters
        for (StringVector::const_iterator sIterV = m_inputClusterListNamesV.begin(), sIterEndV = m_inputClusterListNamesV.end(); sIterV != sIterEndV; ++sIterV){
            const ClusterList *pClusterListV = NULL;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, *sIterV, pClusterListV));
            for(ClusterList::const_iterator cIterV = pClusterListV->begin(), cIterVEnd = pClusterListV->end(); cIterV != cIterVEnd; ++cIterV)VClusters.push_back(*cIterV);
        }
        
        
        // look for U,V,W combinations
            
        for(unsigned int iW=0; iW<WClusters.size();iW++){
            const Cluster* const pWCluster = WClusters[iW];
            float xminW(std::numeric_limits<float>::max()), xmaxW(std::numeric_limits<float>::min());
            if( pWCluster->IsAvailable() && pWCluster->GetNCaloHits()>1   )
                LArClusterHelper::GetClusterSpanX(pWCluster,  xminW, xmaxW);
            for(unsigned int iU=0; iU<UClusters.size();iU++){
                const Cluster* const pUCluster = UClusters[iU];
                float xminU(std::numeric_limits<float>::max()), xmaxU(std::numeric_limits<float>::min());
                if( pUCluster->IsAvailable() && pUCluster->GetNCaloHits()>1  )
                    LArClusterHelper::GetClusterSpanX(pUCluster,  xminU, xmaxU);
                for(unsigned int iV=0; iV<VClusters.size();iV++){
                    const Cluster* const pVCluster = VClusters[iV];
                    float xminV(std::numeric_limits<float>::max()), xmaxV(std::numeric_limits<float>::min());
                    if( pVCluster->IsAvailable() && pVCluster->GetNCaloHits()>1  )
                        LArClusterHelper::GetClusterSpanX(pVCluster,  xminV, xmaxV);
                    const float xmin = std::max((std::max(xminU,xminV)),xminW);
                    const float xmax = std::min((std::min(xmaxU,xmaxV)),xmaxW);
                    // we now have an overlap in time in all three views
                    if(xmax>=xmin){
                          float pseudoChi2(std::numeric_limits<float>::max());
                        this->CompareClusterTriplet(pUCluster,pVCluster,pWCluster,pseudoChi2);
                        float distance(std::numeric_limits<float>::max());
                        if(pseudoChi2<m_chi2For3ViewMatching){
                            ParticleFlowObject *pBestPFO = NULL;
                            this->FindBestCosmicPFO(pUCluster,pVCluster,pWCluster, pBestPFO, distance);
                            if(m_debugPrintingOn){
                                std::cout << "UVW overlap : " << xminU << " - " << xmaxU << " and " << xminV << " - " << xmaxV << " and " << xminW << " - " << xmaxW << " : " << pUCluster->GetNCaloHits() << " " << pVCluster->GetNCaloHits() << " " << pWCluster->GetNCaloHits() << std::endl;

                                std::cout << " d = " << distance << " chi2 = " << pseudoChi2 << std::endl;
                            }
                            if(distance<m_distanceFor3ViewMatching)this->AssociateClusterWithPfo(clusterToPfoMap, pfoAssociatedClusterMap, pUCluster,pVCluster,pWCluster,pBestPFO);
                            if(distance>10 && pseudoChi2 < 5){
                                ParticleFlowObject* pPfo=NULL;
                                clusterToPfoMap.insert( ctopValType(pUCluster,pPfo) );
                                clusterToPfoMap.insert( ctopValType(pVCluster,pPfo) );
                                clusterToPfoMap.insert( ctopValType(pWCluster,pPfo) );
                            }
                            
                            if(m_debugPrintingOn){
                                try {
                                    ClusterList pfoClusterList(pBestPFO->GetClusterList());
                                    ClusterList::const_iterator cIter = pfoClusterList.begin();
                                    const Cluster *const pPfoCluster = *cIter;
                                    const MCParticle *pPfoMCParticle(MCParticleHelper::GetMainMCParticle(pPfoCluster));
                                    const MCParticle *pPfoParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pPfoMCParticle));
                                    const MCParticle *pMCParticleU(MCParticleHelper::GetMainMCParticle(pUCluster));
                                    const MCParticle *pParentMCParticleU(LArMCParticleHelper::GetParentMCParticle(pMCParticleU));
                                    const MCParticle *pMCParticleV(MCParticleHelper::GetMainMCParticle(pVCluster));
                                    const MCParticle *pParentMCParticleV(LArMCParticleHelper::GetParentMCParticle(pMCParticleV));
                                    // const MCParticle *pMCParticleW(MCParticleHelper::GetMainMCParticle(pWCluster)); // --- UNUSED ---
                                    // const MCParticle *pParentMCParticleW(LArMCParticleHelper::GetParentMCParticle(pMCParticleW)); // --- UNUSED ---
                                    if( (pParentMCParticleU==pParentMCParticleV) && (pParentMCParticleU==pParentMCParticleV)){
                                        if(pParentMCParticleU==pPfoParentMCParticle)std::cout << " GOOD MATCH " << std::endl;
                                        if(pParentMCParticleU!=pPfoParentMCParticle)std::cout << " CONSISTENT MATCH - WRONG PFO" << std::endl;
                                    }else{
                                    std::cout << " ERRONEOUS MATCH" << std::endl;
                                    }
                                } catch (StatusCodeException &) {
                                    std::cout << " Failed to find MC particle/parent" << std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }

        
        for(unsigned int iU=0; iU<UClusters.size();iU++){
            const Cluster* const pUCluster = UClusters[iU];
            float xminU(std::numeric_limits<float>::max()), xmaxU(std::numeric_limits<float>::min());
            if( pUCluster->IsAvailable() && pUCluster->GetNCaloHits()>1  )
                LArClusterHelper::GetClusterSpanX(pUCluster,  xminU, xmaxU);
            for(unsigned int iV=0; iV<VClusters.size();iV++){
                const Cluster* const pVCluster = VClusters[iV];
                float xminV(std::numeric_limits<float>::max()), xmaxV(std::numeric_limits<float>::min());
                if( pVCluster->IsAvailable() && pVCluster->GetNCaloHits()>1  )
                    LArClusterHelper::GetClusterSpanX(pVCluster,  xminV, xmaxV);
                for(unsigned int iW=0; iW<WClusters.size();iW++){
                    const Cluster* const pWCluster = WClusters[iW];
                    float xminW(std::numeric_limits<float>::max()), xmaxW(std::numeric_limits<float>::min());
                    if( pWCluster->IsAvailable() && pWCluster->GetNCaloHits()>1   )
                        LArClusterHelper::GetClusterSpanX(pWCluster,  xminW, xmaxW);
                    if(std::min(xmaxU,xmaxV)>std::max(xminU,xminV) && iW==0){
                        float distance(std::numeric_limits<float>::max());
                        ParticleFlowObject *pBestPFO = NULL;
                        this->FindBestCosmicPFO(pUCluster, pVCluster, pBestPFO, distance);
                        if(m_debugPrintingOn)std::cout << "UV overlap : " << xminU << " - " << xmaxU << " and " << xminV << " - " << xmaxV << " : " << pUCluster->GetNCaloHits() << " " << pVCluster->GetNCaloHits() <<  "d = " << distance << std::endl;
                        if(distance<m_distanceFor2ViewMatching)this->AssociateClusterWithPfo(clusterToPfoMap, pfoAssociatedClusterMap,pUCluster,pVCluster,pBestPFO);

                    }
                    if(std::min(xmaxU,xmaxW)>std::max(xminU,xminW) && iV==0){
                        
                        float distance(std::numeric_limits<float>::max());
                        ParticleFlowObject *pBestPFO = NULL;
                        this->FindBestCosmicPFO(pUCluster, pWCluster, pBestPFO, distance);
                        if(m_debugPrintingOn)std::cout << "UW overlap : " << xminU << " - " << xmaxU << " and " << xminW << " - " << xmaxW << " : " << pUCluster->GetNCaloHits() << " " << pWCluster->GetNCaloHits() << " d = " << distance  << std::endl;
                        if(distance<m_distanceFor2ViewMatching)this->AssociateClusterWithPfo(clusterToPfoMap, pfoAssociatedClusterMap, pUCluster,pWCluster,pBestPFO);
                    }
                    if(std::min(xmaxV,xmaxW)>std::max(xminV,xminW) && iU==0){
                       
                        float distance(std::numeric_limits<float>::max());
                        ParticleFlowObject *pBestPFO = NULL;
                        this->FindBestCosmicPFO(pVCluster, pWCluster, pBestPFO, distance);
                        if(m_debugPrintingOn)std::cout << "VW overlap : " << xminV << " - " << xmaxV << " and " << xminW << " - " << xmaxW << " : " << pVCluster->GetNCaloHits() << " " << pWCluster->GetNCaloHits() << " d = " << distance  << std::endl;
                        if(distance<m_distanceFor2ViewMatching)this->AssociateClusterWithPfo(clusterToPfoMap, pfoAssociatedClusterMap,  pVCluster,pWCluster,pBestPFO);
                    }
                }
            }
        }
        
        // now the rest
        for(unsigned int iU=0; iU<UClusters.size();iU++){
            const Cluster* const pUCluster = UClusters[iU];
            float distance(std::numeric_limits<float>::max());
            ParticleFlowObject *pBestPFO = NULL;
            if( pUCluster->IsAvailable() && pUCluster->GetNCaloHits()>0 )this->FindBestCosmicPFO(pUCluster, pBestPFO, distance);
            if(distance<m_distanceFor1ViewMatching){
                this->AssociateClusterWithPfo(clusterToPfoMap, pfoAssociatedClusterMap,pUCluster,pBestPFO);
                if(m_debugPrintingOn)std::cout << " UCluster " << pUCluster->GetNCaloHits() << " d = " << distance << std::endl;
            }
        }
        for(unsigned int iV=0; iV<VClusters.size();iV++){
            const Cluster* const pVCluster = VClusters[iV];
            float distance(std::numeric_limits<float>::max());
            ParticleFlowObject *pBestPFO = NULL;
            if( pVCluster->IsAvailable() && pVCluster->GetNCaloHits()>0 )this->FindBestCosmicPFO(pVCluster, pBestPFO, distance);
            if(distance<m_distanceFor1ViewMatching){
                this->AssociateClusterWithPfo(clusterToPfoMap, pfoAssociatedClusterMap, pVCluster,pBestPFO);
                if(m_debugPrintingOn)std::cout << " VCluster " << pVCluster->GetNCaloHits() << " d = " << distance << std::endl;
            }
        }
        for(unsigned int iW=0; iW<WClusters.size();iW++){
            const Cluster* const pWCluster = WClusters[iW];
            float distance(std::numeric_limits<float>::max());
            ParticleFlowObject *pBestPFO = NULL;
            if( pWCluster->IsAvailable() && pWCluster->GetNCaloHits()>0 )this->FindBestCosmicPFO(pWCluster, pBestPFO, distance);
            if(distance<m_distanceFor1ViewMatching){
                this->AssociateClusterWithPfo(clusterToPfoMap, pfoAssociatedClusterMap, pWCluster,pBestPFO);
                if(m_debugPrintingOn)std::cout << " WCluster " << pWCluster->GetNCaloHits() << " d = " << distance << std::endl;
            }
        }
        
        return STATUS_CODE_SUCCESS;
        
    }
    

    //------------------------------------------------------------------------------------------------------------------------------------------
    
    StatusCode CosmicRayShowerMatchingAlgorithm::CompareClusterTriplet(const Cluster *const pClusterU, const Cluster *const pClusterV,
                                                                   const Cluster *const pClusterW, float &pseudoChi2) const
    {
 
        // find x overlap
        pseudoChi2 = std::numeric_limits<float>::max();
        CartesianVector minimumCoordinates(std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<float>::max());
        CartesianVector maximumCoordinates(std::numeric_limits<float>::min(),std::numeric_limits<float>::min(),std::numeric_limits<float>::min());
        LArClusterHelper::GetClusterSpanXZ(pClusterU,  minimumCoordinates, maximumCoordinates);
        float xminU(minimumCoordinates.GetX()), xmaxU(maximumCoordinates.GetX());
        LArClusterHelper::GetClusterSpanXZ(pClusterV,  minimumCoordinates, maximumCoordinates);
        float xminV(minimumCoordinates.GetX()), xmaxV(maximumCoordinates.GetX());
        LArClusterHelper::GetClusterSpanXZ(pClusterW,  minimumCoordinates, maximumCoordinates);
        float xminW(minimumCoordinates.GetX()), xmaxW(maximumCoordinates.GetX());
        const float xmin = std::max(  (std::max(xminU,xminV)), xminW);
        const float xmax = std::min(  (std::min(xmaxU,xmaxV)), xmaxW);
        if(xmin>=xmax)return STATUS_CODE_FAILURE;
        const float x =(xmin+xmax)/2.;
        const int span =  std::min(    (std::min( (pClusterU)->GetNCaloHits(),  (pClusterV)->GetNCaloHits() )),(pClusterW)->GetNCaloHits());
        float u = this->GetCoordinateAtX(pClusterU,x,xmin,xmax,span);
        float v = this->GetCoordinateAtX(pClusterV,x,xmin,xmax,span);
        float w = this->GetCoordinateAtX(pClusterW,x,xmin,xmax,span);
        const float uv2w(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_V, u, v));
        const float uw2v(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_W, u, w));
        const float vw2u(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_V, TPC_VIEW_W, v, w));
        pseudoChi2 = ((u-vw2u)*(u-vw2u)+(v-uw2v)*(v-uw2v)+(w-uv2w)*(w-uv2w))/3.f;
    
        return STATUS_CODE_SUCCESS;
    }
    
    //------------------------------------------------------------------------------------------------------------------------------------------
    
    StatusCode CosmicRayShowerMatchingAlgorithm::FindBestCosmicPFO(const Cluster *const pClusterU, const Cluster *const pClusterV,
                                                                       const Cluster *const pClusterW, ParticleFlowObject* &pBestPFO, float &distanceToBestPFO) const
    {
        
        distanceToBestPFO = std::numeric_limits<float>::max();
        
        PfoList tempPfoList;
        this->GetAllPfos(tempPfoList);
        
        for (PfoList::const_iterator iter = tempPfoList.begin(), iterEnd = tempPfoList.end(); iter != iterEnd; ++iter)
        {
            ParticleFlowObject *pPfo = *iter;
            float dU(std::numeric_limits<float>::max());
            float dV(std::numeric_limits<float>::max());
            float dW(std::numeric_limits<float>::max());
            const ClusterList &pfoClusterList(pPfo->GetClusterList());
            for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
            {
                const Cluster *const pPfoCluster = *cIter;
                const HitType hitType(LArThreeDHelper::GetClusterHitType(pPfoCluster));
                if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
                {
                    if (TPC_3D == hitType)
                        continue;

                    std::cout << "CosmicRayShowerMatchingAlgorithm: Encountered unexpected hit type " << std::endl;
                    return STATUS_CODE_INVALID_PARAMETER;
                }
                if(hitType==TPC_VIEW_U)dU = std::min( dU, ClusterHelper::GetDistanceToClosestHit(pPfoCluster, pClusterU));
                if(hitType==TPC_VIEW_V)dV = std::min( dV, ClusterHelper::GetDistanceToClosestHit(pPfoCluster, pClusterV));
                if(hitType==TPC_VIEW_W)dW = std::min( dW, ClusterHelper::GetDistanceToClosestHit(pPfoCluster, pClusterW));
            }
            const float distanceToPFO = sqrt(dU*dU + dV*dV + dW*dW);

            
            
            if(distanceToPFO<distanceToBestPFO){
                distanceToBestPFO = distanceToPFO;
                pBestPFO = pPfo;
            }
        }

        return STATUS_CODE_SUCCESS;
    }
    
    //------------------------------------------------------------------------------------------------------------------------------------------
    
    StatusCode CosmicRayShowerMatchingAlgorithm::FindBestCosmicPFO(const Cluster *const pClusterView1, const Cluster *const pClusterView2,
                                                                   ParticleFlowObject* &pBestPFO, float &distanceToBestPFO) const
    {
        
        const HitType view1HitType(LArThreeDHelper::GetClusterHitType(pClusterView1));
        const HitType view2HitType(LArThreeDHelper::GetClusterHitType(pClusterView2));
        
        distanceToBestPFO = std::numeric_limits<float>::max();
 
        PfoList tempPfoList;
        this->GetAllPfos(tempPfoList);
                
        for (PfoList::const_iterator iter = tempPfoList.begin(), iterEnd = tempPfoList.end(); iter != iterEnd; ++iter){
            float d1(std::numeric_limits<float>::max());
            float d2(std::numeric_limits<float>::max());
            ParticleFlowObject *pPfo = *iter;
            const ClusterList &pfoClusterList(pPfo->GetClusterList());
            
            for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
            {
                const Cluster *const pPfoCluster = *cIter;
                const HitType hitType(LArThreeDHelper::GetClusterHitType(pPfoCluster));
                
                if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
                {
                    if (TPC_3D == hitType)
                        continue;

                    std::cout << "CosmicRayShowerMatchingAlgorithm: Encountered unexpected hit type " << std::endl;
                    return STATUS_CODE_INVALID_PARAMETER;
                }
                if(hitType==view1HitType)d1 = std::min( d1, ClusterHelper::GetDistanceToClosestHit(pPfoCluster, pClusterView1));
                if(hitType==view2HitType)d2 = std::min( d2, ClusterHelper::GetDistanceToClosestHit(pPfoCluster, pClusterView2));
            }
            const float distanceToPFO = sqrt(d1*d1+d2*d2);
            if(distanceToPFO<distanceToBestPFO){
                distanceToBestPFO = distanceToPFO;
                pBestPFO = pPfo;
            }
        }
        
        return STATUS_CODE_SUCCESS;
    }
    
    //------------------------------------------------------------------------------------------------------------------------------------------
    
    StatusCode CosmicRayShowerMatchingAlgorithm::FindBestCosmicPFO(const Cluster *const pClusterView1,
                                                                   ParticleFlowObject* &pBestPFO, float &distanceToBestPFO) const
    {
        
        const HitType view1HitType(LArThreeDHelper::GetClusterHitType(pClusterView1));
        
        distanceToBestPFO = std::numeric_limits<float>::max();
        
 
        PfoList tempPfoList;
        this->GetAllPfos(tempPfoList);
               
        
        for (PfoList::const_iterator iter = tempPfoList.begin(), iterEnd = tempPfoList.end(); iter != iterEnd; ++iter){
            float d1(std::numeric_limits<float>::max());
            ParticleFlowObject *pPfo = *iter;
            const ClusterList &pfoClusterList(pPfo->GetClusterList());
            
            for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
            {
                const Cluster *const pPfoCluster = *cIter;
                const HitType hitType(LArThreeDHelper::GetClusterHitType(pPfoCluster));
                
                if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
                {
                    if (TPC_3D == hitType)
                        continue;

                    std::cout << "CosmicRayShowerMatchingAlgorithm: Encountered unexpected hit type " << std::endl;
                    return STATUS_CODE_INVALID_PARAMETER;
                }
                if(hitType==view1HitType)d1 = std::min( d1, ClusterHelper::GetDistanceToClosestHit(pPfoCluster, pClusterView1));
            }
            const float distanceToPFO = d1;
            if(distanceToPFO<distanceToBestPFO){
                distanceToBestPFO = distanceToPFO;
                pBestPFO = pPfo;
            }
        }
        
        return STATUS_CODE_SUCCESS;
    }
    
    
    
    //------------------------------------------------------------------------------------------------------------------------------------------

    float CosmicRayShowerMatchingAlgorithm::GetCoordinateAtX(const Cluster *const pCluster, const float x, const float xmin, const float xmax, const int span) const
    {

        bool fitFailed(false);
    
        CartesianVector fitVector(0.f, 0.f, 0.f);
        const bool useX(true);
    
        // first try sliding fits
        TwoDSlidingFitResult slidingFitResult;
        try {
            LArClusterHelper::LArTwoDSlidingFit(pCluster, span, slidingFitResult);
            slidingFitResult.GetGlobalFitPosition(x, useX, fitVector);
        } catch (StatusCodeException &){fitFailed = true; }
    
    
        // the resort to simple averages if fits have failed
        if(fitFailed)fitVector.SetValues(x, 0.f, LArClusterHelper::GetAverageZ(pCluster,xmin,xmax));
        return fitVector.GetZ();
    
    }
                            
 //----------------------------------------------------------------------------------------------------------------------------
                            
                            
bool CosmicRayShowerMatchingAlgorithm::SortPfosByNHits(const ParticleFlowObject *const pLhs, const ParticleFlowObject *const pRhs)
{
    unsigned int nHitsLhs(0);
    for (ClusterList::const_iterator iter = pLhs->GetClusterList().begin(), iterEnd = pLhs->GetClusterList().end(); iter != iterEnd; ++iter)
        nHitsLhs += (*iter)->GetNCaloHits();

    unsigned int nHitsRhs(0);
    for (ClusterList::const_iterator iter = pRhs->GetClusterList().begin(), iterEnd = pRhs->GetClusterList().end(); iter != iterEnd; ++iter)
        nHitsRhs += (*iter)->GetNCaloHits();

    return (nHitsLhs > nHitsRhs);
}

//--------------------------------------------------------------------------------------------------------------------------------------


    StatusCode CosmicRayShowerMatchingAlgorithm::AssociateClusterWithPfo(ctopmap_t &clusterToPfoMap, ptocmultimap_t &pfoAssociatedClusterMap, const Cluster *const pCluster, const ParticleFlowObject *const pPfo) const

    {

        ctopIter iter = clusterToPfoMap.find(pCluster);
        if(iter!=clusterToPfoMap.end()){
            if((*iter).second == pPfo){
                //std::cout << "Found in same PFO" << std::endl;
            }else{
                //std::cout << "ALREADY IN DIFFERENT PFO" << std::endl;
            }
        }else{

            clusterToPfoMap.insert( ctopValType(pCluster,pPfo) );
            pfoAssociatedClusterMap.insert(ptocValType(pPfo,pCluster));
            
            bool neutrinoCluster(false);
            try {
                const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
                const MCParticle *pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
                if(abs(pParentMCParticle->GetParticleId())==12 || abs(pParentMCParticle->GetParticleId())==14 )neutrinoCluster = true;
            } catch (StatusCodeException &) { }
            
            
            bool neutrinoPfo(false);
            ClusterList pfoClusterList(pPfo->GetClusterList());
            for (ClusterList::const_iterator cIter = pfoClusterList.begin(), cIterEnd = pfoClusterList.end(); cIter != cIterEnd; ++cIter)
            {
                const Cluster *const pPfoCluster = *cIter;
                try {
                    const MCParticle *pPfoMCParticle(MCParticleHelper::GetMainMCParticle(pPfoCluster));
                    const MCParticle *pPfoParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pPfoMCParticle));
                    if(abs(pPfoParentMCParticle->GetParticleId())==12 || abs(pPfoParentMCParticle->GetParticleId())==14 )neutrinoPfo = true;
                } catch (StatusCodeException &) { }
            }
            if(m_debugPrintingOn){
                if(!neutrinoPfo&& !neutrinoCluster)std::cout  << " SUCCESS   Added cosmic cluster to cosmic pfo" << std::endl;
                if(neutrinoPfo&& neutrinoCluster)std::cout  << " SUCCESS   Added neutrino cluster to neutrino pfo" << std::endl;
                if(!neutrinoPfo&& neutrinoCluster)std::cout << " ERROR - - - - - Added neutrino cluster to cosmic pfo < - - - - - - - - - - " << std::endl;
                if(neutrinoPfo&& !neutrinoCluster)std::cout << " ERROR - - - - - Added cosmic cluster to neutrino pfo < - - - - - - - - - - " << std::endl;
            }
        }
        return STATUS_CODE_SUCCESS;
    }
    
    //------------------------------------------------------------------------------------------------------------------------------------------
    
    StatusCode CosmicRayShowerMatchingAlgorithm::AssociateClusterWithPfo(ctopmap_t &clusterToPfoMap, ptocmultimap_t &pfoAssociatedClusterMap, const Cluster *const pC1, const Cluster *const pC2, const ParticleFlowObject *const pPfo) const
    {
        this->AssociateClusterWithPfo(clusterToPfoMap, pfoAssociatedClusterMap,pC1,pPfo);
        this->AssociateClusterWithPfo(clusterToPfoMap, pfoAssociatedClusterMap,pC2,pPfo);
        return STATUS_CODE_SUCCESS;
        
    }
    
    //------------------------------------------------------------------------------------------------------------------------------------------
    
    StatusCode CosmicRayShowerMatchingAlgorithm::AssociateClusterWithPfo(ctopmap_t &clusterToPfoMap, ptocmultimap_t &pfoAssociatedClusterMap, const Cluster *const pC1, const Cluster *const pC2, const Cluster *const pC3, const ParticleFlowObject *const pPfo) const
    {
        
        this->AssociateClusterWithPfo(clusterToPfoMap, pfoAssociatedClusterMap, pC1,pPfo);
        this->AssociateClusterWithPfo(clusterToPfoMap, pfoAssociatedClusterMap, pC2,pPfo);
        this->AssociateClusterWithPfo(clusterToPfoMap, pfoAssociatedClusterMap, pC3,pPfo);
        
        return STATUS_CODE_SUCCESS;
        
    }

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayShowerMatchingAlgorithm::GetAllPfos(PfoList &tempPfoList) const
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCosmicRayPfoListName, pPfoList));

    const PfoList *pOtherPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pOtherPfoList));

    tempPfoList.insert(pPfoList->begin(), pPfoList->end());
    tempPfoList.insert(pOtherPfoList->begin(), pOtherPfoList->end());

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CosmicRayShowerMatchingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "InputCosmicRayPfoListName", m_inputCosmicRayPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "InputPfoListName", m_inputPfoListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesU", m_inputClusterListNamesU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesV", m_inputClusterListNamesV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNamesW", m_inputClusterListNamesW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DebugPrintingOn", m_debugPrintingOn));

    m_visualiseBadMatches = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualiseBadMatches", m_visualiseBadMatches));

    m_visualiseAllMatches = false;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualiseAllMatches", m_visualiseAllMatches));

    m_chi2For3ViewMatching = 100.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Chi2For3ViewMatching", m_chi2For3ViewMatching));

    m_distanceFor3ViewMatching = 10.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceFor3ViewMatching", m_distanceFor3ViewMatching));

    m_distanceFor2ViewMatching = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceFor2ViewMatching", m_distanceFor2ViewMatching));

    m_distanceFor1ViewMatching = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceFor1ViewMatching", m_distanceFor1ViewMatching));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
