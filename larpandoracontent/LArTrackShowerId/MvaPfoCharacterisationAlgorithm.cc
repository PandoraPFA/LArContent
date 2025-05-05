/**
 *  @file   larpandoracontent/LArTrackShowerId/MvaPfoCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the mva pfo characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArTrackShowerId/MvaPfoCharacterisationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

template <typename T>
MvaPfoCharacterisationAlgorithm<T>::MvaPfoCharacterisationAlgorithm() :
    m_persistFeatures(false),
    m_trainingSetMode(false),
    m_testBeamMode(false),
    m_enableProbability(true),
    m_useThreeDInformation(true),
    m_minProbabilityCut(0.5f),
    m_minCaloHitsCut(5),
    m_applyFiducialCut(false),
    m_fiducialMinX(-std::numeric_limits<float>::max()),
    m_fiducialMaxX(std::numeric_limits<float>::max()),
    m_fiducialMinY(-std::numeric_limits<float>::max()),
    m_fiducialMaxY(std::numeric_limits<float>::max()),
    m_fiducialMinZ(-std::numeric_limits<float>::max()),
    m_fiducialMaxZ(std::numeric_limits<float>::max()),
    m_applyReconstructabilityChecks(false),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool MvaPfoCharacterisationAlgorithm<T>::IsClearTrack(const Cluster *const pCluster) const
{
    if (pCluster->GetNCaloHits() < m_minCaloHitsCut)
        return false;

    StringVector featureOrder;
    const LArMvaHelper::MvaFeatureMap featureMap(LArMvaHelper::CalculateFeatures(m_algorithmToolNames, m_featureToolMap, featureOrder, this, pCluster));

    if (m_trainingSetMode)
    {
        bool isTrueTrack(false);

        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
            isTrueTrack = ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId())));
        }
        catch (const StatusCodeException &)
        {
        }

        LArMvaHelper::ProduceTrainingExample(m_trainingOutputFile, isTrueTrack, featureOrder, featureMap);
        return isTrueTrack;
    }

    if (!m_enableProbability)
    {
        return LArMvaHelper::Classify(m_mva, featureOrder, featureMap);
    }
    else
    {
        return (LArMvaHelper::CalculateProbability(m_mva, featureOrder, featureMap) > m_minProbabilityCut);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool MvaPfoCharacterisationAlgorithm<T>::IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const
{
    
    //std::cout << "==============================================" << std::endl;
    //std::cout << "Starting now with IsClearTrack..." << std::endl;

    if (!LArPfoHelper::IsThreeD(pPfo))
    {
        if (m_enableProbability)
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["TrackScore"] = -1.f;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
        }
        return (pPfo->GetParticleId() == MU_MINUS);
    }

    /* 
     *  Standard code!
     */

    // // Charge related features are only calculated using hits in W view
    // ClusterList wClusterList;
    // LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, wClusterList);

    // const PfoCharacterisationFeatureTool::FeatureToolMap &chosenFeatureToolMap(wClusterList.empty() ? m_featureToolMapNoChargeInfo : m_featureToolMapThreeD);
    // const StringVector chosenFeatureToolOrder(wClusterList.empty() ? m_algorithmToolNamesNoChargeInfo : m_algorithmToolNames);
    // StringVector featureOrder;
    // const LArMvaHelper::MvaFeatureMap featureMap(
    //     LArMvaHelper::CalculateFeatures(chosenFeatureToolOrder, chosenFeatureToolMap, featureOrder, this, pPfo));

    //std::cout << "Looking for the Collection plane..." << std::endl;

    /*
     *  Find Collection plane.
     */

    bool IsCollectionEmpty;

    // W (Ind-1)
    bool selectViewW = true;
    ClusterList wClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, wClusterList);
    float minX(9999.), maxX(9999.);
    if (wClusterList.empty()) { selectViewW = false; }
    else { wClusterList.front()->GetClusterSpanX(minX, maxX); }

    // U (Ind-2 or Coll, based on TPC)
    bool selectViewU = true;
    ClusterList uClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_VIEW_U, uClusterList);
    float minX_U(9999.), maxX_U(9999.);
    if (uClusterList.empty()) { selectViewU = false; }
    else { uClusterList.front()->GetClusterSpanX(minX_U, maxX_U); }

    // V (Ind-2 or Coll, based on TPC)
    bool selectViewV = true;
    ClusterList vClusterList;
    LArPfoHelper::GetClusters(pPfo, TPC_VIEW_V, vClusterList);
    float minX_V(9999.), maxX_V(9999.);
    if (vClusterList.empty()) { selectViewV = false; }
    else  { vClusterList.front()->GetClusterSpanX(minX_V, maxX_V); }

    // If no view is available
    if (!selectViewW && !selectViewU && !selectViewV)
        IsCollectionEmpty = true;

    // Create the CaloHitList based on Collection whenever possible                                                                                      
    CaloHitList orderedCaloHitList;

    // Fallback to standard case if there are no U/V views                                                                               
    if (!selectViewU && !selectViewV) 
        this->OrderCaloHitsByDistanceToVertex(wClusterList.front(), orderedCaloHitList);

    // Find the drift span with the available views, and find out whether the particle is cathode-crossing
    if (selectViewU && selectViewV) { minX = std::min(minX_U, minX_V); maxX = std::max(maxX_U, maxX_V); }
    else if(selectViewU && !selectViewV) { minX = minX_U; maxX = maxX_U; }
    else if(!selectViewU && selectViewV) { minX = minX_V; maxX = maxX_V; }
    bool particleCrossingCathode = false;
    if(LocatePointInCryostat(minX) == PositionInCryostat::BelowCathode &&
       LocatePointInCryostat(maxX) == PositionInCryostat::AboveCathode)
        particleCrossingCathode = true;

    //std::cout << "Particle crosses the cathode? " << particleCrossingCathode << std::endl;
    //std::cout << "Particle X span " << minX << "\t" << maxX << std::endl;
    
    // Cathode-crossing particle
    if(!particleCrossingCathode)
    {

        //std::cout << "Particle does not cross the cathode..." << std::endl;

        // TPC 2/3 -> V
        if(LocatePointInCryostat(minX) == PositionInCryostat::AboveCathode && selectViewV) {                                                            
            this->OrderCaloHitsByDistanceToVertex(vClusterList.front(), orderedCaloHitList);
        }

        // TPC 0/1 -> U
        else if(LocatePointInCryostat(maxX) == PositionInCryostat::BelowCathode && selectViewU) {                                          
            this->OrderCaloHitsByDistanceToVertex(uClusterList.front(), orderedCaloHitList);
        }

        // TPC 2/3, but V ill-defined
        else if( LocatePointInCryostat( minX ) == PositionInCryostat::AboveCathode && !selectViewV ) {
            this->OrderCaloHitsByDistanceToVertex(uClusterList.front(),orderedCaloHitList);
        }

        // TPC 0/1, but U ill-defined
        else if( LocatePointInCryostat( maxX ) == PositionInCryostat::BelowCathode && !selectViewU ) {
            this->OrderCaloHitsByDistanceToVertex(vClusterList.front(), orderedCaloHitList);
        }

        //std::cout << "-> Chose U or V." << std::endl;
    }


    // Particles that cross the cathode (need both U and V)
    else {

        //std::cout << "Particle does cross the cathode..." << std::endl;

        if (!selectViewU || !selectViewV) {
            //std::cout << "-> But falling back to Induction-1 either way." << std::endl;
            this->OrderCaloHitsByDistanceToVertex(wClusterList.front(), orderedCaloHitList);
        }

        else {                                                                                                                     
                                                                                                                                                                               
            CaloHitList orderedCaloHitListU;
            this->OrderCaloHitsByDistanceToVertex(uClusterList.front(), orderedCaloHitListU);
                                                                                                                                                                                 
            CaloHitList orderedCaloHitListV;
            this->OrderCaloHitsByDistanceToVertex(vClusterList.front(), orderedCaloHitListV);

            if (!orderedCaloHitListU.empty() && !orderedCaloHitListV.empty()) {                                                   
                
                // get vertex
                const VertexList *pVertexList(nullptr);
                (void) PandoraContentApi::GetCurrentList(*this, pVertexList);
                const float pVertexX = pVertexList->front()->GetPosition().GetX();
                
                // TPC 2/3 -> V, then U
                if(LocatePointInCryostat(pVertexX) == PositionInCryostat::AboveCathode) 
                    this->CombineCaloHitListsToHaveCollection(orderedCaloHitListV, orderedCaloHitListU, orderedCaloHitList);

                // TPC 0/1 -> U, then V
                else if(LocatePointInCryostat(pVertexX) == PositionInCryostat::BelowCathode)
                    this->CombineCaloHitListsToHaveCollection(orderedCaloHitListU, orderedCaloHitListV, orderedCaloHitList);

                // vertex within cathode, cannot enstablish the ordering
                else {
                    this->CombineCaloHitListsToHaveCollection(orderedCaloHitListU, orderedCaloHitListV, orderedCaloHitList);
                }
            } 

            // this means no vertex exists and either list is fine (similarly to what happens for W)                                                                                                                                                               
            else                                                    
                orderedCaloHitList = orderedCaloHitListU;

            //std::cout << "-> Merged U and V views." << std::endl;
        }
      }

    // Determine if the Collection plane hit list is empty
    IsCollectionEmpty = orderedCaloHitList.empty();

    // clear hit list
    orderedCaloHitList.clear();

    //std::cout << "Created the Collection hit list. Is it empty? " << IsCollectionEmpty << std::endl;

    /*
     *  End of new code. After this, change wClusterList.empty() to IsCollectionEmpty.
     */

    // Use the dedicated BDT without charge information if the Collection plane is not usable
    const PfoCharacterisationFeatureTool::FeatureToolMap &chosenFeatureToolMap(IsCollectionEmpty ? m_featureToolMapNoChargeInfo : m_featureToolMapThreeD);
    const StringVector chosenFeatureToolOrder(IsCollectionEmpty ? m_algorithmToolNamesNoChargeInfo : m_algorithmToolNames);
    StringVector featureOrder;
    const LArMvaHelper::MvaFeatureMap featureMap(
        LArMvaHelper::CalculateFeatures(chosenFeatureToolOrder, chosenFeatureToolMap, featureOrder, this, pPfo));

    for (auto const &[featureKey, featureValue] : featureMap)
    {
        (void)featureKey;

        if (!featureValue.IsInitialized())
        {
            if (m_enableProbability)
            {
                object_creation::ParticleFlowObject::Metadata metadata;
                metadata.m_propertiesToAdd["TrackScore"] = -1.f;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
            }
            return (pPfo->GetParticleId() == MU_MINUS);
        }
    }

    if (m_trainingSetMode && m_applyReconstructabilityChecks)
    {
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

        // Mapping target MCParticles -> truth associated Hits
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        if (!m_testBeamMode)
            LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, m_primaryParameters,
                LArMCParticleHelper::IsBeamNeutrinoFinalState, targetMCParticleToHitsMap);
        else
            LArMCParticleHelper::SelectReconstructableMCParticles(
                pMCParticleList, pCaloHitList, m_primaryParameters, LArMCParticleHelper::IsBeamParticle, targetMCParticleToHitsMap);

        LArMCParticleHelper::MCContributionMapVector mcParticlesToGoodHitsMaps({targetMCParticleToHitsMap});

        LArMCParticleHelper::PfoContributionMap pfoToReconstructable2DHitsMap;
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(
            PfoList(1, pPfo), mcParticlesToGoodHitsMaps, pfoToReconstructable2DHitsMap, m_primaryParameters.m_foldBackHierarchy);
        if (pfoToReconstructable2DHitsMap.empty())
            return false;

        LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;
        LArMCParticleHelper::MCParticleToPfoHitSharingMap mcParticleToPfoHitSharingMap;
        LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(
            pfoToReconstructable2DHitsMap, mcParticlesToGoodHitsMaps, pfoToMCParticleHitSharingMap, mcParticleToPfoHitSharingMap);
        if (pfoToMCParticleHitSharingMap.empty())
            return false;

        unsigned int nHitsInBestMCParticleTotal(0);
        unsigned int nHitsSharedWithBestMCParticleTotal(0);
        int bestMCParticlePdgCode(0);
        CartesianVector threeDVertexPosition(0.f, 0.f, 0.f);
        float hitsShower(0), hitsTrack(0);
        const LArMCParticleHelper::MCParticleToSharedHitsVector &mcParticleToSharedHitsVector(pfoToMCParticleHitSharingMap.at(pPfo));

        for (const LArMCParticleHelper::MCParticleCaloHitListPair &mcParticleCaloHitListPair : mcParticleToSharedHitsVector)
        {
            const pandora::MCParticle *const pAssociatedMCParticle(mcParticleCaloHitListPair.first);
            const CaloHitList &allMCHits(targetMCParticleToHitsMap.at(pAssociatedMCParticle));
            const CaloHitList &associatedMCHits(mcParticleCaloHitListPair.second);

            if ((PHOTON == pAssociatedMCParticle->GetParticleId()) || (E_MINUS == std::abs(pAssociatedMCParticle->GetParticleId())))
                hitsShower += associatedMCHits.size();
            else
                hitsTrack += associatedMCHits.size();

            if (associatedMCHits.size() > nHitsSharedWithBestMCParticleTotal)
            {
                nHitsSharedWithBestMCParticleTotal = associatedMCHits.size();
                nHitsInBestMCParticleTotal = allMCHits.size();
                bestMCParticlePdgCode = pAssociatedMCParticle->GetParticleId();
                threeDVertexPosition = pAssociatedMCParticle->GetVertex();
            }
        }

        const float trackShowerHitsRatio((hitsTrack + hitsShower) > 0 ? hitsTrack / (hitsTrack + hitsShower) : 0.f);
        const bool isTrueTrack(trackShowerHitsRatio >= 0.5);

        const int nHitsInPfoTotal(pfoToReconstructable2DHitsMap.at(pPfo).size());
        const float purity((nHitsInPfoTotal > 0) ? nHitsSharedWithBestMCParticleTotal / static_cast<float>(nHitsInPfoTotal) : 0.f);
        const float completeness(
            (nHitsInBestMCParticleTotal > 0) ? nHitsSharedWithBestMCParticleTotal / static_cast<float>(nHitsInBestMCParticleTotal) : 0.f);

        CaloHitList checkHitListW;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, checkHitListW);
        CaloHitList checkHitListU;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, checkHitListU);
        CaloHitList checkHitListV;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, checkHitListV);
        CaloHitList checkHitListAll;
        checkHitListAll.splice(checkHitListAll.end(), checkHitListW);
        checkHitListAll.splice(checkHitListAll.end(), checkHitListU);
        checkHitListAll.splice(checkHitListAll.end(), checkHitListV);

        LArMCParticleHelper::MCRelationMap mcPrimaryMap;
        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);

        LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;
        LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
        LArMCParticleHelper::GetMCParticleToCaloHitMatches(&checkHitListAll, mcPrimaryMap, hitToMCMap, mcToTrueHitListMap);

        unsigned int showerCount(0), allCount(0);
        for (const CaloHit *pHit : checkHitListAll)
        {
            if (hitToMCMap.find(pHit) != hitToMCMap.end())
            {
                const MCParticle *pHitMCParticle(hitToMCMap.at(pHit));
                if ((PHOTON == pHitMCParticle->GetParticleId()) || (E_MINUS == std::abs(pHitMCParticle->GetParticleId())))
                    ++showerCount;
                ++allCount;
            }
        }

        if (allCount == 0)
            return false;
        const float showerProbability(showerCount / static_cast<float>(allCount));
        const bool mischaracterisedPfo((showerProbability < 0.5f && !isTrueTrack) || (showerProbability > 0.5 && isTrueTrack) ? true : false);
        const bool isMainMCParticleSet(bestMCParticlePdgCode != 0);

        if (isMainMCParticleSet)
        {
            if (completeness >= 0.f && purity >= 0.f && !mischaracterisedPfo && (!m_applyFiducialCut || this->PassesFiducialCut(threeDVertexPosition)))
            {
                std::string outputFile(m_trainingOutputFile);
                const std::string end = (IsCollectionEmpty ? "noChargeInfo.txt" : ".txt");
                outputFile.append(end);
                LArMvaHelper::ProduceTrainingExample(outputFile, isTrueTrack, featureOrder, featureMap);
            }
        }

        return isTrueTrack;
    }
    else if (m_trainingSetMode)
    {
        //std::cout << "Entering training set mode..." << std::endl;

        bool isTrueTrack(false);
        bool isMainMCParticleSet(false);

        try
        {
            const MCParticle *const pMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
            isTrueTrack = ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId())));
            isMainMCParticleSet = (pMCParticle->GetParticleId() != 0);
        }
        catch (const StatusCodeException &)
        {
        }

        if (isMainMCParticleSet)
        {
            std::string outputFile(m_trainingOutputFile);
            outputFile.append(IsCollectionEmpty ? "noChargeInfo.txt" : ".txt");
            LArMvaHelper::ProduceTrainingExample(outputFile, isTrueTrack, featureOrder, featureMap);
        }

        return isTrueTrack;
    }

    // If no failures, proceed with MvaPfoCharacterisationAlgorithm classification
    if (!m_enableProbability)
    {
        //std::cout << "Proceeding with classification..." << std::endl;
        return LArMvaHelper::Classify((IsCollectionEmpty ? m_mvaNoChargeInfo : m_mva), featureOrder, featureMap);
        //std::cout << "Classified the particle. " << std::endl;
    }
    else
    {
        //std::cout << "Calculating the track score..." << std::endl;
        const double score(LArMvaHelper::CalculateProbability((IsCollectionEmpty ? m_mvaNoChargeInfo : m_mva), featureOrder, featureMap));
        //std::cout << "-> Calculated the track score: " << score << std::endl;

        object_creation::ParticleFlowObject::Metadata metadata;
        metadata.m_propertiesToAdd["TrackScore"] = score;
        if (m_persistFeatures)
        {
            for (auto const &[name, value] : featureMap)
            {
                metadata.m_propertiesToAdd[name] = value.Get();
            }
        }
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, metadata));
        //std::cout << "-> Updated metadata with features." << std::endl;

        //std::cout << "DONE WITH MVA!" << std::endl;
        //std::cout << "==============================================" << std::endl;
        
        return (m_minProbabilityCut <= score);
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void MvaPfoCharacterisationAlgorithm<T>::OrderCaloHitsByDistanceToVertex(
    // const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, CaloHitList &caloHitList)
    const pandora::Cluster *const pCluster, CaloHitList &caloHitList) const
{
    const VertexList *pVertexList(nullptr);
    (void)PandoraContentApi::GetCurrentList(*this, pVertexList);

    if (!pVertexList || pVertexList->empty())
        return;

    unsigned int nInteractionVertices(0);
    const Vertex *pInteractionVertex(nullptr);

    for (const Vertex *pVertex : *pVertexList)
    {
        if ((pVertex->GetVertexLabel() == VERTEX_INTERACTION) && (pVertex->GetVertexType() == VERTEX_3D))
        {
            ++nInteractionVertices;
            pInteractionVertex = pVertex;
        }
    }
    bool debug = false;
    
    if (pInteractionVertex && (1 == nInteractionVertices))
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), pInteractionVertex->GetPosition(), hitType));
        if (debug)
        {
            if( hitType == TPC_VIEW_U )
                std::cout << "[OrderCaloHitsByDistanceToVertex] This is interaction vertex (view U)" << vertexPosition2D.GetX() << "\t" << vertexPosition2D.GetY() << "\t" << vertexPosition2D.GetZ() << std::endl;
            else if( hitType == TPC_VIEW_V )
                std::cout << "[OrderCaloHitsByDistanceToVertex] This is interaction vertex (view V)" << vertexPosition2D.GetX() << "\t" << vertexPosition2D.GetY() << "\t" << vertexPosition2D.GetZ() << std::endl;
            else
                std::cout << "[OrderCaloHitsByDistanceToVertex] This is interaction vertex (view W)" << vertexPosition2D.GetX() << "\t" << vertexPosition2D.GetY() << "\t" << vertexPosition2D.GetZ() << std::endl;
        }
        CaloHitList clusterCaloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);

        clusterCaloHitList.sort(MvaPfoCharacterisationAlgorithm<T>::VertexComparator(vertexPosition2D));
        caloHitList.insert(caloHitList.end(), clusterCaloHitList.begin(), clusterCaloHitList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void MvaPfoCharacterisationAlgorithm<T>::CombineCaloHitListsToHaveCollection(
    const pandora::CaloHitList &orderedCaloHitList1, const pandora::CaloHitList &orderedCaloHitList2, 
    pandora::CaloHitList &mergedCaloHitList) const
 {
   bool debug = false;
   // identify the interaction vertex                                                                                                                                                                    
   const VertexList *pVertexList(nullptr);
   (void)PandoraContentApi::GetCurrentList(*this, pVertexList);
   const CartesianVector vertexPosition2DU(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_U));
   const CartesianVector vertexPosition2DV(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_V));
   if( debug )
     std::cout << "This is vertex position on view U " << vertexPosition2DU.GetX() << "\t " << vertexPosition2DU.GetY() << "\t " << vertexPosition2DU.GetZ() << std::endl;
   if( debug )
     std::cout << "This is vertex position on view V " << vertexPosition2DV.GetX() << "\t " << vertexPosition2DV.GetY() << "\t " << vertexPosition2DV.GetZ() << std::endl;
   float distance_hit_vertex = 0.;
   float min_distance_hit_vertex_included_1 = 0.;
   float max_distance_hit_vertex_included_1 = 0.;
   float min_distance_hit_vertex_included_2 = 0.;
   float max_distance_hit_vertex_included_2 = 0.;
   bool set_min_1 = false;
   bool set_min_2 = false;
   bool reorderingNeeded = false;
   // loop over first CaloHitList                                                                                                                                                                        
   /*                                                                                                                                                                                                    
    * view U: COLL below cathode, IND2 above cathode                                                                                                                                                     
    * view V: IND2 below cathode, COLL above cathode                                                                                                                                                     
    */
   // if this is from view V take hits ABOVE cathode otherwise below 
   for (CaloHitList::const_iterator hIter = orderedCaloHitList1.begin(); hIter != orderedCaloHitList1.end(); hIter++) {
     const CaloHit *const pCaloHit = *hIter;
     const CartesianVector &hit(pCaloHit->GetPositionVector());
     if( pCaloHit->GetHitType() == TPC_VIEW_U )
       distance_hit_vertex = pow( vertexPosition2DU.GetX() - hit.GetX(), 2 ) + pow( vertexPosition2DU.GetY()- hit.GetY(), 2 ) + pow( vertexPosition2DU.GetZ()- hit.GetZ(), 2 );
     else
       distance_hit_vertex = pow( vertexPosition2DV.GetX() - hit.GetX(), 2 ) + pow( vertexPosition2DV.GetY()- hit.GetY(), 2 ) + pow( vertexPosition2DV.GetZ()- hit.GetZ(), 2 );
     if( debug )
       std::cout << hit.GetX() << "\t" << hit.GetY() << "\t" << hit.GetZ() << " whose distance^2 from vertex is " << distance_hit_vertex << std::endl;
     bool selectHitViewU = pCaloHit->GetHitType() == TPC_VIEW_U && LocatePointInCryostat( hit.GetX() ) == PositionInCryostat::BelowCathode;
     bool selectHitViewV = pCaloHit->GetHitType() == TPC_VIEW_V && LocatePointInCryostat( hit.GetX() ) == PositionInCryostat::AboveCathode;
     bool selectHit = selectHitViewU || selectHitViewV;
     if( selectHit ){
       if( debug )
         std::cout << "This hit should be added to hit list to have COLLECTION plane" << std::endl;
       mergedCaloHitList.push_back( pCaloHit );
       if( !set_min_1 ) {
         min_distance_hit_vertex_included_1 = distance_hit_vertex;
         set_min_1 = true;
       }
       if( distance_hit_vertex > max_distance_hit_vertex_included_1 )
         max_distance_hit_vertex_included_1 = distance_hit_vertex;
     }
   }
   std::cout << "This is (min,max) from the 1st CaloHitList: ( " << min_distance_hit_vertex_included_1 << ", " << max_distance_hit_vertex_included_1 << " )." << std::endl;

   // loop over second CaloHitList                                                                                                                                                                       
   for (CaloHitList::const_iterator hIter = orderedCaloHitList2.begin(); hIter != orderedCaloHitList2.end(); hIter++) {
     const CaloHit *const pCaloHit = *hIter;
     const CartesianVector &hit(pCaloHit->GetPositionVector());
     if( pCaloHit->GetHitType() == TPC_VIEW_U )
       distance_hit_vertex = pow( vertexPosition2DU.GetX() - hit.GetX(), 2 ) + pow( vertexPosition2DU.GetY()- hit.GetY(), 2 ) + pow( vertexPosition2DU.GetZ()- hit.GetZ(), 2 );
     else
       distance_hit_vertex = pow( vertexPosition2DV.GetX() - hit.GetX(), 2 ) + pow( vertexPosition2DV.GetY()- hit.GetY(), 2 ) + pow( vertexPosition2DV.GetZ()- hit.GetZ(), 2 );
     if( debug )
       std::cout << hit.GetX() << "\t" << hit.GetY() << "\t" << hit.GetZ() << " whose distance^2 from vertex is " << distance_hit_vertex << std::endl;
     bool selectHitViewU = pCaloHit->GetHitType() == TPC_VIEW_U && LocatePointInCryostat( hit.GetX() ) == PositionInCryostat::BelowCathode;
     bool selectHitViewV = pCaloHit->GetHitType() == TPC_VIEW_V && LocatePointInCryostat( hit.GetX() ) == PositionInCryostat::AboveCathode;
     bool selectHit = selectHitViewU || selectHitViewV;
     if( selectHit ){
       if( debug )
         std::cout << "This hit should be added to hit list to have COLLECTION plane" << std::endl;
       mergedCaloHitList.push_back( pCaloHit );
       if( !set_min_2 ) {
         min_distance_hit_vertex_included_2 = distance_hit_vertex;
         set_min_2 = true;
       }
       if( distance_hit_vertex > max_distance_hit_vertex_included_2 )
         max_distance_hit_vertex_included_2 = distance_hit_vertex;
     }
   }
   std::cout << "This is (min,max) from the 2nd CaloHitList: ( " << min_distance_hit_vertex_included_2 << ", " << max_distance_hit_vertex_included_2 << " )." << std::endl;

   // reorder if needed                                                                                                                                                                                  
   reorderingNeeded = ( min_distance_hit_vertex_included_2 < max_distance_hit_vertex_included_1 ) ? 1 : 0;
   if( reorderingNeeded ){
     std::cerr << "[ThreeDChargeFeatureTool::CombineCaloHitListsToHaveCollection()] Reordering necessary, proceeding to do it." << std::endl;
     mergedCaloHitList.sort(MvaPfoCharacterisationAlgorithm<T>::DistanceToVertexComparator( vertexPosition2DU, vertexPosition2DV, TPC_VIEW_U, TPC_VIEW_V ) );
   }
   
   if( debug ){
     std::cout << " Listing the final content of the mergedCaloHitList " << std::endl;
     std::string isBeforeCathode = "";
     std::string belongsToView = "";
     for(CaloHitList::const_iterator it = mergedCaloHitList.begin() ; it != mergedCaloHitList.end() ; it++){
       const CaloHit *const pCaloHit = *it;
       const CartesianVector &hit(pCaloHit->GetPositionVector());
       if( pCaloHit->GetHitType() == TPC_VIEW_U )
         distance_hit_vertex = pow( vertexPosition2DU.GetX() - hit.GetX(), 2 ) + pow( vertexPosition2DU.GetY()- hit.GetY(), 2 ) + pow( vertexPosition2DU.GetZ()- hit.GetZ(), 2 );
       else
         distance_hit_vertex = pow( vertexPosition2DV.GetX() - hit.GetX(), 2 ) + pow( vertexPosition2DV.GetY()- hit.GetY(), 2 ) + pow( vertexPosition2DV.GetZ()- hit.GetZ(), 2 );
       isBeforeCathode = (LocatePointInCryostat( hit.GetX() ) == PositionInCryostat::BelowCathode )
         ? ", is before the cathode" : ", is after the cathode";
       belongsToView = ( pCaloHit->GetHitType() == TPC_VIEW_U ) ? " and belongs to view U " : " and belongs to view V";
       if( pCaloHit->GetHitType() == TPC_VIEW_U )
         std::cout << "This hit has distance " << distance_hit_vertex << isBeforeCathode << belongsToView
              << ". Its position is ( " << hit.GetX() << ", " << hit.GetY() << ", " << hit.GetZ() << " )."
              << " Vertex is ( " << vertexPosition2DU.GetX() << ", " <<  vertexPosition2DU.GetY() << ", " << vertexPosition2DU.GetZ() << " )"
              << std::endl;
       else
         std::cout << "This hit has distance " << distance_hit_vertex << isBeforeCathode << belongsToView
              << ". Its position is ( " << hit.GetX() << ", " << hit.GetY() << ", " << hit.GetZ() << " )."
              << " Vertex is ( " << vertexPosition2DV.GetX() << ", " <<  vertexPosition2DV.GetY() << ", " << vertexPosition2DV.GetZ() << " )"
              << std::endl;
     }
   }
 }

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
MvaPfoCharacterisationAlgorithm<T>::VertexComparator::VertexComparator(const CartesianVector vertexPosition2D) : m_neutrinoVertex(vertexPosition2D)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool MvaPfoCharacterisationAlgorithm<T>::VertexComparator::operator()(const CaloHit *const left, const CaloHit *const right) const
{
    const float distanceL((left->GetPositionVector() - m_neutrinoVertex).GetMagnitudeSquared());
    const float distanceR((right->GetPositionVector() - m_neutrinoVertex).GetMagnitudeSquared());
    return distanceL < distanceR;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
MvaPfoCharacterisationAlgorithm<T>::DistanceToVertexComparator::DistanceToVertexComparator(const CartesianVector vertexPosition2D_A, const CartesianVector vertexPosition2D_B,
                                                                                const HitType hitType_A, const HitType hitType_B) :
  m_neutrinoVertex_A(vertexPosition2D_A), m_neutrinoVertex_B(vertexPosition2D_B), m_hitType_A(hitType_A), m_hitType_B(hitType_B)
{
  // FIX ME: Throw error if the hit types passed coincide?                                                                                                                                              
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool MvaPfoCharacterisationAlgorithm<T>::DistanceToVertexComparator::operator()(const CaloHit *const left, const CaloHit *const right) const
{
  // FIX ME: Throw an error here if the hit type is different from any of the two possibilities?                                                                                                         
  float distanceL;
  float distanceR;
  if( left->GetHitType() == m_hitType_A )
    distanceL = ( left->GetPositionVector() - m_neutrinoVertex_A ).GetMagnitudeSquared();
  else
    distanceL = ( left->GetPositionVector() - m_neutrinoVertex_B ).GetMagnitudeSquared();
  if( right->GetHitType() == m_hitType_A )
    distanceR = ( right->GetPositionVector() - m_neutrinoVertex_A ).GetMagnitudeSquared();
  else
    distanceR = ( right->GetPositionVector() - m_neutrinoVertex_B ).GetMagnitudeSquared();
  return distanceL < distanceR;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode MvaPfoCharacterisationAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinPrimaryGoodHits", m_primaryParameters.m_minPrimaryGoodHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinHitsForGoodView", m_primaryParameters.m_minHitsForGoodView));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinPrimaryGoodViews", m_primaryParameters.m_minPrimaryGoodViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SelectInputHits", m_primaryParameters.m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinHitSharingFraction", m_primaryParameters.m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxPhotonPropagation", m_primaryParameters.m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "FoldToPrimaries", m_primaryParameters.m_foldBackHierarchy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PersistFeatures", m_persistFeatures));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingSetMode", m_trainingSetMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCaloHitsCut", m_minCaloHitsCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseThreeDInformation", m_useThreeDInformation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    // ATTN Support legacy XML configurations (note an order of precedence of XML keys exists)
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BdtFileName", m_mvaFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SvmFileName", m_mvaFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaFileName", m_mvaFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BdtName", m_mvaName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SvmName", m_mvaName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaName", m_mvaName));

    if (m_useThreeDInformation)
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
            XmlHelper::ReadValue(xmlHandle, "BdtFileNameNoChargeInfo", m_mvaFileNameNoChargeInfo));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
            XmlHelper::ReadValue(xmlHandle, "SvmFileNameNoChargeInfo", m_mvaFileNameNoChargeInfo));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
            XmlHelper::ReadValue(xmlHandle, "MvaFileNameNoChargeInfo", m_mvaFileNameNoChargeInfo));

        PANDORA_RETURN_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BdtNameNoChargeInfo", m_mvaNameNoChargeInfo));
        PANDORA_RETURN_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SvmNameNoChargeInfo", m_mvaNameNoChargeInfo));
        PANDORA_RETURN_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaNameNoChargeInfo", m_mvaNameNoChargeInfo));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnableProbability", m_enableProbability));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinProbabilityCut", m_minProbabilityCut));

    if (m_trainingSetMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TestBeamMode", m_testBeamMode));
        PANDORA_RETURN_RESULT_IF_AND_IF(
            STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ApplyFiducialCut", m_applyFiducialCut));
        if (m_applyFiducialCut)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutMinX", m_fiducialMinX));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutMaxX", m_fiducialMaxX));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutMinY", m_fiducialMinY));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutMaxY", m_fiducialMaxY));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutMinZ", m_fiducialMinZ));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FiducialCutMaxZ", m_fiducialMaxZ));
        }
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
            XmlHelper::ReadValue(xmlHandle, "ApplyReconstructabilityChecks", m_applyReconstructabilityChecks));
    }
    else
    {
        if (m_mvaFileName.empty() || m_mvaName.empty())
        {
            std::cout << "MvaPfoCharacterisationAlgorithm: MvaFileName and MvaName must be set if in classification mode " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }

        const std::string fullMvaFileName(LArFileHelper::FindFileInPath(m_mvaFileName, m_filePathEnvironmentVariable));
        m_mva.Initialize(fullMvaFileName, m_mvaName);

        if (m_useThreeDInformation)
        {
            if (m_mvaFileNameNoChargeInfo.empty() || m_mvaNameNoChargeInfo.empty())
            {
                std::cout << "MvaPfoCharacterisationAlgorithm: MvaFileNameNoChargeInfo and MvaNameNoChargeInfo must be set if in classification mode for no charge info in 3D mode "
                          << std::endl;
                return STATUS_CODE_INVALID_PARAMETER;
            }
            const std::string fullMvaFileNameNoChargeInfo(LArFileHelper::FindFileInPath(m_mvaFileNameNoChargeInfo, m_filePathEnvironmentVariable));
            m_mvaNoChargeInfo.Initialize(fullMvaFileNameNoChargeInfo, m_mvaNameNoChargeInfo);
        }
    }

    LArMvaHelper::AlgorithmToolMap algorithmToolMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        LArMvaHelper::ProcessAlgorithmToolListToMap(*this, xmlHandle, "FeatureTools", m_algorithmToolNames, algorithmToolMap));

    if (m_useThreeDInformation)
    {
        // and the map for NoChargeInfo
        LArMvaHelper::AlgorithmToolMap algorithmToolMapNoChargeInfo;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            LArMvaHelper::ProcessAlgorithmToolListToMap(
                *this, xmlHandle, "FeatureToolsNoChargeInfo", m_algorithmToolNamesNoChargeInfo, algorithmToolMapNoChargeInfo));

        for (auto const &[pAlgorithmToolName, pAlgorithmTool] : algorithmToolMap)
            PANDORA_RETURN_RESULT_IF(
                STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToMap(pAlgorithmTool, pAlgorithmToolName, m_featureToolMapThreeD));

        for (auto const &[pAlgorithmToolName, pAlgorithmTool] : algorithmToolMapNoChargeInfo)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                LArMvaHelper::AddFeatureToolToMap(pAlgorithmTool, pAlgorithmToolName, m_featureToolMapNoChargeInfo));
    }
    else
    {
        for (auto const &[pAlgorithmToolName, pAlgorithmTool] : algorithmToolMap)
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArMvaHelper::AddFeatureToolToMap(pAlgorithmTool, pAlgorithmToolName, m_featureToolMap));
    }

    return PfoCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool MvaPfoCharacterisationAlgorithm<T>::PassesFiducialCut(const CartesianVector &vertex) const
{
    const float vx(vertex.GetX()), vy(vertex.GetY()), vz(vertex.GetZ());
    return m_fiducialMinX <= vx && vx <= m_fiducialMaxX && m_fiducialMinY <= vy && vy <= m_fiducialMaxY && m_fiducialMinZ <= vz && vz <= m_fiducialMaxZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template class MvaPfoCharacterisationAlgorithm<AdaBoostDecisionTree>;
template class MvaPfoCharacterisationAlgorithm<SupportVectorMachine>;

} // namespace lar_content
