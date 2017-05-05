/**
 *  @file   larpandoracontent/LArTrackShowerId/SVMClusterCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArSVMHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArTrackShowerId/SVMClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.h"

using namespace pandora;

namespace lar_content
{

SVMClusterCharacterisationAlgorithm::SVMClusterCharacterisationAlgorithm() :
    m_produceSamplesMode(true),
    m_isPfoLevel(false),
    m_postBranchAddition(false),
    m_overwriteExistingId(false),
    m_useUnavailableClusters(false),
   // m_minTrackLikeViews(2),
    m_minCaloHitsCut(6)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SVMClusterCharacterisationAlgorithm::Run()
{
    //TODO: learn to indent correctly in codelite... 
  //  bool isTrackLike(false);
    
    if (m_isPfoLevel)
    {
    
       for (const std::string &pfoListName : m_inputPfoListNames)
        {
        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

        if (!pPfoList || pPfoList->empty())
        {
          if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "SVMPfoCharacterisationAlgorithm: unable to find pfo list " << pfoListName << std::endl;

          continue;
        }
        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
          PandoraContentApi::ParticleFlowObject::Metadata pfoMetadata;  
         // unsigned int nTrackLikeViews(0);  
          ClusterList twoDClusterList;
          LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);

          for (const Cluster *const pCluster : twoDClusterList)
            {

              if (pCluster->GetNCaloHits() < m_minCaloHitsCut)
                continue;     

                ClusterFeatureInfoMap clusterFeatureInfoMap;
                this->PopulateClusterFeatureInfoMap(pCluster,clusterFeatureInfoMap);
                SupportVectorMachine::DoubleVector featureList = this->GenerateFeatureList(pCluster, clusterFeatureInfoMap);                                                                                                                                                 

        if(m_produceSamplesMode)
        {
            bool isTrueTrack(false);
            try
              {
                const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
                isTrueTrack = ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId())));
              }
            catch (StatusCodeException &)
              {
              }
              
            SVMHelper::ProduceTrainingExample(m_trainingOutputFile + ".txt", isTrueTrack, featureList);
        }     
       /* else    
        {
            if (this->IsTrackLike(pCluster))
                ++nTrackLikeViews;
        }*/

       }//each cluster
       
      /* if(!m_produceSamplesMode)
       {
       if (nTrackLikeViews >= m_minTrackLikeViews)
            isTrackLike = true;
       else
            isTrackLike = false;
        //TODO: convert into: this->Characterise()
        if(isTrackLike)
        {
            pfoMetadata.m_particleId = MU_MINUS;

            if (m_showerPfoListName == pfoListName)
                showersToTracks.push_back(pPfo);
        }
        else
        {
           pfoMetadata.m_particleId = E_MINUS;

            if (m_trackPfoListName == pfoListName)
                tracksToShowers.push_back(pPfo);  
        }
                   if (pPfo->GetParticleId() != pfoMetadata.m_particleId.Get())
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pPfo, pfoMetadata));

            if (!m_updateClusterIds)
                continue;

            ClusterList twoDClusterList;
            LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);

            for (const Cluster *const pCluster : twoDClusterList)
            {
                if (pCluster->GetParticleId() == pfoMetadata.m_particleId.Get())
                    continue;

                PandoraContentApi::Cluster::Metadata clusterMetadata;
                clusterMetadata.m_particleId = pfoMetadata.m_particleId.Get();
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, clusterMetadata));
            }
       }*/

        }//each pfo
       
   /* if(!m_produceSamplesMode)
    { 
        if (!tracksToShowers.empty())
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_trackPfoListName, m_showerPfoListName, tracksToShowers));

        if (!showersToTracks.empty())
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_showerPfoListName, m_trackPfoListName, showersToTracks));
    }*/
    } //pfo list
    }   // if pfo level 
    else
    {
    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList = NULL;
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "SVMClusterCharacterisationAlgorithm: unable to find cluster list " << clusterListName << std::endl;

            continue;
        }

        for (const Cluster *const pCluster : *pClusterList)
        {
            if (!m_overwriteExistingId && (UNKNOWN_PARTICLE_TYPE != pCluster->GetParticleId()))
                continue;

            if (!m_useUnavailableClusters && !PandoraContentApi::IsAvailable(*this, pCluster))
                continue;

     if (pCluster->GetNCaloHits() < m_minCaloHitsCut)
         continue;
  
 
            //PandoraContentApi::Cluster::Metadata metadata;
	    ClusterFeatureInfoMap clusterFeatureInfoMap;
	    this->PopulateClusterFeatureInfoMap(pCluster,clusterFeatureInfoMap);

	    SupportVectorMachine::DoubleVector featureList = this->GenerateFeatureList(pCluster, clusterFeatureInfoMap);
	

	    
	    bool isTrueTrack(false);  
	    try
	      { 
		const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster)); 
		isTrueTrack = ((PHOTON != pMCParticle->GetParticleId()) && (E_MINUS != std::abs(pMCParticle->GetParticleId()))); 
	      }
	    catch (StatusCodeException &)                                                                      
	      {	
	      }    
          
     SVMHelper::ProduceTrainingExample(m_trainingOutputFile + ".txt", isTrueTrack, featureList);
     
        }
    }
    
    }//else
    

    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------
void SVMClusterCharacterisationAlgorithm::PopulateClusterFeatureInfoMap(const pandora::Cluster *const pCluster, ClusterFeatureInfoMap &clusterFeatureInfoMap) const
{
    // TODO: Reorganize this, use CalculateFeatures, and deal with different numbers of features in the feature vector if some are dropped off the xml file
    
 //   const double hitType(static_cast<double>(LArClusterHelper::GetClusterHitType(pCluster)));
    const double showerFitWidth(SVMHelper::CalculateFeaturesOfType<TrackShowerIdFeatureTool::ShowerFitWidthFeatureTool>(m_featureToolVector, this, pCluster).at(0));
    //std::cout << " showerFitWidth = " << showerFitWidth << std::endl;
  //  const double nHits(SVMHelper::CalculateFeaturesOfType<TrackShowerIdFeatureTool::NHitsFeatureTool>(m_featureToolVector, this, pCluster).at(0));
    //std::cout << " nHits = " << nHits << std::endl;
    //const double nGoodHits(SVMHelper::CalculateFeaturesOfType<TrackShowerIdFeatureTool::NHitsFeatureTool>(m_featureToolVector, this, pCluster).at(1));
//std::cout << " nGoodHits = " << nGoodHits << std::endl;   
    const double showerFitGapLength(SVMHelper::CalculateFeaturesOfType<TrackShowerIdFeatureTool::ShowerFitGapLengthFeatureTool>(m_featureToolVector, this, pCluster).at(0));
//std::cout << " showerFitGapLength = " << showerFitGapLength << std::endl; 
 //  const double straightLineLength(SVMHelper::CalculateFeaturesOfType<TrackShowerIdFeatureTool::StraightLineLengthFeatureTool>(m_featureToolVector, this, pCluster).at(0));
//std::cout << " straightLineLength = " << straightLineLength << std::endl; 
   const double straightLineLengthLarge(SVMHelper::CalculateFeaturesOfType<TrackShowerIdFeatureTool::StraightLineLengthFeatureTool>(m_featureToolVector, this, pCluster).at(1));
//std::cout <<  " straightLineLengthLarge " << std::endl;
    const double diffWithStraigthLine(SVMHelper::CalculateFeaturesOfType<TrackShowerIdFeatureTool::StraightLineLengthFeatureTool>(m_featureToolVector, this, pCluster).at(2));
//std::cout << "diffWithStraigthLine = " << diffWithStraigthLine << std::endl;
    const double widthDirectionX(SVMHelper::CalculateFeaturesOfType<TrackShowerIdFeatureTool::StraightLineLengthFeatureTool>(m_featureToolVector, this, pCluster).at(3));
//std::cout << "widthDirectionX = " << widthDirectionX << std::endl;
  //  const double pointsOfContact(SVMHelper::CalculateFeaturesOfType<TrackShowerIdFeatureTool::PointsOfContactFeatureTool>(m_featureToolVector, this, pCluster).at(0));
//std::cout << "pointsOfContact = " << pointsOfContact << std::endl; 
  const double nNearbyClusters(SVMHelper::CalculateFeaturesOfType<TrackShowerIdFeatureTool::NNearbyClustersFeatureTool>(m_featureToolVector, this, pCluster).at(0));
//std::cout << "nNearbyClusters = " << nNearbyClusters << std::endl; 
  // const double mipEnergy(SVMHelper::CalculateFeaturesOfType<TrackShowerIdFeatureTool::MipEnergyFeatureTool>(m_featureToolVector, this, pCluster).at(0));
 //  std::cout <<  "mipEnergy = " << mipEnergy << std::endl;
    
    //ClusterFeatureInfo clusterFeatureInfo(/*hitType, nHits,*/ nGoodHits, /*straightLineLength,*/ straightLineLengthLarge, showerFitWidth, showerFitGapLength, 
                        //  diffWithStraigthLine, widthDirectionX, nNearbyClusters, mipEnergy);
                                                    
    ClusterFeatureInfo clusterFeatureInfo( straightLineLengthLarge, showerFitWidth/straightLineLengthLarge, showerFitGapLength/straightLineLengthLarge, 
                          diffWithStraigthLine/straightLineLengthLarge, widthDirectionX/straightLineLengthLarge, nNearbyClusters/straightLineLengthLarge);
                          
    clusterFeatureInfoMap.emplace(pCluster, clusterFeatureInfo);
    // TODO - need to have one calculate features and read when there are less features calculated if we dont use all the variables... 
}

//------------------------------------------------------------------------------------------------------------------------------------------
/*SupportVectorMachine::DoubleVector SVMClusterCharacterisationAlgorithm::GenerateFeatureList(const pandora::Cluster *const pCluster, ClusterFeatureInfoMap &clusterFeatureInfoMap)
{
    ClusterFeatureInfo clusterFeatureInfo(clusterFeatureInfoMap.at(pCluster));
    
    SupportVectorMachine::DoubleVector featureVector;
   // ATTN: temporary, first tried with 11 features, but has redundant information 
 
 //  featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_hitType));
 //  featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_nHits));
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_nGoodHits));
  //  featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_straightLineLength));
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_straightLineLengthLarge));  
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_showerFitWidth));
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_showerFitGapLength));
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_diffWithStraigthLine));
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_widthDirectionX));
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_nNearbyClusters));
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_mipEnergy));
    
    return featureVector;
}*/
//------------------------------------------------------------------------------------------------------------------------------------------
SupportVectorMachine::DoubleVector SVMClusterCharacterisationAlgorithm::GenerateFeatureList(const pandora::Cluster *const pCluster, ClusterFeatureInfoMap &clusterFeatureInfoMap)
{
    ClusterFeatureInfo clusterFeatureInfo(clusterFeatureInfoMap.at(pCluster));
    
    SupportVectorMachine::DoubleVector featureVector;
   // ATTN: temporary, first tried with 11 features, but has redundant information 
 
 //  featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_hitType));
 //  featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_nHits));
   // featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_nGoodHits));
  //  featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_straightLineLength));
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_straightLineLengthLarge));  
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_showerFitWidthRatio));
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_showerFitGapLengthRatio));
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_diffWithStraigthLineRatio));
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_widthDirectionXRatio));
    featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_nNearbyClustersRatio));
  //  featureVector.push_back(static_cast<double>(clusterFeatureInfo.m_mipEnergy));
    
    return featureVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SVMClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    //same algorithm used for cluster/pfo characterisation, provide information about which stage we are (different SVMs to be read)
      PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
    "IsPfoLevel", m_isPfoLevel));  
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PostBranchAddition", m_postBranchAddition));
  
    // for cluster/pfo characterisation 
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OverwriteExistingId", m_overwriteExistingId));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseUnavailableClusters", m_useUnavailableClusters));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsCut", m_minCaloHitsCut));
        
  //  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
    //    "MinTrackLikeViews", m_minTrackLikeViews));
    
    //SVM specific options
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ProduceSamplesMode", m_produceSamplesMode));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ParameterInputFile", m_parameterInputFile));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SVMName", m_svmName));
    
    //to read depending on options above
   if(m_isPfoLevel)
    {
        
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
            m_inputPfoListNames.push_back(m_trackPfoListName);

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
            m_inputPfoListNames.push_back(m_showerPfoListName);
   
    }
    else
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
            "InputClusterListNames", m_inputClusterListNames));
    }
    
    if(m_produceSamplesMode)
    {
            
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "TrainingOutputFileName", m_trainingOutputFile));
    }
	else
    {
        if (m_parameterInputFile.empty() || m_svmName.empty())
        {
            std::cout << "SVMClusterCharacterisationAlgorithm: ParameterInputFile and SVMName must be set if in classification mode " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
        
        m_svMachine.Initialize(m_parameterInputFile, m_svmName);
    }								    
   
    //Read tool vector
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureTools", algorithmToolVector));
    
    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SVMHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVector));

    std::cout << " End of Read Settings" << std::endl;
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
