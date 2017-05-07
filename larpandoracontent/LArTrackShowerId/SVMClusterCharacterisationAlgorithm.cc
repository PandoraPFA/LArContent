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
    m_produceSamplesMode(false),
    m_isPfoLevel(false),
    m_postBranchAddition(false),
    m_ratioVariables(true),   
    m_overwriteExistingId(false),
    m_updateClusterIds(true),
    m_useUnavailableClusters(false),
    m_minTrackLikeViews(2),
    m_minCaloHitsCut(6)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SVMClusterCharacterisationAlgorithm::Run()
{
    //TODO: learn to indent correctly in codelite... 
    //TODO: factorize this below...
    if (m_isPfoLevel)
    {
        
        PfoList tracksToShowers, showersToTracks;
    
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
            unsigned int nTrackLikeViews(0);  
          ClusterList twoDClusterList;
          LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);

          for (const Cluster *const pCluster : twoDClusterList)
            {

              if (pCluster->GetNCaloHits() < m_minCaloHitsCut)
                continue;     
                
            SupportVectorMachine::DoubleVector featureList = this->GenerateFeatureVector(pCluster);                                                                                                                                              

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
        else    
        {
            if (this->IsClearTrack(featureList))
                ++nTrackLikeViews;
        }

       }//each cluster in a pfo
       
       if(!m_produceSamplesMode)
       {
           
            bool isTrackLike(false);
            if (nTrackLikeViews >= m_minTrackLikeViews)
                isTrackLike = true;
            else
                isTrackLike = false;
        //TODO: convert below into: this->Characterise(pfo)
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

            for (const Cluster *const pCluster : twoDClusterList)
            {
                if (pCluster->GetParticleId() == pfoMetadata.m_particleId.Get())
                    continue;

                PandoraContentApi::Cluster::Metadata clusterMetadata;
                clusterMetadata.m_particleId = pfoMetadata.m_particleId.Get();
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, clusterMetadata));
            }
       }

        }//each pfo
    
    } //pfo list
    
       
    if(!m_produceSamplesMode)
    { 
        if (!tracksToShowers.empty())
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_trackPfoListName, m_showerPfoListName, tracksToShowers));

        if (!showersToTracks.empty())
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_showerPfoListName, m_trackPfoListName, showersToTracks));
    }
    
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
    
        SupportVectorMachine::DoubleVector featureList = this->GenerateFeatureVector(pCluster);

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
        else
        {
            PandoraContentApi::Cluster::Metadata metadata;
            if (this->IsClearTrack(featureList))
            {
                metadata.m_particleId = MU_MINUS;
            }
            else
            {
                metadata.m_particleId = E_MINUS;
            }

	    if (pCluster->GetParticleId() != metadata.m_particleId.Get())
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));
        }
        }
    }
    
    }//else
    

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
SupportVectorMachine::DoubleVector SVMClusterCharacterisationAlgorithm::GenerateFeatureVector(const pandora::Cluster *const pCluster) const
{
    SupportVectorMachine::DoubleVector featureVector;
  
    featureVector = SVMHelper::CalculateFeatures(m_featureToolVector, this, pCluster);
    
    if(m_ratioVariables)
    {
        SupportVectorMachine::DoubleVector processedFeatureVector(this->ProcessFeatureVector(featureVector));
        return processedFeatureVector;
    }
    return featureVector;
  
}

//------------------------------------------------------------------------------------------------------------------------------------------
SupportVectorMachine::DoubleVector SVMClusterCharacterisationAlgorithm::ProcessFeatureVector(const SupportVectorMachine::DoubleVector &featureVector) const
{
    //This is in case we need to do any work on the variables (like doing ratios)...
    //TODO check that the linear tool has been included - also, it should be the first one, think about how to check that or make it independent...
    SupportVectorMachine::DoubleVector processedFeatureVector;
    double straightLineLength = featureVector.at(0);
    processedFeatureVector.push_back(straightLineLength);
    
    for (unsigned int i=1; i<featureVector.size(); ++i)
    {
        double variable(featureVector.at(i));
        processedFeatureVector.push_back(variable/straightLineLength);
    }
        return processedFeatureVector;
    
}
//------------------------------------------------------------------------------------------------------------------------------------------
bool SVMClusterCharacterisationAlgorithm::IsClearTrack(const SupportVectorMachine::DoubleVector &featureVector) const
{
    //Track is defined as true (signal) in the training, so nothing to change, but potentially could give more instructions here for classification
    return SVMHelper::Classify(m_svMachine, featureVector);
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
        "UpdateClusterIds", m_updateClusterIds));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseUnavailableClusters", m_useUnavailableClusters));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHitsCut", m_minCaloHitsCut));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinTrackLikeViews", m_minTrackLikeViews));
    
    //SVM specific options
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ProduceSamplesMode", m_produceSamplesMode));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SVMFileName", m_svmFileName));
        
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SVMName", m_svmName));
        
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "RatioVariables", m_ratioVariables));  
    
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
        if (m_svmFileName.empty() || m_svmName.empty())
        {
            std::cout << "SVMClusterCharacterisationAlgorithm: SVMFileName and SVMName must be set if in classification mode " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
        
        m_svMachine.Initialize(m_svmFileName, m_svmName);
    }								    
   
    //Read tool vector
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "FeatureTools", algorithmToolVector));
    
    for (AlgorithmTool *const pAlgorithmTool : algorithmToolVector)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, SVMHelper::AddFeatureToolToVector(pAlgorithmTool, m_featureToolVector));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
