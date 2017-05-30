/**
 *  @file   larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.cc
 * 
 *  @brief  Implementation of the track shower id feature fool class
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h" 

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/BranchGrowingAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/ClusterCharacterisationAlgorithm.h"


using namespace pandora;

namespace lar_content
{
  
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------ 

TrackShowerIdFeatureTool::ShowerFitFeatureTool::ShowerFitFeatureTool() :
    m_slidingShowerFitWindow(3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void TrackShowerIdFeatureTool::ShowerFitFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster)
{
  if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
    std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
  
  const float showerFitWidth(this->CalculateShowerFitWidth(pAlgorithm, pCluster));
  
  featureVector.push_back(showerFitWidth);
  
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

float TrackShowerIdFeatureTool::ShowerFitFeatureTool::CalculateShowerFitWidth(const SVMClusterCharacterisationAlgorithm *const pAlgorithm, const pandora::Cluster * const pCluster) const
{
 //this lives in ClusterCharacterisationAlgorithm, could be moved here instead
    float showerFitWidth(-1.f);
   
    showerFitWidth = ClusterCharacterisationAlgorithm::GetShowerFitWidth(pAlgorithm, pCluster, m_slidingShowerFitWindow);  
   //TODO: try fit again if it fails with a small fit window
    
    return showerFitWidth;
  
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackShowerIdFeatureTool::ShowerFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
  
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
						       "SlidingShowerFitWindow", m_slidingShowerFitWindow));
  
  return STATUS_CODE_SUCCESS;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------ 
  
 TrackShowerIdFeatureTool::NHitsFeatureTool::NHitsFeatureTool()
{
}  

//------------------------------------------------------------------------------------------------------------------------------------------
 
void TrackShowerIdFeatureTool::NHitsFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm *const pAlgorithm,
    const pandora::Cluster * const pCluster)
{
  if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
    std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
  
  int nHits(this->CalculateNHits(pCluster));
  int nGoodHits(this->CalculateNGoodHits(pCluster, pAlgorithm));
  float nHitsF(nHits), nGoodHitsF(nGoodHits);
  featureVector.push_back(nHitsF); 
  featureVector.push_back(nGoodHitsF); 
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

int TrackShowerIdFeatureTool::NHitsFeatureTool::CalculateNHits(const pandora::Cluster * const pCluster) const
{
  return pCluster->GetNCaloHits();
}

//------------------------------------------------------------------------------------------------------------------------------------------

int TrackShowerIdFeatureTool::NHitsFeatureTool::CalculateNGoodHits(const pandora::Cluster * const pCluster, const SVMClusterCharacterisationAlgorithm *const pAlgorithm) const
   // TODO: NGoodHits could be used for true track/shower definition only, and it is now available in MCHelper - so this needs to change                   
{
	const MCParticleList *pMCParticleList = nullptr;
	PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, m_mcParticleListName, pMCParticleList));
	// Obtain map: [mc particle -> primary mc particle]                                  
	LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
	LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);
  
	//first get the list of calo hits
    CaloHitList caloHitList;
    pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

    int nGoodHits(0);
    float m_minHitSharingFraction(0.9f);
    
    //then proceed as in eventvalidation algorithm
    for (const CaloHit *const pCaloHit : caloHitList)
      {
        MCParticleVector mcParticleVector;
        for (const auto &mapEntry : pCaloHit->GetMCParticleWeightMap()) mcParticleVector.push_back(mapEntry.first);
		std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

        MCParticleWeightMap primaryWeightMap;

        for (const MCParticle *const pMCParticle : mcParticleVector)
		{
			const float weight(pCaloHit->GetMCParticleWeightMap().at(pMCParticle));
			LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToPrimaryMCMap.find(pMCParticle);

            if (mcToPrimaryMCMap.end() != mcIter)
				primaryWeightMap[mcIter->second] += weight;
		}

		MCParticleVector mcPrimaryVector;
		for (const auto &mapEntry : primaryWeightMap) mcPrimaryVector.push_back(mapEntry.first);
		std::sort(mcPrimaryVector.begin(), mcPrimaryVector.end(), PointerLessThan<MCParticle>());

		const MCParticle *pBestPrimaryParticle(nullptr);
		float bestPrimaryWeight(0.f), primaryWeightSum(0.f);
	
		for (const MCParticle *const pPrimaryMCParticle : mcPrimaryVector)
		{
			const float primaryWeight(primaryWeightMap.at(pPrimaryMCParticle));
			primaryWeightSum += primaryWeight;

			if (primaryWeight > bestPrimaryWeight)
			{
				bestPrimaryWeight = primaryWeight;
				pBestPrimaryParticle = pPrimaryMCParticle;
			}
		}

		if (!pBestPrimaryParticle || (primaryWeightSum < std::numeric_limits<float>::epsilon()) || ((bestPrimaryWeight / primaryWeightSum) < m_minHitSharingFraction))
			continue;

			nGoodHits++;
	}

	return nGoodHits;
        
}   

//------------------------------------------------------------------------------------------------------------------------------------------          
   
StatusCode TrackShowerIdFeatureTool::NHitsFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
  
   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
  
  return STATUS_CODE_SUCCESS;
}
 
  
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------ 
    
TrackShowerIdFeatureTool::ShowerFitGapLengthFeatureTool::ShowerFitGapLengthFeatureTool() :
    m_slidingShowerFitWindow(3)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void TrackShowerIdFeatureTool::ShowerFitGapLengthFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm *const pAlgorithm, 
           const pandora::Cluster * const pCluster)
{
  if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
    std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
  
  float showerFitGapLength(this->CalculateShowerFitGapLength(pCluster));
  featureVector.push_back(showerFitGapLength); 
 
} 
   
//------------------------------------------------------------------------------------------------------------------------------------------

float TrackShowerIdFeatureTool::ShowerFitGapLengthFeatureTool::CalculateShowerFitGapLength(const pandora::Cluster * const pCluster) const
{
	float showerFitGapLength(-1.f);
	CartesianVector globalMinLayerPositionOnAxis(0.f, 0.f, 0.f), globalMaxLayerPositionOnAxis(0.f, 0.f, 0.f);                                   
    try                                                                                                                                                      
	{                                                                                                                                                    
		const TwoDSlidingShowerFitResult showerFitResult(pCluster, m_slidingShowerFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));                            
        const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());                                            
        if (layerFitResultMapS.size() > 1)                                                                                                              
		{                                                                                                                                           
  
			showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.begin()->second.GetL(), 0.f, globalMinLayerPositionOnAxis);       
            showerFitResult.GetShowerFitResult().GetGlobalPosition(layerFitResultMapS.rbegin()->second.GetL(), 0.f, globalMaxLayerPositionOnAxis);       

            float straightLinePathLength=((globalMaxLayerPositionOnAxis - globalMinLayerPositionOnAxis).GetMagnitude());                                 
            if (straightLinePathLength > std::numeric_limits<float>::epsilon())                                                                           
			{                                                                                                                                            
                showerFitGapLength = 0.f;                                                                                                                    
                CartesianVector previousLayerPosition(globalMinLayerPositionOnAxis);                                                                         
                                                                                                                                                       
                for (LayerFitResultMap::const_iterator iterS = layerFitResultMapS.begin(); iterS != layerFitResultMapS.end(); ++iterS)                      
				{                                                                                                                                        
					CartesianVector thisLayerPosition(0.f, 0.f, 0.f);                                                                                    
                    showerFitResult.GetShowerFitResult().GetGlobalPosition(iterS->second.GetL(), 0.f, thisLayerPosition);     
                
                    const float thisGapLength((thisLayerPosition - previousLayerPosition).GetMagnitude());                                              

					if (thisGapLength > showerFitGapLength)                                                                                          
					{                                                                                                                                
                        const float minZ(std::min(thisLayerPosition.GetZ(), previousLayerPosition.GetZ()));                                             
                        const float maxZ(std::max(thisLayerPosition.GetZ(), previousLayerPosition.GetZ()));                                            
                        if ((maxZ - minZ) < std::numeric_limits<float>::epsilon())                                                                        
                          throw StatusCodeException(STATUS_CODE_FAILURE);                                                                                
                        const float gapZ(LArGeometryHelper::CalculateGapDeltaZ(this->GetPandora(), minZ, maxZ, LArClusterHelper::GetClusterHitType(pCluster)));            
                        const float correctedGapLength(thisGapLength * (1.f - gapZ / (maxZ - minZ)));                                                        
                        if (correctedGapLength > showerFitGapLength)                                                                                     
                          showerFitGapLength = correctedGapLength;                                                                                        
					}                                                                                                                                 
				}                                                                                                                                     
			}                                                                                                                                         
                                                                                                                                                         
		}
           
	}
    catch (StatusCodeException &)                                                                                                                           
      {                                                                                                                                                   
      }                         
  
  return showerFitGapLength;
      
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackShowerIdFeatureTool::ShowerFitGapLengthFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
  
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
						       "SlidingShowerFitWindow", m_slidingShowerFitWindow));
  
  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

    TrackShowerIdFeatureTool::LinearFitFeatureTool::LinearFitFeatureTool() :
    m_slidingLinearFitWindow(3),
    m_addDiffWithStraightLine(true),
    m_adddTdLWidth(true),
    m_addMaxFitGapLength(true),
    m_addRMSLinearFit(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------   
 
void TrackShowerIdFeatureTool::LinearFitFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm *const pAlgorithm, 
const pandora::Cluster * const pCluster) 
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
  
  
    float dTdLWidth(-1.f), straightLineLengthLarge(-1.f), diffWithStraightLine(0.f), diffWithStraightLineMean(0.f),diffWithStraightLineSigma(0.f) maxFitGapLength(-std::numeric_limits<float>::max()), rmsSlidingLinearFit(0.f);
      
    this->CalculateVariablesSlidingLinearFit(pCluster,straightLineLengthLarge,diffWithStraightLine,diffWithStraightLineMean,diffWithStraightLineSigma,dTdLWidth,maxFitGapLength, rmsSlidingLinearFit);
					     
	//TODO: try with a larget fit window if fails?
  
    featureVector.push_back(straightLineLengthLarge); 
    if (m_addDiffWithStraightLine)
      {
	featureVector.push_back(diffWithStraightLineMean);                                                       
	featureVector.push_back(diffWithStraightLineSigma);
	//        featureVector.push_back(diffWithStraightLine);
      }
    if (m_adddTdLWidth)
        featureVector.push_back(dTdLWidth);
    if (m_addMaxFitGapLength)
        featureVector.push_back(maxFitGapLength);
    if (m_addRMSLinearFit)                                                                                           
      featureVector.push_back(rmsSlidingLinearFit); 
    
} 

//------------------------------------------------------------------------------------------------------------------------------------------
void TrackShowerIdFeatureTool::LinearFitFeatureTool::CalculateVariablesSlidingLinearFit(const pandora::Cluster * const pCluster, float &straightLineLengthLarge, float &diffWithStraightLine, float &diffWithStraightLineMean, float &diffWithStraightLineSigma, float &dTdLWidth, float &maxFitGapLength, float &rmsSlidingLinearFit) const
{
    
       float dTdLMin(+std::numeric_limits<float>::max()), dTdLMax(-std::numeric_limits<float>::max());                         
       FloatVector dTdLVector, diffWithStraightLineVector; 
	CartesianVector globalMinLayerPositionLarge(0.f,0.f,0.f), globalMaxLayerPositionLarge(0.f,0.f,0.f);

	try
	{
	  diffWithStraightLine = 0.f;                                                                              
	  rmsSlidingLinearFit = 0.f;
		const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult slidingFitResultLarge(pCluster, 10000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        
        straightLineLengthLarge = (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
                                                                                                   
        CartesianVector previousFitPosition(slidingFitResult.GetGlobalMinLayerPosition());                                                                                        
        const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());
        const LayerFitResultMap &layerFitResultMapLarge(slidingFitResultLarge.GetLayerFitResultMap());

        for (const auto &mapEntry : layerFitResultMap)                                                                                                      
        {                                                                                                                                                   
            CartesianVector thisFitPosition(0.f, 0.f, 0.f), thisFitDirection(0.f, 0.f, 0.f), thisFitPositionLarge(0.f, 0.f, 0.f);//LORENA
    
            slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);  
            slidingFitResultLarge.GetGlobalFitPosition(mapEntry.second.GetL(), thisFitPositionLarge);   
            LayerFitResultMap::const_iterator iterLarge = layerFitResultMapLarge.find(slidingFitResultLarge.GetLayer(mapEntry.second.GetL()));
            float diff(std::abs(mapEntry.second.GetFitT()-iterLarge->second.GetFitT()));  
            
            diffWithStraightLine += diff; 
	    diffWithStraightLineVector.push_back(diff);                                                      
	    rmsSlidingLinearFit += mapEntry.second.GetRms();

            const float thisGapLength((thisFitPosition - previousFitPosition).GetMagnitude());
  
            const float minZ(std::min(thisFitPosition.GetZ(), previousFitPosition.GetZ()));                                             
            const float maxZ(std::max(thisFitPosition.GetZ(), previousFitPosition.GetZ())); 
            if ((maxZ - minZ) > std::numeric_limits<float>::epsilon())                                                                        
            {
                const float gapZ(LArGeometryHelper::CalculateGapDeltaZ(this->GetPandora(), minZ, maxZ, LArClusterHelper::GetClusterHitType(pCluster)));						
                const float correctedGapLength(thisGapLength * (1.f - gapZ / (maxZ - minZ))); //check this
                                                                                                        
                if (correctedGapLength > maxFitGapLength)                                                                                     
                    maxFitGapLength = correctedGapLength;     
            }                                                                                              
                                                                                                          
            dTdLMin = std::min(dTdLMin, static_cast<float>(mapEntry.second.GetGradient()));
            dTdLMax = std::max(dTdLMax, static_cast<float>(mapEntry.second.GetGradient()));
	    dTdLVector.push_back(static_cast<float>(mapEntry.second.GetGradient()));

            previousFitPosition = thisFitPosition;     
                                                                        
        }               
        dTdLWidth = dTdLMax - dTdLMin;
	for (unsigned i=0; i < diffWithStraightLineVector.size(); i++)                                           
	  diffWithStraightLineMean += diffWithStraightLineVector[i];                                       
	diffWithStraightLineMean = std::abs(diffWithStraightLineMean/diffWithStraightLineVector.size());         
	for (unsigned i=0; i < diffWithStraightLineVector.size(); i++)                                           
	  diffWithStraightLineSigma += (diffWithStraightLineVector[i]-diffWithStraightLineMean)*(diffWithStraightLineVector[i]-diffWithStraightLineMean);                                                                                  
	diffWithStraightLineSigma = std::sqrt(diffWithStraightLineSigma/diffWithStraightLineVector.size()); 
    }
    catch (const StatusCodeException &)
    {
    }
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackShowerIdFeatureTool::LinearFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
  
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "SlidingLinearFitWindow", m_slidingLinearFitWindow));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "AddDiffWithStraightLine", m_addDiffWithStraightLine));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "AddDTdLWidth", m_adddTdLWidth));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "AddMaxFitGapLength", m_addMaxFitGapLength));      
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,        
		"AddRMSLinearFit", m_addRMSLinearFit));        
  
  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
   
 TrackShowerIdFeatureTool::NNearbyClustersFeatureTool::NNearbyClustersFeatureTool() :
    m_nearbyClusterDistance(2.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------   

void TrackShowerIdFeatureTool::NNearbyClustersFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm *const pAlgorithm, 
const pandora::Cluster * const pCluster) 
{
	if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
		std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
        
  
    int pointsOfContact(this->CalculatePointsOfContact(pCluster,pAlgorithm));
  
    featureVector.push_back(pointsOfContact); 
    
} 

//------------------------------------------------------------------------------------------------------------------------------------------

int TrackShowerIdFeatureTool::NNearbyClustersFeatureTool::CalculatePointsOfContact(const pandora::Cluster * const pCluster, const SVMClusterCharacterisationAlgorithm *const pAlgorithm) const
{
	pandora::ClusterVector candidateClusters;  

	for (const std::string &clusterListName : m_clusterListNames) 
	{      
		const pandora::ClusterList *pClusterList = nullptr;
		PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, clusterListName, pClusterList));
                                                                                                          
		for (const Cluster *const pCandidateCluster : *pClusterList)                                                                                       
		{ 
			if ((LArClusterHelper::GetClusterHitType(pCluster))!=(LArClusterHelper::GetClusterHitType(pCandidateCluster)))
				continue;                                                                                                                                                 
			if ((pCandidateCluster != pCluster) && (pCandidateCluster->GetNCaloHits() > 5)) 
				candidateClusters.push_back(pCandidateCluster);                                                                                                     
		}   
	}
    int nPointsOfContact(0);
    
    for (ClusterVector::iterator iterI = candidateClusters.begin(), iterIEnd = candidateClusters.end(); iterI != iterIEnd; ++iterI)    
    {
		const Cluster *const pCandidateCluster = *iterI;                                                                                                 
		if (NULL == pCandidateCluster)                                                                                                                  
			continue;
        
		if (this->IsClusterNearby(pCluster,pCandidateCluster))
			++nPointsOfContact;
    }
 return nPointsOfContact;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackShowerIdFeatureTool::NNearbyClustersFeatureTool::IsClusterNearby(const pandora::Cluster *const pCluster,const pandora::Cluster *const pCandidateCluster) const
{

	float closestDistance(LArClusterHelper::GetClosestDistance(pCluster, pCandidateCluster));
	if (closestDistance<=m_nearbyClusterDistance)
		return true;
            
	return false;           
}

//------------------------------------------------------------------------------------------------------------------------------------------   
       
StatusCode TrackShowerIdFeatureTool::NNearbyClustersFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
  
   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));
   PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "NearbyClusterDistance", m_nearbyClusterDistance));

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------   

 TrackShowerIdFeatureTool::MipEnergyFeatureTool::MipEnergyFeatureTool() :
    m_mipCorrectionPerHit(1.0f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------   

void TrackShowerIdFeatureTool::MipEnergyFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm *const pAlgorithm, 
           const pandora::Cluster * const pCluster) 
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
        
  
    float mipEnergy(this->CalculateMipEnergy(pCluster));
  
    featureVector.push_back(mipEnergy); 
    
} 

//------------------------------------------------------------------------------------------------------------------------------------------

float TrackShowerIdFeatureTool::MipEnergyFeatureTool::CalculateMipEnergy(const pandora::Cluster * const pCluster) const
{

	float mipEnergy(0.f);
 
	CaloHitList clusterHitList; 
	pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterHitList);                                                           

	for (const CaloHit *const pCaloHit : clusterHitList)                                                                                                                                                                                                 
		mipEnergy += pCaloHit->GetMipEquivalentEnergy()*m_mipCorrectionPerHit;     

	return mipEnergy;                                                               
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackShowerIdFeatureTool::MipEnergyFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "MipCorrectionPerHit", m_mipCorrectionPerHit));  
    return STATUS_CODE_SUCCESS;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------   
 
 TrackShowerIdFeatureTool::VertexDistanceFeatureTool::VertexDistanceFeatureTool() :
 m_addVertexDistance(true)
{
}
//------------------------------------------------------------------------------------------------------------------------------------------   

void TrackShowerIdFeatureTool::VertexDistanceFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm *const pAlgorithm, 
const pandora::Cluster * const pCluster) 
{
	if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
        
  
    float vertexDistance(this->CalculateVertexDistance(pAlgorithm, pCluster));
  
	if (m_addVertexDistance)
		featureVector.push_back(vertexDistance); 
} 

//------------------------------------------------------------------------------------------------------------------------------------------ 

float TrackShowerIdFeatureTool::VertexDistanceFeatureTool::CalculateVertexDistance(const SVMClusterCharacterisationAlgorithm *const pAlgorithm, const pandora::Cluster * const pCluster) const
{

	float vertexDistance(ClusterCharacterisationAlgorithm::GetVertexDistance(pAlgorithm, pCluster));
 
	return vertexDistance;                                                               
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackShowerIdFeatureTool::VertexDistanceFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
  
  //why, why, why???
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "AddVertexDistance", m_addVertexDistance));  
  return STATUS_CODE_SUCCESS;
}


} // namespace lar_content


