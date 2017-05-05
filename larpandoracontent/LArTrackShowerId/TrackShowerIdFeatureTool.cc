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
  
// ***** ShowerFitWidthFeatureTool ***** 
TrackShowerIdFeatureTool::ShowerFitWidthFeatureTool::ShowerFitWidthFeatureTool() :
    m_slidingShowerFitWindow(3)
{
}
//------------------------------------------------------------------------------------------------------------------------------------------
    
void TrackShowerIdFeatureTool::ShowerFitWidthFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster)
{
  if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
    std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
  
  const float showerFitWidth(this->CalculateShowerFitWidth(pCluster));
  
  featureVector.push_back(showerFitWidth);
  
}
  
  
//------------------------------------------------------------------------------------------------------------------------------------------

float TrackShowerIdFeatureTool::ShowerFitWidthFeatureTool::CalculateShowerFitWidth(const pandora::Cluster * const pCluster) const
{
  try
    {
      const TwoDSlidingShowerFitResult showerFitResult(pCluster, m_slidingShowerFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
      const LayerFitResultMap &layerFitResultMapS(showerFitResult.GetShowerFitResult().GetLayerFitResultMap());
      const LayerFitResultMap &layerFitResultMapP(showerFitResult.GetPositiveEdgeFitResult().GetLayerFitResultMap());
      const LayerFitResultMap &layerFitResultMapN(showerFitResult.GetNegativeEdgeFitResult().GetLayerFitResultMap());
      
      if (!layerFitResultMapS.empty())
	{
	  float showerFitWidth(0.f);
	  
	  for (const auto &mapEntryS : layerFitResultMapS)
            {
	      
	      LayerFitResultMap::const_iterator iterP = layerFitResultMapP.find(mapEntryS.first);
	      LayerFitResultMap::const_iterator iterN = layerFitResultMapN.find(mapEntryS.first);
	      
	      if ((layerFitResultMapP.end() != iterP) && (layerFitResultMapN.end() != iterN))
		showerFitWidth += std::fabs(iterP->second.GetFitT() - iterN->second.GetFitT());
            } 
	  
          return showerFitWidth;
	}
    }
  catch (const StatusCodeException &)
    {
    }
  
  return -1.f;
  
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackShowerIdFeatureTool::ShowerFitWidthFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
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
  //check this
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
   // Obtain map: [mc particle -> primary mc particle]                    
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
    TrackShowerIdFeatureTool::StraightLineLengthFeatureTool::StraightLineLengthFeatureTool() :
    m_slidingLinearFitWindow(3),
    m_addDiffWithStraightLine(true),
    m_addWidthDirectionX(true)
{
}

 //------------------------------------------------------------------------------------------------------------------------------------------   
void TrackShowerIdFeatureTool::StraightLineLengthFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm *const pAlgorithm, 
           const pandora::Cluster * const pCluster) 
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
  
  
    float straightLineLength(0.f), straightLineLengthLarge(0.f), diffWithStraigthLine(0.f), widthDirectionX(0.f);
    this->CalculateVariablesSlidingFit(pCluster,straightLineLength, straightLineLengthLarge, diffWithStraigthLine, widthDirectionX);
  
    featureVector.push_back(straightLineLength); 
    featureVector.push_back(straightLineLengthLarge); 
    if(m_addDiffWithStraightLine)
        featureVector.push_back(diffWithStraigthLine);
    if(m_addWidthDirectionX)
        featureVector.push_back(widthDirectionX);
    
} 

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackShowerIdFeatureTool::StraightLineLengthFeatureTool::CalculateVariablesSlidingFit(const pandora::Cluster * const pCluster, float &straightLineLength, 
float &straightLineLengthLarge, float &diffWithStraigthLine, float &widthDirectionX) const
    {
                                        
  float minXDir(+std::numeric_limits<float>::max()), maxXDir(-std::numeric_limits<float>::max());                           
                                                                                                                
  CartesianVector globalMinLayerPosition(0.f,0.f,0.f), globalMaxLayerPosition(0.f,0.f,0.f), globalMinLayerPositionLarge(0.f,0.f,0.f), globalMaxLayerPositionLarge(0.f,0.f,0.f);

  try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        globalMinLayerPosition=(slidingFitResult.GetGlobalMinLayerPosition());                        
        globalMaxLayerPosition=(slidingFitResult.GetGlobalMaxLayerPosition());
        straightLineLength = (slidingFitResult.GetGlobalMaxLayerPosition() - slidingFitResult.GetGlobalMinLayerPosition()).GetMagnitude();

        //the large one - in principle the detector size should be enough
        const TwoDSlidingFitResult slidingFitResultLarge(pCluster, 10000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
	globalMinLayerPositionLarge=(slidingFitResultLarge.GetGlobalMinLayerPosition());                        
	globalMaxLayerPositionLarge=(slidingFitResultLarge.GetGlobalMaxLayerPosition());
	straightLineLengthLarge = (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();

    diffWithStraigthLine = 0.f;                                                                                                                    
	CartesianVector previousFitPosition(globalMinLayerPosition);                                                                                        
	const LayerFitResultMap &layerFitResultMap(slidingFitResult.GetLayerFitResultMap());
    const LayerFitResultMap &layerFitResultMapLarge(slidingFitResultLarge.GetLayerFitResultMap());
	
        for (const auto &mapEntry : layerFitResultMap)                                                                                                      
	  {                                                                                                                                                   
	    CartesianVector thisFitPosition(0.f, 0.f, 0.f), thisFitDirection(0.f, 0.f, 0.f), thisFitPositionLarge(0.f, 0.f, 0.f);//LORENA
	    
	    slidingFitResult.GetGlobalPosition(mapEntry.second.GetL(), mapEntry.second.GetFitT(), thisFitPosition);  
        slidingFitResultLarge.GetGlobalFitPosition(mapEntry.second.GetL(), thisFitPositionLarge);   
	    LayerFitResultMap::const_iterator iterLarge = layerFitResultMapLarge.find(slidingFitResultLarge.GetLayer(mapEntry.second.GetL()));
            float diff(std::abs(mapEntry.second.GetFitT()-iterLarge->second.GetFitT()));  
            diffWithStraigthLine += diff; 
            
                                                        
	    previousFitPosition = thisFitPosition;                                                                                                          
                                                                                                          
	    slidingFitResult.GetGlobalFitDirection(mapEntry.second.GetL(), thisFitDirection);                                   
	    if(thisFitDirection.GetX()<minXDir)                                                                                 
	      minXDir=thisFitDirection.GetX();                                                                                  
	    if(thisFitDirection.GetX()>maxXDir)                                                                                 
	      maxXDir=thisFitDirection.GetX();                                                                                  
	  }               

    }
    catch (const StatusCodeException &)
    {
    }
    
      widthDirectionX=maxXDir-minXDir; 
    }

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackShowerIdFeatureTool::StraightLineLengthFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
  
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "SlidingLinearFitWindow", m_slidingLinearFitWindow));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "AddDiffwithStraigthLine", m_addDiffWithStraightLine));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "AddWidthDirectionX", m_addWidthDirectionX));              
  
  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
 //------------------------------------------------------------------------------------------------------------------------------------------   
 TrackShowerIdFeatureTool::PointsOfContactFeatureTool::PointsOfContactFeatureTool() :
     m_nearbyClusterDistance(2.5f),                                                                                                                           
    m_remoteClusterDistance(10.f),
      m_directionTanAngle(1.732f),
    m_directionApexShift(0.333f)
{
}
 //------------------------------------------------------------------------------------------------------------------------------------------   

void TrackShowerIdFeatureTool::PointsOfContactFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm *const pAlgorithm, 
           const pandora::Cluster * const pCluster) 
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
        
  
    int pointsOfContact(this->CalculatePointsOfContact(pCluster,pAlgorithm));
  
    featureVector.push_back(pointsOfContact); 
    
} 

//------------------------------------------------------------------------------------------------------------------------------------------
int TrackShowerIdFeatureTool::PointsOfContactFeatureTool::CalculatePointsOfContact(const pandora::Cluster * const pCluster, const SVMClusterCharacterisationAlgorithm *const pAlgorithm) const
{
    
    pandora::ClusterVector candidateClusters;  

 for (const std::string &clusterListName : m_clusterListNames) 
 {      
   const pandora::ClusterList *pClusterList = nullptr;
   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, clusterListName, pClusterList));
                                                                                                          
    for (const Cluster *const pCandidateCluster : *pClusterList)                                                                                       
      { 
        if((LArClusterHelper::GetClusterHitType(pCluster))!=(LArClusterHelper::GetClusterHitType(pCandidateCluster)))
            continue;                                                                                                                                                 
        if ((pCandidateCluster != pCluster) && (pCandidateCluster->GetNCaloHits() > 5)) //LORENA: check this
            candidateClusters.push_back(pCandidateCluster);                                                                                                     
      }   
 }
    int nPointsOfContact(0);/*, hitsInEnvelope(0);*/ // think if I want to add hitsInEnvelope...
    
    BranchGrowingAlgorithm::ClusterUsageMap forwardUsageMap, backwardUsageMap; 
	this->FindAssociatedClusters(pCluster, candidateClusters, forwardUsageMap, backwardUsageMap, pAlgorithm);                                                      
	BranchGrowingAlgorithm::SeedAssociationList seedAssociationList; 
	this->IdentifyClusterMerges(ClusterVector(1, pCluster), backwardUsageMap, seedAssociationList);                                                    
	
   // for (const pandora::Cluster *const pBranchCluster : seedAssociationList.at(pCluster))                                                                       
	 //   ++nPointsOfContact; 
     nPointsOfContact=(seedAssociationList.at(pCluster)).size(); // LORENA: check this 
 
/*   TODO    
    //and hits in envelope                                                                                   
	    CaloHitList branchHitList;                                       
	    pBranchCluster->GetOrderedCaloHitList().FillCaloHitList(branchHitList);
	    for (const CaloHit *const pCaloHit : branchHitList)                                                            
	      {
             CartesianVector hitPos(pCaloHit->GetPositionVector().GetX(), 0, pCaloHit->GetPositionVector().GetZ());
             double angle(axisExt.GetOpeningAngle(hitPos-startExt));
             if((std::abs(angle)<std::abs(angleEnvelope))&&(std::abs(pCaloHit->GetPositionVector().GetZ())<std::abs(endExt.GetZ()))&&(std::abs(pCaloHit->GetPositionVector().GetZ())>std::abs(startExt.GetZ())))
		    ++hitsInEnvelope;
          }
   */
 return nPointsOfContact;
}

//ATT: Two next functions are copied from BranchGrowingAlgorithm - use from there in future
//------------------------------------------------------------------------------------------------------------------------------------------
void TrackShowerIdFeatureTool::PointsOfContactFeatureTool::FindAssociatedClusters(const pandora::Cluster *const pParticleSeed, pandora::ClusterVector &candidateClusters, 
     BranchGrowingAlgorithm::ClusterUsageMap &forwardUsageMap, BranchGrowingAlgorithm::ClusterUsageMap &backwardUsageMap, const SVMClusterCharacterisationAlgorithm *const pAlgorithm) const 
{                                                                                                                                                            
  ClusterVector currentSeedAssociations, newSeedAssociations;                                                                                               
  currentSeedAssociations.push_back(pParticleSeed);                                                                                                          
  ClusterList strongClusters, standardClusters, singleClusters;
  
  unsigned int associationOrder(1);                                                                                                                          
  while (!currentSeedAssociations.empty())                                                                                                              
    {	
      for (ClusterVector::iterator iterI = candidateClusters.begin(), iterIEnd = candidateClusters.end(); iterI != iterIEnd; ++iterI)                 
        {	
          const Cluster *const pCandidateCluster = *iterI;                                                                                                 
	  if (NULL == pCandidateCluster)                                                                                                                  
	    continue;      
                                                                 
	  for (ClusterVector::iterator iterJ = currentSeedAssociations.begin(), iterJEnd = currentSeedAssociations.end(); iterJ != iterJEnd; ++iterJ)  
            {	
              const Cluster *const pAssociatedCluster = *iterJ;                                                                                           
              const BranchGrowingAlgorithm::AssociationType associationType(this->AreClustersAssociated(pAssociatedCluster, pCandidateCluster, pAlgorithm)); 
              if (BranchGrowingAlgorithm::NONE == associationType)     
		           continue;

	      // Check we store best association between this seed and candidate                                                 
	      BranchGrowingAlgorithm::Association association(associationOrder, associationType); 
              const BranchGrowingAlgorithm::Association &existingAssociation = forwardUsageMap[pParticleSeed][pCandidateCluster];
	      
	      if (association.GetType() > existingAssociation.GetType())                                                                                       
                {  
                  // If not first association, check strength of previous association in chain                                                         
                  if (pParticleSeed != pAssociatedCluster)                                                                                             
                    association.SetType(std::min(association.GetType(), backwardUsageMap[pAssociatedCluster][pParticleSeed].GetType()));                     
                  
		  forwardUsageMap[pParticleSeed][pCandidateCluster] = association;                                                                    
                  backwardUsageMap[pCandidateCluster][pParticleSeed] = association;                                                                   
                }	

              newSeedAssociations.push_back(pCandidateCluster);                                                                                     
              *iterI = NULL;                                                                                                                            
            }	
	  

        }          
      currentSeedAssociations = newSeedAssociations;                                                                                                           
      newSeedAssociations.clear();                                                                                                                             
      ++associationOrder;
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------
BranchGrowingAlgorithm::AssociationType TrackShowerIdFeatureTool::PointsOfContactFeatureTool::AreClustersAssociated(const pandora::Cluster *const pClusterSeed, const pandora::Cluster *const pCluster, const SVMClusterCharacterisationAlgorithm *const pAlgorithm) const 
{

  const VertexList *pVertexList(nullptr);                                                                                                                    
  PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList));                                                   
  const Vertex *const pVertex(((pVertexList->size() == 1) && (VERTEX_3D == (*(pVertexList->begin()))->GetVertexType())) ? *(pVertexList->begin()) : nullptr);
  // Direction of seed cluster (cache for efficiency)                                                                                                        
  ClusterDirectionMap::const_iterator seedIter = m_clusterDirectionMap.find(pClusterSeed);
  
  if (m_clusterDirectionMap.end() == seedIter)                                                                                                               
    {                                                                                   
      const LArVertexHelper::ClusterDirection direction((nullptr == pVertex) ? LArVertexHelper::DIRECTION_UNKNOWN : LArVertexHelper::GetClusterDirectionInZ(this->GetPandora(), pVertex, pClusterSeed, m_directionTanAngle, m_directionApexShift));                                                             
      
      seedIter = m_clusterDirectionMap.insert(ClusterDirectionMap::value_type(pClusterSeed, direction)).first;                                               
    } 
  const LArVertexHelper::ClusterDirection seedDirection(seedIter->second);                                                                                   
  const bool checkSeedForward(seedDirection != LArVertexHelper::DIRECTION_BACKWARD_IN_Z);                                                                    
  const bool checkSeedBackward(seedDirection != LArVertexHelper::DIRECTION_FORWARD_IN_Z);  
  
  // Direction of candidate cluster (cache for efficiency)                                                                                                   
  ClusterDirectionMap::const_iterator candIter = m_clusterDirectionMap.find(pCluster);              
  if (m_clusterDirectionMap.end() == candIter)                                                                                                               
    {                                                                                                                                                
      const LArVertexHelper::ClusterDirection direction((nullptr == pVertex) ? LArVertexHelper::DIRECTION_UNKNOWN : LArVertexHelper::GetClusterDirectionInZ(this->GetPandora(), pVertex, pCluster, m_directionTanAngle, m_directionApexShift));                                                                 
							
      candIter = m_clusterDirectionMap.insert(ClusterDirectionMap::value_type(pCluster, direction)).first;                                                   
    }                                                                                                     
  const LArVertexHelper::ClusterDirection candidateDirection(candIter->second);                                                                              
  const bool checkCandidateForward(candidateDirection != LArVertexHelper::DIRECTION_BACKWARD_IN_Z);                                                          
  const bool checkCandidateBackward(candidateDirection != LArVertexHelper::DIRECTION_FORWARD_IN_Z); 

  // Calculate distances of association                                                                                                                      
  const float sOuter(LArClusterHelper::GetClosestDistance(pClusterSeed->GetCentroid(pClusterSeed->GetOuterPseudoLayer()), pCluster));                        
  const float cOuter(LArClusterHelper::GetClosestDistance(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()), pClusterSeed));                            
  const float sInner(LArClusterHelper::GetClosestDistance(pClusterSeed->GetCentroid(pClusterSeed->GetInnerPseudoLayer()), pCluster));                        
  const float cInner(LArClusterHelper::GetClosestDistance(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()), pClusterSeed));                            
  
  // Association check 1(a), look for enclosed clusters                                                                                                      
  if ((cOuter < m_nearbyClusterDistance && cInner < m_nearbyClusterDistance) &&                                                                              
      (!checkSeedForward || (sInner > m_nearbyClusterDistance)) &&                                                                                           
      (!checkSeedBackward || (sOuter > m_nearbyClusterDistance)))                                                                                            
    { 
      return BranchGrowingAlgorithm::STRONG;                                                                                                                                         
    }                                                      
  // Association check 1(b), look for overlapping clusters                                                                                                   
  if ((checkSeedForward == checkCandidateForward) && (checkSeedBackward == checkCandidateBackward))                                                          
    {                                                                                                                                                    
      if ((cInner < m_nearbyClusterDistance && sOuter < m_nearbyClusterDistance) &&                                                                          
	  (!checkSeedForward || (sInner > m_nearbyClusterDistance)) &&                                                                                       
	  (!checkSeedBackward || (cOuter > m_nearbyClusterDistance)))                                                                                        
        {                                                                                                                                                
	  return BranchGrowingAlgorithm::STRONG;                                                                                                                                     
        }                                                                                                                                            
      if ((cOuter < m_nearbyClusterDistance && sInner < m_nearbyClusterDistance) &&                                                                          
	  (!checkSeedBackward || (sOuter > m_nearbyClusterDistance)) &&                                                                                      
	  (!checkSeedForward || (cInner > m_nearbyClusterDistance)))                                                                                         
        {                                                                                                                                              
	  return BranchGrowingAlgorithm::STRONG;                                                                                                                                     
        }                                                                                                                                             
    }         

  // Association check 2, look for branching clusters                                                                                                        
  if ((!checkSeedForward || (sInner > m_remoteClusterDistance)) &&                                                                                           
      (!checkSeedBackward || (sOuter > m_remoteClusterDistance)) &&                                                                                          
      ((checkCandidateForward && (cInner < m_nearbyClusterDistance)) || (checkCandidateBackward && (cOuter < m_nearbyClusterDistance))))                     
    {                                                                                                                                                   
      return BranchGrowingAlgorithm::STANDARD;                                                                                                                                       
    }                                                                                                                                            
  // Association check 3, look any distance below threshold                                                                                                  
  if ((sOuter < m_nearbyClusterDistance) || (cOuter < m_nearbyClusterDistance) || (sInner < m_nearbyClusterDistance) || (cInner < m_nearbyClusterDistance))  
    return BranchGrowingAlgorithm::SINGLE_ORDER;                                                                                                                                   
                                                                                                                                                           
  return BranchGrowingAlgorithm::NONE;                                                                                                                                               
} 

//------------------------------------------------------------------------------------------------------------------------------------------
void TrackShowerIdFeatureTool::PointsOfContactFeatureTool::IdentifyClusterMerges(const pandora::ClusterVector &particleSeedVector, 
     const BranchGrowingAlgorithm::ClusterUsageMap &backwardUsageMap, BranchGrowingAlgorithm::SeedAssociationList &seedAssociationList) const
{

  ClusterVector sortedCandidates;                                                                                                                              
  for (const auto &mapEntry : backwardUsageMap) sortedCandidates.push_back(mapEntry.first);                                                                    
  std::sort(sortedCandidates.begin(), sortedCandidates.end(), LArClusterHelper::SortByNHits);

  for (const Cluster *const pCluster : sortedCandidates)                                                                                                      
    { 
      const BranchGrowingAlgorithm::ClusterAssociationMap &particleSeedUsageMap(backwardUsageMap.at(pCluster)); 
      if (particleSeedUsageMap.empty())                                                                                                                       
        throw StatusCodeException(STATUS_CODE_FAILURE);
      ClusterVector sortedSeeds;                                                                                                                               
      for (const auto &mapEntry : particleSeedUsageMap) sortedSeeds.push_back(mapEntry.first);                                                                 
      std::sort(sortedSeeds.begin(), sortedSeeds.end(), LArClusterHelper::SortByNHits);
      
      const Cluster *pBestParticleSeed = NULL;                                                                                                                 
      BranchGrowingAlgorithm::AssociationType bestType(BranchGrowingAlgorithm::NONE); 
      
      unsigned int bestOrder(std::numeric_limits<unsigned int>::max());

      for (const Cluster *const pParticleSeed : sortedSeeds)
	{
	  const BranchGrowingAlgorithm::Association &association(particleSeedUsageMap.at(pParticleSeed));
	  if ((association.GetType() > bestType) || ((association.GetType() == bestType) && (association.GetOrder() < bestOrder)))
	    {
	      if ((BranchGrowingAlgorithm::SINGLE_ORDER == association.GetType()) && (association.GetOrder() > 1)) 
		continue;
	      pBestParticleSeed = pParticleSeed;                                                                                                               
              bestType = association.GetType();                                                                                                                
              bestOrder = association.GetOrder(); 
	    }
	  else if ((association.GetType() == bestType) && (association.GetOrder() == bestOrder))
	    {
	      pBestParticleSeed = NULL; 
	    }
	}
      if (NULL == pBestParticleSeed)
	continue;

      seedAssociationList[pBestParticleSeed].push_back(pCluster);
    }
  for (ClusterVector::const_iterator iter = particleSeedVector.begin(), iterEnd = particleSeedVector.end(); iter != iterEnd; ++iter)
    {
      const Cluster *const pParticleSeed = *iter; 
      if (seedAssociationList.end() == seedAssociationList.find(pParticleSeed)) 
	seedAssociationList[pParticleSeed] = ClusterVector();
    }
}
//------------------------------------------------------------------------------------------------------------------------------------------          
   
StatusCode TrackShowerIdFeatureTool::PointsOfContactFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
  
   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));

  
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "NearbyClusterDistance", m_nearbyClusterDistance));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "RemoteClusterDistance", m_remoteClusterDistance));
  PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "DirectionTanAngle", m_directionTanAngle));
   PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "DirectionApexShift", m_directionApexShift));  
 
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
        if((LArClusterHelper::GetClusterHitType(pCluster))!=(LArClusterHelper::GetClusterHitType(pCandidateCluster)))
            continue;                                                                                                                                                 
        if ((pCandidateCluster != pCluster) && (pCandidateCluster->GetNCaloHits() > 5)) //LORENA: check this
            candidateClusters.push_back(pCandidateCluster);                                                                                                     
      }   
 }
    int nPointsOfContact(0);
    
    for (ClusterVector::iterator iterI = candidateClusters.begin(), iterIEnd = candidateClusters.end(); iterI != iterIEnd; ++iterI)    
    {
      const Cluster *const pCandidateCluster = *iterI;                                                                                                 
	  if (NULL == pCandidateCluster)                                                                                                                  
	    continue;
        
      if(this->IsClusterNearby(pCluster,pCandidateCluster))
        ++nPointsOfContact;
    }
 return nPointsOfContact;
}
//------------------------------------------------------------------------------------------------------------------------------------------
bool TrackShowerIdFeatureTool::NNearbyClustersFeatureTool::IsClusterNearby(const pandora::Cluster *const pCluster,const pandora::Cluster *const pCandidateCluster) const
{
       
    /*pandora::CartesianVector innerCoordinate(0.f,0.f,0.f), outerCoordinate(0.f,0.f,0.f); 
    LArClusterHelper::GetExtremalCoordinates(pCandidateCluster,innerCoordinate,outerCoordinate);
    pandora::CartesianVector minimumCoordinate(0.f,0.f,0.f), maximumCoordinate(0.f,0.f,0.f); 
    LArClusterHelper::GetClusterBoundingBox(pCluster, minimumCoordinate, maximumCoordinate);
    
    if(((innerCoordinate.GetX()<=maximumCoordinate.GetX()) && (innerCoordinate.GetX()>=minimumCoordinate.GetX()) && 
        (innerCoordinate.GetZ()<=maximumCoordinate.GetZ()) && (innerCoordinate.GetZ()>=minimumCoordinate.GetZ())) ||
        ((outerCoordinate.GetX()<=maximumCoordinate.GetX()) && (outerCoordinate.GetX()>=minimumCoordinate.GetX()) && 
        (outerCoordinate.GetZ()<=maximumCoordinate.GetZ()) && (outerCoordinate.GetZ()>=minimumCoordinate.GetZ())))
        return true;*/
        
        
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
//LORENA: why do I need this?
StatusCode TrackShowerIdFeatureTool::MipEnergyFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
      PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
                "MipCorrectionPerHit", m_mipCorrectionPerHit));  
    return STATUS_CODE_SUCCESS;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------
 //------------------------------------------------------------------------------------------------------------------------------------------   
 TrackShowerIdFeatureTool::VertexDistanceFeatureTool::MipEnergyFeatureTool() 
{
}
//------------------------------------------------------------------------------------------------------------------------------------------   

void TrackShowerIdFeatureTool::VertexDistanceFeatureTool::Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm *const pAlgorithm, 
           const pandora::Cluster * const pCluster) 
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;
        
  
    float vertexDistancethis->CalculateVertexDistance(pCluster));
  
    featureVector.push_back(vertexDistance); 
    
} 


//------------------------------------------------------------------------------------------------------------------------------------------
float TrackShowerIdFeatureTool::VertexDistanceFeatureTool::CalculateVertexDistance(const pandora::Cluster * const pCluster) const
{
    
    CartesianVector vertexPosition2D(0.f, 0.f, 0.f);
    //this is copied from ClusterCharacterisation, remove from one 
    float vertexDistance = ClusterCharacterisationAlgorithm::GetVertexDistance(this, pCluster, vertexPosition2D);

  return vertexDistance;                                                               
}

} // namespace lar_content


