/**
 *  @file   RemnantClusteringAlgorithm.cc
 * 
 *  @brief  Implementation of the remnant clustering algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "RemnantClusteringAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode RemnantClusteringAlgorithm::Run()
{
    // Get the non-seed cluster list
    const ClusterList *pNonSeedClusterList = NULL;
    const StatusCode statusCode(PandoraContentApi::GetClusterList(*this, m_nonSeedClusterListName, pNonSeedClusterList));

    if ((STATUS_CODE_SUCCESS != statusCode) && (STATUS_CODE_NOT_INITIALIZED != statusCode))
        return statusCode;

    if (STATUS_CODE_NOT_INITIALIZED == statusCode)
        return STATUS_CODE_SUCCESS;

    // ----- Draw OldClusters ----- 
    // ClusterList drawOldClusters(*pNonSeedClusterList);
    // PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
    // PandoraMonitoringApi::VisualizeClusters(&drawOldClusters, "CBMold", BLUE);
    // PandoraMonitoringApi::ViewEvent();
    // ----------------------------

    // Delete the current list of non-seed clusters to free the underlying hits
    ClusterList clustersToDelete(*pNonSeedClusterList);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::DeleteClusters(*this, clustersToDelete, m_nonSeedClusterListName));

    // Create a temporary cluster list
    const ClusterList *pClusterList = NULL; std::string clusterListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryClusterListAndSetCurrent(*this, pClusterList, clusterListName));

    // Generate a list of available hits for re-clustering
    const CaloHitList* pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentCaloHitList(*this, pCaloHitList));

    // Run a simple clustering algorithm over the available hits
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FindClusters(pCaloHitList));

    // Build the new clusters
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->BuildClusters(pCaloHitList));

    // ----- Draw New Clusters -----
    // ClusterList drawNewClusters(*pClusterList);
    // PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
    // PandoraMonitoringApi::VisualizeClusters(&drawNewClusters, "CBMnew", GREEN);
    // PandoraMonitoringApi::ViewEvent();
    // -----------------------------

    // Copy any new clusters into the old list
    if (!pClusterList->empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveClusterList(*this, clusterListName, m_nonSeedClusterListName));

    // ----- Draw Replaced Clusters ----- 
    // const ClusterList *pNonSeedClusterListAgain = NULL;
    // PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetClusterList(*this, m_nonSeedClusterListName, pNonSeedClusterListAgain));
    // ClusterList drawReplacedClusters(*pNonSeedClusterListAgain);
    // PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
    // PandoraMonitoringApi::VisualizeClusters(&drawReplacedClusters, "CBMreplaced", RED);
    // PandoraMonitoringApi::ViewEvent();
    // ----------------------------------

    return STATUS_CODE_SUCCESS;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RemnantClusteringAlgorithm::FindClusters(const CaloHitList* pCaloHitList)
{
  // Generate an ordered list of hits
  OrderedCaloHitList orderedCaloHitList;

  for (CaloHitList::const_iterator iter = pCaloHitList->begin(), iterEnd = pCaloHitList->end(); iter != iterEnd; ++iter)
  {
    CaloHit *pCaloHit = *iter;
    if ((pCaloHit->GetMipEquivalentEnergy() >= m_minPulseHeight) && PandoraContentApi::IsCaloHitAvailable(*this, pCaloHit))
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, orderedCaloHitList.Add(pCaloHit));
  }

  // Run clustering algorithm
  fClusterID = 0;

  fHitAssociations.clear();
  fClusterAssociations.clear();
  fClusterSizes.clear();
  fClusterMap.clear();

  for( OrderedCaloHitList::const_iterator iterI = orderedCaloHitList.begin(); iterI != orderedCaloHitList.end(); ++iterI ) {
     
    PseudoLayer  LayerI = iterI->first;
    CaloHitList* ListI  = iterI->second;

    unsigned int ilayer = 0;

    for( OrderedCaloHitList::const_iterator iterJ = iterI; ilayer<m_clusterLayers+1 && iterJ != orderedCaloHitList.end(); ++iterJ ) {
      ++ilayer;

      PseudoLayer  LayerJ = iterJ->first;
      CaloHitList* ListJ  = iterJ->second;
         
      if( LayerJ-LayerI>m_clusterLayers ) continue;

      for( CaloHitList::iterator hitI = ListI->begin(); hitI != ListI->end(); ++hitI ) {
        CaloHit* pCaloHitI = *hitI;

        for( CaloHitList::iterator hitJ = ListJ->begin(); hitJ != ListJ->end(); ++hitJ ) {
          CaloHit* pCaloHitJ = *hitJ;
              
          if( IsAssociated( pCaloHitI, pCaloHitJ ) ){
            MakeAssociation( pCaloHitI, pCaloHitJ );
	  }
	}

      }
    }
  }

  // Find cluster sizes
  int firstClusterSize = 1;

  for( OrderedCaloHitList::const_iterator iterI = orderedCaloHitList.begin(); iterI != orderedCaloHitList.end(); ++iterI ) {
    CaloHitList* List  = iterI->second;
         
    for( CaloHitList::iterator hit = List->begin(); hit != List->end(); ++hit ) {
      CaloHit* pCaloHit = *hit;
 
      int clusterID = GetClusterID( pCaloHit );
      if( clusterID<0 ) continue;

      ClusterAssociationMap::iterator iterJ = fClusterSizes.find( clusterID );
      
      if( iterJ==fClusterSizes.end() ){
        fClusterSizes.insert( std::pair<int,int>(clusterID,firstClusterSize) );  
      }
      else{
        ++(iterJ->second);
      } 
    }
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RemnantClusteringAlgorithm::BuildClusters(const CaloHitList* pCaloHitList)
{
  // ----- Debug -----
  // std::string ListName = "No Name";
  //
  // PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentCaloHitListName(*this,ListName));
  // std::cout << " Current Hit List: " << ListName << std::endl;
  //
  // PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentClusterListName(*this,ListName));
  // std::cout << " Current Cluster List: " << ListName << std::endl;
  // ----------------

  OrderedCaloHitList orderedCaloHitList;
  PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, orderedCaloHitList.Add(*pCaloHitList));

  for( OrderedCaloHitList::const_iterator iterI = orderedCaloHitList.begin(); iterI != orderedCaloHitList.end(); ++iterI ) {
    CaloHitList* List  = iterI->second;
         
    for( CaloHitList::iterator hit = List->begin(); hit != List->end(); ++hit ) {
      CaloHit* pCaloHit = *hit;
      Cluster* pCluster = NULL;

      int   clusterID   = GetClusterID( pCaloHit );
      int   clusterSize = 0;

      if( clusterID<0 ) continue;
      
      ClusterAssociationMap::iterator iterJ = fClusterSizes.find( clusterID );      
      
      if( iterJ!=fClusterSizes.end() ){
        clusterSize = iterJ->second;
      }

      if( clusterSize<m_clusterMinSize ) continue;

      ClusterMap::iterator iterK = fClusterMap.find( clusterID );
      
      if( iterK==fClusterMap.end() ){
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, pCaloHit, pCluster));   
        fClusterMap.insert( std::pair<int,Cluster*>(clusterID,pCluster) );  
      }
      else{
        pCluster = iterK->second;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddCaloHitToCluster(*this, pCluster, pCaloHit));
      } 
    }
  }

  return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool RemnantClusteringAlgorithm::IsAssociated( CaloHit* pCaloHitI, CaloHit* pCaloHitJ )
{
  assert( pCaloHitI && pCaloHitJ );

  if( pCaloHitI==pCaloHitJ ) return false;

  const CartesianVector hitSeparation( pCaloHitJ->GetPositionVector()-pCaloHitI->GetPositionVector() );

  if( hitSeparation.GetMagnitudeSquared()<m_clusterRadiusSquared ) return true;
  else return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void RemnantClusteringAlgorithm::MakeAssociation( CaloHit* pCaloHitI, CaloHit* pCaloHitJ )
{
  int clusterI = GetClusterID( pCaloHitI );
  int clusterJ = GetClusterID( pCaloHitJ );

  if( clusterI<0 && clusterJ<0 ){
    ++fClusterID;
    SetClusterID( pCaloHitI, fClusterID );
    SetClusterID( pCaloHitJ, fClusterID );
  }
  else if( clusterI<0 ){
    SetClusterID( pCaloHitI, GetClusterID( pCaloHitJ ) );
  }
  else if( clusterJ<0 ){
    SetClusterID( pCaloHitJ, GetClusterID( pCaloHitI ) );
  }
  else if( clusterI!=clusterJ ){
    ResetClusterID( pCaloHitI, pCaloHitJ );
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int RemnantClusteringAlgorithm::GetClusterID( CaloHit* pCaloHit )
{ 
  int clusterID = -1;
  
  HitAssociationMap::iterator iterI = fHitAssociations.find( pCaloHit );

  if( iterI!=fHitAssociations.end() ){
    clusterID = iterI->second;

    ClusterAssociationMap::iterator iterJ = fClusterAssociations.find( clusterID );

    if( iterJ!=fClusterAssociations.end() ){
      clusterID = std::min(clusterID,iterJ->second);
    }
  }

  return clusterID;
}
 
//------------------------------------------------------------------------------------------------------------------------------------------
   
void RemnantClusteringAlgorithm::SetClusterID( CaloHit* pCaloHit, int clusterID )
{
  HitAssociationMap::iterator iter = fHitAssociations.find( pCaloHit );

  if( iter==fHitAssociations.end() ){
    fHitAssociations.insert( std::pair<CaloHit*,int>(pCaloHit,clusterID) );
  }

  return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void RemnantClusteringAlgorithm::ResetClusterID( CaloHit* pCaloHitI, CaloHit* pCaloHitJ )
{ 
  int clusterI = GetClusterID( pCaloHitI );
  int clusterJ = GetClusterID( pCaloHitJ );
  
  if( clusterI==clusterJ ) return;

  int newClusterI = clusterI;

  ClusterAssociationMap::iterator iterI = fClusterAssociations.find( clusterI );

  if( iterI!=fClusterAssociations.end() ){
    newClusterI = iterI->second;
  }

  int newClusterJ = clusterJ;

  ClusterAssociationMap::iterator iterJ = fClusterAssociations.find( clusterJ );

  if( iterJ!=fClusterAssociations.end() ){
    newClusterJ = iterJ->second;
  }

  int newClusterID = std::min(newClusterI,newClusterJ);

  if( newClusterID<newClusterI ){
    if( newClusterI==clusterI ){
      fClusterAssociations.insert( std::pair<int,int>(clusterI,newClusterID) ); 
    }
    else{
      iterI->second = newClusterID;
    }
  }

  if( newClusterID<newClusterJ ){
    if( newClusterJ==clusterJ ){
      fClusterAssociations.insert( std::pair<int,int>(clusterJ,newClusterID) ); 
    }
    else{
      iterJ->second = newClusterID;
    }
  }

  return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RemnantClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));
 
    m_minPulseHeight = 0.25; // MIPS
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumPulseHeight", m_minPulseHeight));

    m_clusterRadius = 2.5; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterRadius", m_clusterRadius));
    m_clusterRadiusSquared = m_clusterRadius*m_clusterRadius;

    m_clusterLayers = 5; // layers
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterLayers", m_clusterLayers));

    m_clusterMinSize = 4; // hits
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ClusterMinimumSize", m_clusterMinSize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
