/**
 *  @file   LArContent/src/LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.cc
 * 
 *  @brief  Implementation of the bounded cluster merging algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArTwoDReco/LArClusterMopUp/BoundedClusterMergingAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode BoundedClusterMergingAlgorithm::Run()
{
    const ClusterList *pSeedClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_seedClusterListName, pSeedClusterList));

    const ClusterList *pNonSeedClusterList = NULL;
    const StatusCode listStatusCode(PandoraContentApi::GetList(*this, m_nonSeedClusterListName, pNonSeedClusterList));

    if ((STATUS_CODE_SUCCESS != listStatusCode) && (STATUS_CODE_NOT_INITIALIZED != listStatusCode))
        return listStatusCode;

    if (STATUS_CODE_NOT_INITIALIZED == listStatusCode)
        return STATUS_CODE_SUCCESS;

    this->PrepareAssociations();

    for (ClusterList::const_iterator iterI = pSeedClusterList->begin(), iterEndI = pSeedClusterList->end(); iterI != iterEndI; ++iterI)
    {
        Cluster *pSeedCluster = *iterI;

        this->ConstructBoundingBox(pSeedCluster);

        for (ClusterList::const_iterator iterJ = pNonSeedClusterList->begin(), iterEndJ = pNonSeedClusterList->end(); iterJ != iterEndJ; ++iterJ)
        {
            Cluster *pCandidateCluster = *iterJ;

            if (this->PassesVertexVeto(pSeedCluster, pCandidateCluster) && this->PassesLengthVeto(pSeedCluster, pCandidateCluster))
            {
                this->MakeAssociation(pSeedCluster, pCandidateCluster, this->ApplyBoundingBox(pCandidateCluster));
            }
        }
    }

    this->PrepareMerges();

    for (ClusterMergeMap::const_iterator iter = fClusterMergeMap.begin(), iterEnd = fClusterMergeMap.end(); iter != iterEnd; ++iter)
    {
        for (ClusterList::const_iterator dauIter = iter->second.begin(), dauIterEnd = iter->second.end(); dauIter != dauIterEnd; ++dauIter)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, iter->first, *dauIter,
                m_seedClusterListName, m_nonSeedClusterListName));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BoundedClusterMergingAlgorithm::ConstructBoundingBox(const Cluster *seedCluster)
{
    return fBoundingBox.BuildBoundingBox(seedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int BoundedClusterMergingAlgorithm::ApplyBoundingBox(const Cluster *candidateCluster)
{
    return fBoundingBox.NumberOfEnclosedEnds(candidateCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BoundedClusterMergingAlgorithm::PassesLengthVeto(const Cluster *seedCluster, const Cluster *candCluster)
{
    if (seedCluster->GetNCaloHits() >= m_minClusterSize && candCluster->GetNCaloHits()>=m_minClusterSize)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BoundedClusterMergingAlgorithm::PassesVertexVeto(const Cluster *seedCluster, const Cluster *candCluster)
{
    const CartesianVector seedInnerVertex = seedCluster->GetCentroid(seedCluster->GetInnerPseudoLayer());
    const CartesianVector seedOuterVertex = seedCluster->GetCentroid(seedCluster->GetOuterPseudoLayer());

    const CartesianVector candInnerVertex = candCluster->GetCentroid(candCluster->GetInnerPseudoLayer());
    const CartesianVector candOuterVertex = candCluster->GetCentroid(candCluster->GetOuterPseudoLayer());

    if ((seedInnerVertex-candInnerVertex).GetMagnitude() > m_vertexVetoRadius && (seedOuterVertex-candOuterVertex).GetMagnitude() > m_vertexVetoRadius)
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BoundedClusterMergingAlgorithm::PrepareAssociations()
{
    fGoodAssociations.clear();
    fBadAssociations.clear();

    fSingleAssociations.clear();
    fRepeatedSingleAssociations.clear();

    fDoubleAssociations.clear();
    fRepeatedDoubleAssociations.clear();

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BoundedClusterMergingAlgorithm::MakeAssociation(Cluster *pSeedCluster, Cluster *pCandidateCluster, unsigned int numEnds)
{
    if (numEnds == 0)
        return;

    if (numEnds == 1)
    {
        ClusterAssociationMap::const_iterator iter1 = fSingleAssociations.find(pCandidateCluster);

        if (iter1 == fSingleAssociations.end())
        {
            fSingleAssociations.insert(std::pair<Cluster*, Cluster*>(pCandidateCluster, pSeedCluster));
        }
        else
        {
            ClusterAssociationMap::const_iterator iter2 = fRepeatedSingleAssociations.find(pCandidateCluster);

            if (iter2 == fRepeatedSingleAssociations.end())
                fRepeatedSingleAssociations.insert( std::pair<Cluster*,Cluster*>(pCandidateCluster,NULL) ); 
        }
    }

    else if (numEnds==2)
    {
        ClusterAssociationMap::const_iterator iter1 = fDoubleAssociations.find(pCandidateCluster);

        if (iter1 == fDoubleAssociations.end())
        {
            fDoubleAssociations.insert(std::pair<Cluster*, Cluster*>(pCandidateCluster,pSeedCluster));
        }
        else
        {
            ClusterAssociationMap::const_iterator iter2 = fRepeatedDoubleAssociations.find(pCandidateCluster);

            if (iter2 == fRepeatedDoubleAssociations.end())
                fRepeatedDoubleAssociations.insert(std::pair<Cluster*, Cluster*>(pCandidateCluster,NULL));
        }
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BoundedClusterMergingAlgorithm::SetGoodAssociation(Cluster *pCandidateCluster)
{
    ClusterAssociationMap::const_iterator iter = fGoodAssociations.find(pCandidateCluster);

    if (iter == fGoodAssociations.end())
        fGoodAssociations.insert(std::pair<Cluster*, Cluster*>(pCandidateCluster, NULL));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BoundedClusterMergingAlgorithm::SetBadAssociation(Cluster *pCandidateCluster)
{
    ClusterAssociationMap::const_iterator iter = fBadAssociations.find(pCandidateCluster);

    if (iter == fBadAssociations.end())
        fBadAssociations.insert(std::pair<Cluster*, Cluster*>(pCandidateCluster, NULL));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BoundedClusterMergingAlgorithm::PrepareMerges()
{
    fClusterMergeMap.clear();

    for (ClusterAssociationMap::const_iterator iter = fDoubleAssociations.begin(), iterEnd = fDoubleAssociations.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCandidateCluster = iter->first;
        Cluster *pSeedCluster = iter->second;

        ClusterAssociationMap::const_iterator iter2A = fGoodAssociations.find(pCandidateCluster);
        ClusterAssociationMap::const_iterator iter2B = fBadAssociations.find(pCandidateCluster);

        if (iter2A == fGoodAssociations.end() && iter2B == fBadAssociations.end())
        {
            ClusterAssociationMap::const_iterator iter3 = fRepeatedDoubleAssociations.find(pCandidateCluster);

            if (iter3 == fRepeatedDoubleAssociations.end())
            {
                fClusterMergeMap[pSeedCluster].insert(pCandidateCluster);
                SetGoodAssociation(pCandidateCluster); 
            }
            else
            {
                SetBadAssociation(pCandidateCluster); 
            }
        }
    }

    for (ClusterAssociationMap::const_iterator iter = fSingleAssociations.begin(), iterEnd = fSingleAssociations.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCandidateCluster = iter->first;
        Cluster *pSeedCluster = iter->second;

        ClusterAssociationMap::const_iterator iter2A = fGoodAssociations.find(pCandidateCluster);
        ClusterAssociationMap::const_iterator iter2B = fBadAssociations.find(pCandidateCluster);

        if (iter2A == fGoodAssociations.end() && iter2B == fBadAssociations.end())
        {
            ClusterAssociationMap::const_iterator iter3 = fRepeatedSingleAssociations.find( pCandidateCluster );

            if (iter3 == fRepeatedSingleAssociations.end())
            {
                fClusterMergeMap[pSeedCluster].insert(pCandidateCluster);
                SetGoodAssociation(pCandidateCluster); 
            }
            else
            {
                SetBadAssociation(pCandidateCluster); 
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

BoundedCluster::BoundedCluster(const pandora::Cluster *pCluster)
{
  fBoxFlag = NULL;
  fBoxMinX = NULL;
  fBoxMaxX = NULL;
  fBoxMinY = NULL;
  fBoxMaxY = NULL;  
  fBoxLayers = 0; 

  fMinLayer = 0;
  fMaxLayer = 0;
  fNumLayers = 0;

  BuildBoundingBox( pCluster );
}

//------------------------------------------------------------------------------------------------------------------------------------------

BoundedCluster::BoundedCluster(const BoundedCluster& rhs)
{
  fMinLayer  = rhs.fMinLayer;
  fMaxLayer  = rhs.fMaxLayer;
  fNumLayers = rhs.fNumLayers;
  fBoxLayers = rhs.fBoxLayers;
  
  if( fBoxLayers>0 ){
    fBoxFlag = new bool[fBoxLayers];
    fBoxMinX = new float[fBoxLayers];
    fBoxMaxX = new float[fBoxLayers];
    fBoxMinY = new float[fBoxLayers];
    fBoxMaxY = new float[fBoxLayers];  

    for( unsigned int n=0; n<fBoxLayers; ++n ) {
      fBoxFlag[n] = rhs.fBoxFlag[n];
      fBoxMinX[n] = rhs.fBoxMinX[n];
      fBoxMaxX[n] = rhs.fBoxMaxX[n];
      fBoxMinY[n] = rhs.fBoxMinY[n];
      fBoxMaxY[n] = rhs.fBoxMaxY[n];   
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

BoundedCluster::~BoundedCluster()
{
  if( fBoxFlag ) delete [] fBoxFlag;
  if( fBoxMinX ) delete [] fBoxMinX;
  if( fBoxMaxX ) delete [] fBoxMaxX; 
  if( fBoxMinY ) delete [] fBoxMinY;
  if( fBoxMaxY ) delete [] fBoxMaxY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BoundedCluster::BuildBoundingBox(const Cluster *pCluster)
{ 
  fMinLayer = 0;
  fMaxLayer = 0;
  fNumLayers = 0;

  if( pCluster==NULL ) return;

  // set lower and upper layer of bounding box
  const unsigned int m_boundingBoxExtraLayers = 3;
  const unsigned int firstLayer = pCluster->GetInnerPseudoLayer();
  const unsigned int lastLayer  = pCluster->GetOuterPseudoLayer();

  fMinLayer  = firstLayer;
  fMaxLayer  = lastLayer + m_boundingBoxExtraLayers;
 
  if( fMinLayer>=m_boundingBoxExtraLayers ) fMinLayer -= m_boundingBoxExtraLayers;
  else                                      fMinLayer = 0;
   
  fNumLayers = 1+fMaxLayer-fMinLayer;

  // reset internal arrays between lower and upper layers
  if( fNumLayers>fBoxLayers ){
    fBoxLayers = fNumLayers;
      
    if( fBoxFlag ) delete [] fBoxFlag;
    if( fBoxMinX ) delete [] fBoxMinX;
    if( fBoxMaxX ) delete [] fBoxMaxX; 
    if( fBoxMinY ) delete [] fBoxMinY;
    if( fBoxMaxY ) delete [] fBoxMaxY;

    fBoxFlag = new bool[fBoxLayers];
    fBoxMinX = new float[fBoxLayers];
    fBoxMaxX = new float[fBoxLayers];
    fBoxMinY = new float[fBoxLayers];
    fBoxMaxY = new float[fBoxLayers];  
  }

  for( unsigned int ilayer=0; ilayer<fNumLayers; ++ilayer ){
    fBoxFlag[ilayer] = false;
    fBoxMinX[ilayer] = 0.0;
    fBoxMaxX[ilayer] = 0.0;
    fBoxMinY[ilayer] = 0.0;
    fBoxMaxY[ilayer] = 0.0;  
  }

  // generate bounding box
  const OrderedCaloHitList &theHitList( pCluster->GetOrderedCaloHitList() );

  for( OrderedCaloHitList::const_iterator iter1 = theHitList.begin(); iter1 != theHitList.end(); ++iter1 ) {
    const unsigned int thisLayer = iter1->first;
    const CaloHitList* thisList  = iter1->second;

    for( CaloHitList::const_iterator iter2 = thisList->begin(); iter2 != thisList->end(); ++iter2 ) {
      const CaloHit* thisHit = *iter2;

      const float hitX = thisHit->GetPositionVector().GetX();
      const float hitY = thisHit->GetPositionVector().GetY();

      unsigned int minLayer = thisLayer;
      unsigned int maxLayer = thisLayer + m_boundingBoxExtraLayers;

      if( thisLayer>=m_boundingBoxExtraLayers ) minLayer -= m_boundingBoxExtraLayers;
      else                                      minLayer = 0;

      for( unsigned int nLayer=minLayer; nLayer<=maxLayer; ++nLayer ){
        unsigned int ilayer = nLayer-fMinLayer;

        if ((nLayer >= fMinLayer) && (ilayer < fNumLayers))
        {
          if( fBoxFlag[ilayer] == false ){
            fBoxMinX[ilayer] = hitX;
            fBoxMaxX[ilayer] = hitX;
            fBoxMinY[ilayer] = hitY;
            fBoxMinY[ilayer] = hitY;
            fBoxFlag[ilayer] = true;
        }
          else{
                 if( hitX<fBoxMinX[ilayer] ) fBoxMinX[ilayer] = hitX;
            else if( hitX>fBoxMaxX[ilayer] ) fBoxMaxX[ilayer] = hitX;
                 if( hitY<fBoxMinY[ilayer] ) fBoxMinY[ilayer] = hitY;
            else if( hitY>fBoxMaxY[ilayer] ) fBoxMaxY[ilayer] = hitY;
            }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int BoundedCluster::NumberOfEnclosedEnds( const Cluster* pCluster)
{
  // sanity check
  if( pCluster==NULL || fNumLayers==0 ) return 0;

  // check that cluster is enclosed in Z
  PseudoLayer innerLayer = pCluster->GetInnerPseudoLayer();
  PseudoLayer outerLayer = pCluster->GetOuterPseudoLayer();

  if( outerLayer<fMinLayer || innerLayer>fMaxLayer ) return 0;

  // count number of ends inside bounding box
  const float m_boundingBoxWindow = 0.5; // cm
  unsigned int numEnds = 0;
  unsigned int ilayer  = 0;

  // check that inner layer of cluster is enclosed in X,Y
  const CartesianVector innerVertex = pCluster->GetCentroid( innerLayer );
  const float innerX = innerVertex.GetX();
  const float innerY = innerVertex.GetY();
 
  if( innerLayer>=fMinLayer && innerLayer<=fMaxLayer ){
    ilayer = innerLayer-fMinLayer;
    if( fBoxFlag[ilayer]==true
     && innerX>fBoxMinX[ilayer]-m_boundingBoxWindow
     && innerX<fBoxMaxX[ilayer]+m_boundingBoxWindow
     && innerY>fBoxMinY[ilayer]-m_boundingBoxWindow
     && innerY<fBoxMaxY[ilayer]+m_boundingBoxWindow ) ++numEnds;
  }

  // check that outer layer of cluster is enclosed in X,Y
  const CartesianVector outerVertex = pCluster->GetCentroid( outerLayer );
  const float outerX = outerVertex.GetX();
  const float outerY = outerVertex.GetY();
  
  if( outerLayer>=fMinLayer && outerLayer<=fMaxLayer ){
    ilayer = outerLayer-fMinLayer;
    if( fBoxFlag[ilayer]==true 
     && outerX>fBoxMinX[ilayer]-m_boundingBoxWindow
     && outerX<fBoxMaxX[ilayer]+m_boundingBoxWindow
     && outerY>fBoxMinY[ilayer]-m_boundingBoxWindow
     && outerY<fBoxMaxY[ilayer]+m_boundingBoxWindow ) ++numEnds;
  }

  return numEnds;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BoundedClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    m_vertexVetoRadius = 2.5; // cm
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "VertexVetoRadius", m_vertexVetoRadius));

    m_minClusterSize = 4;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "MinClusterSize", m_minClusterSize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
