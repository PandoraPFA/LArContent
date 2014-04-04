/**
 *  @file   LArContent/src/LArMonitoring/NtupleWritingAlgorithm.cc
 * 
 *  @brief  Implementation of the ntuple writing algorithm
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArMonitoring/NtupleWritingAlgorithm.h"

#include "LArObjects/LArPointingCluster.h"

using namespace pandora;

namespace lar
{

StatusCode NtupleWritingAlgorithm::Run()
{ 
    const ClusterList *pSeedClusterList = NULL;
    const ClusterList *pNonSeedClusterList = NULL;

    if ( STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, m_seedClusterListName, pSeedClusterList) )
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_nonSeedClusterListName, pNonSeedClusterList));
    }
    else if ( STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_nonSeedClusterListName, pNonSeedClusterList) )
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pNonSeedClusterList));
    }
    else
    {
        std::cout << " Cannot find an input cluster list " << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);
    }

    // Get view
    HitType viewType(INNER_DETECTOR);
    
    if( NULL != pSeedClusterList && false == pSeedClusterList->empty() )
    {
        Cluster* pCluster = *(pSeedClusterList->begin());
        CaloHitList* caloHitList = (pCluster->GetOrderedCaloHitList()).begin()->second;
        CaloHit* pCaloHit = *(caloHitList->begin());
        viewType = pCaloHit->GetHitType();
    }
    else if( NULL != pNonSeedClusterList && false == pNonSeedClusterList->empty() )
    {
        Cluster* pCluster = *(pNonSeedClusterList->begin());
        CaloHitList* caloHitList = (pCluster->GetOrderedCaloHitList()).begin()->second;
        CaloHit* pCaloHit = *(caloHitList->begin());
        viewType = pCaloHit->GetHitType();
    }
    else 
    {
        return STATUS_CODE_NOT_INITIALIZED;
    }


    if ( INNER_DETECTOR == viewType )  return STATUS_CODE_FAILURE;


    // Book-keeping
    static int eventNumber(0);
    static bool foundU(false), foundV(false), foundW(false);

    bool foundNewEvent(false);

    if ( TPC_VIEW_U == viewType )
    {
        if ( false == foundU ) foundU = true; else foundNewEvent = true;
    }

    if ( TPC_VIEW_V == viewType )
    {
        if ( false == foundV ) foundV = true; else foundNewEvent = true;
    }

    if ( TPC_VIEW_W == viewType )
    {
        if ( false == foundW ) foundW = true; else foundNewEvent = true;
    }

    if ( true == foundNewEvent )
    {
        ++eventNumber; foundU = false; foundV = false; foundW = false;
    }


    // Create ntuple writer
    static NtupleWriter ntupleWriter( m_outputFileName, m_outputTreeName );

    
    // Add the new event
    ntupleWriter.NewEntry( eventNumber, viewType );


    // Add the seed clusters
    if( NULL != pSeedClusterList )
    {
        ClusterVector seedClusterVector;
        for ( ClusterList::const_iterator iter = pSeedClusterList->begin(), iterEnd = pSeedClusterList->end(); iter != iterEnd; ++iter ) 
            seedClusterVector.push_back(*iter);
        std::sort(seedClusterVector.begin(), seedClusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);

        for (ClusterVector::const_iterator iter = seedClusterVector.begin(), iterEnd = seedClusterVector.end(); iter != iterEnd; ++iter)
            ntupleWriter.AddCluster(*iter,true);
    }

    // Add the non-seed clusters
    if( NULL != pNonSeedClusterList )
    {
        ClusterVector nonSeedClusterVector;
        for ( ClusterList::const_iterator iter = pNonSeedClusterList->begin(), iterEnd = pNonSeedClusterList->end(); iter != iterEnd; ++iter ) 
            nonSeedClusterVector.push_back(*iter);
        std::sort(nonSeedClusterVector.begin(), nonSeedClusterVector.end(), LArClusterHelper::SortByNOccupiedLayers);

        for (ClusterVector::const_iterator iter = nonSeedClusterVector.begin(), iterEnd = nonSeedClusterVector.end(); iter != iterEnd; ++iter)
            ntupleWriter.AddCluster(*iter,false);
    }
    


    
    ntupleWriter.WriteEntry();
    
    

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void NtupleWritingAlgorithm::NtupleWriter::AddCluster( Cluster* pCluster, bool isSeed )
{   
    // Write out hit info
    const OrderedCaloHitList& orderedCaloHitList = pCluster->GetOrderedCaloHitList();

    for ( OrderedCaloHitList::const_iterator hitListIter = orderedCaloHitList.begin(), hitListIterEnd = orderedCaloHitList.end(); hitListIter != hitListIterEnd; ++hitListIter )
    { 
        CaloHitList* hitList = hitListIter->second;

        for( CaloHitList::const_iterator hitIter = hitList->begin(), hitIterEnd = hitList->end(); hitIter != hitIterEnd; ++hitIter ) 
        {
            CaloHit* pCaloHit = *hitIter;

            m_hitID.push_back( m_nHits++ ); 
            m_hitClusterID.push_back( m_nClusters ); 
            m_hitPosX.push_back( pCaloHit->GetPositionVector().GetX() );
            m_hitPosY.push_back( pCaloHit->GetPositionVector().GetY() );
            m_hitPosZ.push_back( pCaloHit->GetPositionVector().GetZ() );
            m_hitEnergy.push_back( pCaloHit->GetHadronicEnergy() );
        }
    }

    // Write out cluster info
    if ( isSeed ) ++m_nSeedClusters;

    m_clusterID.push_back( m_nClusters++ );
    m_clusterIsSeed.push_back( isSeed );
    m_clusterEnergy.push_back( pCluster->GetHadronicEnergy() );
    m_clusterLayers.push_back( orderedCaloHitList.size() );

    if ( orderedCaloHitList.size() > 1 )
    {
        const LArPointingCluster &pointingCluster(pCluster);

        // TODO use new vertex functionality to determine direction
        const int directionInZ(0);

        m_clusterDirectionInZ.push_back( directionInZ );
        m_clusterLength.push_back( pointingCluster.GetLength() );

        m_clusterStartPosX.push_back( pointingCluster.GetInnerVertex().GetPosition().GetX() );
        m_clusterStartPosY.push_back( pointingCluster.GetInnerVertex().GetPosition().GetY() );
        m_clusterStartPosZ.push_back( pointingCluster.GetInnerVertex().GetPosition().GetZ() );
        m_clusterStartDirX.push_back( pointingCluster.GetInnerVertex().GetDirection().GetX() );
        m_clusterStartDirY.push_back( pointingCluster.GetInnerVertex().GetDirection().GetY() );
        m_clusterStartDirZ.push_back( pointingCluster.GetInnerVertex().GetDirection().GetZ() );
        m_clusterEndPosX.push_back( pointingCluster.GetOuterVertex().GetPosition().GetX() );
        m_clusterEndPosY.push_back( pointingCluster.GetOuterVertex().GetPosition().GetY() );
        m_clusterEndPosZ.push_back( pointingCluster.GetOuterVertex().GetPosition().GetZ() );
        m_clusterEndDirX.push_back( pointingCluster.GetOuterVertex().GetDirection().GetX() );
        m_clusterEndDirY.push_back( pointingCluster.GetOuterVertex().GetDirection().GetY() );
        m_clusterEndDirZ.push_back( pointingCluster.GetOuterVertex().GetDirection().GetZ() );
    }
    else
    {
        CartesianVector dummyDirection( 0.f, 0.f, 1.f);

        m_clusterDirectionInZ.push_back(2);
        m_clusterLength.push_back( 0.f );

        m_clusterStartPosX.push_back( pCluster->GetCentroid( pCluster->GetInnerPseudoLayer() ).GetX() );
        m_clusterStartPosY.push_back( pCluster->GetCentroid( pCluster->GetInnerPseudoLayer() ).GetY() );
        m_clusterStartPosZ.push_back( pCluster->GetCentroid( pCluster->GetInnerPseudoLayer() ).GetZ() );
        m_clusterStartDirX.push_back( dummyDirection.GetX() );
        m_clusterStartDirY.push_back( dummyDirection.GetY() );
        m_clusterStartDirZ.push_back( dummyDirection.GetZ() );
        m_clusterEndPosX.push_back( pCluster->GetCentroid( pCluster->GetOuterPseudoLayer() ).GetX() );
        m_clusterEndPosY.push_back( pCluster->GetCentroid( pCluster->GetOuterPseudoLayer() ).GetY() );
        m_clusterEndPosZ.push_back( pCluster->GetCentroid( pCluster->GetOuterPseudoLayer() ).GetZ() );
        m_clusterEndDirX.push_back( dummyDirection.GetX() );
        m_clusterEndDirY.push_back( dummyDirection.GetY() );
        m_clusterEndDirZ.push_back( dummyDirection.GetZ() );
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

NtupleWritingAlgorithm::NtupleWriter::NtupleWriter() :
    m_fileName("file.root"), m_treeExistsU(0), m_treeExistsV(0), m_treeExistsW(0),
    m_event(0), m_view(0), m_nHits(0), m_nClusters(0), m_nSeedClusters(0)
{    
    std::string treeName("tree");
 
    m_treeNameU = treeName;  m_treeNameU.append("U");
    m_treeNameV = treeName;  m_treeNameV.append("V");
    m_treeNameW = treeName;  m_treeNameW.append("W");

    PANDORA_MONITORING_API(Create());
}

//------------------------------------------------------------------------------------------------------------------------------------------

NtupleWritingAlgorithm::NtupleWriter::NtupleWriter( std::string fileName, std::string treeName ) :
    m_fileName(fileName), m_treeExistsU(0), m_treeExistsV(0), m_treeExistsW(0),
    m_event(0), m_view(0), m_nHits(0), m_nClusters(0), m_nSeedClusters(0)
{    
    m_treeNameU = treeName;  m_treeNameU.append("U");
    m_treeNameV = treeName;  m_treeNameV.append("V");
    m_treeNameW = treeName;  m_treeNameW.append("W");

    PANDORA_MONITORING_API(Create());
}


//------------------------------------------------------------------------------------------------------------------------------------------
    
NtupleWritingAlgorithm::NtupleWriter::~NtupleWriter()
{
    bool fileExists(false);

    if ( m_treeExistsU )
    {
        if ( false==fileExists )
        {
            PANDORA_MONITORING_API(SaveTree(m_treeNameU.c_str(), m_fileName.c_str(), "RECREATE")); fileExists = true;          
        }
        else PANDORA_MONITORING_API(SaveTree(m_treeNameU.c_str(), m_fileName.c_str(), "UPDATE"));
    }

    if ( m_treeExistsV )
    {
        if ( false==fileExists )
        {
            PANDORA_MONITORING_API(SaveTree(m_treeNameV.c_str(), m_fileName.c_str(), "RECREATE")); fileExists = true;
        }
        else PANDORA_MONITORING_API(SaveTree(m_treeNameV.c_str(), m_fileName.c_str(), "UPDATE"));
    }

    if ( m_treeExistsW )
    {
        if ( false==fileExists )
        {
            PANDORA_MONITORING_API(SaveTree(m_treeNameW.c_str(), m_fileName.c_str(), "RECREATE")); fileExists = true;
        }
        else PANDORA_MONITORING_API(SaveTree(m_treeNameW.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

void NtupleWritingAlgorithm::NtupleWriter::NewEntry( int eventNumber, int viewType )
{
    m_event = eventNumber;
    m_view  = viewType;

    m_nHits = 0;
    m_nClusters = 0;
    m_nSeedClusters = 0;

    m_hitID.clear();
    m_hitClusterID.clear();
    m_hitPosX.clear();
    m_hitPosY.clear();
    m_hitPosZ.clear();
    m_hitEnergy.clear();

    m_clusterID.clear();
    m_clusterIsSeed.clear();
    m_clusterDirectionInZ.clear();
    m_clusterLayers.clear();
    m_clusterStartPosX.clear();
    m_clusterStartPosY.clear();
    m_clusterStartPosZ.clear();
    m_clusterStartDirX.clear();
    m_clusterStartDirY.clear();
    m_clusterStartDirZ.clear();
    m_clusterEndPosX.clear();
    m_clusterEndPosY.clear();
    m_clusterEndPosZ.clear();
    m_clusterEndDirX.clear();
    m_clusterEndDirY.clear();
    m_clusterEndDirZ.clear();
    m_clusterEnergy.clear();
    m_clusterLength.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NtupleWritingAlgorithm::NtupleWriter::WriteEntry()
{    
    std::string treeName = "";

    if ( m_view == TPC_VIEW_U )      { treeName = m_treeNameU; m_treeExistsU = true; }
    else if ( m_view == TPC_VIEW_V ) { treeName = m_treeNameV; m_treeExistsV = true; }
    else if ( m_view == TPC_VIEW_W ) { treeName = m_treeNameW; m_treeExistsW = true; }
    else return;

    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "event", m_event));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "view", m_view));

    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "nHits", m_nHits)); 
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "hitID", &m_hitID)); 
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "hitClusterID", &m_hitClusterID)); 
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "hitPosX", &m_hitPosX));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "hitPosY", &m_hitPosY));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "hitPosZ", &m_hitPosZ));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "hitEnergy", &m_hitEnergy));

    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "nClusters", m_nClusters));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "nSeedClusters", m_nSeedClusters));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterID", &m_clusterID));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterIsSeed", &m_clusterIsSeed));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterDirectionInZ", &m_clusterDirectionInZ));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterLayers", &m_clusterLayers));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterLength", &m_clusterLength));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterEnergy", &m_clusterEnergy));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterStartPosX", &m_clusterStartPosX));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterStartPosY", &m_clusterStartPosY));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterStartPosZ", &m_clusterStartPosZ));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterStartDirX", &m_clusterStartDirX));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterStartDirY", &m_clusterStartDirY));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterStartDirZ", &m_clusterStartDirZ));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterEndPosX", &m_clusterEndPosX));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterEndPosY", &m_clusterEndPosY));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterEndPosZ", &m_clusterEndPosZ));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterEndDirX", &m_clusterEndDirX));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterEndDirY", &m_clusterEndDirY));
    PANDORA_MONITORING_API(SetTreeVariable(treeName.c_str(), "clusterEndDirZ", &m_clusterEndDirZ));

    PANDORA_MONITORING_API(FillTree(treeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NtupleWritingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SeedClusterListName", m_seedClusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NonSeedClusterListName", m_nonSeedClusterListName));

    m_outputFileName = "LArRecoResults.root";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputFileName", m_outputFileName));

    m_outputTreeName = "LArRecoTree";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputTreeName", m_outputTreeName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
