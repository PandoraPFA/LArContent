/**
 *  @file   PandoraPFANew/Framework/include/Pandora/VertexFindingAlgorithm.h
 * 
 *  @brief  Header file for the cluster creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_FINDING_ALGORITHM_H
#define LAR_VERTEX_FINDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArPointingCluster.h"

#include <vector>

namespace lar
{

 

/**
 *  @brief  VertexFindingAlgorithm class
 */
class VertexFindingAlgorithm : public pandora::Algorithm
{
 public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

 private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StatusCode SetTrueVertex();
    pandora::StatusCode SetVertex(const pandora::CartesianVector& eventVertex);
    pandora::StatusCode GetFirstPassVertex( const pandora::ClusterList* const pClusterList, pandora::CartesianVector& firstVertex);
    

    void GetListOfCleanClusters(const pandora::ClusterList* const pClusterList, pandora::ClusterVector &clusterVector);  

    void GetListOfCleanVertexClusters( const LArPointingClusterVertexList& inputList, LArPointingClusterVertexList& outputList );

    void CollectMarkers(  const LArPointingClusterVertexList& inputList, pandora::CartesianPointList& outputList );
    void CollectClusters( const LArPointingClusterVertexList& inputList, pandora::ClusterList& outputList );

    void FindBestConnectedVertex( const LArPointingClusterMap& inputMap, const LArPointingClusterVertexList& inputList, LArPointingClusterVertexList& outputList,
                                  pandora::CartesianVector& outputVertex, float& outputEnergy, pandora::CartesianVector& outputMomentum, float& outputFoM );

    void FindBestDisplacedVertex( const LArPointingClusterMap& inputMap, const LArPointingClusterVertexList& inputList, LArPointingClusterVertexList& outputList, 
                                  pandora::CartesianVector& outputVertex, float& outputEnergy, pandora::CartesianVector& outputMomentum, float& outputFoM );


    void RunFastReconstruction( const pandora::CartesianVector& seedVertexPosition, const pandora::CartesianVector& seedVertexDirection,
                                const LArPointingClusterMap& inputMap, const LArPointingClusterVertexList& inputList, LArPointingClusterVertexList& outputList, 
                                unsigned int& primaryClusters, float& primaryEnergy, float& outputEnergy, pandora::CartesianVector& outputMomentum );

    void GetEnergyAndMomentum( const LArPointingClusterVertexList& clusterList, const pandora::CartesianVector& thisVertex, 
                               float& outputEnergy, pandora::CartesianVector& outputMomentum );

    void GetEnergyAndMomentum( const LArPointingCluster::Vertex& thisCluster, const pandora::CartesianVector& thisVertex, 
                               float& outputEnergy, pandora::CartesianVector& outputMomentum );

    void GetIntersection( const LArPointingCluster::Vertex& firstCluster, const LArPointingCluster::Vertex& secondCluster,
                          pandora::CartesianVector& intersectVertex, pandora::CartesianVector& intersectDirection, bool& isPhysical );

    bool IsStrongAssociation( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );
    bool IsWeakAssociation( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );
    bool IsDaughterDiverging( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );

    bool IsNode( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );
    bool IsEmission( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );
    bool IsStrongSecondary( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );
    bool IsWeakSecondary( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );

    bool IsPrimary( const LArPointingCluster::Vertex& clusterTrajectory, const pandora::CartesianVector& targetVertex, const pandora::CartesianVector& targetDirection );
    bool IsAdjacent( const LArPointingCluster::Vertex& clusterTrajectory, const pandora::CartesianVector& targetVertex );

    bool IsNode( const pandora::CartesianVector& clusterVertex, const pandora::CartesianVector& clusterDirection, const pandora::CartesianVector& targetVertex, const pandora::CartesianVector& targetDirection ) const;
    bool IsEmission( const pandora::CartesianVector& clusterVertex, const pandora::CartesianVector& clusterDirection, const pandora::CartesianVector& targetVertex, const pandora::CartesianVector& targetDirection ) const;

    bool IsNode( const pandora::CartesianVector& clusterVertex, const pandora::CartesianVector& clusterDirection, const pandora::CartesianVector& targetVertex ) const;
    bool IsEmission( const pandora::CartesianVector& clusterVertex, const pandora::CartesianVector& clusterDirection, const pandora::CartesianVector& targetVertex ) const;
    
    bool IsConsistentWithBeamDirection( const pandora::CartesianVector& thisMomentum ) const;

    float GetEnergy( const pandora::Cluster* const pCluster ) const;
    float GetLength( const pandora::Cluster* const pCluster ) const;

    typedef std::map<const pandora::Cluster*, bool> LArPointingClusterVertexVeto;

    std::string     m_vertexName;       ///< 

    bool            m_useTrueVertex;
    bool            m_runBeamMode;



    ///////////////////////////////////////////////////////////////////////////////////
    //

    class MyCluster {
     public:
      int icluster;
      double vz;
      double vx;
      double ez;
      double ex;
      double pz;
      double px;
      double qz;
      double qx;
      double n;
    };

    class MyClusterHit {
     public:
      int icluster;
      int ilayer;
      int ihit;
      double z;
      double x;
      double n1;
      double n2;
    };

    pandora::StatusCode RunOldMethod(const pandora::ClusterList *const pClusterList); 
    pandora::StatusCode PrepareData(const pandora::ClusterList *const pClusterList); 
    pandora::StatusCode FindVertex();

    pandora::CartesianVector CalcNearestHitToVertex();

    double CalcLikelihoodFunction( double x, double z );        // 
    double CalcLikelihoodFunctionMethod1( double x, double z ); // assume two dimensions!! 
    double CalcLikelihoodFunctionMethod2( double x, double z ); // 

    unsigned int nClusters;
    unsigned int nClusterHits;

    std::vector<MyCluster*> vClusters;
    std::vector<MyClusterHit*> vClusterHits;  

    MyCluster*    NewCluster();
    MyClusterHit* NewClusterHit();

    // 
    ///////////////////////////////////////////////////////////////////////////////////
    

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *VertexFindingAlgorithm::Factory::CreateAlgorithm() const
{
    return new VertexFindingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_VERTEX_FINDING_ALGORITHM_H
