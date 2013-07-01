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



    class LArPointingVertex
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  
         */
        LArPointingVertex(const pandora::CartesianVector position, const pandora::CartesianVector direction);

        /**
         *  @brief  Get the vertex position
         * 
         *  @return the vertex position
         */
        const pandora::CartesianVector &GetPosition() const;


        /**
         *  @brief  Get the direction
         * 
         *  @return the direction
         */
        const pandora::CartesianVector &GetDirection() const;


    private:
        pandora::CartesianVector m_position;
        pandora::CartesianVector m_direction;

    };

    typedef std::map<const LArPointingVertex*,float> VertexFigureOfMeritList;

 private:


    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StatusCode SetTrueVertex();
    pandora::StatusCode SetVertex(const pandora::CartesianVector& eventVertex, std::string vertexName);
    
    
    

    


    
    
    void GetListOfCleanClusters(const pandora::ClusterList* const pClusterList, pandora::ClusterVector &clusterVector);  

    void GetListOfCleanPointingClusters(const pandora::ClusterList* const pClusterList, LArPointingClusterMap& pointingClusterMap );

    void GetListOfCleanVertexClusters( const LArPointingClusterMap& pointingClusterMap, LArPointingClusterVertexList& outputList );



    void GetListOfCandidateVertexPositions( const LArPointingClusterVertexList& pointingVertexList, pandora::CartesianPointList& candidateVertexList );



    void CollectMarkers( const LArPointingClusterVertexList& inputList, pandora::CartesianPointList& outputList );

    void CollectClusters( const LArPointingClusterVertexList& inputList, pandora::ClusterList& outputList );







    void CleanUp( VertexFigureOfMeritList& figureOfMeritList );     


    void ProcessSingleView( const LArPointingClusterMap& inputMap, const LArPointingClusterVertexList& inputList,
                            VertexFigureOfMeritList& outputFigureOfMeritList );

    void ProcessSingleVertex( const pandora::ClusterList* const pClusterList, const VertexFigureOfMeritList figureOfMeritList, 
                              pandora::CartesianVector& bestVertex );


    void GetFirstPassVertex( const pandora::ClusterList* const pClusterList, pandora::CartesianVector& firstVertex);

    

    void FindPossibleConnectedVertices( const LArPointingClusterMap& inputMap, const LArPointingClusterVertexList& inputList, 
                                        VertexFigureOfMeritList& outputFigureOfMeritList );

    void FindPossibleDisplacedVertices( const LArPointingClusterMap& inputMap, const LArPointingClusterVertexList& inputList, 
                                        VertexFigureOfMeritList& outputFigureOfMeritList );

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
    

    


    std::string     m_vertexNameU;       ///<
    std::string     m_vertexNameV;       ///<
    std::string     m_vertexNameW;       ///<


    std::string     m_clusterListNameU;  ///< 
    std::string     m_clusterListNameV;  ///< 
    std::string     m_clusterListNameW;  ///< 


    bool            m_useTrueVertex;
    bool            m_runBeamMode;



    ///////////////////////////////////////////////////////////////////////////////////
    //
    // OLD CODE LIVES IN v217 
    //
    ///////////////////////////////////////////////////////////////////////////////////


 
};

//------------------------------------------------------------------------------------------------------------------------------------------

VertexFindingAlgorithm::LArPointingVertex::LArPointingVertex(const pandora::CartesianVector position, const pandora::CartesianVector direction) :
  m_position(position), m_direction(direction)
{
  
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::CartesianVector &VertexFindingAlgorithm::LArPointingVertex::GetPosition() const
{
  return m_position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const pandora::CartesianVector &VertexFindingAlgorithm::LArPointingVertex::GetDirection() const
{
  return m_direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *VertexFindingAlgorithm::Factory::CreateAlgorithm() const
{
    return new VertexFindingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_VERTEX_FINDING_ALGORITHM_H
