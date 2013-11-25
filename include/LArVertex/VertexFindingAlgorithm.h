/**
 *  @file   LArContent/include/LArVertex/VertexFindingAlgorithm.h
 * 
 *  @brief  Header file for the cluster creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_FINDING_ALGORITHM_H
#define LAR_VERTEX_FINDING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArPointingCluster.h"

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



    class LArVertexCandidate
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  
         */
      LArVertexCandidate(const pandora::CartesianVector position, const pandora::CartesianVector direction,
                         const float energy, const float momentum);

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

        /**
         *  @brief  Get the energy
         * 
         *  @return the energy
         */
        const float &GetEnergy() const;

        /**
         *  @brief  Get the momentum
         * 
         *  @return the momentum
         */
        const float &GetMomentum() const;

    private:

        pandora::CartesianVector  m_position;
        pandora::CartesianVector  m_direction;
        float                     m_energy;
        float                     m_momentum;

    };

    typedef std::map<const LArVertexCandidate*,float> VertexFigureOfMeritMap;


    typedef std::vector<const LArVertexCandidate> VertexFigureOfMeritList;


 private:


    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StatusCode SetVertex(const pandora::CartesianVector& eventVertex, std::string vertexName);
    
    
    
    
    void GetListOfCleanClusters( const pandora::ClusterList* const pClusterList, pandora::ClusterVector &clusterVector);  

    void GetListOfCleanPointingClusters( const pandora::ClusterList* const pClusterList, LArPointingClusterMap& pointingClusterMap );

    void GetListOfCleanVertexClusters( const LArPointingClusterMap& pointingClusterMap, LArPointingClusterVertexList& cleanClusterList );




    void GetListOfCandidateVertexPositions( const LArPointingClusterMap& pointingClusterMap, const LArPointingClusterVertexList& cleanClusterList,
                                            pandora::CartesianPointList& candidateVertexList );
 
    void GetListOfCleanVertexPositions( const LArPointingClusterMap& pointingClusterMap, const pandora::CartesianPointList& inputVertexList,
                                        pandora::CartesianPointList& outputVertexList );

    void GetListOfMatchedVertexPositions( const LArPointingClusterMap& pointingClusterMapU, const LArPointingClusterMap& pointingClusterMapV, 
                                          const LArPointingClusterMap& pointingClusterMapW,
                                          const pandora::CartesianPointList& inputListU, const pandora::CartesianPointList& inputListV, 
                                          const pandora::CartesianPointList& inputListW, 
                                          pandora::CartesianPointList& outputListU, pandora::CartesianPointList& outputListV, 
                                          pandora::CartesianPointList& outputListW );

    void GetListOfMatchedVertexPositions2D( const pandora::HitType view1, const pandora::HitType view2, 
                                            const pandora::CartesianPointList& inputList1, const pandora::CartesianPointList& inputList2,
                                            pandora::CartesianPointList& projectedList3 );

    void GetListOfMatchedVertexPositions3D( const pandora::CartesianPointList& inputListU, const pandora::CartesianPointList& inputListV,
                                            const pandora::CartesianPointList& inputListW, pandora::CartesianPointList& outputListU, 
                                            pandora::CartesianPointList& outputListV, pandora::CartesianPointList& outputListW );

    void GetReducedListOfVertexPositions( const pandora::CartesianPointList& inputList, pandora::CartesianPointList& outputList );


    void CollectMarkers( const LArPointingClusterVertexList& inputList, pandora::CartesianPointList& outputList );

    void CollectClusters( const LArPointingClusterVertexList& inputList, pandora::ClusterList& outputList );

    void CollectClusters( const LArPointingClusterMap& inputList, pandora::ClusterList& outputList );




    void ProcessView1D( const LArPointingClusterMap& pointingClusterMap, const pandora::CartesianPointList& inputVertexList,
                            VertexFigureOfMeritMap& outputFigureOfMeritMap );

    void ProcessVertex1D( const VertexFigureOfMeritMap figureOfMeritMap, pandora::CartesianVector& bestVertex );

    void ProcessVertex3D( const VertexFigureOfMeritMap figureOfMeritMapU, const VertexFigureOfMeritMap figureOfMeritMapV, 
                          const VertexFigureOfMeritMap figureOfMeritMapW, pandora::CartesianVector& bestVertexU, 
                          pandora::CartesianVector& bestVertexV, pandora::CartesianVector& bestVertexW );

    void CleanUp( VertexFigureOfMeritMap& figureOfMeritMap );   



   
    void RunFastReconstruction( const LArPointingClusterMap& inputMap,
                                const pandora::CartesianVector& seedPosition, const pandora::CartesianVector& seedDirection, 
                                LArPointingClusterVertexList& outputList, float& outputEnergy, pandora::CartesianVector& outputMomentum );

    void GetEnergyAndMomentum( const LArPointingClusterVertexList& clusterList, const pandora::CartesianVector& thisVertex, 
                               float& outputEnergy, pandora::CartesianVector& outputMomentum );

    void GetEnergyAndMomentum( const LArPointingCluster::Vertex& thisCluster, const pandora::CartesianVector& thisVertex, 
                               float& outputEnergy, pandora::CartesianVector& outputMomentum );

    void GetIntersection( const LArPointingCluster::Vertex& firstCluster, const LArPointingCluster::Vertex& secondCluster,
                          pandora::CartesianVector& intersectVertex, pandora::CartesianVector& intersectDirection, bool& isPhysical );

    void GetIntersection( const LArPointingCluster::Vertex& vertexCluster, const pandora::Cluster* pTargetCluster, 
                          pandora::CartesianVector& intersectPosition, float& intersectDisplacement, bool& isPhysical );

    bool IsStrongAssociation( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );

    bool IsWeakAssociation( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );

    bool IsDaughterDiverging( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );

    bool IsNode( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );

    bool IsEmission( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );

    bool IsStrongSecondary( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );

    bool IsWeakSecondary( const LArPointingCluster::Vertex& parentCluster, const LArPointingCluster::Vertex& daughterCluster );

    bool IsPrimary( const LArPointingCluster::Vertex& clusterTrajectory, const pandora::CartesianVector& targetVertex, const pandora::CartesianVector& targetDirection );

    bool IsStrongSecondary( const LArPointingCluster::Vertex& parentCluster, const pandora::CartesianVector& targetVertex, const pandora::CartesianVector& targetDirection );

    bool IsWeakSecondary( const LArPointingCluster::Vertex& parentCluster, const pandora::CartesianVector& targetVertex, const pandora::CartesianVector& targetDirection );

    bool IsAdjacent( const LArPointingCluster::Vertex& clusterTrajectory, const pandora::CartesianVector& targetVertex );

    bool IsDownstream( const LArPointingCluster::Vertex& clusterTrajectory, const pandora::CartesianVector& targetVertex );

    bool IsNode( const pandora::CartesianVector& clusterVertex, const pandora::CartesianVector& clusterDirection, const pandora::CartesianVector& targetVertex, const pandora::CartesianVector& targetDirection ) const;

    bool IsEmission( const pandora::CartesianVector& clusterVertex, const pandora::CartesianVector& clusterDirection, const pandora::CartesianVector& targetVertex, const pandora::CartesianVector& targetDirection ) const;

    bool IsNode( const pandora::CartesianVector& clusterVertex, const pandora::CartesianVector& clusterDirection, const pandora::CartesianVector& targetVertex ) const;
    bool IsEmission( const pandora::CartesianVector& clusterVertex, const pandora::CartesianVector& clusterDirection, const pandora::CartesianVector& targetVertex ) const;
    
    bool IsConsistentWithBeamDirectionW( const pandora::CartesianVector& thisMomentumW ) const;

    bool IsConsistentWithBeamDirectionUV( const pandora::CartesianVector& thisMomentumU, const pandora::CartesianVector& thisMomentumV ) const;

    

    std::string     m_vertexNameU;       ///<
    std::string     m_vertexNameV;       ///<
    std::string     m_vertexNameW;       ///<
    std::string     m_vertexName3D;      ///<


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

inline VertexFindingAlgorithm::LArVertexCandidate::LArVertexCandidate(const pandora::CartesianVector position, const pandora::CartesianVector direction,
        const float energy, const float momentum) :
    m_position(position),
    m_direction(direction),
    m_energy(energy),
    m_momentum(momentum)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &VertexFindingAlgorithm::LArVertexCandidate::GetPosition() const
{
    return m_position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &VertexFindingAlgorithm::LArVertexCandidate::GetDirection() const
{
    return m_direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const float &VertexFindingAlgorithm::LArVertexCandidate::GetEnergy() const
{
  return m_energy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const float &VertexFindingAlgorithm::LArVertexCandidate::GetMomentum() const
{
  return m_momentum;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *VertexFindingAlgorithm::Factory::CreateAlgorithm() const
{
    return new VertexFindingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_VERTEX_FINDING_ALGORITHM_H
