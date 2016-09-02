/**
 *  @file   LArContent/include/LArVertex/VertexSelectionAlgorithm.h
 * 
 *  @brief  Header file for the vertex selection algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_SELECTION_ALGORITHM_H
#define LAR_VERTEX_SELECTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include <list>

#include "LArVertex/VertexClusteringTool.h"
#include "LArVertex/VertexScoringTool.h"

namespace lar_content
{

template<typename, unsigned int> class KDTreeLinkerAlgo;
template<typename, unsigned int> class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  VertexSelectionAlgorithm::Algorithm class
 */
class VertexSelectionAlgorithm : public pandora::Algorithm
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

    /**
     *  @brief  Default constructor
     */
    VertexSelectionAlgorithm();

private:
    /**
     *  @brief  VertexScore class
     */
    class VertexScore
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pVertex the address of the vertex
         *  @param  score the score
         */
        VertexScore(const pandora::Vertex *const pVertex, const float score);

        /**
         *  @brief  Get the address of the vertex
         * 
         *  @return the address of the vertex
         */
        const pandora::Vertex *GetVertex() const;

        /**
         *  @brief  Get the score
         * 
         *  @return the score
         */
        float GetScore() const;

        /**
         *  @brief  Multiply the score
         * 
         *  @return void
         */
        void MultiplyScore(const float &multiplier);

        /**
         *  @brief  operator<
         * 
         *  @param  rhs the value for comparison
         * 
         *  @return boolean
         */
        bool operator< (const VertexScore &rhs) const;

    private:
        const pandora::Vertex  *m_pVertex;          ///< The address of the vertex
        float                   m_score;            ///< The score
    };
    
    typedef std::vector<VertexScore> VertexScoreList;

    pandora::StatusCode Run();

    /**
     *  @brief  From the top-scoring candidate vertices, select a subset for further investigation
     * 
     *  @param  vertexScoreList the vertex score list
     *  @param  selectedVertexList to receive the selected vertex list
     */
    void SelectTopScoreVertices(std::vector<VertexScoringTool::VertexScoreList> &scoredClusterCollection, pandora::VertexList &selectedVertexList) const;

    /**
     *  @brief  Whether to accept a candidate vertex, based on its spatial position in relation to other selected candidates
     * 
     *  @param  pVertex the address of the vertex
     *  @param  selectedVertexList the selected vertex list
     * 
     *  @return boolean
     */
    bool AcceptVertexLocation(const pandora::Vertex *const pVertex, const pandora::VertexList &selectedVertexList) const;

    /**
     *  @brief  Finds the top 5 vertex candidates from a scored vertex cluster collection and an energy vertex score list
     * 
     *  @param  scoredClusterCollection the scored cluster collection by reference
     *  @param  energyVertexScoreList the energy vertex score list by reference
     *  @param  topNVertexScoreList the output vertex score list containing the top 5 vertices, by reference
     *  
     */
    void FindTopNVertices(std::vector<VertexScoringTool::VertexScoreList> &scoredClusterCollection, VertexScoringTool::VertexScoreList &energyVertexScoreList, VertexScoringTool::VertexScoreList &topNVertexScoreList);

    /**
     *  @brief  Stores the top 5 vertices and outputs this list so it can be propagated to the next algorithm
     * 
     *  @param  topNVertexScoreList the input vertex score list containing the top 5 vertices that is to be stored
     *  
     */
    void StoreTopNInformation(VertexScoringTool::VertexScoreList &topNVertexScoreList);
    
    /**
     *  @brief  Stores all vertex information and outputs this list so it can be propagated to the next algorithm
     * 
     *  @param  pTopologyVertexList the address of the topology vertex list
     *  @param  selectedVertexList the list containing the single selected output vertex
     *  @param  pEnergyVertexList the address of the energy vertex list
     *  
     */
    void StoreTopAllInformation(const pandora::VertexList* pTopologyVertexList, pandora::VertexList selectedVertexList, const pandora::VertexList* pEnergyVertexList);
    
    
    void GetClusters(pandora::ClusterList &clusterListU, pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW);
    const pandora::Cluster* GetLongestCluster(pandora::ClusterList &clusterList);
    void FilterTopNVertices(VertexScoringTool::VertexScoreList &topNVertexScoreList, VertexScoringTool::VertexScoreList &filteredTopNVertexScoreList);
    
    //--------------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    VertexClusteringTool *m_pVertexClusteringTool; ///< Vertex Tool 1
    VertexScoringTool    *m_pVertexScoringTool;    ///< Vertex Tool 2

    std::string     m_outputVertexListName;         ///< The name under which to save the output vertex list of selected vertices

    std::string     m_topologyVertexListName;       ///< The name under which to save the output topology vertex list
    std::string     m_energyVertexListName;         ///< The name under which to save the output energy vertex list

    std::string     m_topNVertexListName;           ///< The name under which to save the output top 5vertex list
    std::string     m_allOtherVertexListName;       ///< The name under which to save the output vertex list containing all created vertices
    
    bool            m_replaceCurrentVertexList;     ///< Whether to replace the current vertex list with the output list

    bool            m_beamMode;                     ///< Whether to run in beam mode, assuming neutrinos travel in positive z-direction
    float           m_nDecayLengthsInZSpan;         ///< The number of score decay lengths to use over the course of the vertex z-span
    float           m_kappa;                        ///< Hit-deweighting offset, of form: weight = 1 / sqrt(distance + kappa), units cm

    bool            m_selectSingleVertex;           ///< Whether to make a final decision and select just one vertex candidate
    unsigned int    m_maxTopScoreSelections;        ///< Max number of top-scoring vertex candidate to select for 

    float           m_minCandidateDisplacement;     ///< Ignore other top-scoring candidates located in close proximity to original
    float           m_minCandidateScoreFraction;    ///< Ignore other top-scoring candidates with score less than a fraction of original

    bool            m_useDetectorGaps;              ///< Whether to account for registered detector gaps in vertex selection
    float           m_gapTolerance;                 ///< The tolerance to use when querying whether a sampling point is in a gap, units cm

    bool            m_isEmptyViewAcceptable;        ///< Whether views entirely empty of hits are classed as 'acceptable' for candidate filtration
    
    bool            m_enableClustering;
    
    bool            m_directionFilter;
    bool            m_beamWeightFilter;
    
    unsigned int    m_nVerticesToSelect;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *VertexSelectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new VertexSelectionAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline VertexSelectionAlgorithm::VertexScore::VertexScore(const pandora::Vertex *const pVertex, const float score) :
    m_pVertex(pVertex),
    m_score(score)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Vertex *VertexSelectionAlgorithm::VertexScore::GetVertex() const
{
    return m_pVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexSelectionAlgorithm::VertexScore::GetScore() const
{
    return m_score;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void VertexSelectionAlgorithm::VertexScore::MultiplyScore(const float &multiplier)
{
    m_score *= multiplier;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool VertexSelectionAlgorithm::VertexScore::operator< (const VertexScore &rhs) const
{
    return (this->GetScore() > rhs.GetScore());
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_SELECTION_ALGORITHM_H
