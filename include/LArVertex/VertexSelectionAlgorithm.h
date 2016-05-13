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
    
    /**
     *  @brief  Default destructor
     */
    ~VertexSelectionAlgorithm();

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
    void SelectTopScoreVertices(VertexScoringTool::VertexScoreList &vertexScoreList, pandora::VertexList &selectedVertexList) const;

    /**
     *  @brief  Whether to accept a candidate vertex, based on its spatial position in relation to other selected candidates
     * 
     *  @param  pVertex the address of the vertex
     *  @param  selectedVertexList the selected vertex list
     * 
     *  @return boolean
     */
    bool AcceptVertexLocation(const pandora::Vertex *const pVertex, const pandora::VertexList &selectedVertexList) const;
    
    //------------------------------------------------------------------------------------------------------------------------------

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    //VertexClusteringTool *m_pVertexClusteringTool; ///< Vertex Tool 1
    VertexScoringTool    *m_pVertexScoringTool;    ///< Vertex Tool 2

    std::string     m_outputVertexListName;         ///< The name under which to save the output vertex list
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
    unsigned int    m_minVertexAcceptableViews;     ///< The minimum number of views in which a candidate must sit on/near a hit or in a gap (or view can be empty)
    
    bool            m_writeToTree;
    std::string     m_treeName;                 ///< Name of output tree
    std::string     m_fileName;                 ///< Name of output file
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

inline bool VertexSelectionAlgorithm::VertexScore::operator< (const VertexScore &rhs) const
{
    return (this->GetScore() > rhs.GetScore());
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_SELECTION_ALGORITHM_H
