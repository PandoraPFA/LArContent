/**
 *  @file   larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h
 * 
 *  @brief  Header file for the vertex selection base algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_SELECTION_BASE_ALGORITHM_H
#define LAR_VERTEX_SELECTION_BASE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

template<typename, unsigned int> class KDTreeLinkerAlgo;
template<typename, unsigned int> class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  VertexSelectionBaseAlgorithm class
 */
class VertexSelectionBaseAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexSelectionBaseAlgorithm();

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

    /**
     *  @brief Beam constants class
     */
    class BeamConstants
    {
    public:
        /**
         *  @brief  Get the min z coordinate
         * 
         *  @return the min z coordinate
         */
        float GetMinZCoordinate() const;

        /**
         *  @brief  Get the decay constant
         * 
         *  @return the decay constant
         */
        float GetDecayConstant() const;

        /**
         *  @brief  Set the beam constants
         * 
         *  @param  minZCoordinate the min z coordinate
         *  @param  decayConstant the decay constant
         */
        void SetConstants(const float minZCoordinate, const float decayConstant);

    private:
        pandora::InputFloat     m_minZCoordinate;   ///< The min z coordinate
        pandora::InputFloat     m_decayConstant;    ///< The decay constant
    };

protected:
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;
    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;

    /**
     *  @brief  Filter the input list of vertices to obtain a reduced number of vertex candidates
     * 
     *  @param  pInputVertexList the address of the input vertex list
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     *  @param  filteredVertices to receive the filtered vertex list
     */
    virtual void FilterVertexList(const pandora::VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
        HitKDTree2D &kdTreeW, pandora::VertexVector &filteredVertices) const;

    /**
     *  @brief  Get the beam score constants for a provided list of candidate vertices
     * 
     *  @param  vertexVector the vertex vector
     *  @param  beamConstants to receive the beam constants
     */
    virtual void GetBeamConstants(const pandora::VertexVector &vertexVector, BeamConstants &beamConstants) const;

    /**
     *  @brief  Get the vertex score list for a provided list of candidate vertices
     * 
     *  @param  vertexVector the vertex vector
     *  @param  beamConstants the beam constants
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     *  @param  vertexScoreList to receive the vertex score list
     */
    virtual void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const = 0;

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Initialize kd trees with details of hits in algorithm-configured cluster lists
     * 
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     */
    void InitializeKDTrees(HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const;

    /**
     *  @brief  Initialize a kd trees with details of hits in a named list
     * 
     *  @param  caloHitListName the calo hit list name
     *  @param  kdTree the kd tree
     */
    void InitializeKDTree(const std::string &caloHitListName, HitKDTree2D &kdTree) const;

    /**
     *  @brief  Whether the vertex lies on a hit in the specified view
     * 
     *  @param  pVertex the address of the vertex
     *  @param  hitType the relevant hit type
     *  @param  kdTree the relevant kd tree
     * 
     *  @return boolean
     */
    bool IsVertexOnHit(const pandora::Vertex *const pVertex, const pandora::HitType hitType, HitKDTree2D &kdTree) const;

    /**
     *  @brief  Whether the vertex lies in a registered gap
     * 
     *  @param  pVertex the address of the vertex
     *  @param  hitType the relevant hit type
     * 
     *  @return boolean
     */
    bool IsVertexInGap(const pandora::Vertex *const pVertex, const pandora::HitType hitType) const;

    /**
     *  @brief  From the top-scoring candidate vertices, select a subset for further investigation
     * 
     *  @param  vertexScoreList the vertex score list
     *  @param  selectedVertexList to receive the selected vertex list
     */
    void SelectTopScoreVertices(VertexScoreList &vertexScoreList, pandora::VertexList &selectedVertexList) const;

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
     *  @brief  Sort vertices by increasing z position
     * 
     *  @param  pLhs address of the lhs vertex
     *  @param  pRhs address of the rhs vertex
     * 
     *  @return whether lhs should precedes rhs
     */
    static bool SortByVertexZPosition(const pandora::Vertex *const pLhs, const pandora::Vertex *const pRhs);

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
// TODO
    bool            m_beamMode;                     ///< Whether to run in beam mode, assuming neutrinos travel in positive z-direction

private:
    std::string     m_inputCaloHitListNameU;        ///< The name of the view U calo hit list
    std::string     m_inputCaloHitListNameV;        ///< The name of the view V calo hit list
    std::string     m_inputCaloHitListNameW;        ///< The name of the view W calo hit list
    std::string     m_outputVertexListName;         ///< The name under which to save the output vertex list

    bool            m_replaceCurrentVertexList;     ///< Whether to replace the current vertex list with the output list

    float           m_nDecayLengthsInZSpan;         ///< The number of score decay lengths to use over the course of the vertex z-span
   
    bool            m_selectSingleVertex;           ///< Whether to make a final decision and select just one vertex candidate
    unsigned int    m_maxTopScoreSelections;        ///< Max number of top-scoring vertex candidate to select for output

    float           m_maxOnHitDisplacement;         ///< Max hit-vertex displacement for declaring vertex to lie on a hit in each view

    float           m_minCandidateDisplacement;     ///< Ignore other top-scoring candidates located in close proximity to original
    float           m_minCandidateScoreFraction;    ///< Ignore other top-scoring candidates with score less than a fraction of original

    bool            m_useDetectorGaps;              ///< Whether to account for registered detector gaps in vertex selection
    float           m_gapTolerance;                 ///< The tolerance to use when querying whether a sampling point is in a gap, units cm

    bool            m_isEmptyViewAcceptable;        ///< Whether views entirely empty of hits are classed as 'acceptable' for candidate filtration
    unsigned int    m_minVertexAcceptableViews;     ///< The minimum number of views in which a candidate must sit on/near a hit or in a gap (or view can be empty)
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline VertexSelectionBaseAlgorithm::VertexScore::VertexScore(const pandora::Vertex *const pVertex, const float score) :
    m_pVertex(pVertex),
    m_score(score)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Vertex *VertexSelectionBaseAlgorithm::VertexScore::GetVertex() const
{
    return m_pVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexSelectionBaseAlgorithm::VertexScore::GetScore() const
{
    return m_score;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool VertexSelectionBaseAlgorithm::VertexScore::operator< (const VertexScore &rhs) const
{
    return (this->GetScore() > rhs.GetScore());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexSelectionBaseAlgorithm::BeamConstants::GetMinZCoordinate() const
{
    return m_minZCoordinate.Get();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexSelectionBaseAlgorithm::BeamConstants::GetDecayConstant() const
{
    return m_decayConstant.Get();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void VertexSelectionBaseAlgorithm::BeamConstants::SetConstants(const float minZCoordinate, const float decayConstant)
{
    m_minZCoordinate = minZCoordinate;
    m_decayConstant = decayConstant;
}

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_SELECTION_BASE_ALGORITHM_H
