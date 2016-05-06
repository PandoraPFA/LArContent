/**
 *  @file   LArContent/LArVertex/VertexSelectionAlgorithm.h
 * 
 *  @brief  Header file for the vertex selection algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_SELECTION_ALGORITHM_H
#define LAR_VERTEX_SELECTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

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

    /**
     *  @brief Kernel estimate class
     */
    class KernelEstimate
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  sigma the width associated with the kernel estimate
         */
        KernelEstimate(const float sigma);

        /**
         *  @brief  Sample the parameterised distribution at a specified x coordinate
         * 
         *  @param  x the position at which to sample
         * 
         *  @return the sample value
         */
        float Sample(const float x) const;

        typedef std::multimap<float, float> ContributionList;   ///< Map from x coord to weight, ATTN avoid map.find, etc. with float key

        /**
         *  @brief  Get the contribution list
         * 
         *  @return the contribution list
         */
        const ContributionList &GetContributionList() const;

        /**
         *  @brief  Get the assigned width
         * 
         *  @return the assigned width
         */
        float GetSigma() const;

        /**
         *  @brief  Add a contribution to the distribution
         * 
         *  @param  x the position
         *  @param  weight the weight
         */
        void AddContribution(const float x, const float weight);

    private:
        ContributionList            m_contributionList;         ///< The contribution list
        const float                 m_sigma;                    ///< The assigned width
    };

    typedef KDTreeLinkerAlgo<const pandora::CaloHit*, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit*, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

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
     *  @brief  Filter the input list of vertices to obtain a reduced number of vertex candidates
     * 
     *  @param  pInputVertexList the address of the input vertex list
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     *  @param  filteredVertexList to receive the filtered vertex score list
     */
    void FilterVertexList(const pandora::VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
        HitKDTree2D &kdTreeW, pandora::VertexList &filteredVertexList) const;

    /**
     *  @brief  Get the beam score constants for a provided list of candidate vertices
     * 
     *  @param  vertexList the vertex list
     *  @param  beamConstants to receive the beam constants
     */
    void GetBeamConstants(const pandora::VertexList &vertexList, BeamConstants &beamConstants) const;

    /**
     *  @brief  Get the vertex score list for a provided list of candidate vertices
     * 
     *  @param  vertexList the vertex list
     *  @param  beamConstants the beam constants
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     *  @param  vertexScoreList to receive the vertex score list
     */
    void GetVertexScoreList(const pandora::VertexList &vertexList, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const;

    /**
     *  @brief  Get the score for a trio of kernel estimations, using fast histogram approach
     * 
     *  @param  kernelEstimateU the kernel estimate for the u view
     *  @param  kernelEstimateV the kernel estimate for the v view
     *  @param  kernelEstimateW the kernel estimate for the w view
     * 
     *  @return the fast score
     */
    float GetFastScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV, const KernelEstimate &kernelEstimateW) const;

    /**
     *  @brief  Get the score for a trio of kernel estimations, using kernel density estimation but with reduced (binned) sampling
     * 
     *  @param  kernelEstimateU the kernel estimate for the u view
     *  @param  kernelEstimateV the kernel estimate for the v view
     *  @param  kernelEstimateW the kernel estimate for the w view
     * 
     *  @return the midway score
     */
    float GetMidwayScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV, const KernelEstimate &kernelEstimateW) const;

    /**
     *  @brief  Get the score for a trio of kernel estimations, using kernel density estimation and full hit-by-hit sampling
     * 
     *  @param  kernelEstimateU the kernel estimate for the u view
     *  @param  kernelEstimateV the kernel estimate for the v view
     *  @param  kernelEstimateW the kernel estimate for the w view
     * 
     *  @return the full score
     */
    float GetFullScore(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV, const KernelEstimate &kernelEstimateW) const;

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
     *  @brief  Use hits in clusters (in the provided kd tree) to fill a provided kernel estimate with hit-vertex relationship information
     * 
     *  @param  pVertex the address of the vertex
     *  @param  hitType the relevant hit type
     *  @param  kdTree the relevant kd tree
     *  @param  kernelEstimate to receive the populated kernel estimate
     */
    void FillKernelEstimate(const pandora::Vertex *const pVertex, const pandora::HitType hitType, HitKDTree2D &kdTree, KernelEstimate &kernelEstimate) const;

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
     *  @brief  Fast estimate of std::atan2 function. Rather coarse (max |error| > 0.01) but should suffice for this use-case.
     * 
     *  @param  y the y coordinate
     *  @param  x the x coordinate
     * 
     *  @return estimate of std::atan2
     */
    float atan2Fast(const float y, const float x) const;

    /**
     *  @brief  Sort vertices by increasing z position
     * 
     *  @param  pLhs address of the lhs vertex
     *  @param  pRhs address of the rhs vertex
     * 
     *  @return whether lhs should precedes rhs
     */
    static bool SortByVertexZPosition(const pandora::Vertex *const pLhs, const pandora::Vertex *const pRhs);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_inputCaloHitListNameU;        ///< The name of the view U calo hit list
    std::string     m_inputCaloHitListNameV;        ///< The name of the view V calo hit list
    std::string     m_inputCaloHitListNameW;        ///< The name of the view W calo hit list
    std::string     m_outputVertexListName;         ///< The name under which to save the output vertex list
    bool            m_replaceCurrentVertexList;     ///< Whether to replace the current vertex list with the output list

    bool            m_fastScoreCheck;               ///< Whether to use the fast histogram based score to selectively avoid calling full or midway scores
    bool            m_fastScoreOnly;                ///< Whether to use the fast histogram based score only
    bool            m_fullScore;                    ///< Whether to use the full kernel density estimation score, as opposed to the midway score

    bool            m_beamMode;                     ///< Whether to run in beam mode, assuming neutrinos travel in positive z-direction
    float           m_nDecayLengthsInZSpan;         ///< The number of score decay lengths to use over the course of the vertex z-span
    float           m_kappa;                        ///< Hit-deweighting offset, of form: weight = 1 / sqrt(distance + kappa), units cm

    bool            m_selectSingleVertex;           ///< Whether to make a final decision and select just one vertex candidate
    unsigned int    m_maxTopScoreSelections;        ///< Max number of top-scoring vertex candidate to select for output

    float           m_kernelEstimateSigma;          ///< The Gaussian width to use for kernel estimation

    float           m_minFastScoreFraction;         ///< Fast score must be at least this fraction of best fast score to calculate full score
    unsigned int    m_fastHistogramNPhiBins;        ///< Number of bins to use for fast score histograms
    float           m_fastHistogramPhiMin;          ///< Min value for fast score histograms
    float           m_fastHistogramPhiMax;          ///< Max value for fast score histograms

    float           m_maxOnHitDisplacement;         ///< Max hit-vertex displacement for declaring vertex to lie on a hit in each view
    float           m_maxHitVertexDisplacement1D;   ///< Max hit-vertex displacement in *any one dimension* for contribution to kernel estimation

    float           m_minCandidateDisplacement;     ///< Ignore other top-scoring candidates located in close proximity to original
    float           m_minCandidateScoreFraction;    ///< Ignore other top-scoring candidates with score less than a fraction of original

    bool            m_enableFolding;                ///< Whether to enable folding of -pi -> +pi phi distribution into 0 -> +pi region only

    bool            m_useDetectorGaps;              ///< Whether to account for registered detector gaps in vertex selection
    float           m_gapTolerance;                 ///< The tolerance to use when querying whether a sampling point is in a gap, units cm

    bool            m_isEmptyViewAcceptable;        ///< Whether views entirely empty of hits are classed as 'acceptable' for candidate filtration
    unsigned int    m_minVertexAcceptableViews;     ///< The minimum number of views in which a candidate must sit on/near a hit or in a gap (or view can be empty)
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
//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexSelectionAlgorithm::BeamConstants::GetMinZCoordinate() const
{
    return m_minZCoordinate.Get();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexSelectionAlgorithm::BeamConstants::GetDecayConstant() const
{
    return m_decayConstant.Get();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void VertexSelectionAlgorithm::BeamConstants::SetConstants(const float minZCoordinate, const float decayConstant)
{
    m_minZCoordinate = minZCoordinate;
    m_decayConstant = decayConstant;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline VertexSelectionAlgorithm::KernelEstimate::KernelEstimate(const float sigma) :
    m_sigma(sigma)
{
    if (m_sigma < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const VertexSelectionAlgorithm::KernelEstimate::ContributionList &VertexSelectionAlgorithm::KernelEstimate::GetContributionList() const
{
    return m_contributionList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexSelectionAlgorithm::KernelEstimate::GetSigma() const
{
    return m_sigma;
}

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_SELECTION_ALGORITHM_H
