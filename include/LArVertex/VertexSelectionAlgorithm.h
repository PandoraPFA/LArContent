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
        const pandora::Vertex  *m_pVertex;      ///< The address of the vertex
        float                   m_score;        ///< The score
    };

    typedef std::vector<VertexScore> VertexScoreList;

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
         *  @brief  Get the sum of the contribution weights
         * 
         *  @return the sum of the contribution weights
         */
        float GetWeightSum() const;

        /**
         *  @brief  Get the bandwidth
         * 
         *  @return the bandwidth
         */
        float GetBandwidth() const;

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
        float                       m_weightSum;                ///< The sum of the contribution weights
        mutable pandora::InputFloat m_bandWidth;                ///< The bandwidth
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
     *  @brief  Get the figure of merit for a candidate vertex, using the provided parameters
     * 
     *  @param  pVertex the address of the vertex
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     * 
     *  @return the figure of merit
     */
    float GetFigureOfMerit(const pandora::Vertex *const pVertex, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const;

    /**
     *  @brief  Get the figure of merit for a trio of kernel estimations
     * 
     *  @param  kernelEstimateU the kernel estimate for the u view
     *  @param  kernelEstimateV the kernel estimate for the v view
     *  @param  kernelEstimateW the kernel estimate for the w view
     * 
     *  @return the figure of merit
     */
    float GetFigureOfMerit(const KernelEstimate &kernelEstimateU, const KernelEstimate &kernelEstimateV, const KernelEstimate &kernelEstimateW) const;

    /**
     *  @brief  Get the figure of merit contribution for a single kernelEstimate
     * 
     *  @param  kernelEstimate the kernel estimate
     * 
     *  @return the figure of merit
     */
    float GetFigureOfMerit(const KernelEstimate &kernelEstimate) const;

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
     *  @param  selectedVertexScoreList to receive the selected vertex score list
     */
    void SelectTopScoreVertices(const VertexScoreList &vertexScoreList, VertexScoreList &selectedVertexScoreList) const;

    /**
     *  @brief  From the top-scoring candidate vertices, select the best using beam-direction assumptions
     * 
     *  @param  vertexScoreList the vertex score list
     *  @param  minZCoordinate the minimum candidate vertex z-coordinate
     *  @param  decayConstant the decay constant for calculating beam-modified vertex scores
     *  @param  selectedVertexScoreList may receive an additional selected vertex, using beam-direction evidence
     */
    void SelectTopScoreBeamVertices(const VertexScoreList &vertexScoreList, const float minZCoordinate, const float decayConstant,
        VertexScoreList &selectedVertexScoreList) const;

    /**
     *  @brief  Whether to accept a candidate vertex, based on its spatial position in relation to other selected candidates
     * 
     *  @param  pVertex the address of the vertex
     *  @param  vertexScoreList the vertex score list
     * 
     *  @return boolean
     */
    bool AcceptVertexLocation(const pandora::Vertex *const pVertex, const VertexScoreList &vertexScoreList) const;

    /**
     *  @brief  Whether to accept a candidate vertex, based on its score in relation to other selected candidates
     * 
     *  @param  score the vertex score
     *  @param  minScoreFraction the minimum acceptable fraction of the current top-score
     *  @param  vertexScoreList the vertex score list
     * 
     *  @return boolean
     */
    bool AcceptVertexScore(const float score, const float minScoreFraction, const VertexScoreList &vertexScoreList) const;

    /**
     *  @brief  Select the final vertex from the list of chosen top-scoring candidates
     * 
     *  @param  vertexScoreList the list of chosen top-scoring vertex candidates
     *  @param  minZCoordinate the minimum candidate vertex z-coordinate
     *  @param  decayConstant the decay constant for calculating beam-modified vertex scores
     *  @param  finalVertexList to receive the list of final vertices
     */
    void SelectFinalVertices(const VertexScoreList &vertexScoreList, const float minZCoordinate, const float decayConstant,
        pandora::VertexList &finalVertexList) const;

    /**
     *  @brief  Fast estimate of std::atan2 function. Rather coarse (max |error| > 0.01) but should suffice for this use-case.
     * 
     *  @param  y the y coordinate
     *  @param  x the x coordinate
     * 
     *  @return estimate of std::atan2
     */
    float atan2Fast(const float y, const float x) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_inputCaloHitListNameU;        ///< The name of the view U calo hit list
    std::string     m_inputCaloHitListNameV;        ///< The name of the view V calo hit list
    std::string     m_inputCaloHitListNameW;        ///< The name of the view W calo hit list
    std::string     m_outputVertexListName;         ///< The name under which to save the output vertex list
    bool            m_replaceCurrentVertexList;     ///< Whether to replace the current vertex list with the output list

    bool            m_beamMode;                     ///< Whether to run in beam mode, assuming neutrinos travel in positive z-direction
    bool            m_selectSingleVertex;           ///< Whether to make a final decision and select just one vertex candidate

    float           m_kernelEstimateSigma;          ///< The Gaussian width to use for kernel estimation

    float           m_maxOnHitDisplacement;         ///< Max hit-vertex displacement for declaring vertex to lie on a hit in each view
    float           m_maxHitVertexDisplacement1D;   ///< Max hit-vertex displacement in *any one dimension* for contribution to kernel estimation

    unsigned int    m_maxTopScoreCandidates;        ///< Max number of top-scoring vertices to examine and put forward for final selection
    unsigned int    m_maxTopScoreSelections;        ///< Max number of top-scoring vertices to select for final investigation
    unsigned int    m_maxBeamTopScoreCandidates;    ///< Max number of top-beam-scoring vertices to examine and put forward for final selection
    unsigned int    m_maxBeamTopScoreSelections;    ///< Max number of top-beam-scoring vertices to select for final investigation

    float           m_minCandidateDisplacement;     ///< Ignore other top-scoring candidates located in close proximity to original
    float           m_minCandidateScoreFraction;    ///< Ignore other top-scoring candidates with score less than a fraction of original
    float           m_minBeamCandidateScoreFraction;///< Ignore other top-beam-scoring candidates with score less than a fraction of original

    float           m_nDecayLengthsInZSpan;         ///< The number of score decay lengths to use over the course of the vertex z-span
    float           m_bestScoreMultiplier;          ///< In beam mode, best vertex must surpass a multiple of current best basic score
    float           m_bestBeamScoreMultiplier;      ///< In beam mode, best vertex must surpass a multiple of current best beam-weighted score
    float           m_mustUseBeamScoreMultiplier;   ///< Use ratio of beam scores to decide when to ensure beam score always overturns basic score
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

inline VertexSelectionAlgorithm::KernelEstimate::KernelEstimate(const float sigma) :
    m_sigma(sigma),
    m_weightSum(0.f)
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

//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexSelectionAlgorithm::KernelEstimate::GetWeightSum() const
{
    return m_weightSum;
}

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_SELECTION_ALGORITHM_H
