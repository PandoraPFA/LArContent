/**
 *  @file   larpandoracontent/LArVertex/HitAngleVertexSelectionAlgorithm.h
 * 
 *  @brief  Header file for the hit angle vertex selection algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_HIT_ANGLE_VERTEX_SELECTION_ALGORITHM_H
#define LAR_HIT_ANGLE_VERTEX_SELECTION_ALGORITHM_H 1

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  HitAngleVertexSelectionAlgorithm class
 */
class HitAngleVertexSelectionAlgorithm : public VertexSelectionBaseAlgorithm
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
    HitAngleVertexSelectionAlgorithm();

private:
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

    void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
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
     *  @brief  Use hits in clusters (in the provided kd tree) to fill a provided kernel estimate with hit-vertex relationship information
     * 
     *  @param  pVertex the address of the vertex
     *  @param  hitType the relevant hit type
     *  @param  kdTree the relevant kd tree
     *  @param  kernelEstimate to receive the populated kernel estimate
     */
    void FillKernelEstimate(const pandora::Vertex *const pVertex, const pandora::HitType hitType, HitKDTree2D &kdTree, KernelEstimate &kernelEstimate) const;

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

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool            m_fastScoreCheck;               ///< Whether to use the fast histogram based score to selectively avoid calling full or midway scores
    bool            m_fastScoreOnly;                ///< Whether to use the fast histogram based score only
    bool            m_fullScore;                    ///< Whether to use the full kernel density estimation score, as opposed to the midway score
 
    float           m_kernelEstimateSigma;          ///< The Gaussian width to use for kernel estimation
    float           m_kappa;                        ///< Hit-deweighting offset, of form: weight = 1 / sqrt(distance + kappa), units cm
    float           m_maxHitVertexDisplacement1D;   ///< Max hit-vertex displacement in *any one dimension* for contribution to kernel estimation

    float           m_minFastScoreFraction;         ///< Fast score must be at least this fraction of best fast score to calculate full score
    unsigned int    m_fastHistogramNPhiBins;        ///< Number of bins to use for fast score histograms
    float           m_fastHistogramPhiMin;          ///< Min value for fast score histograms
    float           m_fastHistogramPhiMax;          ///< Max value for fast score histograms

    bool            m_enableFolding;                ///< Whether to enable folding of -pi -> +pi phi distribution into 0 -> +pi region only
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *HitAngleVertexSelectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new HitAngleVertexSelectionAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline HitAngleVertexSelectionAlgorithm::KernelEstimate::KernelEstimate(const float sigma) :
    m_sigma(sigma)
{
    if (m_sigma < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const HitAngleVertexSelectionAlgorithm::KernelEstimate::ContributionList &HitAngleVertexSelectionAlgorithm::KernelEstimate::GetContributionList() const
{
    return m_contributionList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float HitAngleVertexSelectionAlgorithm::KernelEstimate::GetSigma() const
{
    return m_sigma;
}

} // namespace lar_content

#endif // #ifndef LAR_HIT_ANGLE_VERTEX_SELECTION_ALGORITHM_H
