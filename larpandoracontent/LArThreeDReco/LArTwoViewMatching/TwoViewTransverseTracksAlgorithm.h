/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.h
 *
 *  @brief  Header file for the two view transverse tracks algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_VIEW_TRANSVERSE_TRACKS_ALGORITHM_H
#define LAR_TWO_VIEW_TRANSVERSE_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArMonitoring/TwoViewTransverseTracksValidationTool.h"

#include "larpandoracontent/LArObjects/LArDiscreteProbabilityVector.h"
#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewTrackMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingControl.h"

#include <random>

namespace lar_content
{

class TransverseMatrixTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TwoViewTransverseTracksAlgorithm class
 */
class TwoViewTransverseTracksAlgorithm : public NViewTrackMatchingAlgorithm<TwoViewMatchingControl<TwoViewTransverseOverlapResult> >
{
public:
    typedef NViewTrackMatchingAlgorithm<TwoViewMatchingControl<TwoViewTransverseOverlapResult> > BaseAlgorithm;
    typedef std::set<unsigned int> UIntSet;

    /**
     *  @brief  Default constructor
     */
    TwoViewTransverseTracksAlgorithm();

private:
    void CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        const pandora::Cluster *const);

    /**
     *  @brief  Calculates the two view overlap result
     *
     *  @param  pCluster1 the view 0 cluster
     *  @param  pCluster2 the view 1 cluster
     *  @param  overlapResult the two view overlap result
     */
    pandora::StatusCode CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        TwoViewTransverseOverlapResult &overlapResult);

    /**
     *  @brief  Calculates the number of the sliding windows that contains charge bins that locally match
     *
     *  @param  discreteProbabilityVector1 the view 0 discrete probability vector containing the charge information
     *  @param  pCluster2 the view 1 discrete probability vector containing the charge information
     *  @param  randomNumberGenerator a seeded random number generator
     *
     *  @result the number of locally matching sampling points
     */
    unsigned int CalculateNumberOfLocallyMatchingSamplingPoints(const DiscreteProbabilityVector &discreteProbabilityVector1,
        const DiscreteProbabilityVector &discreteProbabilityVector2, std::mt19937 &randomNumberGenerator);

    /**
     *  @brief  Get the dot product between the cluster's primary axis and the drift axis
     *
     *  @param  pCluster the cluster
     *
     *  @result dot product
     */
    float GetPrimaryAxisDotDriftAxis(const pandora::Cluster *const pCluster);

    void ExamineOverlapContainer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<TransverseMatrixTool*> MatrixToolVector;
    MatrixToolVector            m_algorithmToolVector;                   ///< The algorithm tool vector

    unsigned int                m_nMaxMatrixToolRepeats;                 ///< The maximum number of repeat loops over matrix tools
    unsigned int                m_downsampleFactor;                      ///< The downsampling (hit merging) applied to hits in the overlap region
    unsigned int                m_minSamples;                            ///< The minimum number of samples needed for comparing charges
    unsigned int                m_nPermutations;                         ///< The number of permutations for calculating p-values
    float                       m_localMatchingScoreThreshold;           ///< The minimum score to classify a local region as matching
    float                       m_maxDotProduct;                         ///M The maximum allowed cluster primary qxis Dot drift axis to fill the overlap result
    float                       m_minOverallMatchingScore;               ///< The minimum required global matching score to fill the overlap result
    float                       m_minOverallLocallyMatchedFraction;      ///< The minimum required lcoally matched fraction to fill the overlap result
    std::mt19937                m_randomNumberGenerator;                 ///< The random number generator

    TwoViewTransverseTracksValidationTool* p_ValidationTool;             ///< The validation tool
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TransverseMatrixTool class
 */
class TransverseMatrixTool : public pandora::AlgorithmTool
{
public:
    typedef TwoViewTransverseTracksAlgorithm::MatchingType::MatrixType MatrixType;
    typedef std::vector<MatrixType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapMatrix the overlap matrix
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_TWO_VIEW_TRANSVERSE_TRACKS_ALGORITHM_H
