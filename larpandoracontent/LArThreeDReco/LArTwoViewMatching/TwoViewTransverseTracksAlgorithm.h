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

#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"
#include "larpandoracontent/LArObjects/LArDiscreteCumulativeDistribution.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/NViewTrackMatchingAlgorithm.h"
#include "larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewMatchingControl.h"

#include <unsupported/Eigen/Splines>

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

    /**
     *  @brief  Default constructor
     */
    TwoViewTransverseTracksAlgorithm();
    virtual ~TwoViewTransverseTracksAlgorithm();

    typedef Eigen::Spline<float, 1, 1> Spline1f;
    typedef std::vector<std::pair<float,float> > ChargeProfile;
    typedef std::vector<std::pair<float,float> > ScoreProfile;
private:
    void CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const pandora::Cluster *const);
    void ExamineOverlapContainer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<TransverseMatrixTool*> MatrixToolVector;
    MatrixToolVector            m_algorithmToolVector;      ///< The algorithm tool vector

    unsigned int                m_nMaxMatrixToolRepeats;    ///< The maximum number of repeat loops over matrix tools

    Spline1f CreateSplineFromCumulativeDistribution(const DiscreteCumulativeDistribution &cumulativeDistribution);
    ChargeProfile CreateProfileFromCumulativeDistribution(const DiscreteCumulativeDistribution &cumulativeDistribution);
    float CalculateCorrelationCoefficient(const ChargeProfile &profile1, const ChargeProfile &profile2);
    float CalculateTTestValue(const float x, const float coefficient, const float dof);
    float CalculateTTestPValue(const ChargeProfile &profile1, const ChargeProfile &profile2);
    ScoreProfile SlidingWindowMatchingScore(const size_t &sizeWindowInBins, const ChargeProfile &profile1, const ChargeProfile &profile2, float &fracGoodScore);
    ChargeProfile   GetWindow(size_t &i, const size_t &sizeWindowInBins, const ChargeProfile &profile);
    ScoreProfile SlidingWindowMatchingScore(const size_t &sizeWindowInBins, const DiscreteCumulativeDistribution &disCumulDist1, const DiscreteCumulativeDistribution &disCumulDist2, float &fracGoodScore);
    ChargeProfile   GetWindow(size_t &i, const size_t &sizeWindowInBins, const DiscreteCumulativeDistribution &disCumulDist);
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
