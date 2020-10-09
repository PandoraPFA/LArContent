/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h
 *
 *  @brief  Header file for the three dimensional hit creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_HIT_CREATION_ALGORITHM_H
#define LAR_THREE_D_HIT_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArUtility/RANSAC/AbstractModel.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <vector>

namespace lar_content
{

class HitCreationBaseTool;
class ThreeDSlidingFitResult;
class RANSACHit;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDHitCreationAlgorithm::Algorithm class
 */
class ThreeDHitCreationAlgorithm : public pandora::Algorithm
{
public:

    typedef std::vector<RANSACHit> RANSACHitVector;

    /**
     *  @brief  Trajectory samples record the results of sampling a particles in a particular view
     */
    class TrajectorySample
    {
    public:
        /**
         *  @brief  Constructor
         */
        TrajectorySample(const pandora::CartesianVector &position, const pandora::HitType hitType, const double sigma);

        /**
         *  @brief  Get the sampling position
         *
         *  @return the sampling position
         */
        const pandora::CartesianVector &GetPosition() const;

        /**
         *  @brief  Get the sampling hit type
         *
         *  @return the sampling hit type
         */
        pandora::HitType GetHitType() const;

        /**
         *  @brief  Get the sampling sigma
         *
         *  @return the sampling sigma
         */
        double GetSigma() const;

    private:
        pandora::CartesianVector    m_position;             ///< The sampling position
        pandora::HitType            m_hitType;              ///< The sampling hit type
        double                      m_sigma;                ///< The sampling sigma
    };

    typedef std::vector<TrajectorySample> TrajectorySampleVector;

    /**
     *  @brief  Proto hits are temporary constructs to be used during iterative 3D hit procedure
     */
    class ProtoHit
    {
    public:
        /**
         *  @brief  Constructor to init ProtoHit.
         *
         */
        ProtoHit();

        /**
         *  @brief  Constructor
         *
         *  @param  pParentCaloHit2D the address of the parent 2D calo hit
         */
        ProtoHit(const pandora::CaloHit *const pParentCaloHit2D);

        /**
         *  @brief  Get the address of the parent 2D calo hit
         *
         *  @return the address of the parent 2D calo hit
         */
        const pandora::CaloHit *GetParentCaloHit2D() const;

        /**
         *  @brief  Whether the proto hit has been initialised
         *
         *  @return boolean
         */
        bool IsInitialised() const;

        /**
         *  @brief  Whether the proto hit position is set
         *
         *  @return boolean
         */
        bool IsPositionSet() const;

        /**
         *  @brief  Whether the proto hit was generated using interpolation.
         *
         *  @return boolean
         */
        bool IsInterpolated() const;

        /**
         *  @brief  Get the output 3D position
         *
         *  @return the output 3D position, if set
         *
         *  @throws StatusCodeException
         */
        const pandora::CartesianVector &GetPosition3D() const;

        /**
         *  @brief  Get the chi squared value
         *
         *  @return the chi squared value, if set
         *
         *  @throws StatusCodeException
         */
        double GetChi2() const;

        /**
         *  @brief  Get the number of trajectory samples
         *
         *  @return the number of trajectory samples
         */
        unsigned int GetNTrajectorySamples() const;

        /**
         *  @brief  Get the first trajectory sample
         *
         *  @return the first trajectory sample, if at least one sample is present
         *
         *  @throws StatusCodeException
         */
        const TrajectorySample &GetFirstTrajectorySample() const;

        /**
         *  @brief  Get the last trajectory sample
         *
         *  @return the last trajectory sample, if at least two samples are present
         *
         *  @throws StatusCodeException
         */
        const TrajectorySample &GetLastTrajectorySample() const;

        /**
         *  @brief  Set position 3D
         *
         *  @param  the output 3D position
         *  @param  the output chi squared value
         */
        void SetPosition3D(const pandora::CartesianVector &position3D, const double chi2);

        /**
         *  @brief  Set the protohit as an interpolated hit or not.
         *
         *  @param  if the hit is interpolated
         */
        void SetInterpolated(const bool interpolated);

        /**
         *  @brief  Add a trajectory sample
         *
         *  @param  the trajectory sample
         */
        void AddTrajectorySample(const TrajectorySample &trajectorySample);

        /**
         * @brief  Equality operator for a ProtoHit, compares position and parent hit.
         */
        bool operator==(const ProtoHit &other) const;

    private:
        const pandora::CaloHit     *m_pParentCaloHit2D;         ///< The address of the parent 2D calo hit
        bool                        m_isInitialised;            ///< Whether the ProtoHit has been initialised
        bool                        m_isPositionSet;            ///< Whether the output 3D position has been set
        bool                        m_isInterpolated;           ///< Whether the 3D position was built with interpolation.
        pandora::CartesianVector    m_position3D;               ///< The output 3D position
        double                      m_chi2;                     ///< The output chi squared value
        TrajectorySampleVector      m_trajectorySampleVector;   ///< The trajectory sample vector
    };

    typedef std::vector<ProtoHit> ProtoHitVector;
    typedef std::map<std::string, ProtoHitVector> ProtoHitVectorMap;

    /**
     *  @brief  Default constructor
     */
    ThreeDHitCreationAlgorithm();

    /**
     *  @brief  Get the subset of a provided calo hit vector corresponding to a specified hit type
     *
     *  @param  inputCaloHitVector the input calo hit vector
     *  @param  hitType the hit type to filter upon
     *  @param  outputCaloHitVector to receive the output calo hit vector
     */
    void FilterCaloHitsByType(const pandora::CaloHitVector &inputCaloHitVector, const pandora::HitType hitType,
        pandora::CaloHitVector &outputCaloHitVector) const;

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Get the list of 2D calo hits in a pfo for which 3D hits have and have not been created
     *
     *  @param  pPfo the address of the pfo
     *  @param  protoHitVector the vector of proto hits, describing current state of 3D hit construction
     *  @param  remainingHitVector to receive the vector of 2D calo hits for which 3D hits have not been created
     */
    void SeparateTwoDHits(const pandora::ParticleFlowObject *const pPfo, const ProtoHitVector &protoHitVector,
        pandora::CaloHitVector &remainingHitVector) const;

    /**
     *  @brief  Improve initial 3D hits by fitting proto hits and iteratively creating consisted 3D hit trajectory
     *
     *  @param  protoHitVector the vector of proto hits, describing current state of 3D hit construction
     */
    void IterativeTreatment(ProtoHitVector &protoHitVector) const;

    /**
     *  @brief  Choose between the map of all protoHitVectors, to get the best
     *  and most appropriate set of hits for the current event. Uses RANSAC under the hood.
     *
     *  @param  pPfo The current PFO, so we can interpolate using it later.
     *  @param  protoHitVectorMap The map of all protoHitVectors, mapped from the algorithm that created them.
     *  @param  protoHitVector An empty protoHitVector, to be filled with the current state of the 3D hit construction.
     */
    void ConsolidatedMethod(const pandora::ParticleFlowObject *const pPfo, ProtoHitVectorMap &protoHitVectorMap,
            ProtoHitVector &protoHitVector);

    /**
     *  @brief  Project a ProtoHit into the given view.
     *
     *  @param  hit The ProtoHit to project.
     *  @param  view The view to project the hit into.
     */
    void Project3DHit(const ProtoHit &hit, const pandora::HitType view, ProtoHit &projectedHit);

    /**
     *  @brief  Take the set intersection of two vectors.
     *
     *  @param  first The first ProtoHitVector, to take into the set intersection.
     *  @param  second The second ProtoHitVector, to take into the set intersection.
     *  @param  result The result ProtoHitVector, to store the result of the intersection.
     */
    void GetSetIntersection(RANSACHitVector &first, RANSACHitVector &second, RANSACHitVector &result);

    /**
     *  @brief  Interpolate over the given hits to get a more complete image of
     *          the 3D reconstruction for the given algorithm.
     *
     *  @param  pfo the address of the pfo
     *  @param  protoHitVector The protoHitVector for the current algorithm, to be interpolated over.
     */
    void InterpolationMethod(const pandora::ParticleFlowObject *const pfo, ProtoHitVector &protoHitVector) const;

    /**
     *  @brief  Extract key results from a provided proto hit vector
     *
     *  @param  protoHitVector the proto hit vector
     *  @param  chi2 to receive the sum of the proto hit chi2 values
     *  @param  pointVector to receive a vector of proto hit 3D positions
     */
    void ExtractResults(const ProtoHitVector &protoHitVector, double &chi2, pandora::CartesianPointVector &pointVector) const;

    /**
     *  @brief  Receive a chi2 value indicating consistency of a list of proto hits with a provided 3D sliding fit trajectory
     *
     *  @param  slidingFitResult the 3D sliding fit result
     *  @param  protoHitVector the proto hit vector
     *
     *  @return the chi2 value
     */
    double GetChi2WrtFit(const ThreeDSlidingFitResult &slidingFitResult, const ProtoHitVector &protoHitVector) const;

    /**
     *  @brief  Receive a chi2 value indicating consistency of a list of proto hits with the original, input hit positions
     *
     *  @param  protoHitVector the proto hit vector
     *
     *  @return the chi2 value
     */
    double GetHitMovementChi2(const ProtoHitVector &protoHitVector) const;

    /**
     *  @brief  Refine the 3D hit positions (and chi2) for a list of proto hits, in accordance with a provided 3D sliding fit trajectory
     *
     *  @param  slidingFitResult the 3D sliding fit result
     *  @param  protoHitVector the proto hit vector, non const as proto hit properties will be updated
     */
    void RefineHitPositions(const ThreeDSlidingFitResult &slidingFitResult, ProtoHitVector &protoHitVector) const;

    /**
     *  @brief  Create new three dimensional hits from two dimensional hits
     *
     *  @param  protoHitVector the input proto hit vector
     *  @param  newThreeDHits to receive the addresses of the new three dimensional calo hits
     */
    void CreateThreeDHits(const ProtoHitVector &protoHitVector, pandora::CaloHitList &newThreeDHits) const;

    /**
     *  @brief  Create a new three dimensional hit from a two dimensional hit
     *
     *  @param  protoHit the proto hit containing all required information
     *  @param  pCaloHit3D to receive the address of the new three dimensional calo hit
     */
    void CreateThreeDHit(const ProtoHit &protoHit, const pandora::CaloHit *&pCaloHit3D) const;

    /**
     *  @brief  Check that a new three dimensional position is not unphysical
     *
     *  @param  protoHit the proto hit
     *
     *  @param  boolean
     */
    bool CheckThreeDHit(const ProtoHit &protoHit) const;

    /**
     *  @brief  Add a specified list of three dimensional hits to a cluster in a pfo, creating the new cluster if required
     *
     *  @param  pPfo the address of the pfo
     *  @param  caloHitList the list of three dimensional hits
     */
    void AddThreeDHitsToPfo(const pandora::ParticleFlowObject *const pPfo, const pandora::CaloHitList &caloHitList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<HitCreationBaseTool*> HitCreationToolVector;
    HitCreationToolVector   m_algorithmToolVector;      ///< The algorithm tool vector

    std::string             m_inputPfoListName;         ///< The name of the input pfo list
    std::string             m_outputCaloHitListName;    ///< The name of the output calo hit list
    std::string             m_outputClusterListName;    ///< The name of the output cluster list

    bool                    m_iterateTrackHits;         ///< Whether to enable iterative improvement of 3D hits for track trajectories
    bool                    m_iterateShowerHits;        ///< Whether to enable iterative improvement of 3D hits for showers
    bool                    m_useRANSACMethod;          ///< Whether to use the RANSAC based consolidated method.
    unsigned int            m_slidingFitHalfWindow;     ///< The sliding linear fit half window
    unsigned int            m_nHitRefinementIterations; ///< The maximum number of hit refinement iterations
    double                  m_sigma3DFitMultiplier;     ///< Multiplicative factor: sigmaUVW (same as sigmaHit and sigma2DFit) to sigma3DFit
    double                  m_iterationMaxChi2Ratio;    ///< Max ratio between current and previous chi2 values to cease iterations
    double                  m_interpolationCutOff;      ///< Max distance for a point to be interpolated from.
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline ThreeDHitCreationAlgorithm::TrajectorySample::TrajectorySample(const pandora::CartesianVector &position, const pandora::HitType hitType, const double sigma) :
    m_position(position),
    m_hitType(hitType),
    m_sigma(sigma)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ThreeDHitCreationAlgorithm::TrajectorySample::GetPosition() const
{
    return m_position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::HitType ThreeDHitCreationAlgorithm::TrajectorySample::GetHitType() const
{
    return m_hitType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ThreeDHitCreationAlgorithm::TrajectorySample::GetSigma() const
{
    return m_sigma;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ThreeDHitCreationAlgorithm::ProtoHit::ProtoHit() :
    m_pParentCaloHit2D(nullptr),
    m_isInitialised(false),
    m_isPositionSet(false),
    m_isInterpolated(false),
    m_position3D(0.f, 0.f, 0.f),
    m_chi2(std::numeric_limits<double>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ThreeDHitCreationAlgorithm::ProtoHit::ProtoHit(const pandora::CaloHit *const pParentCaloHit2D) :
    m_pParentCaloHit2D(pParentCaloHit2D),
    m_isInitialised(true),
    m_isPositionSet(false),
    m_isInterpolated(false),
    m_position3D(0.f, 0.f, 0.f),
    m_chi2(std::numeric_limits<double>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHit *ThreeDHitCreationAlgorithm::ProtoHit::GetParentCaloHit2D() const
{
    return m_pParentCaloHit2D; // TODO: I've now basically changed this by allowing nullptr. Should check its usage.
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ThreeDHitCreationAlgorithm::ProtoHit::IsPositionSet() const
{
    return m_isPositionSet;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ThreeDHitCreationAlgorithm::ProtoHit::IsInitialised() const
{
    return m_isInitialised;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ThreeDHitCreationAlgorithm::ProtoHit::IsInterpolated() const
{
    return IsPositionSet() && m_isInterpolated;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int ThreeDHitCreationAlgorithm::ProtoHit::GetNTrajectorySamples() const
{
    return m_trajectorySampleVector.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ThreeDHitCreationAlgorithm::ProtoHit::SetPosition3D(const pandora::CartesianVector &position3D, const double chi2)
{
    m_position3D = position3D;
    m_chi2 = chi2;
    m_isPositionSet = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ThreeDHitCreationAlgorithm::ProtoHit::SetInterpolated(const bool interpolated)
{
    m_isInterpolated = interpolated;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ThreeDHitCreationAlgorithm::ProtoHit::AddTrajectorySample(const TrajectorySample &trajectorySample)
{
    m_trajectorySampleVector.push_back(trajectorySample);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ThreeDHitCreationAlgorithm::ProtoHit::operator==(const ProtoHit &other) const
{
    return this->m_pParentCaloHit2D == other.GetParentCaloHit2D() && this->m_position3D == other.GetPosition3D();
}

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_HIT_CREATION_ALGORITHM_H
