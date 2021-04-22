/**
 *  @file   larpandoracontent/LArControlFlow/TrackDirectionTool.h
 *
 *  @brief  Header file for the track direction finding tool class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_DIRECTION_TOOL_H
#define LAR_TRACK_DIRECTION_TOOL_H 1

#include "Pandora/AlgorithmTool.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArPfoObjects.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include <numeric>

namespace lar_content{

class TrackDirectionTool : public TrackDirectionBaseTool
{

public:
    TrackDirectionTool();
    pandora::StatusCode Initialize();

    float FindDirections( const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo);

    class HitCharge
    {
    public:
        HitCharge(const pandora::CaloHit *caloHit, float longitudinalPosition, float hitWidth, float hitCharge, float uncertainty);
        HitCharge();

	bool GetInTails() const;
        const pandora::CaloHit *GetCaloHit() const;
        float GetLongitudinalPosition() const;
        float GetHitWidth() const;
        float GetCharge() const;
        float GetChargeOverWidth() const;
        float GetUncertainty() const;

        void SetDistanceToNN(float distance);
        float GetDistanceToNN() const;

        void SetForwardsFitCharge(float qFitForwards);
        void SetForwardsSigma(float sigmaForwards);
        void SetForwardsDelta(float forwardsDelta);
        void SetForwardsChiSquared(float forwardsHitChiSquared);

        void SetBackwardsFitCharge(float qFitBackwards);
        void SetBackwardsSigma(float sigmaBackwards);
        void SetBackwardsDelta(float backwardsDelta);
        void SetBackwardsChiSquared(float backwardsHitChiSquared);

        float GetForwardsFitCharge() const;
        float GetForwardsSigma() const;
        float GetForwardsDelta() const;
        float GetForwardsChiSquared() const;

        float GetBackwardsFitCharge() const;
        float GetBackwardsSigma() const;
        float GetBackwardsDelta() const;
        float GetBackwardsChiSquared() const;

    private:

	bool m_intails;                     ///< Is this in the tails of the distribution?

        const pandora::CaloHit *m_caloHit;
        float m_longitudinalPosition;       ///< Get longitudinal position of hit
        float m_hitWidth;                   ///< Get hit width
        float m_charge;                     ///< Get hit charge
	float m_qoverx;                     ///< Get hit charge over width
        float m_uncertainty;                ///< Get uncertainty on hit charge
        float m_distanceToNearestNeighbour; ///< Get distence to nearest hit neighbour

        float m_forwardsFitCharge;          ///< Get the charge fit in the forwards firection
        float m_forwardsSigma;              ///< Get the sigma for the hit charge in the forwards direction
        float m_forwardsDelta;              ///< Get the delta chi squared value for the hit charge in the forwards direction
        float m_forwardsChiSquared;         ///< Get the chi squared value in the forards direction

        float m_backwardsFitCharge;          ///< Get the charge fit in the backwards firection
        float m_backwardsSigma;              ///< Get the sigma for the hit charge in the backwards direction
        float m_backwardsDelta;              ///< Get the delta chi squared value for the hit charge in the backwards direction
        float m_backwardsChiSquared;         ///< Get the chi squared value in the backwards direction
    };

    typedef std::vector<HitCharge> HitChargeVector;

    class LookupTable
    {
    public:
        LookupTable();
        LookupTable(float initialEnergy, float binWidth);
        std::map<int, float> GetMap();
        void SetMap(const std::map<int, float> &map);
        std::map<float, int> GetReverseMap();
        void SetReverseMap(const std::map<float, int> &map);

        float GetInitialEnergy() const;
        void SetInitialEnergy(float initialEnergy);

        float GetBinWidth() const;
        void SetBinWidth(float binWidth);

        void SetMaxRange(float maxRange);
        float GetMaxRange() const;

    private:
        std::map<int, float> m_map;
        std::map<float, int> m_reversemap;
        float m_binWidth;         ///< Energy bin width for lookup table
        float m_initialEnergy;    ///< Initial energy for the lookip table
        float m_maxRange;         ///< Lookup table max range
    };

    class DirectionFitObject
    {
    public:
        DirectionFitObject();
        DirectionFitObject(HitChargeVector hitChargeVector, int numberHits, float meanChargeOverWidth, float forwardsChiSquared,
            float backwardsChiSquared);
        DirectionFitObject(HitChargeVector hitChargeVector, HitChargeVector forwardsRecoHits, HitChargeVector backwardsRecoHits,
            int numberHits, float meanChargeOverWidth, float forwardsChiSquared, float backwardsChiSquared);

        TrackDirectionTool::HitChargeVector GetHitChargeVector();
        TrackDirectionTool::HitChargeVector GetForwardsFitCharges();
        TrackDirectionTool::HitChargeVector GetBackwardsFitCharges();

        void SetForwardsFitCharges(const TrackDirectionTool::HitChargeVector &hitChargeVector);
        void SetBackwardsFitCharges(const TrackDirectionTool::HitChargeVector &hitChargeVector);
        float GetForwardsChiSquared();
        float GetBackwardsChiSquared();
        void SetForwardsChiSquared(float forwardsChiSquared);
        void SetBackwardsChiSquared(float backwardsChiSquared);
        float GetForwardsChiSquaredPerHit();
        float GetBackwardsChiSquaredPerHit();

        void SetNHits(float numberHits);
        int GetNHits();
        int GetDirectionEstimate();

        float GetMinChiSquared();
        float GetMinChiSquaredPerHit();
        float GetDeltaChiSquaredPerHit();
        float GetMeanChargeOverWidth();

        void SetBeginPoint(const pandora::CartesianVector &beginPoint);
        void SetEndPoint(const pandora::CartesianVector &endPoint);
        const pandora::CartesianVector GetBeginPoint();
        const pandora::CartesianVector GetEndPoint();

        void SetProbability(float probability);
        float GetProbability();
        void SetMCDirection(int direction);
        int GetMCDirection();

        void Print();

    private:
        HitChargeVector m_hitChargeVector;     ///< Vector of the hit charge points
        HitChargeVector m_forwardsRecoHits;    ///< Vector of the hit charge points gathered for the forwards hypothesis
        HitChargeVector m_backwardsRecoHits;   ///< Vector of the hit charge points gathered for the backwards hypothesis

        int m_nHits;                           ///< The number of hits

        float m_meanqoverx;                    ///< Mean charge over width value
        float m_forwardsChiSquared;            ///< Chi squared value for the forwards hypothosis
        float m_backwardsChiSquared;           ///< Chi squared value for the backwards hypothosis
        float m_probability;                   ///< The probability of the considered object being a CR

        float m_beginx;                        ///< Beginpoint is defined as the track endpoint with the lowest Z coordinate
        float m_beginy;
        float m_beginz;

        float m_endx;                          ///< Endpoint is defined as the track endpoint with the highest Z coordinate
        float m_endy;
        float m_endz;

        int m_mcDirection;                     ///< The MC direction of a given object
    };

    void WriteLookupTableToTree(LookupTable &lookupTable) const;

    void GetClusterDirection(const pandora::Cluster *const pTargetClusterW, TrackDirectionTool::DirectionFitObject &finalDirectionFitObject);
    TrackDirectionTool::DirectionFitObject GetPfoDirection(const pandora::ParticleFlowObject *const pPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    //-----------------------------------------------------------------------------------------------

    const pandora::Cluster *GetTargetClusterFromPFO(const pandora::ParticleFlowObject *PFO, const LArTrackStateVector &trackStateVector);

    void SetEndPoints(DirectionFitObject &fitResult, const pandora::Cluster *const pCluster);

    void SetEndPoints(DirectionFitObject &fitResult, const LArTrackStateVector &trackStateVector);

    void FillHitChargeVector(const pandora::Cluster *const pCluster, HitChargeVector &hitChargeVector);

    void TrackInnerFilter(HitChargeVector &hitChargeVector, HitChargeVector &filteredHitChargeVector);

    void SetNearestNeighbourValues(HitChargeVector &innerHitChargeVector, const int nNeighboursToConsider);

    void SimpleTrackEndFilter(HitChargeVector &hitChargeVector);

    void SplitHitCollectionBySize(const HitChargeVector &hitChargeVector, float splitPosition, HitChargeVector &smallHitChargeVector,
        HitChargeVector &largeHitChargeVector);

    void SplitHitCollectionByLeftRight(
        const HitChargeVector &hitChargeVector, float splitPosition, HitChargeVector &leftHitChargeVector, HitChargeVector &rightHitChargeVector);

    float GetTrackLength(const HitChargeVector &hitChargeVector) const;

    void FitHitChargeVector(HitChargeVector &hitChargeVector, TrackDirectionTool::DirectionFitObject &fitResult, int numberHitsToConsider = 1000000);

    void ComputeProbability(DirectionFitObject &fitResult);

    void PerformFits(HitChargeVector &hitChargeVector, HitChargeVector &forwardsFitPoints, HitChargeVector &backwardsFitPoints,
        int numberHitsToConsider, float &forwardsChiSquared, float &backwardsChiSquared);

    void GetFitParameters(bool isForwardsFit, LookupTable lookupTable, HitChargeVector binnedHitChargeVector, const float trackLength, const float totalCharge, const float totalHitWidth, const float maxScale, const float particleMass, float &holdp0Value, float &holdp1Value, float &holdp2Value);

    void GetCalorimetricDirection(const pandora::Cluster *pTargetClusterW, DirectionFitObject &finalDirectionFitObject);

    void AddToSlidingFitCache(const pandora::Cluster *const pCluster);

    const TwoDSlidingFitResult &GetCachedSlidingFit(const pandora::Cluster *const pCluster) const;

    static bool SortHitChargeVectorByRL(HitCharge &hitCharge1, HitCharge &hitCharge2);
    static bool SortHitChargeVectorByChargeOverWidth(HitCharge &hitCharge1, HitCharge &hitCharge2);
    static bool SortByDistanceToNN(HitCharge &hitCharge1, HitCharge &hitCharge2);
    void BinHitChargeVector(HitChargeVector &hitChargeVector, HitChargeVector &binnedHitChargeVector);

    float DensityCorrection(const float T, const float M);
    float BetheBloch(float &T, const float M);
    void FillLookupTable(LookupTable &lookupTable, const float M);
    float GetEnergyfromLength(LookupTable &lookupTable, const float &trackLength);
    float GetLengthfromEnergy(LookupTable &lookupTable, const float &currentEnergy);

    //-----------------------------------------------------------------------------------------------

    unsigned int m_slidingFitWindow;               ///< The layer window for the sliding linear fits
    TwoDSlidingFitResultMap m_slidingFitResultMap; ///< The sliding fit result map

    unsigned int m_minClusterCaloHits;             ///< The min number of hits in base cluster selection method
    float m_minClusterLength;                      ///< The min length (squared) in base cluster selection method

    float m_tableInitialEnergy;                    ///< The initial energy used in creating the lookup table
    float m_tableStepSize;                         ///< The step size used in the lookup table

    std::string m_lookupTableFileName;             ///< If using an external file for the lookup table, this is the file name
    std::string m_treeName;                        ///< Name of the tree in lookup table external file

    float m_endpointProtectionRange;               ///< How much of the hit vector to set as the end points
    int m_nNeighboursToConsider;                   ///< How many nearest neighbour values to set in the innerHitChargeVector
    float m_lowerBound;                            ///< Lower bound for SimpleTrackEndFilter
    float m_upperBound;                            ///< Upper bound for SimpleTrackEndFilter
    float m_endFilterMultiplierTrack;              ///< Multiplier for the Track Length in SimpleTrackEndFilter
    float m_endFilterMultiplierCharge;             ///< Multiplier for the charge in SimpleTrackEndFilter
    float m_sigmaFitMultiplier;                    ///< Multiplier for the sigma values in PerformFits

    float m_ADCToElectron;                         ///< Convert ADC to electron energy for caloHitEnergy
    float m_ionEPerElectron;                       ///< Ionisation Energy per electron in Mev for caloHitEnergy
    float m_energyScale;                           ///< Energy sclae factor for caloHitEnergy

    float m_uncertaintyCalibration1;               ///< Parameter 1 for calibratedUncertainty in FillHitChargeVector
    float m_uncertaintyCalibration2;               ///< Parameter 2 for calibratedUncertainty in FillHitChargeVector
    float m_innerHitChargeMultiplier;              ///< Multiplier for the inner HitChargeVector
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitCharge::HitCharge(
    const pandora::CaloHit *caloHit, float longitudinalPosition, float hitWidth, float hitCharge, float uncertainty) :
    m_intails(false),
    m_caloHit(caloHit),
    m_longitudinalPosition(longitudinalPosition),
    m_hitWidth(hitWidth),
    m_charge(hitCharge),
    m_qoverx(hitCharge / hitWidth),
    m_uncertainty(uncertainty),
    m_forwardsFitCharge(0.f),
    m_forwardsSigma(0.f),
    m_forwardsDelta(0.f),
    m_forwardsChiSquared(0.f),
    m_backwardsFitCharge(0.f),
    m_backwardsSigma(0.f),
    m_backwardsDelta(0.f),
    m_backwardsChiSquared(0.f)
{
  if (hitWidth != 0) 
  {
      m_intails = (hitCharge / hitWidth) <= 1.4 ? true : false;
  }
  else 
  {
    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackDirectionTool::HitCharge::GetInTails() const
{
    return m_intails;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHit *TrackDirectionTool::HitCharge::GetCaloHit() const
{
    return m_caloHit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetLongitudinalPosition() const
{
    return m_longitudinalPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetHitWidth() const
{
    return m_hitWidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetCharge() const
{
    return m_charge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetChargeOverWidth() const
{
   return m_qoverx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetUncertainty() const
{
    return m_uncertainty;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::HitCharge::SetDistanceToNN(float distance)
{
   m_distanceToNearestNeighbour = distance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetDistanceToNN() const
{
   return m_distanceToNearestNeighbour;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::HitCharge::SetForwardsFitCharge(float qFitForwards)
{
    m_forwardsFitCharge = qFitForwards;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::HitCharge::SetForwardsSigma(float sigmaForwards)
{
    m_forwardsSigma = sigmaForwards;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::HitCharge::SetForwardsDelta(float forwardsDelta)
{
    m_forwardsDelta = forwardsDelta;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::HitCharge::SetForwardsChiSquared(float forwardsHitChiSquared)
{
    m_forwardsChiSquared = forwardsHitChiSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::HitCharge::SetBackwardsFitCharge(float qFitBackwards) 
{
    m_backwardsFitCharge = qFitBackwards;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::HitCharge::SetBackwardsSigma(float sigmaBackwards)
{
    m_backwardsSigma = sigmaBackwards;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::HitCharge::SetBackwardsDelta(float backwardsDelta)
{
    m_backwardsDelta = backwardsDelta;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::HitCharge::SetBackwardsChiSquared(float backwardsHitChiSquared)
{
    m_backwardsChiSquared = backwardsHitChiSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetForwardsFitCharge() const
{
    return m_forwardsFitCharge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetForwardsSigma() const
{
    return m_forwardsSigma;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetForwardsDelta() const
{
    return m_forwardsDelta;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetForwardsChiSquared() const
{
    return m_forwardsChiSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetBackwardsFitCharge() const
{
    return m_backwardsFitCharge;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetBackwardsSigma() const
{
    return m_backwardsSigma;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetBackwardsDelta() const
{
    return m_backwardsDelta;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetBackwardsChiSquared() const
{
    return m_backwardsChiSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::LookupTable::LookupTable() :
  m_binWidth(0.f),
  m_initialEnergy(0.f), 
  m_maxRange(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

 inline TrackDirectionTool::LookupTable::LookupTable(float initialEnergy, float binWidth) :
  m_binWidth(binWidth),
  m_initialEnergy(initialEnergy), 
  m_maxRange(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::map<int, float> TrackDirectionTool::LookupTable::GetMap()
{
    return m_map;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetMap(const std::map<int, float> &map)
{
    m_map = map;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::map<float, int> TrackDirectionTool::LookupTable::GetReverseMap()
{
    return m_reversemap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetReverseMap(const std::map<float, int> &reverseMap)
{
    m_reversemap = reverseMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::LookupTable::GetInitialEnergy() const
{
    return m_initialEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetInitialEnergy(float initialEnergy)
{
    m_initialEnergy = initialEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::LookupTable::GetBinWidth() const
{
    return m_binWidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetBinWidth(float binWidth)
{
    m_binWidth = binWidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetMaxRange(float maxRange)
{
    m_maxRange = maxRange;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::LookupTable::GetMaxRange() const
{
    return m_maxRange;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::DirectionFitObject::DirectionFitObject() :
  m_nHits(0),
  m_meanqoverx(0.f),
  m_forwardsChiSquared(0.f),
  m_backwardsChiSquared(0.f),
  m_probability(0.5),
  m_mcDirection(-1)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::DirectionFitObject::DirectionFitObject(
    HitChargeVector hitChargeVector, int numberHits, float meanChargeOverWidth, float forwardsChiSquared, float backwardsChiSquared)
{
    m_hitChargeVector = hitChargeVector;
    m_nHits = numberHits;
    m_meanqoverx = meanChargeOverWidth;
    m_forwardsChiSquared = forwardsChiSquared;
    m_backwardsChiSquared = backwardsChiSquared;
    m_probability = 0.5;
    m_mcDirection = -1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::DirectionFitObject::DirectionFitObject(HitChargeVector hitChargeVector, HitChargeVector forwardsRecoHits,
    HitChargeVector backwardsRecoHits, int numberHits, float meanChargeOverWidth, float forwardsChiSquared, float backwardsChiSquared) :
    m_hitChargeVector(hitChargeVector),
    m_forwardsRecoHits(forwardsRecoHits),
    m_backwardsRecoHits(backwardsRecoHits),
    m_nHits(numberHits),
    m_meanqoverx(meanChargeOverWidth),
    m_forwardsChiSquared(forwardsChiSquared),
    m_backwardsChiSquared(backwardsChiSquared)
{
    m_probability = 0.5;
    m_mcDirection = -1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitChargeVector TrackDirectionTool::DirectionFitObject::GetHitChargeVector()
{
    return m_hitChargeVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitChargeVector TrackDirectionTool::DirectionFitObject::GetForwardsFitCharges()
{
    return m_forwardsRecoHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitChargeVector TrackDirectionTool::DirectionFitObject::GetBackwardsFitCharges()
{
    return m_backwardsRecoHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetForwardsFitCharges(const TrackDirectionTool::HitChargeVector &hitChargeVector)
{
    m_forwardsRecoHits = hitChargeVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetBackwardsFitCharges(const TrackDirectionTool::HitChargeVector &hitChargeVector)
{
    m_backwardsRecoHits = hitChargeVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetForwardsChiSquared()
{
    return m_forwardsChiSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetBackwardsChiSquared()
{
    return m_backwardsChiSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetForwardsChiSquared(float forwardsChiSquared)
{
    m_forwardsChiSquared = forwardsChiSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetBackwardsChiSquared(float backwardsChiSquared)
{
    m_backwardsChiSquared = backwardsChiSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetForwardsChiSquaredPerHit()
{
    return m_forwardsChiSquared / m_nHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetBackwardsChiSquaredPerHit()
{
    return m_backwardsChiSquared / m_nHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetNHits(float numberHits)
{
    m_nHits = numberHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int TrackDirectionTool::DirectionFitObject::GetNHits()
{
    return m_hitChargeVector.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int TrackDirectionTool::DirectionFitObject::GetDirectionEstimate()
{
    return (m_forwardsChiSquared <= m_backwardsChiSquared ? 1 : 0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetMinChiSquared()
{
    return std::min(m_forwardsChiSquared, m_backwardsChiSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetMinChiSquaredPerHit()
{
    return (m_nHits != 0 ? std::min(m_forwardsChiSquared, m_backwardsChiSquared) / m_nHits : std::numeric_limits<float>::max());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetDeltaChiSquaredPerHit()
{
    return (m_nHits != 0 ? (m_forwardsChiSquared - m_backwardsChiSquared) / m_nHits : std::numeric_limits<float>::max());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetMeanChargeOverWidth()
{
    return m_meanqoverx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetBeginPoint(const pandora::CartesianVector &beginPoint)
{
    m_beginx = beginPoint.GetX();
    m_beginy = beginPoint.GetY();
    m_beginz = beginPoint.GetZ();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetEndPoint(const pandora::CartesianVector &endPoint)
{
    m_endx = endPoint.GetX();
    m_endy = endPoint.GetY();
    m_endz = endPoint.GetZ();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector TrackDirectionTool::DirectionFitObject::GetBeginPoint()
{
    const pandora::CartesianVector beginPoint(m_beginx, m_beginy, m_beginz);
    return beginPoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector TrackDirectionTool::DirectionFitObject::GetEndPoint()
{
    const pandora::CartesianVector endPoint(m_endx, m_endy, m_endz);
    return endPoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetProbability(float probability)
{
    m_probability = probability;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetProbability()
{
    return m_probability;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetMCDirection(int direction)
{
    m_mcDirection = direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int TrackDirectionTool::DirectionFitObject::GetMCDirection()
{
    return m_mcDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::Print()
{
    std::cout << "Probability: " << m_probability << std::endl;
    std::cout << "Vertex position: (" << m_beginx << ", " << m_beginy << ", " << m_beginz << ")" << std::endl;
    std::cout << "Endpoint position: (" << m_endx << ", " << m_endy << ", " << m_endz << ")" << std::endl;
    std::cout << "Best chi squared per hit: " << std::min(m_forwardsChiSquared, m_backwardsChiSquared) / m_nHits << std::endl;
    std::cout << "Number of hits: " << m_nHits << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_TRACK_DIRECTION_TOOL_H
