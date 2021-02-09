/**
 *  @file   larpandoracontent/LArControlFlow/TrackDirectionTool.h
 *
 *  @brief  Header file for the track direction finding tool AlgorithmTool class.
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_DIRECTION_TOOL_H
#define LAR_TRACK_DIRECTION_TOOL_H 1

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include <numeric>

namespace lar_content
{

class TrackDirectionTool : public TrackDirectionBaseTool
{

friend void GetSplitChiSquared(int &npar, double *gin, double &f, double *par, int flag);

public:

    TrackDirectionTool();
    pandora::StatusCode Initialize();

    void FindDirections(const pandora::ParticleFlowObject *const pPfo, float &downProbability, const MasterAlgorithm *const pAlgorithm);

    class HitCharge
    {
    public:

        HitCharge(const pandora::CaloHit* caloHit, float &longitudinalPosition, float &hitWidth, float &hitCharge, float &uncertainty);
        HitCharge();

        const pandora::CaloHit* GetCaloHit() const;
        float GetLongitudinalPosition() const;
        float GetHitWidth() const;
        float GetCharge() const;
        float GetChargeOverWidth() const;
        float GetUncertainty() const;

        void SetDistanceToNN(float &distance);
        float GetDistanceToNN() const;

        void SetForwardsFitCharge(float &Q_fit_f); 
        void SetForwardsSigma(float &f_sigma);
        void SetForwardsDelta(float &forwardsDelta);
        void SetForwardsChiSquared(float &forwardsHitChisquared);

        void SetBackwardsFitCharge(float &Q_fit_b); 
        void SetBackwardsSigma(float &b_sigma);
        void SetBackwardsDelta(float &backwardsDelta);
        void SetBackwardsChiSquared(float &backwardsHitChisquared);

        float GetForwardsFitCharge(); 
        float GetForwardsSigma();
        float GetForwardsDelta();
        float GetForwardsChiSquared();

        float GetBackwardsFitCharge(); 
        float GetBackwardsSigma();
        float GetBackwardsDelta();
        float GetBackwardsChiSquared();

        bool                                       m_intails;

    private:
        const pandora::CaloHit*                    m_calohit;
        float                                      m_longitudinalposition;    ///<
        float                                      m_hitwidth;                ///<
        float                                      m_charge;                  ///<
        float                                      m_qoverx;                  ///<
        float                                      m_uncertainty;          ///<
        float                                      m_distancetonearestneighbour;          ///<

        float                                      m_forwardsfitcharge;
        float                                      m_forwardssigma;
        float                                      m_forwardsdelta;
        float                                      m_forwardschisquared;

        float                                      m_backwardsfitcharge;
        float                                      m_backwardssigma;
        float                                      m_backwardsdelta;
        float                                      m_backwardschisquared;
    };

    typedef std::vector<HitCharge> HitChargeVector;

    class SplitObject
    {
    public:

        SplitObject();
        SplitObject(int beforeNumerHits, int afterNumberHits, float beforeMinChiSquaredPerHit, float afterMinChiSquaredPerhit, float chiSquaredPerHitChange, float splitPosition);
        SplitObject(int beforeNumberHits, int afterNumberHits, float beforeMinChiSquaredPerHit, float afterMinChiSquaredPerHit, float chiSquaredPerHitChange, float splitPosition, bool splitApplied, float beforeDeltaChiSquaredPerHit); 

        void   SetBeforeNHits(int beforeNumberHits);
        void   SetAfterNHits(int afterNumberHits);
        void   SetBeforeMinChiSquaredPerHit(float beforeMinChiSquaredPerHit);
        void   SetAfterMinChiSquaredPerHit(float afterMinChiSquaredPerHit);
        void   SetMinChiSquaredPerHitChange(float chiSquaredPerHitChange);
        void   SetSplitPosition(float splitPosition);
        void   SetSplitApplied(bool splitApplied);
        void   SetBeforeDeltaChiSquaredPerHit(float beforeDeltaChiSquaredPerHit);

        int    GetBeforeNHits();
        int    GetAfterNHits();
        float  GetBeforeMinChiSquaredPerHit();
        float  GetAfterMinChiSquaredPerHit();
        float  GetMinChiSquaredPerHitChange();
        float  GetSplitPosition();
        bool   GetSplitApplied();
        float  GetBeforeDeltaChiSquaredPerHit();
    

    private:
        int                                         m_beforenhits;
        int                                         m_afternhits;
        float                                       m_beforeminchisquaredperhit;
        float                                       m_afterminchisquaredperhit;
        float                                       m_chisquaredperhitchange;
        float                                       m_splitposition;
        bool                                        m_splitapplied;
        float                                       m_beforedeltachisquaredperhit;        
    };

    class LookupTable
    {
    public:

        LookupTable();

        LookupTable(double &initialEnergy, double &binWidth);

        std::map<int, double> GetMap();

        void SetMap(std::map<int, double> &map);

        std::map<double, int> GetReverseMap();

        void SetReverseMap(std::map<double, int> &map);

        double GetInitialEnergy();

        void SetInitialEnergy(double &initialEnergy);

        double GetBinWidth();

        void SetBinWidth(double &binWidth);

        void SetMaxRange(double &maxRange);

        double GetMaxRange();

    private:
        std::map<int, double>                       m_map;
        std::map<double, int>                       m_reversemap;
        double                                      m_binwidth;               ///<
        double                                      m_initialenergy;
        double                                      m_maxrange;
    };

    class DirectionFitObject
    {
    public:

        DirectionFitObject();
        DirectionFitObject(HitChargeVector &hitChargeVector, int &numberHits, float &meanChargeOverWidth, float &forwardsChiSquared, float &backwardsChiSquared);
        DirectionFitObject(HitChargeVector &hitChargeVector, HitChargeVector &forwardsRecoHits, HitChargeVector &backwardsRecoHits, int &numberHits, float &meanChargeOverWidth, float &forwardsChiSquared, float &backwardsChiSquared);

        TrackDirectionTool::HitChargeVector GetHitChargeVector();
        TrackDirectionTool::HitChargeVector GetForwardsFitCharges();
        TrackDirectionTool::HitChargeVector GetBackwardsFitCharges();

        void SetForwardsFitCharges(TrackDirectionTool::HitChargeVector hitChargeVector);
        void SetBackwardsFitCharges(TrackDirectionTool::HitChargeVector hitChargeVector);
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

        void SetBeginpoint(const pandora::CartesianVector &beginPoint);
        void SetEndpoint(const pandora::CartesianVector &endPoint);
        const pandora::CartesianVector GetBeginpoint();
        const pandora::CartesianVector GetEndpoint();

        void SetProbability(float &probability);
        float GetProbability();
        void SetHypothesis(int hypothesis);
        int GetHypothesis();

        void SetSplitObject(SplitObject splitObject);
        SplitObject GetSplitObject();
        void SetTEFObject(SplitObject tefObject);
        SplitObject GetTEFObject();
        void SetFRObject(SplitObject frObject);
        SplitObject GetFRObject();

        void SetMCDirection(int direction);
        int GetMCDirection();

        void Print();

    private:
        HitChargeVector     m_hitchargevector;
        HitChargeVector     m_forwardsrecohits;
        HitChargeVector     m_backwardsrecohits;

        int                 m_nhits;
        int                 m_hypothesis;

        float               m_meanqoverx;
        float               m_forwardschisquared;
        float               m_backwardschisquared;
        float               m_probability;
    
        float               m_beginx; //Beginpoint is defined as the track endpoint with the lowest Z coordinate
        float               m_beginy;
        float               m_beginz;

        float               m_endx; //Endpoint is defined as the track endpoint with the highest Z coordinate
        float               m_endy;
        float               m_endz;

        int                 m_mcdirection;

        SplitObject         m_splitobject;
        SplitObject         m_tefobject;
        SplitObject         m_frobject;
    };

    class JumpObject
    {
        public:

            JumpObject(float &longitudinalPosition, float &jumpValue);
            JumpObject(float &longitudinalPosition, float &jumpValue, float &openingAngle);

            float GetLongitudinalPosition();

            float GetJumpValue();

            float GetOpeningAngle();

        private:
            float       m_longitudinalposition;
            float       m_jumpvalue;
            float       m_openingangle;

    };


    TrackDirectionTool::DirectionFitObject GetClusterDirection(const pandora::Cluster *const pTargetClusterW);

    TrackDirectionTool::DirectionFitObject GetPfoDirection(const pandora::ParticleFlowObject *const pPfo);

    void WriteLookupTableToTree(LookupTable &lookupTable);

    private:

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    //-----------------------------------------------------------------------------------------------

    unsigned int            m_slidingFitWindow;                 ///< The layer window for the sliding linear fits
    TwoDSlidingFitResultMap m_slidingFitResultMap;              ///< The sliding fit result map

    unsigned int            m_minClusterCaloHits;               ///< The min number of hits in base cluster selection method
    float                   m_minClusterLength;          ///< The min length (squared) in base cluster selection method

    int                     m_targetParticlePDG;
    int                     m_numberTrackEndHits;
    bool                    m_enableFragmentRemoval;
    bool                    m_enableSplitting;

    double                  m_tableInitialEnergy;
    double                  m_tableStepSize;

    bool                    m_writeTable;

    std::string             m_lookupTableFileName;
    std::string             m_probabilityFileName;
    std::string             m_treeName;

    //-----------------------------------------------------------------------------------------------


    const pandora::Cluster* GetTargetClusterFromPFO(const pandora::ParticleFlowObject* PFO, const LArTrackStateVector &trackStateVector);

    void SetEndpoints(DirectionFitObject &fitResult, const pandora::Cluster *const pCluster);

    void SetEndpoints(DirectionFitObject &fitResult, const LArTrackStateVector &trackStateVector);

    void SetMCTruth(DirectionFitObject &fitResult, const pandora::Cluster *const pCluster);

    void FillHitChargeVector(const pandora::Cluster *const pCluster, HitChargeVector &hitChargeVector);

    void TrackInnerFilter(HitChargeVector &hitChargeVector, HitChargeVector &filteredHitChargeVector);

    void SetNearestNeighbourValues(HitChargeVector &innerHitChargeVector, int &nNeighboursToConsider);

    void FragmentRemoval(HitChargeVector &hitChargeVector, HitChargeVector &filteredHitChargeVector, float &splitPosition);

    void SimpleTrackEndFilter(HitChargeVector &hitChargeVector);

    void TrackEndFilter(HitChargeVector &hitChargeVector, DirectionFitObject &directionFitObject);

    void AttemptFragmentRemoval(HitChargeVector &hitChargeVector, std::vector<JumpObject> &jumpsVector, HitChargeVector &filteredHitChargeVector, float &finalSplitPosition);

    void FindLargestJumps(HitChargeVector &hitChargeVector, std::vector<JumpObject> &leftJumps);

    void FindPeakJumps(HitChargeVector &hitChargeVector, std::vector<JumpObject> &peakJumps);

    void FindTrackEndJumps(HitChargeVector &hitChargeVector, std::vector<JumpObject> &trackEndJumps);

    void ParticleSplitting(HitChargeVector &hitChargeVector, DirectionFitObject &fitResult1, DirectionFitObject &fitResult2, bool &splitApplied, SplitObject &splitObject);

    void FindKinkSize(const pandora::Cluster *const pCluster, float &splitPosition, float &kinkSize);

    void CreateCalorimetricSplitHitVector(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions);

    void SplitHitCollectionBySize(HitChargeVector &hitChargeVector, float &splitPosition, HitChargeVector &smallHitChargeVector, HitChargeVector &largeHitChargeVector);

    void SplitHitCollectionByLeftRight(HitChargeVector &hitChargeVector, float &splitPosition, HitChargeVector &leftHitChargeVector, HitChargeVector &rightHitChargeVector);

    void GetTrackLength(HitChargeVector &hitChargeVector, float &trackLength);

    void GetAverageQoverWTrackBody(HitChargeVector &hitChargeVector, float &averageChargeTrackBody);

    void GetQoverWRange(HitChargeVector &hitChargeVector, float &QoverWRange);

    void FindKinkSplit(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions);

    void FindPlateauSplit(HitChargeVector &hitChargeVector, std::vector<JumpObject> &jumpObjects);

    void FindJumpSplit(HitChargeVector &hitChargeVector, std::vector<JumpObject> &jumpObjects);

    void FindBowlSplit(HitChargeVector &hitChargeVector, std::vector<float> &splitPositions);

    void FitHitChargeVector(HitChargeVector &hitChargeVector, TrackDirectionTool::DirectionFitObject &fitResult, int numberHitsToConsider=1000000);

    void FitHitChargeVector(HitChargeVector &hitChargeVector1, HitChargeVector &hitChargeVector2, TrackDirectionTool::DirectionFitObject &fitResult1, TrackDirectionTool::DirectionFitObject &fitResult2, int numberHitsToConsider=1000000);

    void ComputeProbability(DirectionFitObject &fitResult);

    void PerformFits(HitChargeVector &hitChargeVector, HitChargeVector &forwardsFitPoints, HitChargeVector &backwardsFitPoints, int numberHitsToConsider, float &forwardsChiSquared, float &backwardsChiSquared, int &fitStatus1, int &fitStatus2);

    void GetCalorimetricDirection(const pandora::Cluster* pTargetClusterW, DirectionFitObject &finalDirectionFitObject);

    void TestHypothesisOne(DirectionFitObject &directionFitObject);

    void TestHypothesisTwo(DirectionFitObject &directionFitObject);

    void TestHypothesisThree(DirectionFitObject &directionFitObject);

    void AddToSlidingFitCache(const pandora::Cluster *const pCluster);

    const TwoDSlidingFitResult &GetCachedSlidingFit(const pandora::Cluster *const pCluster) const;

    static bool SortHitChargeVectorByRL(HitCharge &hitCharge1, HitCharge &hitCharge2);
    static bool SortHitChargeVectorByChargeOverWidth(HitCharge &hitCharge1, HitCharge &hitCharge2);
    static bool SortByDistanceToNN(HitCharge &hitCharge1, HitCharge &hitCharge2);
    static bool SortJumpVector(JumpObject &jumpObject1, JumpObject &jumpObject2);

    void BinHitChargeVector(HitChargeVector &hitChargeVector, HitChargeVector &binnedHitChargeVector);
    double DensityCorrection(double &T, double &M);
    double BetheBloch(double &T, double &M);
    void FillLookupTable(LookupTable &lookupTable, double M);
    double GetEnergyfromLength(LookupTable &lookupTable, double &trackLength);
    double GetLengthfromEnergy(LookupTable &lookupTable, double &currentEnergy);


};

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitCharge::HitCharge(const pandora::CaloHit* caloHit, float &longitudinalPosition, float &hitWidth, float &hitCharge, float &uncertainty) :
    m_intails((hitCharge/hitWidth) <= 1.4 ? true : false),
    m_calohit(caloHit),
    m_longitudinalposition(longitudinalPosition),
    m_hitwidth(hitWidth),
    m_charge(hitCharge),
    m_qoverx(hitCharge/hitWidth),
    m_uncertainty(uncertainty),
    m_forwardsfitcharge(0.f),
    m_forwardssigma(0.f),
    m_forwardsdelta(0.f),
    m_forwardschisquared(0.f),
    m_backwardsfitcharge(0.f),
    m_backwardssigma(0.f),
    m_backwardsdelta(0.f),
    m_backwardschisquared(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitCharge::HitCharge() :
    m_intails(false),
    m_calohit(NULL),
    m_longitudinalposition(0.f),
    m_hitwidth(0.f),
    m_charge(0.f),
    m_qoverx(0.f),
    m_uncertainty(0.f),
    m_forwardsfitcharge(0.f),
    m_forwardssigma(0.f),
    m_forwardsdelta(0.f),
    m_forwardschisquared(0.f),
    m_backwardsfitcharge(0.f),
    m_backwardssigma(0.f),
    m_backwardsdelta(0.f),
    m_backwardschisquared(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHit* TrackDirectionTool::HitCharge::GetCaloHit() const
{
    return m_calohit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetLongitudinalPosition() const
{
    return m_longitudinalposition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetHitWidth() const
{
    return m_hitwidth;
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

inline void TrackDirectionTool::HitCharge::SetDistanceToNN(float &distance)
{
    m_distancetonearestneighbour = distance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetDistanceToNN() const
{
    return m_distancetonearestneighbour;
}

//------------------------------------------------------------------------------------------------------------------------------------------
inline void TrackDirectionTool::HitCharge::SetForwardsFitCharge(float &Q_fit_f) 
{
    m_forwardsfitcharge = Q_fit_f;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void TrackDirectionTool::HitCharge::SetForwardsSigma(float &f_sigma)
{
    m_forwardssigma = f_sigma;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void TrackDirectionTool::HitCharge::SetForwardsDelta(float &forwardsDelta)
{
    m_forwardsdelta = forwardsDelta;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void TrackDirectionTool::HitCharge::SetForwardsChiSquared(float &forwardsHitChisquared)
{
    m_forwardschisquared = forwardsHitChisquared;
}
//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::HitCharge::SetBackwardsFitCharge(float &Q_fit_b) 
{
    m_backwardsfitcharge = Q_fit_b;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void TrackDirectionTool::HitCharge::SetBackwardsSigma(float &b_sigma)
{
    m_backwardssigma= b_sigma;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void TrackDirectionTool::HitCharge::SetBackwardsDelta(float &backwardsDelta)
{
    m_backwardsdelta = backwardsDelta;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline void TrackDirectionTool::HitCharge::SetBackwardsChiSquared(float &backwardsHitChisquared)
{
    m_backwardschisquared = backwardsHitChisquared;
}
//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetForwardsFitCharge() 
{
    return m_forwardsfitcharge;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float TrackDirectionTool::HitCharge::GetForwardsSigma()
{
    return m_forwardssigma;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float TrackDirectionTool::HitCharge::GetForwardsDelta()
{
    return m_forwardsdelta;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float TrackDirectionTool::HitCharge::GetForwardsChiSquared()
{
    return m_forwardschisquared;
}
//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::HitCharge::GetBackwardsFitCharge() 
{
    return m_backwardsfitcharge;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float TrackDirectionTool::HitCharge::GetBackwardsSigma()
{
    return m_backwardssigma;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float TrackDirectionTool::HitCharge::GetBackwardsDelta()
{
    return m_backwardsdelta;
}
//------------------------------------------------------------------------------------------------------------------------------------------
inline float TrackDirectionTool::HitCharge::GetBackwardsChiSquared()
{
    return m_backwardschisquared;
}
//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::SplitObject::SplitObject() :
    m_beforenhits(0),
    m_afternhits(0),
    m_beforeminchisquaredperhit(0.f),
    m_afterminchisquaredperhit(0.f),
    m_chisquaredperhitchange(0.f),
    m_splitposition(0.f),
    m_splitapplied(false),
    m_beforedeltachisquaredperhit(0.f)        
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::SplitObject::SplitObject(int beforeNumberHits, int afterNumberHits, float beforeMinChiSquaredPerHit, float afterMinChiSquaredPerHit, float chiSquaredPerHitChange, float splitPosition) :
    m_beforenhits(beforeNumberHits),
    m_afternhits(afterNumberHits),
    m_beforeminchisquaredperhit(beforeMinChiSquaredPerHit),
    m_afterminchisquaredperhit(afterMinChiSquaredPerHit),
    m_chisquaredperhitchange(chiSquaredPerHitChange),
    m_splitposition(splitPosition),
    m_splitapplied(false),
    m_beforedeltachisquaredperhit(0.f)        
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::SplitObject::SplitObject(int beforeNumberHits, int afterNumberHits, float beforeMinChiSquaredPerHit, float afterMinChiSquaredPerHit, float chiSquaredPerHitChange, float splitPosition, bool splitApplied, float beforeDeltaChiSquaredPerHit) :
    m_beforenhits(beforeNumberHits),
    m_afternhits(afterNumberHits),
    m_beforeminchisquaredperhit(beforeMinChiSquaredPerHit),
    m_afterminchisquaredperhit(afterMinChiSquaredPerHit),
    m_chisquaredperhitchange(chiSquaredPerHitChange),
    m_splitposition(splitPosition),
    m_splitapplied(splitApplied),
    m_beforedeltachisquaredperhit(beforeDeltaChiSquaredPerHit)        
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void     TrackDirectionTool::SplitObject::SetBeforeNHits(int beforeNumberHits)
{
    m_beforenhits = beforeNumberHits;
}

inline void     TrackDirectionTool::SplitObject::SetAfterNHits(int afterNumberHits)
{
    m_afternhits = afterNumberHits;
}

inline void   TrackDirectionTool::SplitObject::SetBeforeMinChiSquaredPerHit(float beforeMinChiSquaredPerHit)
{
    m_beforeminchisquaredperhit = beforeMinChiSquaredPerHit;
}

inline void   TrackDirectionTool::SplitObject::SetAfterMinChiSquaredPerHit(float afterMinChiSquaredPerHit)
{
    m_afterminchisquaredperhit = afterMinChiSquaredPerHit;
}

inline void   TrackDirectionTool::SplitObject::SetMinChiSquaredPerHitChange(float chiSquaredPerHitChange)
{
    m_chisquaredperhitchange = chiSquaredPerHitChange;
}

inline void   TrackDirectionTool::SplitObject::SetSplitPosition(float splitPosition)
{
    m_splitposition = splitPosition;
}

inline void   TrackDirectionTool::SplitObject::SetSplitApplied(bool splitApplied)
{
    m_splitapplied = splitApplied;
}

inline void   TrackDirectionTool::SplitObject::SetBeforeDeltaChiSquaredPerHit(float beforeDeltaChiSquaredPerHit)
{
    m_beforedeltachisquaredperhit = beforeDeltaChiSquaredPerHit;
}

inline int     TrackDirectionTool::SplitObject::GetBeforeNHits()
{
    return m_beforenhits;
}

inline int     TrackDirectionTool::SplitObject::GetAfterNHits()
{
    return m_afternhits;
}

inline float   TrackDirectionTool::SplitObject::GetBeforeMinChiSquaredPerHit()
{
    return m_beforeminchisquaredperhit;
}

inline float   TrackDirectionTool::SplitObject::GetAfterMinChiSquaredPerHit()
{
    return m_afterminchisquaredperhit;
}

inline float   TrackDirectionTool::SplitObject::GetMinChiSquaredPerHitChange()
{
    return m_chisquaredperhitchange;
}

inline float   TrackDirectionTool::SplitObject::GetSplitPosition()
{
    return m_splitposition;
}

inline bool TrackDirectionTool::SplitObject::GetSplitApplied()
{
    return m_splitapplied;
}

inline float TrackDirectionTool::SplitObject::GetBeforeDeltaChiSquaredPerHit()
{
    return m_beforedeltachisquaredperhit;
}

inline TrackDirectionTool::LookupTable::LookupTable()
{
    std::map<int, double> emptyMap;
    std::map<double, int> emptyReverseMap;

    m_map = emptyMap;
    m_reversemap = emptyReverseMap;
    m_binwidth = 0.f;
    m_initialenergy = 0.f,
    m_maxrange = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::LookupTable::LookupTable(double &initialEnergy, double &binWidth)
{
    std::map<int, double> emptyMap;
    std::map<double, int> emptyReverseMap;

    m_map = emptyMap;
    m_reversemap = emptyReverseMap;
    m_binwidth = binWidth;
    m_initialenergy = initialEnergy,
    m_maxrange = 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::map<int, double> TrackDirectionTool::LookupTable::GetMap()
{
    return m_map;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetMap(std::map<int, double> &map)
{
    m_map = map;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::map<double, int> TrackDirectionTool::LookupTable::GetReverseMap()
{
    return m_reversemap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetReverseMap(std::map<double, int> &reverseMap)
{
    m_reversemap = reverseMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TrackDirectionTool::LookupTable::GetInitialEnergy()
{
    return m_initialenergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetInitialEnergy(double &initialEnergy)
{
    m_initialenergy = initialEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TrackDirectionTool::LookupTable::GetBinWidth()
{
    return m_binwidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetBinWidth(double &binWidth)
{
    m_binwidth = binWidth;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::LookupTable::SetMaxRange(double &maxRange)
{
    m_maxrange = maxRange;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double TrackDirectionTool::LookupTable::GetMaxRange()
{
    return m_maxrange;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::DirectionFitObject::DirectionFitObject()
{
    HitChargeVector emptyVector;
    m_hitchargevector = (emptyVector);
    m_forwardsrecohits = (emptyVector);
    m_backwardsrecohits = (emptyVector);
    m_nhits = 0;
    m_hypothesis = 0;
    m_meanqoverx = 0.f;
    m_forwardschisquared = 0.f;
    m_backwardschisquared = 0.f;
    m_probability = 0.5;
    m_mcdirection = -1;
    SplitObject splitObject;
    m_splitobject = splitObject;
    m_tefobject = splitObject;
    m_frobject = splitObject;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::DirectionFitObject::DirectionFitObject(HitChargeVector &hitChargeVector, int &numberHits, float &meanChargeOverWidth, float &forwardsChiSquared, float &backwardsChiSquared)
{
    HitChargeVector emptyVector;
    m_hitchargevector = hitChargeVector;
    m_forwardsrecohits = emptyVector;
    m_backwardsrecohits = emptyVector;
    m_nhits = numberHits;
    m_hypothesis = 0;
    m_meanqoverx = meanChargeOverWidth;
    m_forwardschisquared = forwardsChiSquared;
    m_backwardschisquared = backwardsChiSquared;
    m_probability = 0.5;
    m_mcdirection = -1;
    SplitObject splitObject;
    m_splitobject = splitObject;
    m_tefobject = splitObject;
    m_frobject = splitObject;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::DirectionFitObject::DirectionFitObject(HitChargeVector &hitChargeVector, HitChargeVector &forwardsRecoHits, HitChargeVector &backwardsRecoHits, int &numberHits, float &meanChargeOverWidth, float &forwardsChiSquared, float &backwardsChiSquared) :
    m_hitchargevector(hitChargeVector),
    m_forwardsrecohits(forwardsRecoHits),
    m_backwardsrecohits(backwardsRecoHits),
    m_nhits(numberHits),
    m_meanqoverx(meanChargeOverWidth),
    m_forwardschisquared(forwardsChiSquared),
    m_backwardschisquared(backwardsChiSquared)
{
    m_hypothesis = 0;
    m_probability = 0.5;
    m_mcdirection = -1;
    SplitObject splitObject;
    m_splitobject = splitObject;
    m_tefobject = splitObject;
    m_frobject = splitObject;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitChargeVector TrackDirectionTool::DirectionFitObject::GetHitChargeVector()
{
    return m_hitchargevector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitChargeVector TrackDirectionTool::DirectionFitObject::GetForwardsFitCharges()
{
    return m_forwardsrecohits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::HitChargeVector TrackDirectionTool::DirectionFitObject::GetBackwardsFitCharges()
{
    return m_backwardsrecohits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetForwardsFitCharges(TrackDirectionTool::HitChargeVector hitChargeVector)
{
    m_forwardsrecohits = hitChargeVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetBackwardsFitCharges(TrackDirectionTool::HitChargeVector hitChargeVector)
{
    m_backwardsrecohits = hitChargeVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetForwardsChiSquared()
{
    return m_forwardschisquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetBackwardsChiSquared()
{
    return m_backwardschisquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetForwardsChiSquared(float forwardsChiSquared)
{
    m_forwardschisquared = forwardsChiSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetBackwardsChiSquared(float backwardsChiSquared)
{
    m_backwardschisquared = backwardsChiSquared;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetForwardsChiSquaredPerHit()
{
    return m_forwardschisquared/m_nhits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetBackwardsChiSquaredPerHit()
{
    return m_backwardschisquared/m_nhits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetNHits(float numberHits)
{
    m_nhits = numberHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int TrackDirectionTool::DirectionFitObject::GetNHits()
{
    return m_hitchargevector.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int TrackDirectionTool::DirectionFitObject::GetDirectionEstimate()
{
    return (m_forwardschisquared <= m_backwardschisquared ? 1 : 0);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetMinChiSquared()
{
    return std::min(m_forwardschisquared, m_backwardschisquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetMinChiSquaredPerHit()
{
    return (m_nhits != 0 ? std::min(m_forwardschisquared, m_backwardschisquared)/m_nhits : std::numeric_limits<float>::max());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetDeltaChiSquaredPerHit()
{
    return (m_nhits != 0 ? (m_forwardschisquared - m_backwardschisquared)/m_nhits : std::numeric_limits<float>::max());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetMeanChargeOverWidth()
{
    return m_meanqoverx;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetBeginpoint(const pandora::CartesianVector &beginPoint)
{
    m_beginx = beginPoint.GetX();
    m_beginy = beginPoint.GetY();
    m_beginz = beginPoint.GetZ();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetEndpoint(const pandora::CartesianVector &endPoint)
{
    m_endx = endPoint.GetX();
    m_endy = endPoint.GetY();
    m_endz = endPoint.GetZ();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector TrackDirectionTool::DirectionFitObject::GetBeginpoint()
{
    const pandora::CartesianVector beginPoint(m_beginx, m_beginy, m_beginz);
    return beginPoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector TrackDirectionTool::DirectionFitObject::GetEndpoint()
{
    const pandora::CartesianVector endPoint(m_endx, m_endy, m_endz);
    return endPoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetProbability(float &probability)
{
    m_probability = probability;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::DirectionFitObject::GetProbability()
{
    return m_probability;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetHypothesis(int hypothesis)
{
    m_hypothesis = hypothesis;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int TrackDirectionTool::DirectionFitObject::GetHypothesis()
{
    return m_hypothesis;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetSplitObject(TrackDirectionTool::SplitObject splitObject)
{
    m_splitobject = splitObject;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::SplitObject TrackDirectionTool::DirectionFitObject::GetSplitObject()
{
    return m_splitobject;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetTEFObject(TrackDirectionTool::SplitObject tefObject)
{
    m_tefobject = tefObject;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::SplitObject TrackDirectionTool::DirectionFitObject::GetTEFObject()
{
    return m_tefobject;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetFRObject(TrackDirectionTool::SplitObject frObject)
{
    m_frobject = frObject;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::SplitObject TrackDirectionTool::DirectionFitObject::GetFRObject()
{
    return m_frobject;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::SetMCDirection(int direction)
{
    m_mcdirection = direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int TrackDirectionTool::DirectionFitObject::GetMCDirection()
{
    return m_mcdirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void TrackDirectionTool::DirectionFitObject::Print()
{
    std::cout << "Probability: " << m_probability << std::endl;
    std::cout << "Vertex position: (" << m_beginx << ", " << m_beginy << ", " << m_beginz << ")" << std::endl;
    std::cout << "Endpoint position: (" << m_endx << ", " << m_endy << ", " << m_endz << ")" << std::endl;
    std::cout << "Best chi squared per hit: " << std::min(m_forwardschisquared, m_backwardschisquared)/m_nhits << std::endl;
    std::cout << "Number of hits: " << m_nhits << std::endl;
    std::cout << "Hypothesis: " << m_hypothesis << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::JumpObject::JumpObject(float &longitudinalPosition, float &jumpValue) :
    m_longitudinalposition(longitudinalPosition),
    m_jumpvalue(jumpValue),
    m_openingangle(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackDirectionTool::JumpObject::JumpObject(float &longitudinalPosition, float &jumpValue, float &openingAngle) :
    m_longitudinalposition(longitudinalPosition),
    m_jumpvalue(jumpValue),
    m_openingangle(openingAngle)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::JumpObject::GetLongitudinalPosition()
{
    return m_longitudinalposition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::JumpObject::GetJumpValue()
{
    return m_jumpvalue;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackDirectionTool::JumpObject::GetOpeningAngle()
{
    return m_openingangle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_TRACK_DIRECTION_TOOL_H
