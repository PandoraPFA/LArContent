/**
 *  @file   larpandoracontent/LArControlFlow/TestBeamCosmicRayTaggingTool.h
 *
 *  @brief  Header file for the cosmic-ray tagging tool class.
 *
 *  $Log: $
 */
#ifndef LAR_TEST_BEAM_COSMIC_RAY_TAGGING_TOOL_H
#define LAR_TEST_BEAM_COSMIC_RAY_TAGGING_TOOL_H 1

#include "larpandoracontent/LArControlFlow/CosmicRayTaggingBaseTool.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  TestBeamCosmicRayTaggingTool class
 */
class TestBeamCosmicRayTaggingTool : public CosmicRayTaggingBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TestBeamCosmicRayTaggingTool();

    pandora::StatusCode Initialize();
    void FindAmbiguousPfos(const pandora::PfoList &parentCosmicRayPfos, pandora::PfoList &ambiguousPfos, const MasterAlgorithm *const pAlgorithm);

private:
    /**
     *  @brief  Class to encapsulate the logic required determine if a Pfo should or shouldn't be tagged as a cosmic ray
     */
    class CRCandidate
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pandora the relevant pandora instance
         *  @param  pPfo the address of the candidate pfo
         *  @param  slice the slice id
         */
        CRCandidate(const pandora::Pandora &pandora, const pandora::ParticleFlowObject *const pPfo, const unsigned int sliceId);

        const pandora::ParticleFlowObject *const m_pPfo; ///< Address of the candidate Pfo
        unsigned int m_sliceId;                          ///< Slice ID
        bool m_canFit;                                   ///< If there are a sufficient number of 3D hits to perform a fitting
        pandora::CartesianVector m_endPoint1;            ///< First fitted end point in 3D
        pandora::CartesianVector m_endPoint2;            ///< Second fitted end point in 3D
        double m_length;                                 ///< Straight line length of the linear fit
        double m_curvature;                              ///< Measure of the curvature of the track
        double m_theta;                                  ///< Direction made with vertical

    private:
        /**
         *  @brief  Calculate all variables which require a fit
         *
         *  @param  slidingFitResult the three dimensional sliding fit result
         */
        void CalculateFitVariables(const ThreeDSlidingFitResult &slidingFitResult);
    };

    typedef std::list<CRCandidate> CRCandidateList;

    /**
     *  @brief  Get the 3D calo hit cluster associated with a given Pfo, and check if it has sufficient hits
     *
     *  @param  pPfo input Pfo
     *  @param  pCluster3D to receive the address of the 3D cluster
     *
     *  @return whether the Pfo has sufficient hits?
     */
    bool GetValid3DCluster(const pandora::ParticleFlowObject *const pPfo, const pandora::Cluster *&pCluster3D) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject *, pandora::PfoList> PfoToPfoListMap;

    /**
     *  @brief  Get mapping between Pfos that are associated with it other by pointing
     *
     *  @param  parentCosmicRayPfos input list of Pfos
     *  @param  pfoAssociationsMap to receive the output mapping between associated Pfos
     */
    void GetPfoAssociations(const pandora::PfoList &parentCosmicRayPfos, PfoToPfoListMap &pfoAssociationMap) const;

    /**
     *  @brief  Check whethe two Pfo endpoints are associated by distance of closest approach
     *
     *  @param  endPoint1 position vector of an endpoint of Pfo 1
     *  @param  endDir1 direction vector of an endpoint of Pfo 1
     *  @param  endPoint2 position vector of an endpoint of Pfo 2
     *  @param  endDir2 direction vector of an endpoint of Pfos
     *
     *  @return whether the Pfos are associated
     */
    bool CheckAssociation(const pandora::CartesianVector &endPoint1, const pandora::CartesianVector &endDir1,
        const pandora::CartesianVector &endPoint2, const pandora::CartesianVector &endDir2) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject *, unsigned int> PfoToSliceIdMap;

    /**
     *  @brief  Break the event up into slices of associated Pfos
     *
     *  @param  parentCosmicRayPfos input list of Pfos
     *  @param  pfoAssociationMap mapping between Pfos and other associated Pfos
     *  @param  pfoToSliceIdMap to receive the mapping between Pfos and their slice ID
     */
    void SliceEvent(const pandora::PfoList &parentCosmicRayPfos, const PfoToPfoListMap &pfoAssociationMap, PfoToSliceIdMap &pfoToSliceIdMap) const;

    /**
     *  @brief  Fill a slice iteratively using Pfo associations
     *
     *  @param  pPfo Pfo to add to the slice
     *  @param  pfoAssociationMap mapping between Pfos and other associated Pfos
     *  @param  slice the slice to add Pfos to
     */
    void FillSlice(const pandora::ParticleFlowObject *const pPfo, const PfoToPfoListMap &pfoAssociationMap, pandora::PfoList &slice) const;

    /**
     *  @brief  Make a list of CRCandidates
     *
     *  @param  parentCosmicRayPfos input list of Pfos
     *  @param  pfoToSliceIdMap input mapping between Pfos and their slice id
     *  @param  candidates to receive the output list of CRCandidates
     */
    void GetCRCandidates(const pandora::PfoList &parentCosmicRayPfos, const PfoToSliceIdMap &pfoToSliceIdMap, CRCandidateList &candidates) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject *, bool> PfoToBoolMap;

    /**
     *  @brief  Check if each candidate is "in time"
     *
     *  @param  candidates input list of candidates
     *  @param  pfoToInTimeMap output mapping between candidates Pfos and if they are in time
     */
    void CheckIfOutOfTime(const CRCandidateList &candidates, PfoToBoolMap &pfoToInTimeMap) const;

    /**
     *  @brief  Check if each candidate is "top to bottom"
     *
     *  @param  candidates input list of candidates
     *  @param  pfoToIsCosmicRayMap output mapping between candidates Pfos and if they are top to bottom
     */
    void CheckIfTopToBottom(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsCosmicRayMap) const;

    /**
     *  @brief  Check if each candidate enters from the top of the detector
     *
     *  @param  candidates input list of candidates
     *  @param  pfoToIsCosmicRayMap output mapping between candidates Pfos and if they are top entering
     */
    void CheckIfTopEntering(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsCosmicRayMap) const;

    /**
     *  @brief  Check if each candidate has its highest value in one of the vetoed TPC volumes
     *
     *  @param  candidates input list of candidates
     *  @param  pfoToIsCosmicRayMap output mapping between candidates Pfos and if they are in the vetoed TPCs
     */
    void CheckIfInVetoedTPC(const CRCandidateList &candidates, PfoToBoolMap &pfoToIsCosmicRayMap) const;

    typedef std::set<unsigned int> UIntSet;
    typedef std::unordered_map<int, bool> IntBoolMap;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::pair<const ThreeDSlidingFitResult, const ThreeDSlidingFitResult> SlidingFitPair;
    typedef std::unordered_map<const pandora::ParticleFlowObject *, SlidingFitPair> PfoToSlidingFitsMap;
    typedef std::vector<pandora::PfoList> SliceList;

    float m_angularUncertainty;    ///< The uncertainty in degrees for the angle of a Pfo
    float m_positionalUncertainty; ///< The uncertainty in cm for the position of Pfo endpoint in 3D
    float m_maxAssociationDist; ///< The maximum distance from endpoint to point of closest approach, typically a multiple of LAr radiation length

    unsigned int m_minimumHits; ///< The minimum number of hits for a Pfo to be considered

    float m_inTimeMargin; ///< The maximum distance outside of the physical detector volume that a Pfo may be to still be considered in time
    float m_inTimeMaxX0;  ///< The maximum pfo x0 (determined from shifted vertex) to allow pfo to still be considered in time
    float m_marginY;      ///< The minimum distance from a detector Y-face for a Pfo to be associated

    bool m_tagTopEntering; ///< Whether to tag all top entering particles as cosmic rays
    bool m_tagTopToBottom; ///< Whether to tag all top-to-bottom particles as cosmic rays
    bool m_tagOutOfTime;   ///< Whether to tag all out-of-time particles as cosmic rays

    bool m_tagInVetoedTPCs;                 ///< Whether to tag all particles with their highest position in vetoed TPCs as cosmic rays
    std::vector<unsigned int> m_vetoedTPCs; ///< List of vetoed TPCs for tagging cosmic rays

    float m_face_Xa; ///< Anode      X face
    float m_face_Xc; ///< Cathode    X face
    float m_face_Yb; ///< Bottom     Y face
    float m_face_Yt; ///< Top        Y face
    float m_face_Zu; ///< Upstream   Z face
    float m_face_Zd; ///< Downstream Z face
};

} // namespace lar_content

#endif // #ifndef LAR_TETS_BEAM_COSMIC_RAY_TAGGING_TOOL_H
