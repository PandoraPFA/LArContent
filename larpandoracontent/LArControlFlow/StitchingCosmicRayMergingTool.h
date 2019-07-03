/**
 *  @file   LArContent/include/LArControlFlow/StitchingCosmicRayMergingTool.h
 *
 *  @brief  Header file for the stitching pfo merging tool class.
 *
 *  $Log: $
 */
#ifndef LAR_STITCHING_COSMIC_RAY_MERGING_TOOL_H
#define LAR_STITCHING_COSMIC_RAY_MERGING_TOOL_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  StitchingCosmicRayMergingTool class
 */
class StitchingCosmicRayMergingTool : public StitchingBaseTool
{
public:
    /**
     *  @brief  Default constructor
     */
    StitchingCosmicRayMergingTool();

    void Run(const MasterAlgorithm *const pAlgorithm, const pandora::PfoList *const pMultiPfoList, PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map);

    /**
     *  @brief  PfoAssociation class
     */
    class PfoAssociation
    {
    public:
        /**
         *  @brief  Vertex enumeration
         */
        enum VertexType
        {
            UNDEFINED = 0,
            INNER     = 1,
            OUTER     = 2
        };

        /**
         *  @brief  Constructor
         *
         *  @param  parent
         *  @param  daughter
         *  @param  fom
         */
        PfoAssociation(const VertexType parent, const VertexType daughter, const float fom);

        /**
         *  @brief  Get parent
         *
         *  @return the parent
         */
        VertexType GetParent() const;

        /**
         *  @brief  Get daughter
         *
         *  @return the daughter
         */
        VertexType GetDaughter() const;

        /**
         *  @brief  Get figure of merit
         *
         *  @return the figure of merit
         */
        float GetFigureOfMerit() const;

    private:
        VertexType      m_parent;           ///<
        VertexType      m_daughter;         ///<
        float           m_fom;              ///<
    };

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Select primary Pfos from the input list of Pfos
     *
     *  @param  pInputPfoList the input list of Pfos
     *  @param  pfoToLArTPCMap the input mapping between Pfos and tpc
     *  @param  outputPfoList the output list of Pfos
     */
    void SelectPrimaryPfos(const pandora::PfoList *pInputPfoList, const PfoToLArTPCMap &pfoToLArTPCMap, pandora::PfoList &outputPfoList) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject*, LArPointingCluster> ThreeDPointingClusterMap;

    /**
     *  @brief  Build a 3D pointing cluster for each Pfo
     *
     *  @param  inputPfoList the input list of Pfos
     *  @param  pfoToLArTPCMap the input mapping between Pfos and tpc
     *  @param  pointingClusterMap  the mapping between Pfos and their corresponding 3D pointing clusters
     */
    void BuildPointingClusterMaps(const pandora::PfoList &inputPfoList, const PfoToLArTPCMap &pfoToLArTPCMap, ThreeDPointingClusterMap &pointingClusterMap) const;

    typedef std::unordered_map<const pandora::LArTPC*, pandora::PfoList> LArTPCToPfoMap;

    /**
     *  @brief  Build a list of Pfos for each tpc
     *
     *  @param  inputPfoList the input list of Pfos
     *  @param  pfoToLArTPCMap the input mapping between Pfos and tpc
     *  @param  larTPCToPfoMap the output mapping between tpc and Pfos
     */
    void BuildTPCMaps(const pandora::PfoList &inputPfoList, const PfoToLArTPCMap &pfoToLArTPCMap, LArTPCToPfoMap &larTPCToPfoMap) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject*, PfoAssociation> PfoAssociationMap;
    typedef std::unordered_map<const pandora::ParticleFlowObject*, PfoAssociationMap> PfoAssociationMatrix;

    /**
     *  @brief  Create associations between Pfos using 3D pointing clusters
     *
     *  @param  larTPCToPfoMap the input mapping between tpc and Pfos
     *  @param  pointingClusterMap the mapping between Pfos and their corresponding 3D pointing clusters
     *  @param  pfoAssociationMatrix the output matrix of associations between Pfos
     */
    void CreatePfoMatches(const LArTPCToPfoMap &larTPCToPfoMap, const ThreeDPointingClusterMap &pointingClusterMap,
        PfoAssociationMatrix &pfoAssociationMatrix) const;

    /**
     *  @brief  Create associations between Pfos using 3D pointing clusters
     *
     *  @param  larTPC1 the tpc description for the first Pfo
     *  @param  larTPC2 the tpc description for the second Pfo
     *  @param  pPfo1 the first Pfo
     *  @param  pPfo2 the second Pfo
     *  @param  pointingClusterMap the mapping between Pfos and their corresponding 3D pointing clusters
     *  @param  pfoAssociationMatrix the output matrix of associations between Pfos
     */
    void CreatePfoMatches(const pandora::LArTPC &larTPC1, const pandora::LArTPC &larTPC2, const pandora::ParticleFlowObject *const pPfo1,
        const pandora::ParticleFlowObject *const pPfo2, const ThreeDPointingClusterMap &pointingClusterMap, PfoAssociationMatrix &pfoAssociationMatrix) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject*, pandora::PfoList> PfoMergeMap;

    /**
     *  @brief  Select the best associations between Pfos; create a mapping between associated Pfos, handling any ambiguities
     *
     *  @param  pfoAssociationMatrix the input list of all associations between Pfos
     *  @param  pfoSelectedMatches the output list of good associations between pfos (candidates for merging)
     */
    void SelectPfoMatches(const PfoAssociationMatrix &pfoAssociationMatrix, PfoMergeMap &pfoSelectedMatches) const;

    /**
     *  @brief  Create an initial map of Pfo merges to be made
     *
     *  @param  pfoMatches the input list of associated Pfos
     *  @param  pfoMerges the output mapping between Pfos to enlarge and Pfos to delete
     */
    void SelectPfoMerges(const PfoMergeMap &pfoMatches, PfoMergeMap &pfoMerges) const;

    /**
     *  @brief  Collect up associations between Pfos
     *
     *  @param  pSeedPfo the seed Pfo (will be enlarged after merging)
     *  @param  pCurrentPfo the target Pfo (will be deleted after merging)
     *  @param  pfoMerges the list of associations between Pfos
     *  @param  vetoSet the set of associated Pfos that have already been considered
     *  @param  associatedList the output list of associated Pfos
     */
    void CollectAssociatedPfos(const pandora::ParticleFlowObject *const pSeedPfo, const pandora::ParticleFlowObject *const pCurrentPfo,
        const PfoMergeMap &pfoMerges, const pandora::PfoSet &vetoSet, pandora::PfoList &associatedList) const;

    /**
     *  @brief  Identify the vertex Pfo and then re-order the map of merges so that the vertex Pfo will be enlarged
     *
     *  @param  pfoToLArTPCMap the mapping between pfos and tpc
     *  @param  pointingClusterMap the mapping between Pfos and their corresponding 3D pointing clusters
     *  @param  inputPfoMerges the input map of Pfo merges
     *  @param  outputPfoMerges the re-ordered map of Pfo merges
     */
    void OrderPfoMerges(const PfoToLArTPCMap &pfoToLArTPCMap, const ThreeDPointingClusterMap &pointingClusterMap,
        const PfoMergeMap &inputPfoMerges, PfoMergeMap &outputPfoMerges) const;

    typedef std::pair<const pandora::LArTPC*, const pandora::LArTPC*> LArTPCPair;
    typedef std::map<const pandora::ParticleFlowObject*, LArPointingCluster::Vertex> PfoToPointingVertexMap;

    /**
     *  @brief  Apply X0 corrections, and then stitch together Pfos
     *
     *  @param  pAlgorithm the address of the parent stitching algorithm
     *  @param  pointingClusterMap  the mapping between Pfos and their corresponding 3D pointing clusters
     *  @param  pfoMerges the input map of Pfo merges
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     *  @param  stitchedPfosToX0Map a map of cosmic-ray pfos that have been stitched between lar tpcs to the X0 shift
     */
    void StitchPfos(const MasterAlgorithm *const pAlgorithm, const ThreeDPointingClusterMap &pointingClusterMap,
        const PfoMergeMap &pfoMerges, PfoToLArTPCMap &pfoToLArTPCMap, PfoToFloatMap &stitchedPfosToX0Map) const;

    /**
     * @brief Reduce the original pfoVector to one of size 2 if its greater than that
     *
     * @param pfoVector vector of pfos being stitched
     * @param reducedPfoVector the reduced vector of pfos
     * @param pPfoToEnlarge the pfo we are enlarging
     * @param pfoToLArTPCMap the pfo to lar tpc map
     *
     * @return particleflow object to enlarge
     */
    const pandora::ParticleFlowObject *ReduceToLongestStitch(pandora::PfoVector &pfoVector, pandora::PfoVector &reducedPfoVector,
        const pandora::ParticleFlowObject *const pPfoToEnlarge, const PfoToLArTPCMap &pfoToLArTPCMap) const;

    // TODO - need doxygen info
    pandora::PfoVector SelectLongestStitch(pandora::PfoVector &pfoVector, const PfoToLArTPCMap &pfoToLArTPCMap) const;
    const pandora::ParticleFlowObject *GetClosestPfo(const pandora::ParticleFlowObject *const pPfoToEnlarge, const pandora::PfoVector &pfoVector) const;

    unsigned int NumberOfTwoDHits(const pandora::ParticleFlowObject *pfo) const;

    /**
     *  @brief  Find the pair of LArTPCs that contain the pfos being stitched
     *
     *  @param  pfoVector vector of pfos being stitched
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     *  @param  stitchedLArTPCs the pair of LArTPCs containing the pfos being stitched
     */
    void FindStitchedLArTPCs(const pandora::PfoVector &pfoVector, const PfoToLArTPCMap &pfoToLArTPCMap, LArTPCPair &stitchedLArTPCs) const;

    /**
     *  @brief  Calculate x0 shift for a group of associated Pfos
     *
     *  @param  pfoToLArTPCMap the mapping between pfos and tpc
     *  @param  pointingClusterMap the mapping between Pfos and their corresponding 3D pointing clusters
     *  @param  pfoVector the vector of parent Pfos to stitch together
     *  @param  x0 the output x0 value
     *  @param  pfoToPointingVertexMap map of pfo to pointing vertex used in stitching
     */
    void CalculateX0(const PfoToLArTPCMap &pfoToLArTPCMap, const ThreeDPointingClusterMap &pointingClusterMap,
        const pandora::PfoVector &pfoVector, float &x0, PfoToPointingVertexMap &pfoToPointingVertexMap) const;

    bool            m_useXcoordinate;
    bool            m_alwaysApplyT0Calculation;
    int             m_halfWindowLayers;
    float           m_minLengthSquared;
    float           m_minCosRelativeAngle;
    float           m_maxLongitudinalDisplacementX;
    float           m_maxTransverseDisplacement;
    float           m_relaxCosRelativeAngle;
    float           m_relaxTransverseDisplacement;
    unsigned int    m_minNCaloHits3D;
};

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_COSMIC_RAY_MERGING_TOOL_H
