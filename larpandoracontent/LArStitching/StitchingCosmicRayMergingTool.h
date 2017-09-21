/**
 *  @file   LArContent/include/LArStitching/StitchingCosmicRayMergingTool.h
 *
 *  @brief  Header file for the stitching pfo merging tool class.
 *
 *  $Log: $
 */
#ifndef LAR_STITCHING_COSMIC_RAY_MERGING_TOOL_H
#define LAR_STITCHING_COSMIC_RAY_MERGING_TOOL_H 1

#include "larpandoracontent/LArStitching/MultiPandoraApi.h"
#include "larpandoracontent/LArStitching/StitchingAlgorithm.h"

#include "larpandoracontent/LArObjects/LArPointingCluster.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  StitchingCosmicRayMergingTool class
 */
class StitchingCosmicRayMergingTool : public StitchingTool
{
public:
    /**
     *  @brief  Default constructor
     */
    StitchingCosmicRayMergingTool();

    void Run(const StitchingAlgorithm *const pAlgorithm, StitchingAlgorithm::StitchingInfo &stitchingInfo);

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
     *  @param  outputPfoList the output list of Pfos
     */
    void SelectPrimaryPfos(const pandora::PfoList *pInputPfoList, pandora::PfoList &outputPfoList) const;

    typedef std::unordered_map<const pandora::ParticleFlowObject*, LArPointingCluster> ThreeDPointingClusterMap;

    /**
     *  @brief  Build a 3D pointing cluster for each Pfo
     *
     *  @param  inputPfoList  the input list of Pfos
     *  @param  pointingClusterMap  the mapping between Pfos and their corresponding 3D pointing clusters
     */
    void BuildPointingClusterMaps(const pandora::PfoList &inputPfoList, ThreeDPointingClusterMap &pointingClusterMap) const;

    typedef StitchingAlgorithm::PfoToLArTPCMap PfoToLArTPCMap;
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

    /**
     *  @brief  Apply X0 corrections, and then stitch together Pfos
     *
     *  @param  pAlgorithm the address of the parent stitching algorithm
     *  @param  pointingClusterMap  the mapping between Pfos and their corresponding 3D pointing clusters
     *  @param  pfoMerges the input map of Pfo merges
     *  @param  stitchingInfo the stitching information block
     */
    void StitchPfos(const StitchingAlgorithm *const pAlgorithm, const ThreeDPointingClusterMap &pointingClusterMap,
        const PfoMergeMap &pfoMerges, StitchingAlgorithm::StitchingInfo &stitchingInfo) const;

    /**
     *  @brief  Calculate x0 shift for a group of associated Pfos
     *
     *  @param  pfoToLArTPCMap the mapping between pfos and tpc
     *  @param  pointingClusterMap the mapping between Pfos and their corresponding 3D pointing clusters
     *  @param  pfoVector the vector of parent Pfos to stitch together
     *  @param  x0 the output x0 value
     */
    void CalculateX0(const PfoToLArTPCMap &pfoToLArTPCMap, const ThreeDPointingClusterMap &pointingClusterMap,
        const pandora::PfoVector &pfoVector, float &x0) const;

    bool  m_useXcoordinate;
    int   m_halfWindowLayers;
    float m_minLengthSquared;
    float m_minCosRelativeAngle;
    float m_maxLongitudinalDisplacementX;
    float m_maxTransverseDisplacement;
    float m_relaxCosRelativeAngle;
    float m_relaxTransverseDisplacement;
};

} // namespace lar_content

#endif // #ifndef LAR_STITCHING_COSMIC_RAY_MERGING_TOOL_H
