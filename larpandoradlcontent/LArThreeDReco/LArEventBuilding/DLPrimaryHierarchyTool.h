/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLPrimaryHierarchyTool.h
 *
 *  @brief  Header file for the DL primary hierarchy tool
 *
 *  $Log: $
 */
#ifndef LAR_DL_PRIMARY_HIERARCHY_TOOL_H
#define LAR_DL_PRIMARY_HIERARCHY_TOOL_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLBaseHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"

namespace lar_dl_content
{
/**
 *   @brief  DLPrimaryHierarchyTool to apply the primary hierarchy DL networks
 */
class DLPrimaryHierarchyTool : public DLBaseHierarchyTool
{
public:
    struct DLPrimaryNetworkParams // all floats because they'll be normalised
    {
        float m_nSpacepoints = -999.f;             ///< the number of 3D spacepoints
        float m_nuSeparation = -999.f;             ///< the sep. between the POI and the nu vertex
        float m_vertexRegionNHits = -999.f;        ///< the number of 3D hits 'close to' the POI
        float m_vertexRegionNParticles = -999.f;   ///< the number of different particles 'close to' the POI
        float m_dca = -999.f;                      ///< the distance of closest approach from the POI to the nu vertex
        float m_connectionExtrapDistance = -999.f; ///< the sep. between the POI and the DCA point
        float m_isPOIClosestToNu = -999.f;         ///< whether the POI is that closest to the nu vertex
        float m_parentConnectionDistance = -999.f; ///< the DCA to the most likely parent pfo endpoint
        float m_childConnectionDistance = -999.f;  ///< the sep. between the POI and the extension point for m_parentConnectionDistance

        /**
         *  @brief  Add the orientation independent primary network parameters
         *          to the model input tensor
         *
         *  @param  insertIndex the index from which to begin the insert
         *  @param  modelInput the torch vector to which to add
         */
        void AddCommonParamsToInput(int &insertIndex, LArDLHelper::TorchInput &modelInput) const;

        /**
         *  @brief  Add the orientation dependent primary network parameters
         *          to the model input tensor
         *
         *  @param  insertIndex the index from which to begin the insert
         *  @param  modelInput the torch vector to which to add
         */
        void AddOrientationParamsToInput(int &insertIndex, LArDLHelper::TorchInput &modelInput) const;
    };

    struct NormalisationLimits
    {
        float m_nSpacepointsMin = 0.f;                ///< the m_nSpacepoints minimum range boundary
        float m_nSpacepointsMax = 2000.f;             ///< the m_nSpacepoints maximum range boundary
        float m_nuSeparationMin = -50.f;              ///< the m_nuSeparation minimum range boundary
        float m_nuSeparationMax = 500.f;              ///< the m_nuSeparation maximum range boundary
        float m_vertexRegionNHitsMin = -10.f;         ///< the m_vertexRegionNHits minimum range boundary
        float m_vertexRegionNHitsMax = 100.f;         ///< the m_vertexRegionNHits maximum range boundary
        float m_vertexRegionNParticlesMin = -1.f;     ///< the m_vertexRegionNParticles minimum range boundary
        float m_vertexRegionNParticlesMax = 8.f;      ///< the m_vertexRegionNParticles maximum range boundary
        float m_dcaMin = -60.f;                       ///< the m_dca minimum range boundary
        float m_dcaMax = 600.f;                       ///< the m_dca maximum range boundary
        float m_connectionExtrapDistanceMin = -700.f; ///< the m_connectionExtrapDistance minimum range boundary
        float m_connectionExtrapDistanceMax = 500.f;  ///< the m_connectionExtrapDistance maximum range boundary
        float m_parentConnectionDistanceMin = -150.f; ///< the m_parentConnectionDistance minimum range boundary
        float m_parentConnectionDistanceMax = 150.f;  ///< the m_parentConnectionDistance maximum range boundary
        float m_childConnectionDistanceMin = -30.f;   ///< the m_childConnectionDistance minimum range boundary
        float m_childConnectionDistanceMax = 300.f;   ///< the m_childConnectionDistance maximum range boundary
    };

    pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pNeutrinoPfo,
        const HierarchyPfoVector &trackPfos, const HierarchyPfo &hierarchyPfo, std::vector<DLPrimaryNetworkParams> &networkParamVector,
        float &primaryScore);

    /**
     *  @brief  Default constructor
     */
    DLPrimaryHierarchyTool();

    /**
     *  @brief  Calculate the variables describing the extension of a child particle to a given parent 
     *          (m_parentConnectionDistance, m_childConnectionDistance)
     *
     *  @param  parentPoint the extremal point of the parent
     *  @param  childPoint the extremal point of the child
     *  @param  parentConnectionDistance the DCA to the parent vertex
     *  @param  childConnectionDistance the extension of the child to the DCA point
     */
    void CalculateConnectionDistances(const ExtremalPoint &parentPoint, const ExtremalPoint &childPoint, float &parentConnectionDistance,
        float &childConnectionDistance) const;

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Function to calculate the DLPrimaryNetworkParams
     *
     *  @param  pAlgorithm a pointer to the pandora algorithm
     *  @param  hierarchyPfo the input hierarchy pfo object
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  trackPfos the HierachyPfoVector of track-like particles
     *  @param  useUpstream whether the POI is the endpoint closest to the nu vertex
     *  @param  primaryNetworkParams the primary network parameters to fill
     *
     *  @return a StatusCode to signify whether the network variables could be correctly calculated
     */
    pandora::StatusCode CalculateNetworkVariables(const pandora::Algorithm *const pAlgorithm, const HierarchyPfo &hierarchyPfo,
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoVector &trackPfos, const bool useUpstream,
        DLPrimaryNetworkParams &primaryNetworkParams) const;

    /**
     *  @brief  Set the vertex region DLPrimaryNetworkParams params
     *          (m_vertexRegionNHits, m_vertexRegionNParticles)
     *
     *  @param  pAlgorithm a pointer to the pandora algorithm
     *  @param  pPfo a pointer to the input pfo
     *  @param  particleVertex the position of the pfo vertex
     *  @param  primaryNetworkParams the primary network parameters
     */
    void SetVertexRegionParams(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo,
        const pandora::CartesianVector &particleVertex, DLPrimaryNetworkParams &primaryNetworkParams) const;

    /**
     *  @brief  Set the connection region DLPrimaryNetworkParams params
     *          (m_dca, m_connectionExtrapDistance)
     *
     *  @param  particlePoint the extremal point of the particle
     *  @param  nuVertex the neutrino vertex
     *  @param  primaryNetworkParams the primary network parameters
     */
    void SetConnectionParams(
        const ExtremalPoint &particlePoint, const pandora::CartesianVector &nuVertex, DLPrimaryNetworkParams &primaryNetworkParams) const;

    /**
     *  @brief  Set the event context DLPrimaryNetworkParams params
     *          (m_parentConnectionDistance, m_childConnectionDistance)
     *
     *  @param  pPfo a pointer to the input pfo
     *  @param  particlePoint the extremal point of the particle
     *  @param  nuVertex the neutrino vertex
     *  @param  trackPfos the HierachyPfoVector of track-like particles
     *  @param  primaryNetworkParams the primary network parameters
     */
    void SetContextParams(const pandora::ParticleFlowObject *const pPfo, const ExtremalPoint &particlePoint,
        const pandora::CartesianVector &nuVertex, const HierarchyPfoVector &trackPfos, DLPrimaryNetworkParams &primaryNetworkParams) const;

    /**
     *  @brief  Shift and normalise the primary network parameters
     *
     *  @param  primaryNetworkParam the input primary network parameters
     */
    void NormaliseNetworkParams(DLPrimaryNetworkParams &primaryNetworkParams) const;

    /**
     *  @brief  Apply the primary track classification network
     *
     *  @param  edgeParamsUp the primary network parameters associated with the upstream endpoint
     *  @param  edgeParamsDown the primary network parameters associated with the downstream endpoint
     *
     *  @return the network score
     */
    float ClassifyTrack(const DLPrimaryNetworkParams &edgeParamsUp, const DLPrimaryNetworkParams &edgeParamsDown);

    /**
     *  @brief  Apply the primary track orientation edge network - determine whether an edge is
     *          signal (with correct orientation), signal (with the wrong orientation) or background
     *
     *  @param  edgeParams the primary network parameters associated with the edge to classify
     *  @param  otherEdgeParams the primary network parameters associated with the other edge
     *
     *  @return the float vector of (signal, wrong orientation, background) scores
     */
    pandora::FloatVector ClassifyTrackEdge(const DLPrimaryNetworkParams &edgeParams, const DLPrimaryNetworkParams &otherEdgeParams);

    /**
     *  @brief  Apply the primary shower classification network
     *
     *  @param  primaryNetworkParams the primary network parameters associated with the shower endpoint
     *
     *  @return the network score
     */
    float ClassifyShower(const DLPrimaryNetworkParams &primaryNetworkParams);

    // For model
    bool m_trainingMode;                                    ///< whether to run the tool in training mode
    std::string m_primaryTrackBranchModelName;              ///< the name of the primary track edge model file
    std::string m_primaryTrackClassifierModelName;          ///< the name of the primary track classification model file
    std::string m_primaryShowerClassifierModelName;         ///< the name of the primary shower classification model file
    LArDLHelper::TorchModel m_primaryTrackBranchModel;      ///< the primary track edge model
    LArDLHelper::TorchModel m_primaryTrackClassifierModel;  ///< the primary track classification model
    LArDLHelper::TorchModel m_primaryShowerClassifierModel; ///< the primary shower classification model
    // For tool
    float m_extrapolationStepSize; ///< the step size used to trace back a child particle's path
    // For normalisation
    bool m_normalise;                 ///< whether to normalise the network parameters
    NormalisationLimits m_normLimits; ///< struct of the normalisation limits
};

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_PRIMARY_HIERARCHY_TOOL_H
