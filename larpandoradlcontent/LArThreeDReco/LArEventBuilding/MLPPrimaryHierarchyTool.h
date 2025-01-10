/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPPrimaryHierarchyTool.h
 *
 *  @brief  Header file for the MLP primary hierarchy tool
 *
 *  $Log: $
 */
#ifndef LAR_MLP_PRIMARY_HIERARCHY_TOOL_H
#define LAR_MLP_PRIMARY_HIERARCHY_TOOL_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPBaseHierarchyTool.h"

namespace lar_dl_content
{
/**
 *   @brief  MLPPrimaryHierarchyTool to apply the primary hierarchy MLP networks
 */
class MLPPrimaryHierarchyTool : public MLPBaseHierarchyTool
{
public:
    struct MLPPrimaryNetworkParams  // all floats because they'll be normalised
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
         *  @brief  Print the value of the MLPPrimaryNetworkParams 
         */
        void Print() const;

        /**
         *  @brief  Return the vector of orientation independent primary network parameters
         *
         *  @return The vector of orientation independent primary network parameters
         */
        pandora::FloatVector GetCommonParamsForModel() const;

        /**
         *  @brief  Return the vector of orientation dependent primary network parameters
         *
         *  @return The vector of orientation dependent primary network parameters
         */
        pandora::FloatVector GetOrientationParamsForModel() const;
    };

    pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pNeutrinoPfo, 
        const HierarchyPfoMap &trackPfos, const HierarchyPfo &hierarchyPfo, std::vector<MLPPrimaryNetworkParams> &networkParamVector, 
        float &primaryScore);

    /**
     *  @brief  Default constructor
     */
    MLPPrimaryHierarchyTool();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Function to calculate the MLPPrimaryNetworkParams
     *
     *  @param  pAlgorithm a pointer to the pandora algorithm
     *  @param  hierarchyPfo the input hierarchy pfo object
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  trackPfos the <pfo -> HierarchyPfo> map for the track-like particles
     *  @param  useUpstream whether the POI is the endpoint closest to the nu vertex
     *  @param  primaryNetworkParams the primary network parameters to fill
     *
     *  @return a StatusCode to signify whether the network variables could be correctly calculated
     */
    pandora::StatusCode CalculateNetworkVariables(const pandora::Algorithm *const pAlgorithm, const HierarchyPfo &hierarchyPfo, 
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoMap &trackPfos, const bool useUpstream, 
        MLPPrimaryNetworkParams &primaryNetworkParams) const;
  
    /**
     *  @brief  Set the vertex region MLPPrimaryNetworkParams params
     *          (m_vertexRegionNHits, m_vertexRegionNParticles)
     *
     *  @param  pAlgorithm a pointer to the pandora algorithm
     *  @param  pPfo a pointer to the input pfo
     *  @param  particleVertex the position of the pfo vertex
     *  @param  primaryNetworkParams the primary network parameters
     */
    void SetVertexRegionParams(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pPfo, 
        const pandora::CartesianVector &particleVertex, MLPPrimaryNetworkParams &primaryNetworkParams) const;

    /**
     *  @brief  Set the connection region MLPPrimaryNetworkParams params
     *          (m_dca, m_connectionExtrapDistance)
     *
     *  @param  particleVertex the position of the pfo vertex
     *  @param  particleDirection the direction at the pfo vertex
     *  @param  nuVertex the neutrino vertex
     *  @param  primaryNetworkParams the primary network parameters
     */
    void SetConnectionParams(const pandora::CartesianVector &particleVertex, const pandora::CartesianVector &particleDirection, 
        const pandora::CartesianVector &nuVertex, MLPPrimaryNetworkParams &primaryNetworkParams) const;

    /**
     *  @brief  Set the event context MLPPrimaryNetworkParams params
     *          (m_parentConnectionDistance, m_childConnectionDistance)
     *
     *  @param  pPfo a pointer to the input pfo
     *  @param  particleVertex the position of the pfo vertex
     *  @param  particleDirection the direction at the pfo vertex
     *  @param  nuVertex the neutrino vertex
     *  @param  trackPfos the <HierarchyPfo -> pfo> map for the track-like particles
     *  @param  primaryNetworkParams the primary network parameters
     */
    void SetContextParams(const pandora::ParticleFlowObject *const pPfo, const pandora::CartesianVector &particleVertex, 
        const pandora::CartesianVector &particleDirection, const pandora::CartesianVector &nuVertex, 
        const HierarchyPfoMap &trackPfos, MLPPrimaryNetworkParams &primaryNetworkParams) const;

    /**
     *  @brief  Calculate the variables describing the extension of a child particle to a given parent 
     *          (m_parentConnectionDistance, m_childConnectionDistance)
     *
     *  @param  parentVertex the position of the vertex of the parent pfo
     *  @param  parentDirection the direction at the parent vertex
     *  @param  childVertex the position of the vertex of the child pfo
     *  @param  childDirection the direction at the child vertex
     *  @param  parentConnectionDistance the DCA to the parent vertex
     *  @param  childConnectionDistance the extension of the child to the DCA point
     */
    void CalculateConnectionDistances(const pandora::CartesianVector &parentVertex, const pandora::CartesianVector &parentDirection, 
        const pandora::CartesianVector &childVertex, const pandora::CartesianVector &childDirection, 
        float &parentConnectionDistance, float &childConnectionDistance) const;

    /**
     *  @brief  Shift and normalise the primary network parameters
     *
     *  @param  primaryNetworkParam the input primary network parameters
     */
    void NormaliseNetworkParams(MLPPrimaryNetworkParams &primaryNetworkParams) const;

    /**
     *  @brief  Apply the primary track classification network
     *
     *  @param  edgeParamsUp the primary network parameters associated with the upstream endpoint
     *  @param  edgeParamsDown the primary network parameters associated with the downstream endpoint
     *
     *  @return the network score
     */
    float ClassifyTrack(const MLPPrimaryNetworkParams &edgeParamsUp, const MLPPrimaryNetworkParams &edgeParamsDown);

    /**
     *  @brief  Apply the primary track orientation edge network - determine whether an edge is
     *          signal (with correct orientation), signal (with the wrong orientation) or background
     *
     *  @param  edgeParams the primary network parameters associated with the edge to classify
     *  @param  otherEdgeParams the primary network parameters associated with the other edge
     *
     *  @return the float vector of (signal, wrong orientation, background) scores
     */
    pandora::FloatVector ClassifyTrackEdge(const MLPPrimaryNetworkParams &edgeParams, 
        const MLPPrimaryNetworkParams &otherEdgeParams);

    /**
     *  @brief  Apply the primary shower classification network
     *
     *  @param  primaryNetworkParams the primary network parameters associated with the shower endpoint
     *
     *  @return the network score
     */
    float ClassifyShower(const MLPPrimaryNetworkParams &primaryNetworkParams);

    // For model
    std::string m_primaryTrackBranchModelName;              ///< the name of the primary track edge model file
    std::string m_primaryTrackClassifierModelName;          ///< the name of the primary track classification model file 
    std::string m_primaryShowerClassifierModelName;         ///< the name of the primary shower classification model file 
    LArDLHelper::TorchModel m_primaryTrackBranchModel;      ///< the primary track edge model
    LArDLHelper::TorchModel m_primaryTrackClassifierModel;  ///< the primary track classification model
    LArDLHelper::TorchModel m_primaryShowerClassifierModel; ///< the primary shower classification model
    // For tool
    float m_extrapolationStepSize;       ///< the step size used to trace back a child particle's path
    // For normalisation
    bool m_normalise;                    ///< whether to normalise the network parameters
    float m_nSpacepointsMin;             ///< the m_nSpacepoints minimum range boundary 
    float m_nSpacepointsMax;             ///< the m_nSpacepoints maximum range boundary
    float m_nuSeparationMin;             ///< the m_nuSeparation minimum range boundary
    float m_nuSeparationMax;             ///< the m_nuSeparation maximum range boundary
    float m_vertexRegionNHitsMin;        ///< the m_vertexRegionNHits minimum range boundary
    float m_vertexRegionNHitsMax;        ///< the m_vertexRegionNHits maximum range boundary
    float m_vertexRegionNParticlesMin;   ///< the m_vertexRegionNParticles minimum range boundary
    float m_vertexRegionNParticlesMax;   ///< the m_vertexRegionNParticles maximum range boundary
    float m_dcaMin;                      ///< the m_dca minimum range boundary
    float m_dcaMax;                      ///< the m_dca maximum range boundary
    float m_connectionExtrapDistanceMin; ///< the m_connectionExtrapDistance minimum range boundary
    float m_connectionExtrapDistanceMax; ///< the m_connectionExtrapDistance maximum range boundary
    float m_parentConnectionDistanceMin; ///< the m_parentConnectionDistance minimum range boundary
    float m_parentConnectionDistanceMax; ///< the m_parentConnectionDistance maximum range boundary
    float m_childConnectionDistanceMin;  ///< the m_childConnectionDistance minimum range boundary
    float m_childConnectionDistanceMax;  ///< the m_childConnectionDistance maximum range boundary
};

//------------------------------------------------------------------------------------------------------------------------------------------

 inline void MLPPrimaryHierarchyTool::MLPPrimaryNetworkParams::Print() const
{
    std::cout << "IsPOIClosestToNu: " << m_isPOIClosestToNu << std::endl;
    std::cout << "NSpacepoints: " << m_nSpacepoints << std::endl;
    std::cout << "NuVertexSep: " << m_nuSeparation << std::endl;
    std::cout << "StartRegionNHits: " << m_vertexRegionNHits << std::endl;
    std::cout << "StartRegionNParticles: " << m_vertexRegionNParticles << std::endl;
    std::cout << "DCA: " << m_dca << std::endl;
    std::cout << "ExtrapDistance: " << m_connectionExtrapDistance << std::endl;
    std::cout << "parentConnectionDistance: " << m_parentConnectionDistance << std::endl;
    std::cout << "childConnectionDistance: " << m_childConnectionDistance << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::FloatVector MLPPrimaryHierarchyTool::MLPPrimaryNetworkParams::GetCommonParamsForModel() const
{
    return { m_nSpacepoints };
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::FloatVector MLPPrimaryHierarchyTool::MLPPrimaryNetworkParams::GetOrientationParamsForModel() const
{
    return { m_nuSeparation, m_vertexRegionNHits, m_vertexRegionNParticles, m_dca, m_connectionExtrapDistance, 
            m_isPOIClosestToNu, m_parentConnectionDistance, m_childConnectionDistance };
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_PRIMARY_HIERARCHY_TOOL_H
