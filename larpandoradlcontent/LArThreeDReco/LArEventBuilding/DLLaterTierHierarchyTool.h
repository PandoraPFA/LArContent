/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLLaterTierHierarchyTool.h
 *
 *  @brief  Header file for the DL later hierarchy tool
 *
 *  $Log: $
 */
#ifndef LAR_DL_LATER_TIER_HIERARCHY_TOOL_H
#define LAR_DL_LATER_TIER_HIERARCHY_TOOL_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLBaseHierarchyTool.h"

namespace lar_dl_content
{

/**
 *   @brief  DLLaterTierHierarchyTool to apply the later-tier hierarchy DL networks
 */
class DLLaterTierHierarchyTool : public DLBaseHierarchyTool
{
public:
    struct DLLaterTierNetworkParams  // all floats because they'll be normalised
    {
        float m_parentTrackScore = -999.f;          ///< the track/shower score of the parent pfo
        float m_childTrackScore = -999.f;           ///< the track/shower score of the child pfo
        float m_parentNSpacepoints = -999.f;        ///< the number of 3D hits in the parent pfo
        float m_childNSpacepoints = -999.f;         ///< the number of 3D hits in the child pfo
        float m_separation3D = -999.f;              ///< the smallest 3D distance between the parent and child 3D hits 
        float m_parentNuVertexSep = -999.f;         ///< the separation between the neutrino vertex and assumed parent start point 
        float m_childNuVertexSep = -999.f;          ///< the separation between the neutrino vertex and assumed child start point
        float m_parentEndRegionNHits = -999.f;      ///< the number of 3D hits 'close to' the parent POI
        float m_parentEndRegionNParticles = -999.f; ///< the number of different particles 'close to' the parent POI
        float m_parentEndRegionRToWall = -999.f;    ///< the smallest parent POI to detector boundary separation
        float m_vertexSeparation = -999.f;          ///< the separation between the parent and child POIs
        float m_doesChildConnect = -999.f;          ///< whether the backwards trace of the child's path connects to the parent
        pandora::CartesianVector m_connectionPoint =
	    pandora::CartesianVector(-999.f, -999.f, -999.f); ///< the connection point of the child on the parent  
        pandora::CartesianVector m_connectionDirection =
	    pandora::CartesianVector(-999.f, -999.f, -999.f); ///< the parent's direction at the connection point
        float m_overshootStartDCA = -999.f;       ///< the DCA of the child to the parent startpoint (set if the child doesn't connect) 
        float m_overshootStartL = -999.f;         ///< the extension distance of the child to the DCA point (set if the child doesn't connect)
        float m_overshootEndDCA = -999.f;         ///< the DCA of the child to the parent endpoint (set if the child doesn't connect)
        float m_overshootEndL = -999.f;           ///< the extension distance of the child to the DCA point (set if the child doesn't connect) 
        float m_childCPDCA = -999.f;              ///< the DCA of the child to the connection point (set if the child connects)
        float m_childCPExtrapDistance = -999.f;   ///< the extension distance of the child to the DCA point (set if the child connects)
        float m_childCPLRatio = -999.f;           ///< the ratio of the parent length at the connection point and the full parent length (set if the child connects)
        float m_parentCPNUpstreamHits = -999.f;   ///< the number of 3D parent hits upstream of the connection point (set if the child connects)
        float m_parentCPNDownstreamHits = -999.f; ///< the number of 3D parent hits downstream of the connection point (set if the child connects)
        float m_parentCPNHitRatio = -999.f;       ///< the ratio of the number of downstream hits and upstream 3D hits (set if the child connects)
        float m_parentCPEigenvalueRatio = -999.f; ///< the ratio of the second eigenvalues of the downstream and upstream 3D hit groups (set if the child connects)
        float m_parentCPOpeningAngle = -999.f;    ///< the opening angle between the first eigenvectors of the downstream and upstream 3D hit groups (set if the child connects)
        float m_parentIsPOIClosestToNu = -999.f;  ///< whether the parent POI is that closest to the neutrino vertex
        float m_childIsPOIClosestToNu = -999.f;   ///< whether the child POI is that closest to the neutrino vertex

        /**
         *  @brief  Return the vector of orientation independent later tier network parameters
         *
         *  @return The vector of orientation independent later tier network parameters
         */
        pandora::FloatVector GetCommonParamsForModel() const;

        /**
         *  @brief  Return the vector of orientation dependent later tier network parameters
         *
         *  @return The vector of orientation dependent later tier network parameters
         */      
        pandora::FloatVector GetOrientationParamsForModel() const;
    };

    pandora::StatusCode Run(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pNeutrinoPfo, 
        const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo, std::vector<DLLaterTierNetworkParams> &networkParamVector,
        float &laterTierScore);
  
    /**
     *  @brief  Default constructor
     */
    DLLaterTierHierarchyTool();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Whether the shower POI is the endpoint closest to the neutrino vertex
     *
     *  @param  parentHierarchyPfo the HierarchyPfo object of the parent pfo
     *  @param  childHierarchyPfo the HierarchyPfo object of the child pfo
     *
     *  @return whether the shower POI is in the upstream position
     */          
    bool IsShowerVertexUpstream(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo) const;

    /**
     *  @brief  Function to calculate the DLLaterTierNetworkParams
     *
     *  @param  pAlgorithm a pointer to the pandora algorithm
     *  @param  parentHierarchyPfo the HierarchyPfo object of the parent pfo
     *  @param  childHierarchyPfo the HierarchyPfo object of the child pfo
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  useUpstreamForParent whether the parent POI is that closest to the nu vertex
     *  @param  useUpstreamForChild whether the child POI is that closest to the nu vertex
     *  @param  laterTierNetworkParams the later tier network parameters to fill
     *
     *  @return a StatusCode to signify whether the network variables could be correctly calculated
     */
    pandora::StatusCode CalculateNetworkVariables(const pandora::Algorithm *const pAlgorithm, 
        const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo,
        const pandora::ParticleFlowObject *const pNeutrinoPfo, const bool useUpstreamForParent, const bool useUpstreamForChild, 
        DLLaterTierNetworkParams &laterTierNetworkParams) const;  

    /**
     *  @brief  Get the track/shower score DLLaterTierNetworkParams
     *          (m_parentTrackScore, m_childTrackScore)
     *
     *  @param  parentHierarchyPfo the HierarchyPfo object of the parent pfo
     *  @param  childHierarchyPfo the HierarchyPfo object of the child pfo     
     *
     *  @return the (parentTrackScore, childTrackScore) pair
     */
    std::pair<float, float> GetTrackScoreParams(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo) const;

    /**
     *  @brief  Get the number of 3D hit DLLaterTierNetworkParams
     *          (m_parentNSpacepoints, m_childNSpacepoints)
     *
     *  @param  parentHierarchyPfo the HierarchyPfo object of the parent pfo
     *  @param  childHierarchyPfo the HierarchyPfo object of the child pfo     
     *
     *  @return the (parentNSpacepoints, childNSpacepoints) pair
     */
    std::pair<float, float> GetNSpacepointsParams(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo) const;

    /**
     *  @brief  Get the m_separation3D DLLaterTierNetworkParam
     *
     *  @param  parentHierarchyPfo the HierarchyPfo object of the parent pfo
     *  @param  childHierarchyPfo the HierarchyPfo object of the child pfo     
     *
     *  @return the smallest 3D distance between the parent and child 3D hits
     */  
    float GetSeparation3D(const HierarchyPfo &parentHierarchyPfo, const HierarchyPfo &childHierarchyPfo) const;

   /**
    *  @brief  Set the orientation independent DLLaterTierNetworkParams
    *          (m_parentTrackScore, m_childTrackScore, m_parentNSpacepoints, m_childNSpacepoints, m_separation3D)
    *
    *  @param  trackScoreParams the track score (parent, child) pair
    *  @param  nSpacepointsParams the number of 3D hit (parent, child) pair
    *  @param  separation3D the smallest 3D distance between the parent and child 3D hits
    *  @param  laterTierNetworkParams the later tier network parameters to fill
    */  
    void SetCommonParams(const std::pair<float, float> &trackScoreParams, const std::pair<float, float> &nSpacepointsParams, 
        const float separation3D, DLLaterTierNetworkParams &laterTierNetworkParams) const;

    /**
     *  @brief  Set the vertex related DLLaterTierNetworkParams
     *          (m_parentNuVertexSep, m_childNuVertexSep, m_vertexSeparation)
     *
     *  @param  nuVertex the neutrino vertex position
     *  @param  parentStart the assumed parent pfo startpoint
     *  @param  parentEnd the assumed parent pfo endpoint
     *  @param  childStart the assumed child pfo startpoint
     *  @param  laterTierNetworkParams the later tier network parameters to fill
     */  
    void SetVertexParams(const pandora::CartesianVector &nuVertex, const pandora::CartesianVector &parentStart, 
        const pandora::CartesianVector &parentEnd, const pandora::CartesianVector &childStart, 
        DLLaterTierNetworkParams &laterTierNetworkParams) const;

    /**
     *  @brief  Set the parent endpoint region DLLaterTierNetworkParams
     *          (m_parentEndRegionNHits, m_parentEndRegionNParticles, m_parentEndRegionRToWall)
     *
     *  @param  pAlgorithm a pointer to the pandora algorithm 
     *  @param  pParentPfo a pointer to the parent pfo
     *  @param  parentEnd the assumed parent pfo endpoint
     *  @param  laterTierNetworkParams the later tier network parameters to fill
     */  
    void SetEndRegionParams(const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pParentPfo, 
        const pandora::CartesianVector &parentEnd, DLLaterTierNetworkParams &laterTierNetworkParams) const;

    /**
     *  @brief  Set the m_parentEndRegionRToWall DLLaterTierNetworkParam
     *
     *  @param  parentEnd the assumed parent pfo endpoint
     *  @param  laterTierNetworkParams the later tier network parameters to fill  
     */
    void SetEndRegionRToWall(const pandora::CartesianVector &parentEnd, DLLaterTierNetworkParams &laterTierNetworkParams) const;

    /**
     *  @brief  Set the child connection DLLaterTierNetworkParams
     *          (m_doesChildConnect, m_connectionPoint, m_connectionDirection, 
     *           m_childCPDCA, m_childCPExtrapDistance, m_childCPLRatio)
     *
     *  @param  parentHierarchyPfo the HierarchyPfo object of the parent pfo
     *  @param  parentStart the assumed parent pfo startpoint
     *  @param  childStart the assumed child pfo startpoint
     *  @param  childStartDirection the direction at the assumed child startpoint
     *  @param  laterTierNetworkParams the later tier network parameters to fill  
     */  
    void SetConnectionParams(const HierarchyPfo &parentHierarchyPfo,
        const pandora::CartesianVector &parentStart, const pandora::CartesianVector &childStart, 
        const pandora::CartesianVector &childStartDirection, DLLaterTierNetworkParams &laterTierNetworkParams) const;

    /**
     *  @brief  Project the (parentPos - childStart) vector onto the child direction axis 
     *          to obtain an extrapolation point
     *
     *  @param  parentPosition the parent position
     *  @param  childStart the assumed child pfo startpoint
     *  @param  childStartDirection the direction at the child startpoint
     *
     *  @return the (extrapolated point, whether the extrpolation was backwards wrt the child direction) pair
     */
    std::pair<pandora::CartesianVector, bool> ExtrapolateChildToParent(const pandora::CartesianVector &parentPosition, 
        const pandora::CartesianVector &childStart, const pandora::CartesianVector &childStartDirection) const;

    /**
     *  @brief  Return whether an input position connects to a line defined by two endpoints
     * 
     *  @param  boundary1 one endpoint
     *  @param  boundary2 the other endpoint
     *  @param  testPoint the input position
     *
     *  @return whether the input position connects to the line defined by the two input endpoints
     */  
    bool DoesConnect(const pandora::CartesianVector &boundary1, const pandora::CartesianVector &boundary2,
        const pandora::CartesianVector &testPoint) const;

    /**
     *  @brief  Set the overshoot DLLaterTierNetworkParams
     *          (m_overshootStartDCA, m_overshootStartL, m_overshootEndDCA, m_overshootEndL)
     *
     *  @param  parentStart the assumed parent pfo startpoint
     *  @param  parentStartDirection the direction at the parent startpoint
     *  @param  parentEnd the assumed parent pfo endpoint
     *  @param  parentEndDirection the direction at the parent endpoint
     *  @param  childStart the assumed child pfo startpoint
     *  @param  childStartDirection the direction at the child startpoint
     *  @param  laterTierNetworkParams the later tier network parameters to fill  
     */
    void SetOvershootParams(const pandora::CartesianVector &parentStart, const pandora::CartesianVector &parentStartDirection, 
        const pandora::CartesianVector &parentEnd, const pandora::CartesianVector &parentEndDirection, const pandora::CartesianVector &childStart, 
        const pandora::CartesianVector &childStartDirection, DLLaterTierNetworkParams &laterTierNetworkParams) const;

    /**
     *  @brief  Set the parent connection point DLLaterTierNetworkParams
     *          (m_parentCPNUpstreamHits, m_parentCPNDownstreamHits, m_parentCPNHitRatio, 
     *           m_parentCPEigenvalueRatio, m_parentCPOpeningAngle)
     *
     *  @param  parentHierarchyPfo the HierarchyPfo object of the parent pfo
     *  @param  laterTierNetworkParams the later tier network parameters to fill  
     */  
    void SetParentConnectionPointVars(const HierarchyPfo &parentHierarchyPfo, DLLaterTierNetworkParams &laterTierNetworkParams) const;

    /**
     *  @brief  Shift and normalise the later tier network parameters
     *
     *  @param  laterTierNetworkParam the input later tier network parameters
     */  
    void NormaliseNetworkParams(DLLaterTierNetworkParams &laterTierNetworkParams) const;

    /**
     *  @brief  Apply the track-track parent-child classification network
     *
     *  @param  edgeParamsUpUp the DLLaterTierNetworkParams for the parent upstream POI and child upstream POI
     *  @param  edgeParamsUpDown the DLLaterTierNetworkParams for the parent upstream POI and child downstream POI
     *  @param  edgeParamsDownUp the DLLaterTierNetworkParams for the parent downstream POI and child upstream POI
     *  @param  edgeParamsDownDown the DLLaterTierNetworkParams for the parent downstream POI and child downstream POI
     *
     *  @return the network score (how likely to be signal)
     */  
    float ClassifyTrackTrack(const DLLaterTierNetworkParams &edgeParamsUpUp, const DLLaterTierNetworkParams &edgeParamsUpDown, 
        const DLLaterTierNetworkParams &edgeParamsDownUp, const DLLaterTierNetworkParams &edgeParamsDownDown);

    /**
     *  @brief  Apply the track-track orientation edge network - determine whether an edge is
     *          signal (with correct orientation), signal (with the wrong orientation) or background
     *
     *  @param  edgeParams the DLLaterTierNetworkParams of the edge to classify
     *  @param  otherEdgeParams1 the DLLaterTierNetworkParams of another edge
     *  @param  otherEdgeParams2 the DLLaterTierNetworkParams of another another edge
     *  @param  otherEdgeParams3 the DLLaterTierNetworkParams of another another another edge
     *
     *  @return the float vector of (signal, wrong orientation, background) scores
     */
    pandora::FloatVector ClassifyTrackTrackEdge(const DLLaterTierNetworkParams &edgeParams, 
        const DLLaterTierNetworkParams &otherEdgeParams1, const DLLaterTierNetworkParams &otherEdgeParams2, 
        const DLLaterTierNetworkParams &otherEdgeParams3);

    /**
     *  @brief  Apply the track-shower parent-child classification network
     *
     *  @param  edgeParamsUp the DLLaterTierNetworkParams for the parent upstream POI
     *  @param  edgeParamsDown the DLLaterTierNetworkParams for the parent upstream POI
     *
     *  @return the network score (how likely to be signal)
     */    
    float ClassifyTrackShower(const DLLaterTierNetworkParams &edgeParamsUp, const DLLaterTierNetworkParams &edgeParamsDown);

    /**
     *  @brief  Apply the track-shower orientation edge network - determine whether an edge is
     *          signal (with correct orientation), signal (with the wrong orientation) or background
     *
     *  @param  edgeParams the DLLaterTierNetworkParams of the edge to classify
     *  @param  otherEdgeParams the DLLaterTierNetworkParams of another edge
     *
     *  @return the float vector of (signal, wrong orientation, background) scores
     */  
    pandora::FloatVector ClassifyTrackShowerEdge(const DLLaterTierNetworkParams &edgeParams, 
        const DLLaterTierNetworkParams &otherEdgeParams);

    // For model
    std::string m_trackTrackBranchModelName;              ///< the name of the track-track edge model file
    std::string m_trackTrackClassifierModelName;          ///< the name of the track-track classification model file
    std::string m_trackShowerBranchModelName;             ///< the name of the track-shower edge model file
    std::string m_trackShowerClassifierModelName;         ///< the name of the track-shower classification model file
    LArDLHelper::TorchModel m_trackTrackBranchModel;      ///< the track-track edge model
    LArDLHelper::TorchModel m_trackTrackClassifierModel;  ///< the track-track classification model
    LArDLHelper::TorchModel m_trackShowerBranchModel;     ///< the track-shower edge model
    LArDLHelper::TorchModel m_trackShowerClassifierModel; ///< the track-shower classification model
    // For tool
    bool m_trainingMode;        ///< whether to run the tool in training mode
    float m_trajectoryStepSize; ///< the size of the steps taken to trace the parent trajectory
    float m_connectionBuffer;   ///< the 'is child connected?' threshold
    float m_searchRegion;       ///< the dimensions of the box used to obtain the up/downstream 3D hit groups
    // For normalisation
    bool m_normalise;                     ///< whether to normalise the network parameters
    float m_trackScoreMin;                ///< the m_parent(child)TrackScore minimum range boundary
    float m_trackScoreMax;                ///< the m_parent(child)TrackScore maximum range boundary
    float m_nSpacepointsMin;              ///< the m_parent(child)NSpacepoints minimum range boundary
    float m_nSpacepointsMax;              ///< the m_parent(child)NSpacepoints maximum range boundary
    float m_separation3DMin;              ///< the m_separation3D minimum range boundary
    float m_separation3DMax;              ///< the m_separation3D maximum range boundary
    float m_nuVertexSepMin;               ///< the m_parent(child)NuVertexSep minimum range boundary
    float m_nuVertexSepMax;               ///< the m_parent(child)NuVertexSep maximum range boundary
    float m_parentEndRegionNHitsMin;      ///< the m_parentEndRegionNHits minimum range boundary
    float m_parentEndRegionNHitsMax;      ///< the m_parentEndRegionNHits maximum range boundary
    float m_parentEndRegionNParticlesMin; ///< the m_parentEndRegionNParticles minimum range boundary
    float m_parentEndRegionNParticlesMax; ///< the m_parentEndRegionNParticles maximum range boundary
    float m_parentEndRegionRToWallMin;    ///< the m_parentEndRegionRToWall minimum range boundary
    float m_parentEndRegionRToWallMax;    ///< the m_parentEndRegionRToWall maximum range boundary
    float m_vertexSepMin;                 ///< the m_vertexSeparation minimum range boundary
    float m_vertexSepMax;                 ///< the m_vertexSeparation maximum range boundary
    float m_doesChildConnectMin;          ///< the m_doesChildConnect minimum range boundary
    float m_doesChildConnectMax;          ///< the m_doesChildConnect maximum range boundary
    float m_overshootDCAMin;              ///< the m_overshootStart(End)DCA minimum range boundary
    float m_overshootDCAMax;              ///< the m_overshootStart(End)DCA maximum range boundary
    float m_overshootLMin;                ///< the m_overshootStart(End)L minimum range boundary
    float m_overshootLMax;                ///< the m_overshootStart(End)L maximum range boundary
    float m_childCPDCAMin;                ///< the m_childCPDCA minimum range boundary
    float m_childCPDCAMax;                ///< the m_childCPDCA maximum range boundary
    float m_childCPExtrapDistanceMin;     ///< the m_childCPExtrapDistance minimum range boundary
    float m_childCPExtrapDistanceMax;     ///< the m_childCPExtrapDistance maximum range boundary
    float m_childCPLRatioMin;             ///< the m_childCPLRatio minimum range boundary
    float m_childCPLRatioMax;             ///< the m_childCPLRatio maximum range boundary
    float m_parentCPNHitsMin;             ///< the m_parentCPNUp(Down)streamHits minimum range boundary
    float m_parentCPNHitsMax;             ///< the m_parentCPNUp(Down)streamHits maximum range boundary
    float m_parentCPNHitRatioMin;         ///< the m_parentCPNHitRatio minimum range boundary
    float m_parentCPNHitRatioMax;         ///< the m_parentCPNHitRatio maximum range boundary
    float m_parentCPEigenvalueRatioMin;   ///< the m_parentCPEigenvalueRatio minimum range boundary
    float m_parentCPEigenvalueRatioMax;   ///< the m_parentCPEigenvalueRatio maximum range boundary
    float m_parentCPOpeningAngleMin;      ///< the m_parentCPOpeningAngle minimum range boundary
    float m_parentCPOpeningAngleMax;      ///< the m_parentCPOpeningAngle maximum range boundary
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::FloatVector DLLaterTierHierarchyTool::DLLaterTierNetworkParams::GetCommonParamsForModel() const
{
    return {m_parentTrackScore, m_childTrackScore, m_parentNSpacepoints, m_childNSpacepoints, m_separation3D};
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::FloatVector DLLaterTierHierarchyTool::DLLaterTierNetworkParams::GetOrientationParamsForModel() const
{
    return {m_parentNuVertexSep, m_childNuVertexSep, m_parentEndRegionNHits, m_parentEndRegionNParticles, m_parentEndRegionRToWall, 
            m_vertexSeparation, m_doesChildConnect, m_overshootStartDCA, m_overshootStartL, m_overshootEndDCA, m_overshootEndL, 
            m_childCPDCA, m_childCPExtrapDistance, m_childCPLRatio, m_parentCPNUpstreamHits, m_parentCPNDownstreamHits, 
            m_parentCPNHitRatio, m_parentCPEigenvalueRatio, m_parentCPOpeningAngle, m_parentIsPOIClosestToNu, m_childIsPOIClosestToNu};
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_LATER_TIER_HIERARCHY_TOOL_H
