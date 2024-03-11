/**
 *  @file   larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h
 *
 *  @brief  Header file for the vertex selection base algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_SELECTION_BASE_ALGORITHM_H
#define LAR_VERTEX_SELECTION_BASE_ALGORITHM_H 1

#include "Objects/Vertex.h"
#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArObjects/LArSupportVectorMachine.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

template <typename, unsigned int>
class KDTreeLinkerAlgo;
template <typename, unsigned int>
class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  VertexSelectionBaseAlgorithm class
 */
class VertexSelectionBaseAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexSelectionBaseAlgorithm();

    /**
     *  @brief  VertexScore class
     */
    class VertexScore
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pVertex the address of the vertex
         *  @param  score the score
         */
        VertexScore(const pandora::Vertex *const pVertex, const float score);

        /**
         *  @brief  Get the address of the vertex
         *
         *  @return the address of the vertex
         */
        const pandora::Vertex *GetVertex() const;

        /**
         *  @brief  Get the score
         *
         *  @return the score
         */
        float GetScore() const;

        /**
         *  @brief  operator<
         *
         *  @param  rhs the value for comparison
         *
         *  @return boolean
         */
        bool operator<(const VertexScore &rhs) const;

    private:
        const pandora::Vertex *m_pVertex; ///< The address of the vertex
        float m_score;                    ///< The score
    };

    typedef std::vector<VertexScore> VertexScoreList;

    /**
     *  @brief Beam constants class
     */
    class BeamConstants
    {
    public:
        /**
         *  @brief  Get the min z coordinate
         *
         *  @return the min z coordinate
         */
        float GetMinZCoordinate() const;

        /**
         *  @brief  Get the decay constant
         *
         *  @return the decay constant
         */
        float GetDecayConstant() const;

        /**
         *  @brief  Set the beam constants
         *
         *  @param  minZCoordinate the min z coordinate
         *  @param  decayConstant the decay constant
         */
        void SetConstants(const float minZCoordinate, const float decayConstant);

    private:
        pandora::InputFloat m_minZCoordinate; ///< The min z coordinate
        pandora::InputFloat m_decayConstant;  ///< The decay constant
    };

    /**
     *  @brief Sliding fit data class.
     */
    class SlidingFitData
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  pCluster pointer to the cluster
         *  @param  slidingFitWindow the sliding fit window
         *  @param  slidingFitPitch the sliding fit pitch
         */
        SlidingFitData(const pandora::Cluster *const pCluster, const int slidingFitWindow, const float slidingFitPitch);

        /**
         *  @brief  Get the min layer direction
         *
         *  @return the min layer direction
         */
        const pandora::CartesianVector &GetMinLayerDirection() const;

        /**
         *  @brief  Get the max layer direction
         *
         *  @return the max layer direction
         */
        const pandora::CartesianVector &GetMaxLayerDirection() const;

        /**
         *  @brief  Get the min layer position
         *
         *  @return the min layer position
         */
        const pandora::CartesianVector &GetMinLayerPosition() const;

        /**
         *  @brief  Get the max layer position
         *
         *  @return the max layer position
         */
        const pandora::CartesianVector &GetMaxLayerPosition() const;

        /**
         *  @brief  Get a pointer to the corresponding cluster
         *
         *  @return pointer to the corresponding cluster
         */
        const pandora::Cluster *GetCluster() const;

    private:
        pandora::CartesianVector m_minLayerDirection; ///< The direction of the fit at the min layer
        pandora::CartesianVector m_maxLayerDirection; ///< The direction of the fit at the min layer
        pandora::CartesianVector m_minLayerPosition;  ///< The position of the fit at the max layer
        pandora::CartesianVector m_maxLayerPosition;  ///< The position of the fit at the max layer
        const pandora::Cluster *m_pCluster;           ///< Pointer to the corresponding cluster
    };

    typedef std::vector<SlidingFitData> SlidingFitDataList;

    /**
     *  @brief Shower cluster class
     */
    class ShowerCluster
    {
    public:
        /**
         *  @brief  Constructor
         *
         *  @param  clusterList the list of clusters
         *  @param  slidingFitWindow the sliding fit window
         *  @param  slidingFitPitch the sliding fit pitch
         */
        ShowerCluster(const pandora::ClusterList &clusterList, const int slidingFitWindow, const float slidingFitPitch);

        /**
         *  @brief  Get the cluster list
         *
         *  @return the cluster list
         */
        const pandora::ClusterList &GetClusters() const;

        /**
         *  @brief  Get the 2D sliding linear fit
         *
         *  @return the fit
         */
        const TwoDSlidingFitResult &GetFit() const;

        /**
         *  @brief  Get the coordinate vector for a cluster list
         *
         *  @param  clusterList the cluster list
         *
         *  @return the coordinate vector
         */
        pandora::CartesianPointVector GetClusterListCoordinateVector(const pandora::ClusterList &clusterList) const;

    private:
        pandora::ClusterList m_clusterList;               ///< The list of clusters
        pandora::CartesianPointVector m_coordinateVector; ///< The coordinate vector
        TwoDSlidingFitResult m_twoDSlidingFitResult;      ///< The fit to the hits of the cluster list
    };

    typedef std::vector<ShowerCluster> ShowerClusterList;

    typedef KDTreeNodeInfoT<const pandora::CaloHit *, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;
    typedef KDTreeLinkerAlgo<const pandora::CaloHit *, 2> HitKDTree2D;

    typedef std::map<pandora::HitType, const pandora::ClusterList &> ClusterListMap;    ///< Map array of cluster lists for passing to tools
    typedef std::map<pandora::HitType, const SlidingFitDataList> SlidingFitDataListMap; ///< Map of sliding fit data lists for passing to tools
    typedef std::map<pandora::HitType, const ShowerClusterList> ShowerClusterListMap; ///< Map of shower cluster lists for passing to tools
    typedef std::map<pandora::HitType, const std::reference_wrapper<HitKDTree2D>> KDTreeMap; ///< Map array of hit kd trees for passing to tools

    typedef MvaFeatureTool<const VertexSelectionBaseAlgorithm *const, const pandora::Vertex *const, const SlidingFitDataListMap &,
        const ClusterListMap &, const KDTreeMap &, const ShowerClusterListMap &, const float, float &>
        VertexFeatureTool; ///< The base type for the vertex feature tools

protected:
    /**
     *  @brief  Filter the input list of vertices to obtain a reduced number of vertex candidates
     *
     *  @param  pInputVertexList the address of the input vertex list
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     *  @param  filteredVertices to receive the filtered vertex list
     */
    virtual void FilterVertexList(const pandora::VertexList *const pInputVertexList, HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV,
        HitKDTree2D &kdTreeW, pandora::VertexVector &filteredVertices) const;

    /**
     *  @brief  Get the beam score constants for a provided list of candidate vertices
     *
     *  @param  vertexVector the vertex vector
     *  @param  beamConstants to receive the beam constants
     */
    virtual void GetBeamConstants(const pandora::VertexVector &vertexVector, BeamConstants &beamConstants) const;

    /**
     *  @brief  Get the vertex score list for a provided list of candidate vertices
     *
     *  @param  vertexVector the vertex vector
     *  @param  beamConstants the beam constants
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     *  @param  vertexScoreList to receive the vertex score list
     */
    virtual void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const = 0;

    /**
     *  @brief  Get the cluster lists
     *
     *  @param  inputClusterListNames the input cluster list names
     *  @param  clusterListU the U-view cluster list to populate
     *  @param  clusterListV the V-view cluster list to populate
     *  @param  clusterListW the W-view cluster list to populate
     */
    void GetClusterLists(const pandora::StringVector &inputClusterListNames, pandora::ClusterList &clusterListU,
        pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW) const;

    /**
     *  @brief  Calculate the cluster sliding fits
     *
     *  @param  inputClusterList the input cluster list
     *  @param  minClusterCaloHits the minimum number of cluster calo hits
     *  @param  slidingFitWindow the sliding fit window
     *  @param  slidingFitDataList the list of sliding fits to fill
     */
    void CalculateClusterSlidingFits(const pandora::ClusterList &inputClusterList, const unsigned int minClusterCaloHits,
        const unsigned int slidingFitWindow, SlidingFitDataList &slidingFitDataList) const;

    /**
     *  @brief  Get the beam deweighting score for a vertex
     *
     *  @param  beamConstants the beam constants
     *  @param  pVertex address of the vertex
     *
     *  @return the score
     */
    float GetBeamDeweightingScore(const BeamConstants &beamConstants, const pandora::Vertex *const pVertex) const;

    /**
     *  @brief  Whether algorithm is running in beam mode, assuming neutrinos travel in positive z-direction
     *
     *  @return boolean
     */
    bool IsBeamModeOn() const;

    /**
     *  @brief  Calculate the energy of a vertex candidate by summing values from all three planes
     *
     *  @param  pVertex the address of the vertex
     *  @param  kdTreeMap the map of 2D hit kd trees
     *
     *  @return the summed vertex energy
     */
    float GetVertexEnergy(const pandora::Vertex *const pVertex, const KDTreeMap &kdTreeMap) const;

    /**
     *  @brief  Finds the energy of the nearest hit to the vertex candidate in this view
     *
     *  @param  pVertex the address of the vertex
     *  @param  hitType the relevant hit type
     *  @param  kdTree the kd tree of 2D hits
     *
     *  @return the energy of the nearest hit
     */
    float VertexHitEnergy(const pandora::Vertex *const pVertex, const pandora::HitType hitType, HitKDTree2D &kdTree) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Initialize kd trees with details of hits in algorithm-configured cluster lists
     *
     *  @param  kdTreeU the kd tree for u hits
     *  @param  kdTreeV the kd tree for v hits
     *  @param  kdTreeW the kd tree for w hits
     */
    void InitializeKDTrees(HitKDTree2D &kdTreeU, HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW) const;

    /**
     *  @brief  Whether the vertex lies on a hit in the specified view
     *
     *  @param  pVertex the address of the vertex
     *  @param  hitType the relevant hit type
     *  @param  kdTree the relevant kd tree
     *
     *  @return boolean
     */
    bool IsVertexOnHit(const pandora::Vertex *const pVertex, const pandora::HitType hitType, HitKDTree2D &kdTree) const;

    /**
     *  @brief  Whether the vertex lies in a registered gap
     *
     *  @param  pVertex the address of the vertex
     *  @param  hitType the relevant hit type
     *
     *  @return boolean
     */
    bool IsVertexInGap(const pandora::Vertex *const pVertex, const pandora::HitType hitType) const;

    /**
     *  @brief  From the top-scoring candidate vertices, select a subset for further investigation
     *
     *  @param  vertexScoreList the vertex score list
     *  @param  selectedVertexList to receive the selected vertex list
     */
    void SelectTopScoreVertices(VertexScoreList &vertexScoreList, pandora::VertexList &selectedVertexList) const;

    /**
     *  @brief  Whether to accept a candidate vertex, based on its spatial position in relation to other selected candidates
     *
     *  @param  pVertex the address of the vertex
     *  @param  selectedVertexList the selected vertex list
     *
     *  @return boolean
     */
    bool AcceptVertexLocation(const pandora::Vertex *const pVertex, const pandora::VertexList &selectedVertexList) const;

    /**
     *  @brief  Sort vertices by increasing z position
     *
     *  @param  pLhs address of the lhs vertex
     *  @param  pRhs address of the rhs vertex
     *
     *  @return whether lhs should precedes rhs
     */
    static bool SortByVertexZPosition(const pandora::Vertex *const pLhs, const pandora::Vertex *const pRhs);

private:
    pandora::StringVector m_inputCaloHitListNames; ///< The list of calo hit list names
    std::string m_outputVertexListName;            ///< The name under which to save the output vertex list

    bool m_replaceCurrentVertexList; ///< Whether to replace the current vertex list with the output list

    bool m_beamMode;              ///< Whether to run in beam mode, assuming neutrinos travel in positive z-direction
    float m_nDecayLengthsInZSpan; ///< The number of score decay lengths to use over the course of the vertex z-span

    bool m_selectSingleVertex;            ///< Whether to make a final decision and select just one vertex candidate
    unsigned int m_maxTopScoreSelections; ///< Max number of top-scoring vertex candidate to select for output

    float m_maxOnHitDisplacement; ///< Max hit-vertex displacement for declaring vertex to lie on a hit in each view

    float m_minCandidateDisplacement;  ///< Ignore other top-scoring candidates located in close proximity to original
    float m_minCandidateScoreFraction; ///< Ignore other top-scoring candidates with score less than a fraction of original

    bool m_useDetectorGaps; ///< Whether to account for registered detector gaps in vertex selection
    float m_gapTolerance;   ///< The tolerance to use when querying whether a sampling point is in a gap, units cm

    bool m_isEmptyViewAcceptable; ///< Whether views entirely empty of hits are classed as 'acceptable' for candidate filtration
    unsigned int m_minVertexAcceptableViews; ///< The minimum number of views in which a candidate must sit on/near a hit or in a gap (or view can be empty)
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexSelectionBaseAlgorithm::GetBeamDeweightingScore(const BeamConstants &beamConstants, const pandora::Vertex *const pVertex) const
{
    const float vertexMinZ(std::max(pVertex->GetPosition().GetZ(), beamConstants.GetMinZCoordinate()));
    return (beamConstants.GetMinZCoordinate() - vertexMinZ) * beamConstants.GetDecayConstant();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool VertexSelectionBaseAlgorithm::IsBeamModeOn() const
{
    return m_beamMode;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline VertexSelectionBaseAlgorithm::VertexScore::VertexScore(const pandora::Vertex *const pVertex, const float score) :
    m_pVertex(pVertex),
    m_score(score)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Vertex *VertexSelectionBaseAlgorithm::VertexScore::GetVertex() const
{
    return m_pVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexSelectionBaseAlgorithm::VertexScore::GetScore() const
{
    return m_score;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool VertexSelectionBaseAlgorithm::VertexScore::operator<(const VertexScore &rhs) const
{
    return (this->GetScore() > rhs.GetScore());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexSelectionBaseAlgorithm::BeamConstants::GetMinZCoordinate() const
{
    return m_minZCoordinate.Get();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexSelectionBaseAlgorithm::BeamConstants::GetDecayConstant() const
{
    return m_decayConstant.Get();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void VertexSelectionBaseAlgorithm::BeamConstants::SetConstants(const float minZCoordinate, const float decayConstant)
{
    m_minZCoordinate = minZCoordinate;
    m_decayConstant = decayConstant;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &VertexSelectionBaseAlgorithm::SlidingFitData::GetMinLayerDirection() const
{
    return m_minLayerDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &VertexSelectionBaseAlgorithm::SlidingFitData::GetMaxLayerDirection() const
{
    return m_maxLayerDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &VertexSelectionBaseAlgorithm::SlidingFitData::GetMinLayerPosition() const
{
    return m_minLayerPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &VertexSelectionBaseAlgorithm::SlidingFitData::GetMaxLayerPosition() const
{
    return m_maxLayerPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *VertexSelectionBaseAlgorithm::SlidingFitData::GetCluster() const
{
    return m_pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ClusterList &VertexSelectionBaseAlgorithm::ShowerCluster::GetClusters() const
{
    return m_clusterList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TwoDSlidingFitResult &VertexSelectionBaseAlgorithm::ShowerCluster::GetFit() const
{
    return m_twoDSlidingFitResult;
}

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_SELECTION_BASE_ALGORITHM_H
