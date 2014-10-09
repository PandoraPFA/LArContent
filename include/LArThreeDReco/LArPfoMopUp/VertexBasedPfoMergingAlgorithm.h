/**
 *  @file   LArContent/include/LArThreeDReco/LArPfoMopUp/VertexBasedPfoMergingAlgorithm.h
 * 
 *  @brief  Header file for the vertex based pfo merging algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_BASED_PFO_MERGING_ALGORITHM_H
#define LAR_VERTEX_BASED_PFO_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  VertexBasedPfoMergingAlgorithm::Algorithm class
 */
class VertexBasedPfoMergingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    VertexBasedPfoMergingAlgorithm();

private:
    /**
     *  @brief  ClusterAssociation class
     */
    class ClusterAssociation
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pVertexCluster the address of the vertex cluster
         *  @param  pDaughterCluster the address of the daughter cluster
         *  @param  boundedFraction the fraction of daughter hits bounded by the cone defined by the vertex cluster
         */
        ClusterAssociation(pandora::Cluster *pVertexCluster, pandora::Cluster *pDaughterCluster, const float boundedFraction);

        /**
         *  @brief  Get the address of the vertex cluster
         *
         *  @return the address of the vertex cluster
         */
        pandora::Cluster *GetVertexCluster() const;

        /**
         *  @brief  Get the address of the daughter cluster
         *
         *  @return the address of the daughter cluster
         */
        pandora::Cluster *GetDaughterCluster() const;

        /**
         *  @brief  Get the fraction of daughter hits bounded by the cone defined by the vertex cluster
         *
         *  @return the bounded fraction
         */
        float GetBoundedFraction() const;

    private:
        pandora::Cluster   *m_pVertexCluster;           ///< The address of the vertex cluster
        pandora::Cluster   *m_pDaughterCluster;         ///< The address of the daughter cluster
        float               m_boundedFraction;          ///< The fraction of daughter hits bounded by the cone defined by the vertex cluster
    };

    /**
     *  @brief  PfoAssociation class
     */
    class PfoAssociation
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pVertexPfo the address of the vertex pfo
         *  @param  pDaughterPfo the address of the daughter candidate pfo
         *  @param  clusterAssociationU the cluster association in the u view
         *  @param  clusterAssociationV the cluster association in the v view
         *  @param  clusterAssociationW the cluster association in the w view
         */
        PfoAssociation(pandora::Pfo *pVertexPfo, pandora::Pfo *pDaughterPfo, const ClusterAssociation &clusterAssociationU,
            const ClusterAssociation &clusterAssociationV, const ClusterAssociation &clusterAssociationW);

        /**
         *  @brief  Get the address of the vertex-associated pfo
         *
         *  @return the address of the vertex-associated pfo
         */
        pandora::Pfo *GetVertexPfo() const;

        /**
         *  @brief  Get the address of the non-vertex-associated candidate daughter pfo
         *
         *  @return the address of the non-vertex-associated candidate daughter pfo
         */
        pandora::Pfo *GetDaughterPfo() const;

        /**
         *  @brief  Get the mean bounded fraction, averaging over the u, v and w views
         *
         *  @return the mean bounded fraction, averaging over the u, v and w views
         */
        float GetMeanBoundedFraction() const;

        /**
         *  @brief  Get the maximum bounded fraction from the u, v and w views
         *
         *  @return the maximum bounded fraction from the u, v and w views
         */
        float GetMaxBoundedFraction() const;

        /**
         *  @brief  Get the minimum bounded fraction from the u, v and w views
         *
         *  @return the minimum bounded fraction from the u, v and w views
         */
        float GetMinBoundedFraction() const;

        /**
         *  @brief  Get the cluster association in the u view
         *
         *  @return the cluster association in the u view
         */
        const ClusterAssociation &GetClusterAssociationU() const;

        /**
         *  @brief  Get the cluster association in the v view
         *
         *  @return the cluster association in the v view
         */
        const ClusterAssociation &GetClusterAssociationV() const;

        /**
         *  @brief  Get the cluster association in the w view
         *
         *  @return the cluster association in the w view
         */
        const ClusterAssociation &GetClusterAssociationW() const;

        /**
         *  @brief  operator<
         * 
         *  @param  rhs the pfo association object for comparison
         * 
         *  @return boolean
         */
        bool operator< (const PfoAssociation &rhs) const;

    private:
        pandora::Pfo       *m_pVertexPfo;               ///< The address of the vertex-associated pfo
        pandora::Pfo       *m_pDaughterPfo;             ///< The address of the non-vertex-associated candidate daughter pfo

        float               m_meanBoundedFraction;      ///< The mean bounded fraction, averaging over the u, v and w views
        float               m_maxBoundedFraction;       ///< The maximum bounded fraction from the u, v and w views
        float               m_minBoundedFraction;       ///< The minimum bounded fraction from the u, v and w views

        ClusterAssociation  m_clusterAssociationU;      ///< The cluster association in the u view
        ClusterAssociation  m_clusterAssociationV;      ///< The cluster association in the v view
        ClusterAssociation  m_clusterAssociationW;      ///< The cluster association in the w view
    };

    typedef std::vector<PfoAssociation> PfoAssociationList;

    /**
     *  @brief  ConeParameters class
     */
    class ConeParameters
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pCluster address of the cluster
         *  @param  vertexPosition2D the event 2D vertex position
         *  @param  coneAngleCentile the cone angle centile
         */
        ConeParameters(const pandora::Cluster *const pCluster, const pandora::CartesianVector &vertexPosition2D, const float coneAngleCentile);

        /**
         *  @brief  Get the fraction of hits in a candidate daughter cluster bounded by the cone
         * 
         *  @param  pDaughterCluster the address of the daughter cluster
         *  @param  coneLengthMultiplier cnsider hits as bound if inside cone with projected distance less than N times cone length
         * 
         *  @return the bounded fraction
         */
        float GetBoundedFraction(const pandora::Cluster *const pDaughterCluster, const float coneLengthMultiplier) const;

    //private: // TODO put back in!
        /**
         *  @brief  Get the cone direction estimate, with apex fixed at the 2d vertex position
         * 
         *  @return the direction estimate
         */
        pandora::CartesianVector GetDirectionEstimate() const;

        /**
         *  @brief  Get the cone length
         * 
         *  @return rhe cone length
         */
        float GetConeLength() const;

        /**
         *  @brief  Get the cone cos half angle estimate
         * 
         *  @param  coneAngleCentile the cone angle centile
         * 
         *  @return the cone cos half angle estimate
         */
        float GetCosHalfAngleEstimate(const float coneAngleCentile) const;

        const pandora::Cluster     *m_pCluster;             ///< The parent cluster
        pandora::CartesianVector    m_apex;                 ///< The cone apex
        pandora::CartesianVector    m_direction;            ///< The cone direction
        float                       m_coneLength;           ///< The cone length
        float                       m_coneCosHalfAngle;     ///< The cone cos half angle
    };

    pandora::StatusCode Run();

    /**
     *  @brief  Get the list of input pfos and divide them into vertex-associated and non-vertex-associated lists
     * 
     *  @param  pVertex the address of the 3d vertex
     *  @param  vertexPfos to receive the list of vertex-associated pfos
     *  @param  nonVertexPfos to receive the list of nonvertex-associated pfos
     */
    void GetInputPfos(const pandora::Vertex *const pVertex, pandora::PfoList &vertexPfos, pandora::PfoList &nonVertexPfos) const;

    /**
     *  @brief  Get the list of associations between vertex-associated pfos and non-vertex-associated pfos
     * 
     *  @param  pVertex the address of the 3d vertex
     *  @param  vertexPfos the list of vertex-associated pfos
     *  @param  nonVertexPfos the list of nonvertex-associated pfos
     *  @param  pfoAssociationList to receive the pfo association list
     */
    void GetPfoAssociations(const pandora::Vertex *const pVertex, const pandora::PfoList &vertexPfos, const pandora::PfoList &nonVertexPfos,
        PfoAssociationList &pfoAssociationList) const;

    /**
     *  @brief  Process the list of pfo associations, merging the best-matching pfo
     * 
     *  @param  pfoAssociationList the pfo association list
     * 
     *  @return whether a pfo merge was made
     */
    bool ProcessPfoAssociations(const PfoAssociationList &pfoAssociationList) const;

    /**
     *  @brief  Whether a specified pfo is associated with a specified vertex
     * 
     *  @param  pPfo the address of the pfo
     *  @param  pVertex the address of the 3d vertex
     * 
     *  @return boolean
     */
    bool IsVertexAssociated(const pandora::Pfo *const pPfo, const pandora::Vertex *const pVertex) const;

    /**
     *  @brief  Get pfo association details between a vertex-associated pfo and a non-vertex associated daughter candidate pfo
     * 
     *  @param  pVertex the address of the 3d vertex
     *  @param  pVertexPfo the address of the vertex-associated pfo
     *  @param  pDaughterPfo the address of the non-vertex-associated pfo
     * 
     *  @return the pfo association details
     */
    PfoAssociation GetPfoAssociation(const pandora::Vertex *const pVertex, pandora::Pfo *const pVertexPfo, pandora::Pfo *const pDaughterPfo) const;

    /**
     *  @brief  Get cluster association details between a vertex-associated cluster and a non-vertex associated daughter candidate cluster
     * 
     *  @param  vertexPosition2D the projected vertex position
     *  @param  pVertexCluster the address of the vertex-associated cluster
     *  @param  pDaughterCluster the address of the non-vertex-associated cluster
     * 
     *  @return the cluster association details
     */
    ClusterAssociation GetClusterAssociation(const pandora::CartesianVector &vertexPosition2D, pandora::Cluster *const pVertexCluster,
        pandora::Cluster *const pDaughterCluster) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::set<pandora::HitType> HitTypeSet;
    typedef std::map<pandora::HitType, ClusterAssociation> HitTypeToAssociationMap;

    std::string     m_trackPfoListName;                 ///< The input track pfo list name
    std::string     m_showerPfoListName;                ///< The input shower pfo list name

    float           m_minVertexLongitudinalDistance;    ///< Vertex association check: min longitudinal distance cut
    float           m_maxVertexTransverseDistance;      ///< Vertex association check: max transverse distance cut
    unsigned int    m_minVertexAssociatedHitTypes;      ///< The min number of vertex associated hit types for a vertex associated pfo
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *VertexBasedPfoMergingAlgorithm::Factory::CreateAlgorithm() const
{
    return new VertexBasedPfoMergingAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *VertexBasedPfoMergingAlgorithm::ClusterAssociation::GetVertexCluster() const
{
    return m_pVertexCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *VertexBasedPfoMergingAlgorithm::ClusterAssociation::GetDaughterCluster() const
{
    return m_pDaughterCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float VertexBasedPfoMergingAlgorithm::ClusterAssociation::GetBoundedFraction() const
{
    return m_boundedFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Pfo *VertexBasedPfoMergingAlgorithm::PfoAssociation::GetVertexPfo() const
{
    return m_pVertexPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Pfo *VertexBasedPfoMergingAlgorithm::PfoAssociation::GetDaughterPfo() const
{
    return m_pDaughterPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const VertexBasedPfoMergingAlgorithm::ClusterAssociation &VertexBasedPfoMergingAlgorithm::PfoAssociation::GetClusterAssociationU() const
{
    return m_clusterAssociationU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const VertexBasedPfoMergingAlgorithm::ClusterAssociation &VertexBasedPfoMergingAlgorithm::PfoAssociation::GetClusterAssociationV() const
{
    return m_clusterAssociationV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const VertexBasedPfoMergingAlgorithm::ClusterAssociation &VertexBasedPfoMergingAlgorithm::PfoAssociation::GetClusterAssociationW() const
{
    return m_clusterAssociationW;
}

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_BASED_PFO_MERGING_ALGORITHM_H
