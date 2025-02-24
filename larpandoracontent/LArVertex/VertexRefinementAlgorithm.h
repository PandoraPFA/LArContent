/**
 *  @file   larpandoracontent/LArVertex/VertexRefinementAlgorithm.h
 *
 *  @brief  Header file for the vertex refinement algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_REFINEMENT_ALGORITHM_H
#define LAR_VERTEX_REFINEMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  VertexRefinementAlgorithm class
 */
class VertexRefinementAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexRefinementAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Get the input cluster lists
     *
     *  @param  inputClusterListNames the input cluster list names
     *  @param  clusterListU the U-view cluster list to populate
     *  @param  clusterListV the V-view cluster list to populate
     *  @param  clusterListW the W-view cluster list to populate
     */
    void GetClusterLists(const pandora::StringVector &inputClusterListNames, pandora::ClusterList &clusterListU,
        pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW) const;

    /**
     *  @brief  Perform the refinement proceduce on a list of vertices
     *
     *  @param  pVertexList address of the vertex list
     *  @param  clusterListU the list of U-view clusters
     *  @param  clusterListV the list of V-view clusters
     *  @param  clusterListW the list of W-view clusters
     */
    void RefineVertices(const pandora::VertexList *const pVertexList, const pandora::ClusterList &clusterListU,
        const pandora::ClusterList &clusterListV, const pandora::ClusterList &clusterListW) const;

    /**
     *  @brief  Refine the position of a two dimensional projection of a vertex using the clusters in that view
     *
     *  @param  clusterList the list of two dimensional clusters
     *  @param  originalVtxPos the original vertex position projected into two dimensions
     *
     *  @return the new refined position
     */
    pandora::CartesianVector RefineVertexTwoD(const pandora::ClusterList &clusterList, const pandora::CartesianVector &originalVtxPos) const;

    /**
     *  @brief  Calculate the best fit point of a set of lines using a matrix equation
     *
     *  @param  intercepts the vector of the defining points of the lines
     *  @param  directions the vector of line directions
     *  @param  weights the vector of weights for each line
     *  @param  bestFitPoint the resulting best fit point
     */
    void GetBestFitPoint(const pandora::CartesianPointVector &intercepts, const pandora::CartesianPointVector &directions,
        const pandora::FloatVector &weights, pandora::CartesianVector &bestFitPoint) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_inputClusterListNames; ///< The list of input cluster list names
    std::string m_inputVertexListName;             ///< The initial vertex list
    std::string m_outputVertexListName;            ///< The refined vertex list to be outputted

    float m_chiSquaredCut;         ///< The maximum chi2 value a refined vertex can have to be kept
    float m_distanceCut;           ///< The maximum distance a refined vertex can be from the original position to be kept
    unsigned int m_minimumHitsCut; ///< The minimum size of a cluster to be used in refinement
    float m_twoDDistanceCut;       ///< The maximum distance a cluster can be from the original position to be used in refinement
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_REFINEMENT_ALGORITHM_H
