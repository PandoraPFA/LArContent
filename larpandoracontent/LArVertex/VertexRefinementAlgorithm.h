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

#include <Eigen/Dense>

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
     *  @brief  Perform the refinement proceduce on a list of vertices
     *
     *  @param  pVertexList address of the vertex list
     *  @param  caloHitList the list of calo hits (all views) to use in refining the vertex
     */
    void RefineVertices(const pandora::VertexList &vertexList, const pandora::CaloHitList &caloHitList) const;

    /**
     *  @brief  Refine the position of a two dimensional projection of a vertex using the calo hits in that view
     *
     *  @param  caloHitList the list of calo hits in a given view
     *  @param  seedVertex the provisional vertex to refine
     *
     *  @return the new refined position
     */
    pandora::CartesianVector RefineVertexTwoD(const pandora::CaloHitList &caloHitList, const pandora::CartesianVector &seedVertex) const;

    /**
     *  @brief  Retrieve the hits within m_hitRadii of a given centroid.
     *
     *  @param  hitVector the vector describing input hit positions
     *  @param  centroid the centre of the region to be searched
     *  @param  nearbyHitList the output calo hit list
     */
    void GetNearbyHits(const pandora::CaloHitVector &hitVector, const pandora::CartesianVector &centroid, pandora::CaloHitList &nearbyHitMatrix) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName; ///< The name of the input calo hit list
    std::string m_inputVertexListName; ///< The initial vertex list
    std::string m_outputVertexListName; ///< The refined vertex list to be outputted
    std::string m_primaryVertexListName; ///< The primary vertex list

    float m_hitRadii; ///< The search radius to collect hits around the vertex
    bool m_vetoPrimaryRegion; ///< Avoid refinements in the vicinity of the primary vertex
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_REFINEMENT_ALGORITHM_H
