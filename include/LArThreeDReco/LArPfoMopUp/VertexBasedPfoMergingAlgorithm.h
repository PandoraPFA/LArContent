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
     *  @brief  PfoAssociation class
     */
    class PfoAssociation
    {
    public:
    private:
        
    };

    typedef std::vector<PfoAssociation> PfoAssociationList;

    pandora::StatusCode Run();

    /**
     *  @brief  Get the list of input pfos and divide them into vertex-associated and non-vertex-associated lists
     * 
     *  @param  vertexPfos to receive the list of vertex-associated pfos
     *  @param  nonVertexPfos to receive the list of nonvertex-associated pfos
     */
    void GetInputPfos(pandora::PfoList &vertexPfos, pandora::PfoList &nonVertexPfos) const;

    /**
     *  @brief  Get the list of associations between vertex-associated pfos and non-vertex-associated pfos
     * 
     *  @param  vertexPfos the list of vertex-associated pfos
     *  @param  nonVertexPfos the list of nonvertex-associated pfos
     *  @param  pfoAssociationList to receive the pfo association list
     */
    void GetPfoAssociations(const pandora::PfoList &vertexPfos, const pandora::PfoList &nonVertexPfos, PfoAssociationList &pfoAssociationList) const;

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

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::set<pandora::HitType> HitTypeSet;      ///< The hit type set typedef

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

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_BASED_PFO_MERGING_ALGORITHM_H
