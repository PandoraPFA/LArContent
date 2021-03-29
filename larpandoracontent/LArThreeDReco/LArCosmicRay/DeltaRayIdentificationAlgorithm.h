/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayIdentificationAlgorithm.h
 *
 *  @brief  Header file for the delta ray identification algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_IDENTIFICATION_ALGORITHM_H
#define LAR_DELTA_RAY_IDENTIFICATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  DeltaRayIdentificationAlgorithm class
 */
class DeltaRayIdentificationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DeltaRayIdentificationAlgorithm();

private:
    pandora::StatusCode Run();

    typedef std::unordered_map<const pandora::ParticleFlowObject *, const pandora::ParticleFlowObject *> PfoAssociationMap;

    /**
     *  @brief Get the vector of Pfos, given the input list name
     *
     *  @param inputPfoListName the input Pfo list name
     *  @param outputPfoVector the output vector of Pfos
     */
    void GetPfos(const std::string &inputPfoListName, pandora::PfoVector &outputPfoVector) const;

    /**
     *  @brief Build parent/daughter associations between PFOs
     *
     *  @param inputPfos the input vector of current parent Pfos
     *  @param outputPfos the input vector of current daughter Pfos
     *  @param pfoAssociationMap the output map of parent/daughter associations
     */
    void BuildAssociationMap(const pandora::PfoVector &inputPfos, const pandora::PfoVector &outputPfos, PfoAssociationMap &pfoAssociationMap) const;

    /**
     *  @brief Determine if a given pair of Pfos have a parent/daughter association
     *
     *  @param pDaughterPfo the input daughter Pfo
     *  @param pParentPfo the input parent Pfo
     *  @param displacement the average displacement between parent and daughter
     *
     *  @return boolean
     */
    bool IsAssociated(const pandora::ParticleFlowObject *const pDaughterPfo, const pandora::ParticleFlowObject *const pParentPfo, float &displacement) const;

    /**
     *  @brief Calculate 2D separation between two Pfos
     *
     *  @param pDaughterPfo the input daughter Pfo
     *  @param pParentPfo the input parent Pfo
     *
     *  @return average displacement between parent and daughter
     */
    float GetTwoDSeparation(const pandora::ParticleFlowObject *const pDaughterPfo, const pandora::ParticleFlowObject *const pParentPfo) const;

    /**
     *  @brief Calculate 2D separation between two Pfos
     *
     *  @param pPfo the input daughter Pfo
     *  @param hitType the hit type
     *  @param vertexVector the vector of possible vertex positions
     */
    void GetTwoDVertices(const pandora::ParticleFlowObject *const pPfo, const pandora::HitType &hitType, pandora::CartesianPointVector &vertexVector) const;

    /**
     *  @brief Calculate closest 2D separation between a set of vertices and a set of clusters
     *
     *  @param vertexVector the vector of possible vertex positions
     *  @param clusterList the list of clusters
     */
    float GetClosestDistance(const pandora::CartesianPointVector &vertexVector, const pandora::ClusterList &clusterList) const;

    /**
     *  @brief Build the parent/daughter links from the map of parent/daughter associations
     *
     *  @param pfoAssociationMap the map of parent/daughter associations
     *  @param outputPfoList the output list of daughter Pfos
     */
    void BuildParentDaughterLinks(const PfoAssociationMap &pfoAssociationMap, pandora::PfoList &outputPfoList) const;

    /**
     *  @brief For a given daughter, follow the parent/daughter links to find the overall parent
     *
     *  @param pfoAssociationMap the map of parent/daughter associations
     *  @param pPfo the daughter Pfo
     *
     *  @return the parent Pfo
     */
    const pandora::ParticleFlowObject *GetParent(const PfoAssociationMap &pfoAssociationMap, const pandora::ParticleFlowObject *const pPfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_parentPfoListName;   ///< The parent pfo list name
    std::string m_daughterPfoListName; ///< The daughter pfo list name

    float m_distanceForMatching;      ///< Maximum allowed distance of delta ray from parent cosmic ray
    float m_minParentLengthSquared;   ///< Minimum allowed length of parent cosmic ray
    float m_maxDaughterLengthSquared; ///< Maximum allowed length of daughter delta ray
};

} // namespace lar_content

#endif // #ifndef LAR_DELTA_RAY_IDENTIFICATION_ALGORITHM_H
