/**
 *  @file   LArContent/include/LArThreeDReco/LArCosmicRay/DeltaRayIdentificationAlgorithm.h
 *
 *  @brief  Header file for the delta ray identification algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_IDENTIFICATION_ALGORITHM_H
#define LAR_DELTA_RAY_IDENTIFICATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  DeltaRayIdentificationAlgorithm class
 */
class DeltaRayIdentificationAlgorithm : public pandora::Algorithm
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

private:
    pandora::StatusCode Run();

    typedef std::map<const pandora::ParticleFlowObject*, const pandora::ParticleFlowObject*> PfoAssociationMap;

    /**
     *  @brief Get the list of Pfos, given the input list name
     *
     *  @param inputPfoListName the input Pfo list name
     *  @param pfoList the output list of Pfos
     */
    void GetPfos(const std::string inputPfoListName, pandora::PfoList &pfoList) const;

    /**
     *  @brief Build parent/daughter associations between PFOs
     *
     *  @param inputPfos the input list of current parent Pfos
     *  @param outputPfos the input list of current daughter Pfos
     *  @param pfoAssociationMap the output map of parent/daughter associations
     */
    void BuildAssociationMap(const pandora::PfoList &inputPfos, const pandora::PfoList &outputPfos,
        PfoAssociationMap &pfoAssociationMap) const;

    /**
     *  @brief Determine if a given pair of Pfos have a parent/daughter association
     *
     *  @param pDaughterPfo the input daughter Pfo
     *  @param pParentPfo the input parent Pfo
     *  @param displacementSquared the average squared displacement between the parent and daughter
     *
     *  @return boolean
     */
    bool IsAssociated(const pandora::ParticleFlowObject *const pDaughterPfo, const pandora::ParticleFlowObject *const pParentPfo,
        float &displacementSquared) const;

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
    pandora::ParticleFlowObject *GetParent(const PfoAssociationMap &pfoAssociationMap, const pandora::ParticleFlowObject *const pPfo) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_inputPfoListName;          ///< The parent pfo list name
    std::string     m_outputPfoListName;         ///< The daughter pfo list name
    float           m_maxDisplacementSquared;    ///< Maximum allowed distance of delta ray from parent cosmic ray
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *DeltaRayIdentificationAlgorithm::Factory::CreateAlgorithm() const
{
    return new DeltaRayIdentificationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_DELTA_RAY_IDENTIFICATION_ALGORITHM_H
