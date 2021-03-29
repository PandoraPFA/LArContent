/**
 *  @file   larpandoracontent/LArThreeDReco/LArEventBuilding/NeutrinoCreationAlgorithm.h
 *
 *  @brief  Header file for the neutrino creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_CREATION_ALGORITHM_H
#define LAR_NEUTRINO_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoCreationAlgorithm class
 */
class NeutrinoCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    NeutrinoCreationAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Force creation of a single neutrino, with no vertex, regardless of number of input vertices
     */
    pandora::StatusCode ForceSingleEmptyNeutrino() const;

    /**
     *  @brief  Fill provided pfo parameters with default/dummy values for later refinement
     *
     *  @param  pfoParameters the pfo parameters
     */
    void FillDefaultNeutrinoParameters(PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_vertexListName;      ///< The name of the neutrino vertex list
    std::string m_neutrinoPfoListName; ///< The name of the neutrino pfo list

    bool m_forceSingleEmptyNeutrino; ///< Whether to force creation of a single neutrino, with no vertex, regardless of number of input vertices
};

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_CREATION_ALGORITHM_H
