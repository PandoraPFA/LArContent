/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/NeutrinoVertexCreationAlgorithm.h
 *
 *  @brief  Header file for the neutrino vertex creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_VERTEX_CREATION_ALGORITHM_H
#define LAR_NEUTRINO_VERTEX_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoVertexCreationAlgorithm class
 */
class NeutrinoVertexCreationAlgorithm : public pandora::Algorithm
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
    NeutrinoVertexCreationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the list of selected vertices
     *
     *  @param  vertexList to received the list of selected vertices
     */
    void GetInputVertexList(pandora::VertexList &vertexList) const;

    /**
     *  @brief  Get the list of daughter pfos to be added to the new neutrino pfo
     *
     *  @param  pfoList to receive the list of daughter pfos
     */
    void GetDaughterPfoList(pandora::PfoList &pfoList) const;

    std::string             m_inputVertexListName;        ///< The name of the input vertex list
    std::string             m_outputVertexListName;       ///< The name of the output vertex list
    pandora::StringVector   m_daughterPfoListNames;       ///< The list of daughter pfo list names
    bool                    m_replaceCurrentVertexList;   ///< Whether to replace the current vertex list with the output list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *NeutrinoVertexCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new NeutrinoVertexCreationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_VERTEX_CREATION_ALGORITHM_H
