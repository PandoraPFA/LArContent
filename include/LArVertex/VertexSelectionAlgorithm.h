/**
 *  @file   LArContent/include/LArUtility/VertexSelectionAlgorithm.h
 * 
 *  @brief  Header file for the vertex selection algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_SELECTION_ALGORITHM_H
#define LAR_VERTEX_SELECTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  VertexSelectionAlgorithm::Algorithm class
 */
class VertexSelectionAlgorithm : public pandora::Algorithm
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
    VertexSelectionAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_outputVertexListName;         ///< The name under which to save the output vertex list
    bool            m_replaceCurrentVertexList;     ///< Whether to replace the current vertex list with the output list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *VertexSelectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new VertexSelectionAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_SELECTION_ALGORITHM_H
