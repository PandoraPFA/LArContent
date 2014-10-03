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

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *VertexSelectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new VertexSelectionAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_SELECTION_ALGORITHM_H
