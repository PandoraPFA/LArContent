/**
 *  @file   LArContent/include/LArUtility/ThreeDPreparationAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional preparation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_PREPARATION_ALGORITHM_H
#define LAR_THREE_D_PREPARATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDPreparationAlgorithm::Algorithm class
 */
class ThreeDPreparationAlgorithm : public pandora::Algorithm
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

    std::string     m_vertexName;       ///< The vertex name to set as the current vertex
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDPreparationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDPreparationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_PREPARATION_ALGORITHM_H
