/**
 *  @file   LArContent/include/LArCheating/CheatingNeutrinoDaughterVerticesAlgorithm.h
 * 
 *  @brief  Header file for the cheating neutrino daughter vertices algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_NEUTRINO_DAUGHTER_VERTICES_ALGORITHM_H
#define LAR_CHEATING_NEUTRINO_DAUGHTER_VERTICES_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingNeutrinoDaughterVerticesAlgorithm::Algorithm class
 */
class CheatingNeutrinoDaughterVerticesAlgorithm : public pandora::Algorithm
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

inline pandora::Algorithm *CheatingNeutrinoDaughterVerticesAlgorithm::Factory::CreateAlgorithm() const
{
    return new CheatingNeutrinoDaughterVerticesAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_NEUTRINO_DAUGHTER_VERTICES_ALGORITHM_H
