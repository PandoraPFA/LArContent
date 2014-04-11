/**
 *  @file   LArContent/include/LArCheating/CheatingPfoCreationAlgorithm.h
 * 
 *  @brief  Header file for the cheating cluster creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_PFO_CREATION_ALGORITHM_H
#define LAR_CHEATING_PFO_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  CheatingPfoCreationAlgorithm class
 */
class CheatingPfoCreationAlgorithm : public pandora::Algorithm
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

inline pandora::Algorithm *CheatingPfoCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CheatingPfoCreationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_CHEATING_PFO_CREATION_ALGORITHM_H
