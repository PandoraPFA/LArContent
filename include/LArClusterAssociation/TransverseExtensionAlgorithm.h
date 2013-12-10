/**
 *  @file   LArContent/include/LArClusterAssociation/TransverseExtensionAlgorithm.h
 * 
 *  @brief  Header file for the hello world algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_TRANSVERSE_EXTENSION_ALGORITHM_H
#define LAR_TRANSVERSE_EXTENSION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  TransverseExtensionAlgorithm class
 */
class TransverseExtensionAlgorithm : public pandora::Algorithm
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

inline pandora::Algorithm *TransverseExtensionAlgorithm::Factory::CreateAlgorithm() const
{
    return new TransverseExtensionAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_TRANSVERSE_EXTENSION_ALGORITHM_H
