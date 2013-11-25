/**
 *  @file   LArContent/include/LArUtility/HelloWorldAlgorithm.h
 * 
 *  @brief  Header file for the hello world algorithm
 * 
 *  $Log: $
 */
#ifndef LAR_HELLO_WORLD_ALGORITHM_H
#define LAR_HELLO_WORLD_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  HelloWorldAlgorithm class
 */
class HelloWorldAlgorithm : public pandora::Algorithm
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

inline pandora::Algorithm *HelloWorldAlgorithm::Factory::CreateAlgorithm() const
{
    return new HelloWorldAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_HELLO_WORLD_ALGORITHM_H
