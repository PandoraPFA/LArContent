/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional hit creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_HIT_CREATION_ALGORITHM_H
#define LAR_THREE_D_HIT_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDHitCreationAlgorithm::Algorithm class
 */
class ThreeDHitCreationAlgorithm : public pandora::Algorithm
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

    std::string     m_inputPfoListName;         ///< The name of the input pfo list
    std::string     m_outputCaloHitListName;    ///< The name of the output calo hit list
    std::string     m_outputClusterListName;    ///< The name of the output cluster list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDHitCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDHitCreationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_HIT_CREATION_ALGORITHM_H
