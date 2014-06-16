/**
 *  @file   LArContent/include/LArThreeDReco/LArCosmicRay/DeltaRayIdentificationAlgorithm.h
 * 
 *  @brief  Header file for the delta ray identification algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_IDENTIFICATION_ALGORITHM_H
#define LAR_DELTA_RAY_IDENTIFICATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  DeltaRayIdentificationAlgorithm class
 */
class DeltaRayIdentificationAlgorithm : public pandora::Algorithm
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

    std::string     m_inputPfoListName;             ///< The input pfo list name
    std::string     m_outputPfoListName;            ///< The output pfo list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *DeltaRayIdentificationAlgorithm::Factory::CreateAlgorithm() const
{
    return new DeltaRayIdentificationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_DELTA_RAY_IDENTIFICATION_ALGORITHM_H
