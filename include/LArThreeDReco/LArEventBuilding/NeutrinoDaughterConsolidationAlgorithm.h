/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/NeutrinoDaughterConsolidationAlgorithm.h
 *
 *  @brief  Header file for the neutrino daughter consolidation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_DAUGHTER_CONSOLIDATION_ALGORITHM_H
#define LAR_NEUTRINO_DAUGHTER_CONSOLIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoDaughterConsolidationAlgorithm class
 */
class NeutrinoDaughterConsolidationAlgorithm : public pandora::Algorithm
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
    NeutrinoDaughterConsolidationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *NeutrinoDaughterConsolidationAlgorithm::Factory::CreateAlgorithm() const
{
    return new NeutrinoDaughterConsolidationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_DAUGHTER_CONSOLIDATION_ALGORITHM_H
