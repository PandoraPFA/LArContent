/**
 *  @file   LArContent/include/LArCheating/CheatingClusterCharacterisationAlg.h
 * 
 *  @brief  Header file for the cluster characterisation cheater class.
 * 
 *  $Log: $
 */
#ifndef LAR_CHEATING_CLUSTER_CHARACTERISATION_ALG_H
#define LAR_CHEATING_CLUSTER_CHARACTERISATION_ALG_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  CheatingClusterCharacterisationAlg class
 */
class CheatingClusterCharacterisationAlg : public pandora::Algorithm
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

    std::string     m_inputClusterListName;     ///< The name of the input cluster list. If not specified, will access current list.
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CheatingClusterCharacterisationAlg::Factory::CreateAlgorithm() const
{
    return new CheatingClusterCharacterisationAlg();
}

} // namespace lar

#endif // #ifndef LAR_CHEATING_CLUSTER_CHARACTERISATION_ALG_H
