/**
 *  @file   LArContent/include/LArUtility/TwoDPreparationAlgorithm.h
 * 
 *  @brief  Header file for the two dimensional preparation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_TWO_D_PREPARATION_ALGORITHM_H
#define LAR_TWO_D_PREPARATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  TwoDPreparationAlgorithm::Algorithm class
 */
class TwoDPreparationAlgorithm : public pandora::Algorithm
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

    std::string     m_clusterListName;  ///< The cluster list name to set as the current cluster list
    std::string     m_caloHitListName;  ///< The calo hit list name to set as the current calo hit list
    std::string     m_vertexName;       ///< The vertex name to set as the current vertex
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *TwoDPreparationAlgorithm::Factory::CreateAlgorithm() const
{
    return new TwoDPreparationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_TWO_D_PREPARATION_ALGORITHM_H
