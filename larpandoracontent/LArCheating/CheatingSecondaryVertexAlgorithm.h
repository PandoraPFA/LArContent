/**
 *  @file   larpandoracontent/LArCheating/CheatingSecondaryVertexAlgorithm.h
 *
 *  @brief  Header file for the cheating vertex creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_SECONDARY_VERTEX_ALGORITHM_H
#define LAR_CHEATING_SECONDARY_VERTEX_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingSecondaryVertexAlgorithm::Algorithm class
 */
class CheatingSecondaryVertexAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingSecondaryVertexAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputCaloHitListName; ///< The name of the input calo hit list
    std::string m_outputVertexListName; ///< The name under which to save the output vertex list
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_SECONDARY_VERTEX_ALGORITHM_H
