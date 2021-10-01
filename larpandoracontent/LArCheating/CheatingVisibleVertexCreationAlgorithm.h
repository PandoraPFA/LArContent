/**
 *  @file   larpandoracontent/LArCheating/CheatingVisibleVertexCreationAlgorithm.h
 *
 *  @brief  Header file for the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_VISIBLE_VERTEX_CREATION_ALGORITHM_H
#define LAR_CHEATING_VISIBLE_VERTEX_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CheatingVisibleVertexCreationAlgorithm class
 */
class CheatingVisibleVertexCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingVisibleVertexCreationAlgorithm();

    void SetVertex(const pandora::ParticleFlowObject *const pPfo) const;

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_inputPfoListNames;
    std::string m_vertexListName;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_PFO_CREATION_ALGORITHM_H
