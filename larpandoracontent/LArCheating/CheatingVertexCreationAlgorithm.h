/**
 *  @file   larpandoracontent/LArCheating/CheatingVertexCreationAlgorithm.h
 *
 *  @brief  Header file for the cheating vertex creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_VERTEX_CREATION_ALGORITHM_H
#define LAR_CHEATING_VERTEX_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CheatingVertexCreationAlgorithm::Algorithm class
 */
class CheatingVertexCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingVertexCreationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Functions for coordinate transformations to match SCE maps
     *
     *  @param xPosition the initial position to transformate, in each coordinate separately
     *
     *  @return the transformed position
     */
    float TransformX(float xPosition) const;
    float TransformY(float yPosition) const;
    float TransformZ(float zPosition) const;

    std::string     m_outputVertexListName;         ///< The name under which to save the output vertex list
    std::string     m_pathToSCEFile;                ///< The path to the file with the SCE map to use for corrections
    bool            m_replaceCurrentVertexList;     ///< Whether to replace the current vertex list with the output list
    float           m_vertexXCorrection;            ///< The vertex x correction, added to reported mc neutrino endpoint x value, in cm
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_VERTEX_CREATION_ALGORITHM_H
