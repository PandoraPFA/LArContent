/**
 *  @file   larpandoracontent/LArMonitoring/VertexAssessmentAlgorithm.h
 *
 *  @brief  Header file for the vertex assessment algorithm class
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_ASSESSMENT_ALGORITHM_H
#define LAR_VERTEX_ASSESSMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief VertexAssessmentAlgorithm class
 */
class VertexAssessmentAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexAssessmentAlgorithm();

    /**
     *  @brief  Default destructor
     */
    ~VertexAssessmentAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Visualize the 2D calo hits and projected 2D reconstructed and (if available) true vertices
     */
    void Visualize() const;

    /**
     *  @brief  Retrieve the map from reconstructable MC particles to their hits
     *
     *  @param  mcToHitsMap the output map from MC particles to hits
     */
    pandora::StatusCode GetMCToHitsMap(LArMCParticleHelper::MCContributionMap &mcToHitsMap) const;

    bool m_visualize;   ///< Whether or not to visualize the vertices
};

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_ASSESSMENT_ALGORITHM_H
