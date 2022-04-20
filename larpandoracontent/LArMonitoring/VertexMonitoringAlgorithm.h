/**
 *  @file   larpandoracontent/LArMonitoring/VertexMonitoringAlgorithm.h
 *
 *  @brief  Header file for the particle visualisation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_MONITORING_ALGORITHM_H
#define LAR_VERTEX_MONITORING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  VertexMonitoringAlgorithm class
 */
class VertexMonitoringAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    VertexMonitoringAlgorithm();

    virtual ~VertexMonitoringAlgorithm();

private:
    pandora::StatusCode AssessVertices() const;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_visualise;               // Whether to produce visual monitoring output
    bool m_writeFile;               // Whether to produce ROOT output file
    std::string m_filename;         // The filename of the ROOT output file
    std::string m_treename;         // The name of the ROOT tree
    float m_transparencyThresholdE; ///< Cell energy for which transparency is saturated (0%, fully opaque)
    float m_energyScaleThresholdE;  ///< Cell energy for which color is at top end of continous color palette
    float m_scalingFactor;          ///< TEve works with [cm], Pandora usually works with [mm] (but LArContent went with cm too)
};

} // namespace lar_content

#endif // LAR_VERTEX_MONITORING_ALGORITHM_H
