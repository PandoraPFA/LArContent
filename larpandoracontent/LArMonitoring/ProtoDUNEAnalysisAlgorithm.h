/**
 *  @file   larpandoracontent/LArMonitoring/ProtoDUNEAnalysisAlgorithm.h
 *
 *  @brief  Header file for the ProtoDUNE data analysis algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_PROTODUNE_ANALYSIS_ALGORITHM_H
#define LAR_PROTODUNE_ANALYSIS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ProtoDUNEAnalysisAlgorithm class
 */
class ProtoDUNEAnalysisAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ProtoDUNEAnalysisAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~ProtoDUNEAnalysisAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_mcParticleListName;           ///< Name of input MC particle list
    std::string             m_pfoListName;                  ///< Name of input Pfo list

    std::string             m_treeName;                     ///< Name of output tree
    std::string             m_fileName;                     ///< Name of output file

    bool                    m_visualDisplay;                ///< Visual display

    int                     m_eventNumber;                  ///< The event number
};

} // namespace lar_content

#endif //LAR_PROTODUNE_ANALYSIS_ALGORITHM_H
