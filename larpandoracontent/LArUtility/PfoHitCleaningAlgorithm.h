/**
 *  @file   larpandoracontent/LArUtility/PfoHitCleaningAlgorithm.h
 *
 *  @brief  Header file for the pfo hit cleaning algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PFO_HIT_CLEANING_ALGORITHM_H
#define LAR_PFO_HIT_CLEANING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  PfoHitCleaningAlgorithm class
 */
class PfoHitCleaningAlgorithm : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_pfoListNames;     ///< The list of pfo list names
    pandora::StringVector m_clusterListNames; ///< The list of cluster list names
};

} // namespace lar_content

#endif // #ifndef LAR_PFO_HIT_CLEANING_ALGORITHM_H
