/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/UnattachedDeltaRaysAlgorithm.h
 *
 *  @brief  Header file for the unattached delta rays algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_UNATTACHED_DELTA_RAYS_ALGORITHM_H
#define LAR_UNATTACHED_DELTA_RAYS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  UnattachedDeltaRaysAlgorithm class
 */
class UnattachedDeltaRaysAlgorithm : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_pfoListName; ///< The pfo list name
};

} // namespace lar_content

#endif // #ifndef LAR_UNATTACHED_DELTA_RAYS_ALGORITHM_H
