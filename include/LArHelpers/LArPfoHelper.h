/**
 *  @file   LArContent/include/LArHelpers/LArPfoHelper.h
 *
 *  @brief  Header file for the cluster helper class.
 *
 *  $Log: $
 */
#ifndef LAR_PFO_HELPER_H
#define LAR_PFO_HELPER_H 1

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"

namespace lar
{

/**
 *  @brief  LArPfoHelper class
 */
class LArPfoHelper
{
public:


 
    /**
     *  @brief  Sort pfos by number of constituent hits
     * 
     *  @param  pLhs address of first pfo
     *  @param  pRhs address of second pfo
     */
    static bool SortByNHits(const pandora::ParticleFlowObject *const pLhs, const pandora::ParticleFlowObject *const pRhs);

    /**
     *  @brief  Read the vertex helper settings
     *
     *  @param  xmlHandle the relevant xml handle
     */
    static pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
 
};

} // namespace lar

#endif // #ifndef LAR_PFO_HELPER_H
