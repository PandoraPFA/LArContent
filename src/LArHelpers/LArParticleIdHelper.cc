/**
 *  @file   LArContent/src/LArHelpers/LArParticleIdHelper.cc
 * 
 *  @brief  Implementation of the lar particle id class.
 * 
 *  $Log: $
 */

#include "Helpers/XmlHelper.h"

#include "Objects/Cluster.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArParticleIdHelper.h"

namespace lar
{

using namespace pandora;

bool LArParticleIdHelper::LArEmShowerId(const Cluster *const /*pCluster*/)
{
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArParticleIdHelper::LArPhotonId(const Cluster *const /*pCluster*/)
{
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArParticleIdHelper::LArElectronId(const Cluster *const /*pCluster*/)
{
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArParticleIdHelper::LArMuonId(const Cluster *const pCluster)
{
    if (LArClusterHelper::GetLayerOccupancy(pCluster) < 0.75f)
        return false;

    if (LArClusterHelper::LArTrackWidth(pCluster) > 0.5f)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArParticleIdHelper::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar
