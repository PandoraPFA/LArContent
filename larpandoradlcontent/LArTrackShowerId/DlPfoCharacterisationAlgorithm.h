/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlPfoCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the cut based pfo characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_DL_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/PfoCharacterisationBaseAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DlPfoCharacterisationAlgorithm class
 */
class DlPfoCharacterisationAlgorithm : public lar_content::PfoCharacterisationBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlPfoCharacterisationAlgorithm();

private:
    bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    bool IsClearTrack(const pandora::ParticleFlowObject *const pPfo) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_PFO_CHARACTERISATION_ALGORITHM_H
