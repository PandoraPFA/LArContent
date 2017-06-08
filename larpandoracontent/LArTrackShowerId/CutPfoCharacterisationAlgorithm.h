/**
 *  @file   larpandoracontent/LArTrackShowerId/CutPfoCharacterisationAlgorithm.h
 * 
 *  @brief  Header file for the cut based pfo characterisation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_PFO_CHARACTERISATION_ALGORITHM_H
#define LAR_PFO_CHARACTERISATION_ALGORITHM_H 1

#include "larpandoracontent/LArTrackShowerId/PfoCharacterisationBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  CutPfoCharacterisationAlgorithm class
 */
class CutPfoCharacterisationAlgorithm : public PfoCharacterisationBaseAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CutPfoCharacterisationAlgorithm();

private:
    bool IsClearTrack(const pandora::Cluster *const pCluster) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool                    m_postBranchAddition;           ///< Whether to use configuration for shower clusters post branch addition
    unsigned int            m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    unsigned int            m_slidingShowerFitWindow;       ///< The layer window for the sliding shower fits
    float                   m_maxShowerLengthCut;           ///< The maximum cluster length to qualify as a shower
    float                   m_dTdLWidthRatioCut;            ///< The maximum ratio of transverse fit gradient width to straight line length to qualify as a track
    float                   m_vertexDistanceRatioCut;       ///< The maximum ratio of vertex separation to straight line length to qualify as a track
    float                   m_showerWidthRatioCut;          ///< The maximum ratio of shower fit width to straight line length to qualify as a track
};

} // namespace lar_content

#endif // #ifndef LAR_PFO_CHARACTERISATION_ALGORITHM_H
