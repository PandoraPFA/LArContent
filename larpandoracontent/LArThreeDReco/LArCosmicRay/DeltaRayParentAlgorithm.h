/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayParentAlgorithm.h
 *
 *  @brief  Header file for the delta ray parent class.
 *
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_PARENT_ALGORITHM_H
#define LAR_DELTA_RAY_PARENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  DeltaRayParentAlgorithm class
 */
class DeltaRayParentAlgorithm : public pandora::Algorithm
{
public:
    typedef std::map<const pandora::ParticleFlowObject *, float> PfoLengthMap;

    /**
     *  @brief  Default constructor
     */
    DeltaRayParentAlgorithm();

    pandora::StatusCode Run();

private:
    /**
     *  @brief  Initialise the delta ray pfo length map
     *
     *  @param  muonPfoList the list of all cosmic ray pfos
     *  @param  deltaRayPfoList the list of all delta ray pfos
     *  @param  pfoLengthMap the output mapping of pfos to their 2D length
     */
    void InitialisePfoLengthMap(const pandora::PfoList *const muonPfoList, const pandora::PfoList *const deltaRayPfoList, PfoLengthMap &pfoLengthMap) const;

    /**
     *  @brief  Identify the parent pfo of a given delta ray pfo (can be either a cosmic ray or delta ray pfo)
     *
     *  @param  pfoLengthMap the mapping of pfos to their 2D length
     *  @param  pPfo the address of the input delta ray pfo
     *  @param  pParentPfo the output address of the parent pfo
     */
    void FindParentPfo(const PfoLengthMap &pfoLengthMap, const pandora::ParticleFlowObject *const pPfo, const pandora::ParticleFlowObject *&pParentPfo) const;

    /**
     *  @brief  Get distance between two Pfos using 2D clusters
     *
     *  @param  pPfo the address of the first Pfo
     *  @param  pPfo the address of the second Pfo
     *  @param  to recieve the output pfo separation
     *
     *  @return  whether the pfo separation could be calculated
     */
    pandora::StatusCode GetTwoDSeparation(
        const pandora::ParticleFlowObject *const pPfo1, const pandora::ParticleFlowObject *const pPfo2, float &separation) const;

    /**
     *  @brief  Apply parent-child link (if parent is a cosmic ray create parent-child link else merge the delta ray cluster into parent delta ray pfo)
     *
     *  @param  muonPfoList the list of all cosmic ray pfos
     *  @param  deltaRayPfoList the list of all delta ray pfos
     *  @param  pPfo the address of the input delta ray pfo
     *  @param  pParentPfo the address of the parent pfo
     *  @param  pfoLengthMap the mapping of pfos to their 2D length
     */
    void AssignToParentPfo(const pandora::PfoList *const muonPfoList, const pandora::PfoList *const deltaRayPfoList,
        const pandora::ParticleFlowObject *const pPfo, const pandora::ParticleFlowObject *const pParentPfo, PfoLengthMap &pfoLengthMap) const;

    /**
     *  @brief  Update the pfo length map after a parent-child delta ray merge
     *
     *  @param  pfosToRemove the list of pfos to remove from the map
     *  @param  pPfoToAdd the address of the pfo to add to the map
     *  @param  pfoLengthMap the mapping of pfos to their 2D length
     */
    void UpdatePfoLengthMap(const pandora::PfoList &pfosToRemove, const pandora::ParticleFlowObject *const pPfoToAdd, PfoLengthMap &pfoLengthMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_muonPfoListName;     ///< The list of reconstructed muon pfos
    std::string m_deltaRayPfoListName; ///< The list of reconstructed delta ray pfos
    float m_distanceForMatching;       ///< The maximum separation of a delta ray pfo from its parent
};

} // namespace lar_content

#endif // #ifndef LAR_DELTA_RAY_PARENT_ALGORITHM_H
