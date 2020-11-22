/**
 *  @file   larpandoracontent/LArThreeDReco/LArCosmicRay/DeltaRayParentAlgorithm.h
 *
 *  @brief  Header file for the delta ray parent class
 *
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_PARENT_ALGORITHM_H
#define LAR_DELTA_RAY_PARENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

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
    void InitialisePfoLengthMap(const pandora::PfoList *const muonPfoList, const pandora::PfoList *const deltaRayPfoList, PfoLengthMap &pfoLengthMap) const;

    void FindParentPfo(const PfoLengthMap &pfoLengthMap, const pandora::ParticleFlowObject *const pPfo, const pandora::ParticleFlowObject *&pParentPfo) const;

    void AssignToParentPfo(const pandora::PfoList *const muonPfoList, const pandora::PfoList *const deltaRayPfoList, const pandora::ParticleFlowObject *const pPfo,
        const pandora::ParticleFlowObject *const pParentPfo, PfoLengthMap &pfoLengthMap) const;

    void UpdatePfoLengthMap(const pandora::PfoList &pfosToRemove, const pandora::ParticleFlowObject *const pPfoToAdd, PfoLengthMap &pfoLengthMap) const;
    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_muonPfoListName;
    std::string m_deltaRayPfoListName;
    float m_distanceForMatching;
};
    
} // namespace lar_content

#endif // #ifndef LAR_THREE_VIEW_SHOWERS_ALGORITHM_H
