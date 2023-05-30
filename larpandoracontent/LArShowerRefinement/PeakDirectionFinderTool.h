/**
 *  @file   larpandoracontent/LArShowerRefinement/PeakDirectionFinderTool.h
 *
 *  @brief  Header file for the peak direction finder tool class.
  *
 *  $Log: $
 */
#ifndef LAR_PEAK_DIRECTION_FINDER_TOOL_H
#define LAR_PEAK_DIRECTION_FINDER_TOOL_H 1

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

namespace lar_content
{

class PeakDirectionFinderTool : public pandora::AlgorithmTool
{
public:
    PeakDirectionFinderTool();

    pandora::StatusCode Run(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D, const pandora::CaloHitList *const pViewHitList, 
        const pandora::HitType hitType, pandora::CartesianPointVector &peakDirectionVector);

private:
    typedef std::map<int, float> AngularDecompositionMap;

    void CollectHitsWithinROI(const pandora::CaloHitList &showerHitList, const pandora::CaloHitList *const pViewHitList, 
        const pandora::CartesianVector &nuVertex2D, pandora::CaloHitList &viewROIHits) const;

    void GetAngularExtrema(const pandora::CaloHitList &showerHitList, const pandora::CartesianVector &nuVertex2D, 
        float &lowestTheta, float &highestTheta) const;

    void CollectHitsWithinExtrema(const pandora::CaloHitList *const pViewHitList, const pandora::CartesianVector &nuVertex2D, 
        const float lowestTheta, const float highestTheta, pandora::CaloHitList &viewROIHits) const;

    void FillAngularDecompositionMap(const pandora::CaloHitList &viewShowerHitList, const pandora::CartesianVector &nuVertex2D, 
        AngularDecompositionMap &angularDecompositionMap) const;

    void SmoothAngularDecompositionMap(AngularDecompositionMap &angularDecompositionMap) const;

    void RetrievePeakDirections(const AngularDecompositionMap &angularDecompositionMap, pandora::CartesianPointVector &peakDirectionVector) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_pathwaySearchRegion;
    float m_theta0XZBinSize;
    int m_smoothingWindow;
    bool m_eventMode;
};

//------------------------------------------------------------------------------------------------------------------------------------------
} // namespace lar_content

#endif // #ifndef LAR_PEAK_DIRECTION_FINDER_TOOL_H
