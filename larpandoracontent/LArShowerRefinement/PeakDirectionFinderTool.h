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
    /**
     *  @brief  Default constructor
     */
    PeakDirectionFinderTool();

    pandora::StatusCode Run(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D,
        const pandora::CaloHitList *const pViewHitList, const pandora::HitType hitType, pandora::CartesianPointVector &peakDirectionVector);

private:
    typedef std::map<int, float> AngularDecompositionMap;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Collect the 2D hits within a region of interest 
     *          (m_ambiguousParticleMode ? hits not in the shower : hits within the initial shower cone [originating from the nu vertex])
     *
     *  @param  showerHitList the 2D shower hit list
     *  @param  pViewHitList the event 2D hits list
     *  @param  nuVertex2D the 2D neutrino vertex
     *  @param  viewROIHits the region of interest 2D hit list
     */
    void CollectHitsWithinROI(const pandora::CaloHitList &showerHitList, const pandora::CaloHitList *const pViewHitList,
        const pandora::CartesianVector &nuVertex2D, pandora::CaloHitList &viewROIHits) const;

    /**
     *  @brief  Determine the angle (from the +ve drift-axis) of the shower cone boundaries (originating from the nu vertex) 
     *
     *  @param  showerHitList the 2D shower hit list
     *  @param  nuVertex2D the 2D neutrino vertex
     *  @param  lowestTheta the lower angle boundary
     *  @param  highestTheta the higher angle boundary
     */
    void GetAngularExtrema(const pandora::CaloHitList &showerHitList, const pandora::CartesianVector &nuVertex2D, float &lowestTheta,
        float &highestTheta) const;

    /**
     *  @brief  Collect the hits that lie within the initial shower cone (originating from the nu vertex)
     *
     *  @param  pViewHitList the event 2D hits list  
     *  @param  nuVertex2D the 2D neutrino vertex
     *  @param  lowestTheta the lower angle (from the +ve drift-axis) boundary
     *  @param  highestTheta the higher angle (from the +ve drift-axis) boundary
     *  @param  viewROIHits the region of interest 2D hit list  
     */
    void CollectHitsWithinExtrema(const pandora::CaloHitList *const pViewHitList, const pandora::CartesianVector &nuVertex2D,
        const float lowestTheta, const float highestTheta, pandora::CaloHitList &viewROIHits) const;

    /**
     *  @brief  Determine the angular distribution of the ROI hits
     *
     *  @param  viewShowerHitList the 2D shower hit list
     *  @param  nuVertex2D the 2D neutrino vertex
     *  @param  angularDecompositionMap the output [angle from drift-axis -> weight] map
     */
    void FillAngularDecompositionMap(const pandora::CaloHitList &viewShowerHitList, const pandora::CartesianVector &nuVertex2D,
        AngularDecompositionMap &angularDecompositionMap) const;

    /**
     *  @brief  Smooth the ROI angular angular distribution
     *
     *  @param  angularDecompositionMap the [angle from drift-axis -> weight] map
     */
    void SmoothAngularDecompositionMap(AngularDecompositionMap &angularDecompositionMap) const;

    /**
     *  @brief  Obtain a vector of directions from the angular distribution peaks
     *
     *  @param  angularDecompositionMap the [angle from drift-axis -> weight] map
     *  @param  peakDirectionVector the output peak direction vector 
     */
    void RetrievePeakDirections(const AngularDecompositionMap &angularDecompositionMap, pandora::CartesianPointVector &peakDirectionVector) const;

    float m_pathwaySearchRegion; ///< The initial shower cone distance
    float m_theta0XZBinSize;     ///< The angular distribution bin size
    int m_smoothingWindow;       ///< On each side, the number of neighbouring bins with which each bin is averaged
    bool m_ambiguousParticleMode;  ///< Whether to find the initial pathway direction of the shower or of the other event particles
};

//------------------------------------------------------------------------------------------------------------------------------------------
} // namespace lar_content

#endif // #ifndef LAR_PEAK_DIRECTION_FINDER_TOOL_H
