/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/UndershootTracksTool.h
 * 
 *  @brief  Header file for the undershoot tracks tool class.
 * 
 *  $Log: $
 */
#ifndef UNDERSHOOT_TRACKS_TOOL_H
#define UNDERSHOOT_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDKinkBaseTool.h"

namespace lar
{

/**
 *  @brief  UndershootTracksTool class
 */
class UndershootTracksTool : public ThreeDKinkBaseTool
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm tool
     */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

    /**
     *  @brief  Default constructor
     */
    UndershootTracksTool();

private:
    /**
     *  @brief  Particle class
     */
    class Particle
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  elementA the tensor element A
         *  @param  elementB the tensor element B
         */
        Particle(const TensorType::Element &elementA, const TensorType::Element &elementB);

        pandora::Cluster           *m_pClusterA;            ///< Address of non-shared cluster in element A
        pandora::Cluster           *m_pClusterB;            ///< Address of non-shared cluster in element B
        pandora::Cluster           *m_pCommonCluster1;      ///< Address of the common cluster in view 1
        pandora::Cluster           *m_pCommonCluster2;      ///< Address of the common cluster in view 2
    };

    void GetIteratorListModifications(ThreeDTransverseTracksAlgorithm *pAlgorithm, const IteratorList &iteratorList, ModificationList &modificationList) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool    m_splitMode;                        ///< Whether to run in cluster splitting mode, as opposed to cluster merging mode
    float   m_minLongitudinalImpactParameter;   ///< The minimum longitudinal impact parameter for connecting broken clusters
    float   m_maxTransverseImpactParameter;     ///< The maximum transverse impact parameter for connecting broken clusters
    float   m_minImpactParameterCosTheta;       ///< The minimum cos theta (angle between vertex directions) for connecting broken clusters
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *UndershootTracksTool::Factory::CreateAlgorithmTool() const
{
    return new UndershootTracksTool();
}

} // namespace lar

#endif // #ifndef UNDERSHOOT_TRACKS_TOOL_H
