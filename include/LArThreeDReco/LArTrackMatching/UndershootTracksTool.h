/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/UndershootTracksTool.h
 * 
 *  @brief  Header file for the undershoot tracks tool class.
 * 
 *  $Log: $
 */
#ifndef UNDERSHOOT_TRACKS_TOOL_H
#define UNDERSHOOT_TRACKS_TOOL_H 1

#include "LArThreeDReco/LArTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar
{

/**
 *  @brief  UndershootTracksTool class
 */
class UndershootTracksTool : public TensorManipulationTool
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

private:
    /**
     *  @brief  ParticleComponent class
     */
    class ParticleComponent
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pClusterU the u cluster
         *  @param  pClusterV the v cluster
         *  @param  pClusterW the w cluster
         *  @param  trackOverlapResult the track overlap result
         */
        ParticleComponent(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW, const TrackOverlapResult &trackOverlapResult);

        /**
         *  @brief  Get the u cluster
         * 
         *  @return the u cluster
         */
        pandora::Cluster *GetClusterU() const;

        /**
         *  @brief  Get the v cluster
         * 
         *  @return the v cluster
         */
        pandora::Cluster *GetClusterV() const;

        /**
         *  @brief  Get the w cluster
         * 
         *  @return the w cluster
         */
        pandora::Cluster *GetClusterW() const;

        /**
         *  @brief  Get the track overlap result
         * 
         *  @return the track overlap result
         */
        const TrackOverlapResult &GetTrackOverlapResult() const;

    private:
        pandora::Cluster   *m_pClusterU;            ///< 
        pandora::Cluster   *m_pClusterV;            ///< 
        pandora::Cluster   *m_pClusterW;            ///< 
        TrackOverlapResult  m_trackOverlapResult;   ///< 
    };

    typedef std::vector<ParticleComponent> ParticleComponentList;

    bool Run(const SlidingFitResultMap &slidingFitResultMap, TrackOverlapTensor &overlapTensor, ProtoParticleVector &protoParticleVector);

    /**
     *  @brief  Build proto particle, starting with provided component and picking up any matched components in the overlap tensor
     * 
     *  @param  firstComponent the first particle component
     *  @param  protoParticle to receive the populated proto particle
     */
    void BuildProtoParticle(const ParticleComponent &firstComponent, ProtoParticle &protoParticle) const;

    /**
     *  @brief  Whether two particle components match, representing the same particle
     * 
     *  @param  firstComponent the first particle component
     *  @param  secondComponent the second particle component
     * 
     *  @return boolean
     */
    bool IsParticleMatch(const ParticleComponent &firstComponent, const ParticleComponent &secondComponent) const;

    /**
     *  @brief  Whether two clusters might match and represent the same particle
     * 
     *  @param  pFirstCluster the address of the first cluster
     *  @param  pSecondCluster the address of the second cluster
     * 
     *  @return boolean
     */
    bool IsPossibleMatch(pandora::Cluster *const pFirstCluster, pandora::Cluster *const pSecondCluster) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::AlgorithmTool *UndershootTracksTool::Factory::CreateAlgorithmTool() const
{
    return new UndershootTracksTool();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline UndershootTracksTool::ParticleComponent::ParticleComponent(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV,
        pandora::Cluster *pClusterW, const TrackOverlapResult &trackOverlapResult) :
    m_pClusterU(pClusterU),
    m_pClusterV(pClusterV),
    m_pClusterW(pClusterW),
    m_trackOverlapResult(trackOverlapResult)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *UndershootTracksTool::ParticleComponent::GetClusterU() const
{
    return m_pClusterU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *UndershootTracksTool::ParticleComponent::GetClusterV() const
{
    return m_pClusterV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Cluster *UndershootTracksTool::ParticleComponent::GetClusterW() const
{
    return m_pClusterW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TrackOverlapResult &UndershootTracksTool::ParticleComponent::GetTrackOverlapResult() const
{
    return m_trackOverlapResult;
}

} // namespace lar

#endif // #ifndef UNDERSHOOT_TRACKS_TOOL_H
