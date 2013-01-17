/**
 *  @file   ThreeDParticleCreationAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional particle creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_PARTICLE_CREATION_ALGORITHM_H
#define LAR_THREE_D_PARTICLE_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDParticleCreationAlgorithm class
 */
class ThreeDParticleCreationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  ClusterOverlapInfo class
     */
    class ClusterOverlapInfo
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pCluster1 address of cluster u
         *  @param  pCluster2 address of cluster v
         *  @param  pCluster3 address of cluster w
         */
        ClusterOverlapInfo(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);

        /**
         *  @brief  Get chi2 calculated using cluster centroids u,v->w
         * 
         *  @return chi2 calculated using cluster centroids u,v->w
         */
        float GetChi2UVW() const;

        /**
         *  @brief  Get chi2 calculated using cluster centroids u,w->v
         * 
         *  @return chi2 calculated using cluster centroids u,w->v
         */
        float GetChi2UWV() const;

        /**
         *  @brief  Get chi2 calculated using cluster centroids v,w->u
         * 
         *  @return chi2 calculated using cluster centroids v,w->u
         */
        float GetChi2VWU() const;

        /**
         *  @brief  Get the x overlap
         * 
         *  @return the x overlap
         */
        float GetXOverlap() const;

    private:
        float           m_chi2UVW;                  ///< Chi2 calculated using cluster centroids u,v->w
        float           m_chi2UWV;                  ///< Chi2 calculated using cluster centroids u,w->v
        float           m_chi2VWU;                  ///< Chi2 calculated using cluster centroids v,w->u
        float           m_xOverlap;                 ///< The x overlap
    };

    typedef std::map<pandora::Cluster*, ClusterOverlapInfo> ClusterOverlapList;
    typedef std::map<pandora::Cluster*, ClusterOverlapList> ClusterOverlapMatrix;

    /**
     *  @brief  ClusterOverlapTensor class
     */
    class ClusterOverlapTensor
    {
    public:
        /**
         *  @brief  Get the cluster overlap matrix for a specified cluster
         * 
         *  @param  pClusterU address of cluster u
         * 
         *  @return the cluster overlap matrix
         */
        const ClusterOverlapMatrix &GetClusterOverlapMatrix(pandora::Cluster *pClusterU) const;

        /**
         *  @brief  Get the cluster overlap list for a specified pair of clusters
         * 
         *  @param  pClusterU address of cluster u
         *  @param  pClusterV address of cluster v
         * 
         *  @return the cluster overlap list
         */
        const ClusterOverlapList &GetClusterOverlapList(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV) const;

        /**
         *  @brief  Get the cluster overlap info for a specified trio of clusters
         * 
         *  @param  pClusterU address of cluster u
         *  @param  pClusterV address of cluster v
         *  @param  pClusterW address of cluster w
         * 
         *  @return the cluster overlap info
         */
        const ClusterOverlapInfo &GetClusterOverlapInfo(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW) const;

        /**
         *  @brief  Calculate and set cluster overlap info
         * 
         *  @param  pClusterU address of cluster u
         *  @param  pClusterV address of cluster v
         *  @param  pClusterW address of cluster w
         */
        void SetClusterOverlapInfo(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);

    private:
        typedef std::map<pandora::Cluster*, ClusterOverlapMatrix> OverlapTensor;

        OverlapTensor   m_overlapTensor;                ///< The cluster overlap tensor
    };

    std::string         m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string         m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string         m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string         m_outputPfoListName;            ///< The output pfo list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDParticleCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDParticleCreationAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDParticleCreationAlgorithm::ClusterOverlapInfo::GetChi2UVW() const
{
    return m_chi2UVW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDParticleCreationAlgorithm::ClusterOverlapInfo::GetChi2UWV() const
{
    return m_chi2UWV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDParticleCreationAlgorithm::ClusterOverlapInfo::GetChi2VWU() const
{
    return m_chi2VWU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ThreeDParticleCreationAlgorithm::ClusterOverlapInfo::GetXOverlap() const
{
    return m_xOverlap;
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_PARTICLE_CREATION_ALGORITHM_H
