/**
 *  @file   LArContent/include/LArThreeDSeed/ThreeDBaseAlgorithm.h
 * 
 *  @brief  Header file for the three dimension algorithm base class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_BASE_ALGORITHM_H
#define LAR_THREE_D_BASE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDBaseAlgorithm class
 */
class ThreeDBaseAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDBaseAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~ThreeDBaseAlgorithm();

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Select a subset of input clusters for processing in this algorithm
     */
    virtual void SelectInputClusters();

    /**
     *  @brief  Initialize cluster overlap tensor in derived algorithm
     */
    virtual void InitializeTensor() = 0;

    /**
     *  @brief  Calculate cluster overlap result and store in tensor
     * 
     *  @param  pClusterU address of U view cluster
     *  @param  pClusterV address of V view cluster
     *  @param  pClusterW address of W view cluster
     */
    virtual void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW) = 0;

    /**
     *  @brief  Examine contents of tensor, collect together best-matching 2D particles and modify clusters as required
     * 
     *  @return whether to proceed with another round of tensor update and examination
     */
    virtual bool ExamineTensor() = 0;

    /**
     *  @brief  Create particles using findings from recent algorithm processing
     */
    virtual void CreateThreeDParticles();

    /**
     *  @brief  Update tensor to reflect changes made during recent processing
     */
    virtual void UpdateTensor();

    /**
     *  @brief  Tidy member variables in derived class
     */
    virtual void TidyUp();

    /**
     *  @brief  OverlapTensor class
     */
    template <typename T>
    class OverlapTensor
    {
    public:
        typedef T OverlapResult;

        /**
         *  @brief  Get the overlap result for a specified trio of clusters
         * 
         *  @param  pClusterU address of cluster u
         *  @param  pClusterV address of cluster v
         *  @param  pClusterW address of cluster w
         * 
         *  @return the address of the overlap result
         */
        const OverlapResult &GetOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW) const;

        typedef std::map<pandora::Cluster*, T> OverlapList;

        /**
         *  @brief  Get the  overlap list for a specified pair of clusters
         * 
         *  @param  pClusterU address of cluster u
         *  @param  pClusterV address of cluster v
         * 
         *  @return the cluster overlap list
         */
        const OverlapList &GetOverlapList(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV) const;

        typedef std::map<pandora::Cluster*, OverlapList> OverlapMatrix;

        /**
         *  @brief  Get the cluster overlap matrix for a specified cluster
         * 
         *  @param  pClusterU address of cluster u
         * 
         *  @return the cluster overlap matrix
         */
        const OverlapMatrix &GetOverlapMatrix(pandora::Cluster *pClusterU) const;

        /**
         *  @brief  Get the cluster list U
         * 
         *  @return the cluster list U
         */
        const pandora::ClusterList &GetClusterListU() const;

        /**
         *  @brief  Get the cluster list V
         * 
         *  @return the cluster list V
         */
        const pandora::ClusterList &GetClusterListV() const;

        /**
         *  @brief  Get the cluster list W
         * 
         *  @return the cluster list W
         */
        const pandora::ClusterList &GetClusterListW() const;

        /**
         *  @brief  Set overlap result
         * 
         *  @param  pClusterU address of cluster u
         *  @param  pClusterV address of cluster v
         *  @param  pClusterW address of cluster w
         *  @param  overlapResult the overlap result
         */
        void SetOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW, const OverlapResult &overlapResult);

        /**
         *  @brief  Remove entries from tensor corresponding to specified cluster
         * 
         *  @param  pCluster address of the cluster
         */
        void RemoveCluster(pandora::Cluster *pCluster);

        /**
         *  @brief  Clear overlap tensor
         */
        void Clear();

    private:
        typedef std::map<pandora::Cluster*, OverlapMatrix> TheTensor;

        TheTensor               m_overlapTensor;                ///< The overlap tensor
        pandora::ClusterList    m_clusterListU;                 ///< The cluster list U
        pandora::ClusterList    m_clusterListV;                 ///< The cluster list V
        pandora::ClusterList    m_clusterListW;                 ///< The cluster list W
    };

    /**
     *  @brief  ProtoParticle
     */
    class ProtoParticle
    {
    public:
        pandora::ClusterVector  m_clusterVectorU;               ///< List of 2D U clusters in a 3D proto particle
        pandora::ClusterVector  m_clusterVectorV;               ///< List of 2D V clusters in a 3D proto particle
        pandora::ClusterVector  m_clusterVectorW;               ///< List of 2D W clusters in a 3D proto particle
    };

    typedef std::vector<ProtoParticle> ProtoParticleVector;
    ProtoParticleVector         m_protoParticleVector;          ///< The proto particle vector

    const pandora::ClusterList *m_pInputClusterListU;           ///< Address of the input cluster list U
    const pandora::ClusterList *m_pInputClusterListV;           ///< Address of the input cluster list V
    const pandora::ClusterList *m_pInputClusterListW;           ///< Address of the input cluster list W

    pandora::ClusterVector      m_clusterVectorU;               ///< The selected modified and sorted cluster vector U
    pandora::ClusterVector      m_clusterVectorV;               ///< The selected modified and sorted cluster vector V
    pandora::ClusterVector      m_clusterVectorW;               ///< The selected modified and sorted cluster vector W

private:
    pandora::StatusCode Run();

    std::string                 m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string                 m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string                 m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string                 m_outputPfoListName;            ///< The output pfo list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const typename ThreeDBaseAlgorithm::OverlapTensor<T>::OverlapResult &ThreeDBaseAlgorithm::OverlapTensor<T>::GetOverlapResult(
    pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW) const
{
    const OverlapList &overlapList(this->GetOverlapList(pClusterU, pClusterV));
    typename OverlapList::const_iterator iter = overlapList.find(pClusterW);

    if (overlapList.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const typename ThreeDBaseAlgorithm::OverlapTensor<T>::OverlapList &ThreeDBaseAlgorithm::OverlapTensor<T>::GetOverlapList(
    pandora::Cluster *pClusterU, pandora::Cluster *pClusterV) const
{
    const OverlapMatrix &overlapMatrix(this->GetOverlapMatrix(pClusterU));
    typename OverlapMatrix::const_iterator iter = overlapMatrix.find(pClusterV);

    if (overlapMatrix.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const typename ThreeDBaseAlgorithm::OverlapTensor<T>::OverlapMatrix &ThreeDBaseAlgorithm::OverlapTensor<T>::GetOverlapMatrix(
    pandora::Cluster *pClusterU) const
{
    typename TheTensor::const_iterator iter = m_overlapTensor.find(pClusterU);

    if (m_overlapTensor.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const pandora::ClusterList &ThreeDBaseAlgorithm::OverlapTensor<T>::GetClusterListU() const
{
    return m_clusterListU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const pandora::ClusterList &ThreeDBaseAlgorithm::OverlapTensor<T>::GetClusterListV() const
{
    return m_clusterListV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const pandora::ClusterList &ThreeDBaseAlgorithm::OverlapTensor<T>::GetClusterListW() const
{
    return m_clusterListW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void ThreeDBaseAlgorithm::OverlapTensor<T>::SetOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV,
    pandora::Cluster *pClusterW, const OverlapResult &overlapResult)
{
    OverlapList &overlapList = m_overlapTensor[pClusterU][pClusterV];
    typename OverlapList::const_iterator iter = overlapList.find(pClusterW);

    if (overlapList.end() != iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!overlapList.insert(typename OverlapList::value_type(pClusterW, overlapResult)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    m_clusterListU.insert(pClusterU);
    m_clusterListV.insert(pClusterV);
    m_clusterListW.insert(pClusterW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void ThreeDBaseAlgorithm::OverlapTensor<T>::RemoveCluster(pandora::Cluster *pCluster)
{
    if (m_clusterListU.erase(pCluster) > 0)
    {
        typename TheTensor::iterator iter = m_overlapTensor.find(pCluster);

        if (m_overlapTensor.end() != iter)
            m_overlapTensor.erase(iter);
    }

    if (m_clusterListV.erase(pCluster) > 0)
    {
        for (typename TheTensor::iterator iterU = m_overlapTensor.begin(), iterUEnd = m_overlapTensor.end(); iterU != iterUEnd; ++iterU)
        {
            typename OverlapMatrix::iterator iter = iterU->second.find(pCluster);

            if (iterU->second.end() != iter)
                iterU->second.erase(iter);
        }
    }

    if (m_clusterListW.erase(pCluster) > 0)
    {
        for (typename TheTensor::iterator iterU = m_overlapTensor.begin(), iterUEnd = m_overlapTensor.end(); iterU != iterUEnd; ++iterU)
        {
            for (typename OverlapMatrix::iterator iterV = iterU->second.begin(), iterVEnd = iterU->second.end(); iterV != iterVEnd; ++iterV)
            {
                typename OverlapList::iterator iter = iterV->second.find(pCluster);

                if (iterV->second.end() != iter)
                    iterV->second.erase(iter);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void ThreeDBaseAlgorithm::OverlapTensor<T>::Clear()
{
    m_overlapTensor.clear();
    m_clusterListU.clear();
    m_clusterListV.clear();
    m_clusterListW.clear();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_BASE_ALGORITHM_H
