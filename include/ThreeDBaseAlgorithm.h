/**
 *  @file   ThreeDBaseAlgorithm.h
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
    /**
     *  @brief  Select a subset of input clusters for processing in this algorithm
     */
    virtual void SelectInputClusters() = 0;

    /**
     *  @brief  Modify selected input clusters as desired, maybe breaking up track-like sections
     */
    virtual void ModifyInputClusters() = 0;

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
     *  @brief  Update tensor to reflect changes made during recent processing
     */
    virtual void UpdateTensor() = 0;

    /**
     *  @brief  Tidy member variables in derived class
     */
    virtual void TidyUp() = 0;

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
        const OverlapResult *GetOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW) const;

        typedef std::map<pandora::Cluster*, T*> OverlapList;

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
         *  @brief  Set overlap result
         * 
         *  @param  pClusterU address of cluster u
         *  @param  pClusterV address of cluster v
         *  @param  pClusterW address of cluster w
         *  @param  pOverlapResult address of the overlap result
         */
        void SetOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW, OverlapResult *pOverlapResult);

    private:
        typedef std::map<pandora::Cluster*, OverlapMatrix> TheTensor;

        TheTensor               m_overlapTensor;                ///< The overlap tensor
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
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Create particles using findings from recent algorithm processing
     */
    void CreateThreeDParticles();

    /**
     *  @brief  Tidy member variables in base class
     */
    void TidyMemberVariables();

    std::string                 m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string                 m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string                 m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string                 m_outputPfoListName;            ///< The output pfo list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const typename ThreeDBaseAlgorithm::OverlapTensor<T>::OverlapResult *ThreeDBaseAlgorithm::OverlapTensor<T>::GetOverlapResult(
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
const typename ThreeDBaseAlgorithm::OverlapTensor<T>::OverlapList &ThreeDBaseAlgorithm::OverlapTensor<T>::GetOverlapList(
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
const typename ThreeDBaseAlgorithm::OverlapTensor<T>::OverlapMatrix &ThreeDBaseAlgorithm::OverlapTensor<T>::GetOverlapMatrix(
    pandora::Cluster *pClusterU) const
{
    typename TheTensor::const_iterator iter = m_overlapTensor.find(pClusterU);

    if (m_overlapTensor.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void ThreeDBaseAlgorithm::OverlapTensor<T>::SetOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV,
    pandora::Cluster *pClusterW, OverlapResult *pOverlapResult)
{
    OverlapList &overlapList = m_overlapTensor[pClusterU][pClusterV];
    typename OverlapList::const_iterator iter = overlapList.find(pClusterW);

    if (overlapList.end() != iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_ALREADY_PRESENT);

    if (!overlapList.insert(OverlapList::value_type(pClusterW, pOverlapResult)).second)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_BASE_ALGORITHM_H
