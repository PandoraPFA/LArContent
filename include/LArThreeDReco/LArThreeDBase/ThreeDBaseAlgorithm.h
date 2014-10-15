/**
 *  @file   LArContent/include/LArThreeDReco/LArThreeDBase/ThreeDBaseAlgorithm.h
 *
 *  @brief  Header file for the three dimension algorithm base class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_BASE_ALGORITHM_H
#define LAR_THREE_D_BASE_ALGORITHM_H 1

#include "Api/PandoraContentApi.h"

#include "Pandora/Algorithm.h"

#include "LArObjects/LArOverlapTensor.h"

namespace lar_content
{

/**
 *  @brief  ProtoParticle
 */
class ProtoParticle
{
public:
    pandora::ClusterList    m_clusterListU;                 ///< List of 2D U clusters in a 3D proto particle
    pandora::ClusterList    m_clusterListV;                 ///< List of 2D V clusters in a 3D proto particle
    pandora::ClusterList    m_clusterListW;                 ///< List of 2D W clusters in a 3D proto particle
};

typedef std::vector<ProtoParticle> ProtoParticleVector;
typedef std::map<pandora::Cluster*, pandora::ClusterList> ClusterMergeMap;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDBaseAlgorithm class
 */
template<typename T>
class ThreeDBaseAlgorithm : public pandora::Algorithm
{
public:
    typedef OverlapTensor<T> TensorType;

    /**
     *  @brief  Default constructor
     */
    ThreeDBaseAlgorithm();

    /**
     *  @brief  Destructor
     */
    virtual ~ThreeDBaseAlgorithm();

    /**
     *  @brief  Create particles using findings from recent algorithm processing
     *
     *  @param  protoParticleVector the proto particle vector
     *
     *  @param  whether particles were created
     */
    virtual bool CreateThreeDParticles(const ProtoParticleVector &protoParticleVector);

    /**
     *  @brief  Calculate Pfo properties from proto particle
     *
     *  @param  protoParticle the input proto particle
     *  @param  pfoParameters the output pfo parameters
     */
    virtual void SetPfoParameters(const ProtoParticle &protoParticle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const = 0;

    /**
     *  @brief  Merge clusters together
     *
     *  @param  clusterMergeMap the cluster merge map
     *
     *  @return whether changes to the tensor have been made
     */
    virtual bool MakeClusterMerges(const ClusterMergeMap &clusterMergeMap);

    /**
     *  @brief  Update to reflect a cluster merge
     *
     *  @param  pEnlargedCluster address of the enlarged cluster
     *  @param  pDeletedCluster address of the deleted cluster
     */
    virtual void UpdateUponMerge(pandora::Cluster *const pEnlargedCluster, pandora::Cluster *const pDeletedCluster);

    /**
     *  @brief  Update to reflect a cluster split
     *
     *  @param  pSplitCluster1 address of the first cluster fragment
     *  @param  pSplitCluster2 address of the second cluster fragment
     *  @param  pDeletedCluster address of the deleted cluster
     */
    virtual void UpdateUponSplit(pandora::Cluster *const pSplitCluster1, pandora::Cluster *const pSplitCluster2,
        pandora::Cluster *const pDeletedCluster);

    /**
     *  @brief  Update to reflect addition of a new cluster to the problem space
     *
     *  @param  pNewCluster address of the new cluster
     */
    virtual void UpdateForNewCluster(pandora::Cluster *const pNewCluster);

    /**
     *  @brief  Update to reflect cluster deletion
     *
     *  @param  pDeletedCluster address of the deleted cluster
     */
    virtual void UpdateUponDeletion(pandora::Cluster *const pDeletedCluster);

    /**
     *  @brief  Update tensor to remove all elements that have been added to pfos and so are unavailable
     */
    virtual void RemoveUnavailableTensorElements();

    /**
     *  @brief  Get the input u cluster list
     */
    const pandora::ClusterList &GetInputClusterListU() const;

    /**
     *  @brief  Get the input v cluster list
     */
    const pandora::ClusterList &GetInputClusterListV() const;

    /**
     *  @brief  Get the input w cluster list
     */
    const pandora::ClusterList &GetInputClusterListW() const;

    /**
     *  @brief  Get the selected u cluster list
     */
    const pandora::ClusterList &GetSelectedClusterListU() const;

    /**
     *  @brief  Get the selected v cluster list
     */
    const pandora::ClusterList &GetSelectedClusterListV() const;

    /**
     *  @brief  Get the selected w cluster list
     */
    const pandora::ClusterList &GetSelectedClusterListW() const;

    /**
     *  @brief  Get the name of the u cluster list
     */
    const std::string &GetClusterListNameU() const;

    /**
     *  @brief  Get the name of the v cluster list
     */
    const std::string &GetClusterListNameV() const;

    /**
     *  @brief  Get the name of the w cluster list
     */
    const std::string &GetClusterListNameW() const;

    /**
     *  @brief  Select a subset of input clusters for processing in this algorithm
     *
     *  @param  pInputClusterList address of an input cluster list
     *  @param  selectedClusterList to receive the selected cluster list
     */
    virtual void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const = 0;

protected:
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Select a subset of input clusters for processing in this algorithm
     */
    virtual void SelectAllInputClusters();

    /**
     *  @brief  Perform any preparatory steps required, e.g. caching expensive fit results for clusters
     */
    virtual void PreparationStep();

    /**
     *  @brief  Main loop over cluster combinations in order to populate the tensor. Responsible for calling CalculateOverlapResult.
     */
    virtual void PerformMainLoop();

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
     */
    virtual void ExamineTensor() = 0;

    /**
     *  @brief  Tidy member variables in derived class
     */
    virtual void TidyUp();

    const pandora::ClusterList *m_pInputClusterListU;           ///< Address of the input cluster list U
    const pandora::ClusterList *m_pInputClusterListV;           ///< Address of the input cluster list V
    const pandora::ClusterList *m_pInputClusterListW;           ///< Address of the input cluster list W

    pandora::ClusterList        m_clusterListU;                 ///< The selected modified cluster list U
    pandora::ClusterList        m_clusterListV;                 ///< The selected modified cluster list V
    pandora::ClusterList        m_clusterListW;                 ///< The selected modified cluster list W

    TensorType                  m_overlapTensor;                ///< The overlap tensor

private:
    pandora::StatusCode Run();

    std::string                 m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string                 m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string                 m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string                 m_outputPfoListName;            ///< The output pfo list name
};

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline const pandora::ClusterList &ThreeDBaseAlgorithm<T>::GetInputClusterListU() const
{
    if (NULL == m_pInputClusterListU)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return (*m_pInputClusterListU);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline const pandora::ClusterList &ThreeDBaseAlgorithm<T>::GetInputClusterListV() const
{
    if (NULL == m_pInputClusterListV)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return (*m_pInputClusterListV);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline const pandora::ClusterList &ThreeDBaseAlgorithm<T>::GetInputClusterListW() const
{
    if (NULL == m_pInputClusterListW)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return (*m_pInputClusterListW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline const pandora::ClusterList &ThreeDBaseAlgorithm<T>::GetSelectedClusterListU() const
{
    return m_clusterListU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline const pandora::ClusterList &ThreeDBaseAlgorithm<T>::GetSelectedClusterListV() const
{
    return m_clusterListV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline const pandora::ClusterList &ThreeDBaseAlgorithm<T>::GetSelectedClusterListW() const
{
    return m_clusterListW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline const std::string &ThreeDBaseAlgorithm<T>::GetClusterListNameU() const
{
    return m_inputClusterListNameU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline const std::string &ThreeDBaseAlgorithm<T>::GetClusterListNameV() const
{
    return m_inputClusterListNameV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline const std::string &ThreeDBaseAlgorithm<T>::GetClusterListNameW() const
{
    return m_inputClusterListNameW;
}

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_BASE_ALGORITHM_H
