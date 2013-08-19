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

#include "LArObjects/LArOverlapTensor.h"

namespace lar
{

/**
 *  @brief  ThreeDBaseAlgorithm class
 */
template<typename T>
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
     *  @brief  Perform any preparatory steps required, e.g. caching expensive fit results for clusters
     */
    virtual void PreparationStep();

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
    ProtoParticleVector         m_protoParticleVector;          ///< The proto particle vector

    const pandora::ClusterList *m_pInputClusterListU;           ///< Address of the input cluster list U
    const pandora::ClusterList *m_pInputClusterListV;           ///< Address of the input cluster list V
    const pandora::ClusterList *m_pInputClusterListW;           ///< Address of the input cluster list W

    pandora::ClusterVector      m_clusterVectorU;               ///< The selected modified and sorted cluster vector U
    pandora::ClusterVector      m_clusterVectorV;               ///< The selected modified and sorted cluster vector V
    pandora::ClusterVector      m_clusterVectorW;               ///< The selected modified and sorted cluster vector W

    OverlapTensor<T>            m_overlapTensor;                ///< The overlap tensor

private:
    pandora::StatusCode Run();

    std::string                 m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string                 m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string                 m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string                 m_outputPfoListName;            ///< The output pfo list name
};

} // namespace lar

#endif // #ifndef LAR_THREE_D_BASE_ALGORITHM_H
