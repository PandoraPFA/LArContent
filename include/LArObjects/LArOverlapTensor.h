/**
 *  @file   LArContent/include/LArObjects/LArOverlapTensor.h
 * 
 *  @brief  Header file for the lar overlap tensor class.
 * 
 *  $Log: $
 */
#ifndef LAR_OVERLAP_TENSOR_H
#define LAR_OVERLAP_TENSOR_H 1

#include "Pandora/PandoraInternal.h"

#include <map>
#include <vector>

namespace lar_content
{

/**
 *  @brief  OverlapTensor class
 */
template <typename T>
class OverlapTensor
{
public:
    typedef T OverlapResult;

    /**
     *  @brief  Element class
     */
    class Element
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pClusterU the address of the u cluster
         *  @param  pClusterV the address of the v cluster
         *  @param  pClusterW the address of the w cluster
         *  @param  overlapResult the overlap result
         */
        Element(pandora::Cluster *const pClusterU, pandora::Cluster *const pClusterV, pandora::Cluster *const pClusterW, const OverlapResult &overlapResult);

        /**
         *  @brief  Get the address of the u cluster
         * 
         *  @return address of the u cluster
         */
        pandora::Cluster *GetClusterU() const;

        /**
         *  @brief  Get the address of the v cluster
         * 
         *  @return address of the v cluster
         */
        pandora::Cluster *GetClusterV() const;

        /**
         *  @brief  Get the address of the w cluster
         * 
         *  @return address of the w cluster
         */
        pandora::Cluster *GetClusterW() const;

        /**
         *  @brief  Get the overlap result
         * 
         *  @return the overlap result
         */
        const OverlapResult &GetOverlapResult() const;

    private:
        pandora::Cluster   *m_pClusterU;                    ///< The address of the u cluster
        pandora::Cluster   *m_pClusterV;                    ///< The address of the v cluster
        pandora::Cluster   *m_pClusterW;                    ///< The address of the w cluster
        OverlapResult       m_overlapResult;                ///< The overlap result
    };

    typedef std::vector<Element> ElementList;

    /**
     *  @brief  Get unambiguous elements
     * 
     *  @param  ignoreUnavailable whether to ignore unavailable clusters
     *  @param  elementList to receive the unambiguous element list
     */
    void GetUnambiguousElements(const bool ignoreUnavailable, ElementList &elementList) const;

    /**
     *  @brief  Default ambiguity function, checking that only one U, V and W cluster is found
     * 
     *  @param  clusterListU cluster list U
     *  @param  clusterListV cluster list V
     *  @param  clusterListW cluster list W
     *  @param  pClusterU to receive the address of the unambiguous U cluster
     *  @param  pClusterV to receive the address of the unambiguous V cluster
     *  @param  pClusterW to receive the address of the unambiguous W cluster
     * 
     *  @return boolean
     */
    bool DefaultAmbiguityFunction(const pandora::ClusterList &clusterListU, const pandora::ClusterList &clusterListV, const pandora::ClusterList &clusterListW,
        pandora::Cluster *&pClusterU, pandora::Cluster *&pClusterV, pandora::Cluster *&pClusterW) const;

    /**
     *  @brief  Get the number of connections for a specified cluster
     * 
     *  @param  pCluster address of a cluster
     *  @param  ignoreUnavailable whether to ignore unavailable clusters
     *  @param  nU to receive the number of u connections
     *  @param  nV to receive the number of v connections
     *  @param  nW to receive the number of w connections
     */
    void GetNConnections(pandora::Cluster *const pCluster, const bool ignoreUnavailable, unsigned int &nU, unsigned int &nV, unsigned int &nW) const;

    /**
     *  @brief  Get a list of elements connected to a specified cluster
     * 
     *  @param  pCluster address of a cluster
     *  @param  ignoreUnavailable whether to ignore unavailable clusters
     *  @param  elementList to receive the connected element list
     */
    void GetConnectedElements(pandora::Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList) const;

    /**
     *  @brief  Get a list of elements connected to a specified cluster
     * 
     *  @param  pCluster address of a cluster
     *  @param  ignoreUnavailable whether to ignore unavailable clusters
     *  @param  elementList to receive the connected element list
     *  @param  nU to receive the number of u connections
     *  @param  nV to receive the number of v connections
     *  @param  nW to receive the number of w connections
     */
    void GetConnectedElements(pandora::Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList, unsigned int &nU,
        unsigned int &nV, unsigned int &nW) const;

    typedef std::map<pandora::Cluster*, pandora::ClusterList> ClusterNavigationMap;
    typedef std::map<pandora::Cluster*, OverlapResult> OverlapList;
    typedef std::map<pandora::Cluster*, OverlapList> OverlapMatrix;
    typedef std::map<pandora::Cluster*, OverlapMatrix> TheTensor;

    typedef typename TheTensor::const_iterator const_iterator;

    /**
     *  @brief  Returns an iterator referring to the first element in the overlap tensor
     */
    const_iterator begin() const;

    /**
     *  @brief  Returns an iterator referring to the past-the-end element in the overlap tensor
     */
    const_iterator end() const;

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

    /**
     *  @brief  Get the  overlap list for a specified pair of clusters
     * 
     *  @param  pClusterU address of cluster u
     *  @param  pClusterV address of cluster v
     * 
     *  @return the cluster overlap list
     */
    const OverlapList &GetOverlapList(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV) const;

    /**
     *  @brief  Get the cluster overlap matrix for a specified cluster
     * 
     *  @param  pClusterU address of cluster u
     * 
     *  @return the cluster overlap matrix
     */
    const OverlapMatrix &GetOverlapMatrix(pandora::Cluster *pClusterU) const;

    /**
     *  @brief  Get the cluster navigation map U->V
     * 
     *  @return the cluster navigation map U->V
     */
    const ClusterNavigationMap &GetClusterNavigationMapUV() const;

    /**
     *  @brief  Get the cluster navigation map V->W
     * 
     *  @return the cluster navigation map V->W
     */
    const ClusterNavigationMap &GetClusterNavigationMapVW() const;

    /**
     *  @brief  Get the cluster navigation map W->U
     * 
     *  @return the cluster navigation map W->U
     */
    const ClusterNavigationMap &GetClusterNavigationMapWU() const;

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
     *  @brief  SetReplace an existing overlap result
     * 
     *  @param  pClusterU address of cluster u
     *  @param  pClusterV address of cluster v
     *  @param  pClusterW address of cluster w
     *  @param  overlapResult the overlap result
     */
    void ReplaceOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW, const OverlapResult &overlapResult);

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
    /**
     *  @brief  Get elements connected to a specified cluster
     * 
     *  @param  pCluster address of the cluster
     *  @param  elementList the element list
     *  @param  clusterListU connected u clusters
     *  @param  clusterListV connected v clusters
     *  @param  clusterListW connected w clusters
     */
    void GetConnectedElements(pandora::Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList,
        pandora::ClusterList &clusterListU, pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW) const;

    /**
     *  @brief  Explore connections associated with a given cluster
     * 
     *  @param  pCluster address of the cluster
     *  @param  clusterListU connected u clusters
     *  @param  clusterListV connected v clusters
     *  @param  clusterListW connected w clusters
     */
    void ExploreConnections(pandora::Cluster *const pCluster, const bool ignoreUnavailable, pandora::ClusterList &clusterListU,
        pandora::ClusterList &clusterListV, pandora::ClusterList &clusterListW) const;

    TheTensor               m_overlapTensor;                ///< The overlap tensor
    ClusterNavigationMap    m_clusterNavigationMapUV;       ///< The cluster navigation map U->V
    ClusterNavigationMap    m_clusterNavigationMapVW;       ///< The cluster navigation map V->W
    ClusterNavigationMap    m_clusterNavigationMapWU;       ///< The cluster navigation map W->U
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void OverlapTensor<T>::GetNConnections(pandora::Cluster *const pCluster, const bool ignoreUnavailable, unsigned int &nU, unsigned int &nV,
    unsigned int &nW) const
{
    ElementList elementList;
    this->GetConnectedElements(pCluster, ignoreUnavailable, elementList, nU, nV, nW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void OverlapTensor<T>::GetConnectedElements(pandora::Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList) const
{
    unsigned int nU(0), nV(0), nW(0);
    this->GetConnectedElements(pCluster, ignoreUnavailable, elementList, nU, nV, nW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline typename OverlapTensor<T>::const_iterator OverlapTensor<T>::begin() const
{
    return m_overlapTensor.begin();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline typename OverlapTensor<T>::const_iterator OverlapTensor<T>::end() const
{
    return m_overlapTensor.end();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const typename OverlapTensor<T>::OverlapResult &OverlapTensor<T>::GetOverlapResult(
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
inline const typename OverlapTensor<T>::OverlapList &OverlapTensor<T>::GetOverlapList(
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
inline const typename OverlapTensor<T>::OverlapMatrix &OverlapTensor<T>::GetOverlapMatrix(
    pandora::Cluster *pClusterU) const
{
    typename TheTensor::const_iterator iter = m_overlapTensor.find(pClusterU);

    if (m_overlapTensor.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const typename OverlapTensor<T>::ClusterNavigationMap &OverlapTensor<T>::GetClusterNavigationMapUV() const
{
    return m_clusterNavigationMapUV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const typename OverlapTensor<T>::ClusterNavigationMap &OverlapTensor<T>::GetClusterNavigationMapVW() const
{
    return m_clusterNavigationMapVW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const typename OverlapTensor<T>::ClusterNavigationMap &OverlapTensor<T>::GetClusterNavigationMapWU() const
{
    return m_clusterNavigationMapWU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void OverlapTensor<T>::Clear()
{
    m_overlapTensor.clear();
    m_clusterNavigationMapUV.clear();
    m_clusterNavigationMapVW.clear();
    m_clusterNavigationMapWU.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
OverlapTensor<T>::Element::Element(pandora::Cluster *const pClusterU, pandora::Cluster *const pClusterV, pandora::Cluster *const pClusterW, const OverlapResult &overlapResult) :
    m_pClusterU(pClusterU),
    m_pClusterV(pClusterV),
    m_pClusterW(pClusterW),
    m_overlapResult(overlapResult)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
pandora::Cluster *OverlapTensor<T>::Element::GetClusterU() const
{
    return m_pClusterU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
pandora::Cluster *OverlapTensor<T>::Element::GetClusterV() const
{
    return m_pClusterV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
pandora::Cluster *OverlapTensor<T>::Element::GetClusterW() const
{
    return m_pClusterW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const typename OverlapTensor<T>::OverlapResult &OverlapTensor<T>::Element::GetOverlapResult() const
{
    return m_overlapResult;
}

} // namespace lar_content

#endif // #ifndef LAR_OVERLAP_TENSOR_H
