/**
 *  @file   larpandoracontent/LArObjects/LArOverlapMatrix.h
 *
 *  @brief  Header file for the lar overlap matrix class.
 *
 *  $Log: $
 */
#ifndef LAR_OVERLAP_MATRIX_H
#define LAR_OVERLAP_MATRIX_H 1

#include "Pandora/PandoraInternal.h"

#include <unordered_map>
#include <vector>

namespace lar_content
{

/**
 *  @brief  OverlapMatrix class
 */
template <typename T>
class OverlapMatrix
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
         *  @param  pCluster1 the address of cluster 1
         *  @param  pCluster2 the address of cluster 2
         *  @param  overlapResult the overlap result
         */
        Element(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const OverlapResult &overlapResult);

        /**
         *  @brief  Get the address of cluster 1
         *
         *  @return address of cluster 1
         */
        const pandora::Cluster *GetCluster1() const;

        /**
         *  @brief  Get the address of cluster 2
         *
         *  @return address of cluster 2
         */
        const pandora::Cluster *GetCluster2() const;

        /**
         *  @brief  Get the overlap result
         *
         *  @return the overlap result
         */
        const OverlapResult &GetOverlapResult() const;

        /**
         *  @brief  Element less than operator
         *
         *  @param  rhs the element for comparison
         */
        bool operator<(const Element &rhs) const;

    private:
        const pandora::Cluster *m_pCluster1; ///< The address of cluster 1
        const pandora::Cluster *m_pCluster2; ///< The address of cluster 2
        OverlapResult m_overlapResult;       ///< The overlap result
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
     *  @brief  Default ambiguity function, checking that only one cluster from view 1 and view 2 is found
     *
     *  @param  clusterList1 cluster list 1
     *  @param  clusterList2 cluster list 2
     *  @param  pCluster1 to receive the address of the unambiguous cluster 1
     *  @param  pCluster2 to receive the address of the unambiguous cluster 2
     *
     *  @return boolean
     */
    bool DefaultAmbiguityFunction(const pandora::ClusterList &clusterList1, const pandora::ClusterList &clusterList2,
        const pandora::Cluster *&pCluster1, const pandora::Cluster *&pCluster2) const;

    /**
     *  @brief  Get the number of connections for a specified cluster
     *
     *  @param  pCluster address of a cluster
     *  @param  ignoreUnavailable whether to ignore unavailable clusters
     *  @param  n1 to receive the number of view 1 connections
     *  @param  n2 to receive the number of view 2 connections
     */
    void GetNConnections(const pandora::Cluster *const pCluster, const bool ignoreUnavailable, unsigned int &n1, unsigned int &n2) const;

    /**
     *  @brief  Get a list of elements connected to a specified cluster
     *
     *  @param  pCluster address of a cluster
     *  @param  ignoreUnavailable whether to ignore unavailable clusters
     *  @param  elementList to receive the connected element list
     */
    void GetConnectedElements(const pandora::Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList) const;

    /**
     *  @brief  Get a list of elements connected to a specified cluster
     *
     *  @param  pCluster address of a cluster
     *  @param  ignoreUnavailable whether to ignore unavailable clusters
     *  @param  elementList to receive the connected element list
     *  @param  n1 to receive the number of view 1 connections
     *  @param  n2 to receive the number of view 2 connections
     */
    void GetConnectedElements(const pandora::Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList,
        unsigned int &n1, unsigned int &n2) const;

    typedef std::unordered_map<const pandora::Cluster *, pandora::ClusterList> ClusterNavigationMap;
    typedef std::unordered_map<const pandora::Cluster *, OverlapResult> OverlapList;
    typedef std::unordered_map<const pandora::Cluster *, OverlapList> TheMatrix;

    typedef typename TheMatrix::const_iterator const_iterator;

    /**
     *  @brief  Returns an iterator referring to the first element in the overlap matrix
     */
    const_iterator begin() const;

    /**
     *  @brief  Returns an iterator referring to the past-the-end element in the overlap matrix
     */
    const_iterator end() const;

    /**
     *  @brief  Get a sorted vector of key clusters (view 1 clusters with current implementation)
     *
     *  @param sortedKeyClusters to receive the sorted vector of key clusters
     */
    void GetSortedKeyClusters(pandora::ClusterVector &sortedKeyClusters) const;

    /**
     *  @brief  Get the overlap result for a specified pair of clusters
     *
     *  @param  pCluster1 address of cluster 1
     *  @param  pCluster2 address of cluster 2
     *
     *  @return the address of the overlap result
     */
    const OverlapResult &GetOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const;

    /**
     *  @brief  Get the  overlap list for a specified cluster
     *
     *  @param  pCluster1 address of cluster 1
     *
     *  @return the cluster overlap list
     */
    const OverlapList &GetOverlapList(const pandora::Cluster *const pCluster1) const;

    /**
     *  @brief  Get the cluster navigation map 1->2
     *
     *  @return the cluster navigation map 1->2
     */
    const ClusterNavigationMap &GetClusterNavigationMap12() const;

    /**
     *  @brief  Get the cluster navigation map 2->1
     *
     *  @return the cluster navigation map 2->1
     */
    const ClusterNavigationMap &GetClusterNavigationMap21() const;

    /**
     *  @brief  Set overlap result
     *
     *  @param  pCluster1 address of cluster 1
     *  @param  pCluster2 address of cluster 2
     *  @param  overlapResult the overlap result
     */
    void SetOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const OverlapResult &overlapResult);

    /**
     *  @brief  SetReplace an existing overlap result
     *
     *  @param  pCluster1 address of cluster 1
     *  @param  pCluster2 address of cluster 2
     *  @param  overlapResult the overlap result
     */
    void ReplaceOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const OverlapResult &overlapResult);

    /**
     *  @brief  Remove entries from matrix corresponding to specified cluster
     *
     *  @param  pCluster address of the cluster
     */
    void RemoveCluster(const pandora::Cluster *const pCluster);

    /**
     *  @brief  Clear overlap matrix
     */
    void Clear();

private:
    /**
     *  @brief  Get elements connected to a specified cluster
     *
     *  @param  pCluster address of the cluster
     *  @param  elementList the element list
     *  @param  clusterList1 connected view 1 clusters
     *  @param  clusterList2 connected view 2 clusters
     */
    void GetConnectedElements(const pandora::Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList,
        pandora::ClusterList &clusterList1, pandora::ClusterList &clusterList2) const;

    /**
     *  @brief  Explore connections associated with a given cluster
     *
     *  @param  pCluster address of the cluster
     *  @param  clusterList1 connected view 1 clusters
     *  @param  clusterList2 connected view 2 clusters
     */
    void ExploreConnections(const pandora::Cluster *const pCluster, const bool ignoreUnavailable, pandora::ClusterList &clusterList1,
        pandora::ClusterList &clusterList2) const;

    TheMatrix m_overlapMatrix;                     ///< The overlap matrix
    ClusterNavigationMap m_clusterNavigationMap12; ///< The cluster navigation map 1->2
    ClusterNavigationMap m_clusterNavigationMap21; ///< The cluster navigation map 2->1
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void OverlapMatrix<T>::GetNConnections(const pandora::Cluster *const pCluster, const bool ignoreUnavailable, unsigned int &n1, unsigned int &n2) const
{
    ElementList elementList;
    this->GetConnectedElements(pCluster, ignoreUnavailable, elementList, n1, n2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void OverlapMatrix<T>::GetConnectedElements(const pandora::Cluster *const pCluster, const bool ignoreUnavailable, ElementList &elementList) const
{
    unsigned int n1(0), n2(0);
    this->GetConnectedElements(pCluster, ignoreUnavailable, elementList, n1, n2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline typename OverlapMatrix<T>::const_iterator OverlapMatrix<T>::begin() const
{
    return m_overlapMatrix.begin();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline typename OverlapMatrix<T>::const_iterator OverlapMatrix<T>::end() const
{
    return m_overlapMatrix.end();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const typename OverlapMatrix<T>::OverlapResult &OverlapMatrix<T>::GetOverlapResult(
    const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2) const
{
    const OverlapList &overlapList(this->GetOverlapList(pCluster1));
    typename OverlapList::const_iterator iter = overlapList.find(pCluster2);

    if (overlapList.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const typename OverlapMatrix<T>::OverlapList &OverlapMatrix<T>::GetOverlapList(const pandora::Cluster *const pCluster1) const
{
    typename TheMatrix::const_iterator iter = m_overlapMatrix.find(pCluster1);

    if (m_overlapMatrix.end() == iter)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const typename OverlapMatrix<T>::ClusterNavigationMap &OverlapMatrix<T>::GetClusterNavigationMap12() const
{
    return m_clusterNavigationMap12;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const typename OverlapMatrix<T>::ClusterNavigationMap &OverlapMatrix<T>::GetClusterNavigationMap21() const
{
    return m_clusterNavigationMap21;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void OverlapMatrix<T>::Clear()
{
    m_overlapMatrix.clear();
    m_clusterNavigationMap12.clear();
    m_clusterNavigationMap21.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline OverlapMatrix<T>::Element::Element(
    const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const OverlapResult &overlapResult) :
    m_pCluster1(pCluster1), m_pCluster2(pCluster2), m_overlapResult(overlapResult)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const pandora::Cluster *OverlapMatrix<T>::Element::GetCluster1() const
{
    return m_pCluster1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const pandora::Cluster *OverlapMatrix<T>::Element::GetCluster2() const
{
    return m_pCluster2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const typename OverlapMatrix<T>::OverlapResult &OverlapMatrix<T>::Element::GetOverlapResult() const
{
    return m_overlapResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool OverlapMatrix<T>::Element::operator<(const Element &rhs) const
{
    if (this == &rhs)
        return false;

    return (this->GetOverlapResult() < rhs.GetOverlapResult());
}

} // namespace lar_content

#endif // #ifndef LAR_OVERLAP_MATRIX_H
