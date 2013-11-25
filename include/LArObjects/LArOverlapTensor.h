/**
 *  @file   LArContent/include/LArObjects/LArOverlapTensor.h
 * 
 *  @brief  Header file for the lar overlap tensor class.
 * 
 *  $Log: $
 */
#ifndef LAR_OVERLAP_TENSOR_H
#define LAR_OVERLAP_TENSOR_H 1

namespace lar
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

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TrackOverlapResult class
 */
class TrackOverlapResult
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackOverlapResult();

    /**
     *  @brief  Constructor
     * 
     *  @param  nMatchedSamplingPoints
     *  @param  nSamplingPoints
     *  @param  chi2
     */
    TrackOverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints, const float chi2);

    /**
     *  @brief  Copy constructor
     * 
     *  @param  rhs
     */
    TrackOverlapResult(const TrackOverlapResult &rhs);

    /**
     *  @brief  Whether the track overlap result has been initialized
     *
     *  @return boolean
     */
    bool IsInitialized() const;

    /**
     *  @brief  Get the number of matched sampling points
     *
     *  @return the number of matched sampling points
     */
    unsigned int GetNMatchedSamplingPoints() const;

    /**
     *  @brief  Get the number of sampling points
     *
     *  @return the number of sampling points
     */
    unsigned int GetNSamplingPoints() const;

    /**
     *  @brief  Get the fraction of sampling points resulting in a match
     *
     *  @return the fraction of sampling points resulting in a match
     */
    float GetMatchedFraction() const;

    /**
     *  @brief  Get the absolute chi2 value
     *
     *  @return the absolute chi2 value
     */
    float GetChi2() const;

    /**
     *  @brief  Get the chi2 per samping point value
     *
     *  @return the chi2 per samping point value
     */
    float GetReducedChi2() const;

    /**
     *  @brief  Track overlap result less than operator
     * 
     *  @param  rhs the track overlap result for comparison
     */
    bool operator<(const TrackOverlapResult &rhs) const;

    /**
     *  @brief  Track overlap result greater than operator
     * 
     *  @param  rhs the track overlap result for comparison
     */
    bool operator>(const TrackOverlapResult &rhs) const;

    /**
     *  @brief  Track overlap result assigment operator
     * 
     *  @param  rhs the track overlap result to assign
     */
    TrackOverlapResult &operator=(const TrackOverlapResult &rhs);

private:
    bool            m_isInitialized;                ///< Whether the track overlap result has been initialized
    unsigned int    m_nMatchedSamplingPoints;       ///< The number of matched sampling points
    unsigned int    m_nSamplingPoints;              ///< The number of sampling points
    float           m_matchedFraction;              ///< The fraction of sampling points resulting in a match
    float           m_chi2;                         ///< The absolute chi2 value
    float           m_reducedChi2;                  ///< The chi2 per samping point value
};

/**
 *  @brief  Track overlap result + operator
 * 
 *  @param  lhs the first track overlap result to add
 *  @param  rhs the second track overlap result to add
 */
TrackOverlapResult operator+(const TrackOverlapResult &lhs, const TrackOverlapResult &rhs);

//------------------------------------------------------------------------------------------------------------------------------------------
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
inline const pandora::ClusterList &OverlapTensor<T>::GetClusterListU() const
{
    return m_clusterListU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const pandora::ClusterList &OverlapTensor<T>::GetClusterListV() const
{
    return m_clusterListV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const pandora::ClusterList &OverlapTensor<T>::GetClusterListW() const
{
    return m_clusterListW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void OverlapTensor<T>::SetOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV,
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
inline void OverlapTensor<T>::RemoveCluster(pandora::Cluster *pCluster)
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
inline void OverlapTensor<T>::Clear()
{
    m_overlapTensor.clear();
    m_clusterListU.clear();
    m_clusterListV.clear();
    m_clusterListW.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackOverlapResult::TrackOverlapResult() :
    m_isInitialized(false),
    m_nMatchedSamplingPoints(0),
    m_nSamplingPoints(0),
    m_matchedFraction(0.f),
    m_chi2(0.f),
    m_reducedChi2(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackOverlapResult::TrackOverlapResult(const unsigned int nMatchedSamplingPoints, const unsigned int nSamplingPoints, const float chi2) :
    m_isInitialized(true),
    m_nMatchedSamplingPoints(nMatchedSamplingPoints),
    m_nSamplingPoints(nSamplingPoints),
    m_matchedFraction(0.f),
    m_chi2(chi2),
    m_reducedChi2(0.f)
{
    if ((0 == m_nSamplingPoints) || (m_nMatchedSamplingPoints > m_nSamplingPoints))
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    m_matchedFraction = static_cast<float>(m_nMatchedSamplingPoints) / static_cast<float>(m_nSamplingPoints);
    m_reducedChi2 = m_chi2 / static_cast<float>(m_nSamplingPoints);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackOverlapResult::TrackOverlapResult(const TrackOverlapResult &rhs) :
    m_isInitialized(rhs.m_isInitialized),
    m_nMatchedSamplingPoints(rhs.m_nMatchedSamplingPoints),
    m_nSamplingPoints(rhs.m_nSamplingPoints),
    m_matchedFraction(rhs.m_matchedFraction),
    m_chi2(rhs.m_chi2),
    m_reducedChi2(rhs.m_reducedChi2)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackOverlapResult::IsInitialized() const
{
    return m_isInitialized;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int TrackOverlapResult::GetNMatchedSamplingPoints() const
{
    return m_nMatchedSamplingPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int TrackOverlapResult::GetNSamplingPoints() const
{
    return m_nSamplingPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackOverlapResult::GetMatchedFraction() const
{
    return m_matchedFraction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackOverlapResult::GetChi2() const
{
    return m_chi2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float TrackOverlapResult::GetReducedChi2() const
{
    return m_reducedChi2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackOverlapResult::operator<(const TrackOverlapResult &rhs) const
{
    if (m_nMatchedSamplingPoints != rhs.m_nMatchedSamplingPoints)
        return (m_nMatchedSamplingPoints < rhs.m_nMatchedSamplingPoints);

    return (m_reducedChi2 > rhs.m_reducedChi2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool TrackOverlapResult::operator>(const TrackOverlapResult &rhs) const
{
    if (this == &rhs)
        return false;

    return !(*this < rhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackOverlapResult &TrackOverlapResult::operator=(const TrackOverlapResult &rhs)
{
    if (this == &rhs)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    m_isInitialized = rhs.m_isInitialized;
    m_nMatchedSamplingPoints = rhs.m_nMatchedSamplingPoints;
    m_nSamplingPoints = rhs.m_nSamplingPoints;
    m_matchedFraction = rhs.m_matchedFraction;
    m_chi2 = rhs.m_chi2;
    m_reducedChi2 = rhs.m_reducedChi2;

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline TrackOverlapResult operator+(const TrackOverlapResult &lhs, const TrackOverlapResult &rhs)
{
    if (!lhs.IsInitialized() && !rhs.IsInitialized())
        return TrackOverlapResult();

    return TrackOverlapResult(lhs.GetNMatchedSamplingPoints() + rhs.GetNMatchedSamplingPoints(),
        lhs.GetNSamplingPoints() + rhs.GetNSamplingPoints(), lhs.GetChi2() + rhs.GetChi2());
}

} // namespace lar

#endif // #ifndef LAR_OVERLAP_TENSOR_H
