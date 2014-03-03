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
     *  @param  elementList to receive the unambiguous element list
     */
    void GetUnambiguousElements(ElementList &elementList) const;

    /**
     *  @brief  Get the number of connections for a specified cluster
     * 
     *  @param  pCluster address of a cluster
     *  @param  nU to receive the number of u connections
     *  @param  nV to receive the number of v connections
     *  @param  nW to receive the number of w connections
     */
    void GetNConnections(const pandora::Cluster *const pCluster, unsigned int &nU, unsigned int &nV, unsigned int &nW) const;

    /**
     *  @brief  Get a list of elements connected to a specified cluster
     * 
     *  @param  pCluster address of a cluster
     *  @param  elementList to receive the connected element list
     */
    void GetConnectedElements(const pandora::Cluster *const pCluster, ElementList &elementList) const;

    /**
     *  @brief  Get a list of elements connected to a specified cluster
     * 
     *  @param  pCluster address of a cluster
     *  @param  elementList to receive the connected element list
     *  @param  nU to receive the number of u connections
     *  @param  nV to receive the number of v connections
     *  @param  nW to receive the number of w connections
     */
    void GetConnectedElements(const pandora::Cluster *const pCluster, ElementList &elementList, unsigned int &nU, unsigned int &nV, unsigned int &nW) const;

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
    TheTensor               m_overlapTensor;                ///< The overlap tensor
    ClusterNavigationMap    m_clusterNavigationMapUV;       ///< The cluster navigation map U->V
    ClusterNavigationMap    m_clusterNavigationMapVW;       ///< The cluster navigation map V->W
    ClusterNavigationMap    m_clusterNavigationMapWU;       ///< The cluster navigation map W->U
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

typedef std::vector<TrackOverlapResult> TrackOverlapResultVector;

//------------------------------------------------------------------------------------------------------------------------------------------

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
inline void OverlapTensor<T>::GetNConnections(const pandora::Cluster *const pCluster, unsigned int &nU, unsigned int &nV, unsigned int &nW) const
{
    ElementList elementList;
    this->GetConnectedElements(pCluster, elementList, nU, nV, nW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline void OverlapTensor<T>::GetConnectedElements(const pandora::Cluster *const pCluster, ElementList &elementList) const
{
    unsigned int nU(0), nV(0), nW(0);
    this->GetConnectedElements(pCluster, elementList, nU, nV, nW);
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
