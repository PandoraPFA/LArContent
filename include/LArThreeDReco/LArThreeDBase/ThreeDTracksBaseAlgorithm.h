/**
 *  @file   LArContent/include/LArThreeDReco/LArThreeDBase/ThreeDTracksBaseAlgorithm.h
 * 
 *  @brief  Header file for the three dimensional tracks algorithm base class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_TRACKS_BASE_ALGORITHM_H
#define LAR_THREE_D_TRACKS_BASE_ALGORITHM_H 1

#include "LArObjects/LArTwoDSlidingFitResult.h"

#include "LArThreeDReco/LArThreeDBase/ThreeDBaseAlgorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDTransverseTracksAlgorithm class
 */
template<typename T>
class ThreeDTracksBaseAlgorithm : public ThreeDBaseAlgorithm<T>
{
public:
    /**
     *  @brief  Get a sliding fit result from the algorithm cache
     * 
     *  @param  pCluster address of the relevant cluster
     */
    const TwoDSlidingFitResult &GetCachedSlidingFitResult(pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Get the layer window for the sliding linear fits
     * 
     *  @return the layer window for the sliding linear fits
     */
    unsigned int GetSlidingFitWindow() const;

    virtual void UpdateForNewCluster(pandora::Cluster *const pNewCluster);
    virtual void UpdateUponDeletion(pandora::Cluster *const pDeletedCluster);
    virtual void PreparationStep();
    virtual void TidyUp();

protected:
    virtual pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int                m_slidingFitWindow;         ///< The layer window for the sliding linear fits
    TwoDSlidingFitResultMap     m_slidingFitResultMap;      ///< The sliding fit result map
};

//------------------------------------------------------------------------------------------------------------------------------------------

template<typename T>
inline unsigned int ThreeDTracksBaseAlgorithm<T>::GetSlidingFitWindow() const
{
    return m_slidingFitWindow;
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_TRACKS_BASE_ALGORITHM_H
