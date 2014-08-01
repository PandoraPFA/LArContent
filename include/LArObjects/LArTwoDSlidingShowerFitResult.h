/**
 *  @file   LArContent/include/LArObjects/LArTwoDSlidingShowerFitResult.h
 *
 *  @brief  Header file for the lar two dimensional sliding shower fit result class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_D_SLIDING_SHOWER_FIT_RESULT_H
#define LAR_TWO_D_SLIDING_SHOWER_FIT_RESULT_H 1

#include "Api/PandoraApi.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

namespace lar
{

/**
 *  @brief  ShowerEdge enum
 */
enum ShowerEdge
{
    POSITIVE_SHOWER_EDGE,
    NEGATIVE_SHOWER_EDGE
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SlidingShowerFitResult class
 */
class TwoDSlidingShowerFitResult
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pCluster address of the candidate shower cluster
     *  @param  slidingFitWindow the sliding fit window
     */
    TwoDSlidingShowerFitResult(const pandora::Cluster *const pCluster, const unsigned int slidingFitWindow);

    /**
     *  @brief  Get the sliding fit result for the full shower cluster
     * 
     *  @return the sliding fit result for the full shower cluster
     */
    const TwoDSlidingFitResult &GetShowerFitResult() const;

    /**
     *  @brief  Get the sliding fit result for the negative shower edge
     * 
     *  @return the sliding fit result for the negative shower edge
     */
    const TwoDSlidingFitResult &GetNegativeEdgeFitResult() const;

    /**
     *  @brief  Get the sliding fit result for the positive shower edge
     * 
     *  @return the sliding fit result for the positive shower edge
     */
    const TwoDSlidingFitResult &GetPositiveEdgeFitResult() const;

    /**
     *  @brief  Get the most appropriate shower edges at a given x coordinate
     * 
     *  @param  x the x coordinate
     *  @param  edgePositions to receive the list of intersections of the shower fit at the given x coordinate
     */
    void GetShowerEdges(const float x, pandora::FloatVector &edgePositions) const;

private:
    TwoDSlidingFitResult    m_showerFitResult;              ///< The sliding fit result for the full shower cluster
    TwoDSlidingFitResult    m_negativeEdgeFitResult;        ///< The sliding fit result for the negative shower edge
    TwoDSlidingFitResult    m_positiveEdgeFitResult;        ///< The sliding fit result for the positive shower edge
};

typedef std::map<pandora::Cluster*, TwoDSlidingShowerFitResult> TwoDSlidingShowerFitResultMap;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ShowerExtent
 */
class ShowerExtent
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  xCoordinate the x coordinate
     *  @param  highEdgeZ the shower high edge z coordinate
     *  @param  lowEdgeZ the shower low edge z coordinate
     */
    ShowerExtent(const float xCoordinate, const float highEdgeZ, const float lowEdgeZ);

    /**
     *  @param  Get the x coordinate
     * 
     *  @return the x coordinate
     */
    float GetXCoordinate() const;

    /**
     *  @param  Get the shower high edge z coordinate
     * 
     *  @return the shower high edge z coordinate
     */
    float GetHighEdgeZ() const;

    /**
     *  @param  Get the shower low edge z coordinate
     * 
     *  @return the shower low edge z coordinate
     */
    float GetLowEdgeZ() const;

private:
    float   m_xCoordinate;      ///< The x coordinate
    float   m_highEdgeZ;        ///< The shower high edge z coordinate
    float   m_lowEdgeZ;         ///< The shower low edge z coordinate
};

typedef std::map<int, ShowerExtent> ShowerPositionMap;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const TwoDSlidingFitResult &TwoDSlidingShowerFitResult::GetShowerFitResult() const
{
    return m_showerFitResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TwoDSlidingFitResult &TwoDSlidingShowerFitResult::GetNegativeEdgeFitResult() const
{
    return m_negativeEdgeFitResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const TwoDSlidingFitResult &TwoDSlidingShowerFitResult::GetPositiveEdgeFitResult() const
{
    return m_positiveEdgeFitResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ShowerExtent::ShowerExtent(const float xCoordinate, const float edge1, const float edge2) :
    m_xCoordinate(xCoordinate),
    m_highEdgeZ(std::max(edge1, edge2)),
    m_lowEdgeZ(std::min(edge1, edge2))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ShowerExtent::GetXCoordinate() const
{
    return m_xCoordinate;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ShowerExtent::GetHighEdgeZ() const
{
    return m_highEdgeZ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float ShowerExtent::GetLowEdgeZ() const
{
    return m_lowEdgeZ;
}

} // namespace lar

#endif // #ifndef LAR_TWO_D_SLIDING_SHOWER_FIT_RESULT_H
