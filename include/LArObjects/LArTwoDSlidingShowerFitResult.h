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

namespace lar_content
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
     *  @param  slidingFitLayerPitch the sliding fit z pitch, units cm
     *  @param  showerEdgeMultiplier artificially tune width of shower envelope so as to make it more/less inclusive
     */
    TwoDSlidingShowerFitResult(const pandora::Cluster *const pCluster, const unsigned int slidingFitWindow, const float slidingFitLayerPitch,
        const float showerEdgeMultiplier = 1.f);

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
     *  @param  widenIfAmbiguity whether to widen the shower edges in cases of ambiguities (i.e. be generous)
     *  @param  edgePositions to receive the list of intersections of the shower fit at the given x coordinate
     */
    void GetShowerEdges(const float x, const bool widenIfAmbiguity, pandora::FloatVector &edgePositions) const;

private:
    /**
     *  @brief  Perform two dimensional sliding fit to shower edge, using specified primary axis
     *
     *  @param  fullShowerFit the result of fitting the full shower
     *  @param  showerEdge the shower edge
     *  @param  showerEdgeMultiplier artificially tune width of shower envelope so as to make it more/less inclusive
     * 
     *  @return the shower edge fit result
     */
    static TwoDSlidingFitResult LArTwoDShowerEdgeFit(const TwoDSlidingFitResult &fullShowerFit, const ShowerEdge showerEdge,
        const float showerEdgeMultiplier);

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

} // namespace lar_content

#endif // #ifndef LAR_TWO_D_SLIDING_SHOWER_FIT_RESULT_H
