/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h
 *
 *  @brief  Header file for the HierarchyPfo class.
 *
 *  $Log: $
 */
#ifndef LAR_HIERARCHY_PFO_H
#define LAR_HIERARCHY_PFO_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

using namespace lar_content;

namespace lar_dl_content
{

/**
  *  @brief  ExtremalPoint class
  */
class ExtremalPoint
{
public:
    /**
     *  @brief  Default constructor
     */
    ExtremalPoint();

    /**
     *  @brief  Get the position
     *
     *  @return the position
     */
    const pandora::CartesianVector &GetPosition() const;

    /**
     *  @brief  Get the direction at the extremal point
     *
     *  @return the direction
     */
    const pandora::CartesianVector &GetDirection() const;

    /**
     *  @brief  Set the the extremal point's position and direction
     *
     *  @param  position the position
     *  @param  direction the direction
     */
    void Set(const pandora::CartesianVector &position, const pandora::CartesianVector &direction);

    /**
     *  @brief  Whether extremal point object has been set
     *
     *  @return whether the extremal point object has been set
     */
    bool IsSet() const;

    /**
     *  @brief  Assignment operator
     */
    ExtremalPoint &operator=(const ExtremalPoint &rhs);

private:
    bool m_isSet;                         ///< whether the extremal point object has been set
    pandora::CartesianVector m_position;  ///< the extremal point position
    pandora::CartesianVector m_direction; ///< the extremal point direction (pointing into the particle)
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
  *  @brief  HierarchyPfo class
  */
class HierarchyPfo
{
public:
    /**
     *  @brief  Default constructor
     */
    HierarchyPfo();

    /**
     *  @brief  Constructor
     *
     *  @param  pPfo pointer to the input pfo
     *  @param  threeDSlidingFitResult 3D sliding fit of the input pfo
     *  @param  upstreamPoint the extremal point closest to the neutrino vertex
     *  @param  downstreamPoint the extremal point furthest from the neutrino vertex
     */
    HierarchyPfo(const pandora::ParticleFlowObject *pPfo, const ThreeDSlidingFitResult &threeDSlidingFitResult,
        const ExtremalPoint &upstreamPoint, const ExtremalPoint &downstreamPoint);

    /**
     *  @brief  Get the pfo
     *
     *  @return a pointer to the pfo
     */
    const pandora::ParticleFlowObject *GetPfo() const;

    /**
     *  @brief  Set the pfo
     *
     *  @param  pPfo a pointer to the pfo
     */
    void SetPfo(const pandora::ParticleFlowObject *pPfo);

    /**
     *  @brief  Get the pfo's 3D sliding fit result
     *
     *  @return the 3D sliding fit result
     */
    const ThreeDSlidingFitResult &GetSlidingFitResult() const;

    /**
     *  @brief  Get the best matched parent pfo
     *
     *  @return a pointer to the best matched parent pfo
     */
    const pandora::ParticleFlowObject *GetPredictedParentPfo() const;

    /**
     *  @brief  Set the best matched parent pfo
     *
     *  @param  pPredictedParentPfo a pointer to the best matched parent pfo
     */
    void SetPredictedParentPfo(const pandora::ParticleFlowObject *pPredictedParentPfo);

    /**
     *  @brief  Get the parent pfo
     *
     *  @return a pointer to the parent pfo
     */
    const pandora::ParticleFlowObject *GetParentPfo() const;

    /**
     *  @brief  Set the parent pfo
     *
     *  @param  pParentPfo a pointer to the parent pfo
     */
    void SetParentPfo(const pandora::ParticleFlowObject *pParentPfo);

    /**
     *  @brief  Get the vector of child pfos
     *
     *  @return the vector of pointers to child pfos
     */
    const pandora::PfoVector &GetChildPfoVector() const;

    /**
     *  @brief  Add a child pfo to the child pfo vector
     *
     *  @param  pChildPfo the pointer to the child pfo to add
     */
    void AddChildPfo(const pandora::ParticleFlowObject *const pChildPfo);

    /**
     *  @brief  Get the upstream extremal point
     *
     *  @return the upstream extremal point
     */
    const ExtremalPoint &GetUpstreamPoint() const;

    /**
     *  @brief  Set the upstream extremal point
     *
     *  @param  upstreamPoint the upstream extremal point
     */
    void SetUpstreamPoint(const ExtremalPoint &upstreamPoint);

    /**
     *  @brief  Get the downstream extremal point
     *
     *  @return the downstream extremal point
     */
    const ExtremalPoint &GetDownstreamPoint() const;

    /**
     *  @brief  Set the downstream extremal point
     *
     *  @param  downstreamPoint the downstream extremal point
     */
    void SetDownstreamPoint(const ExtremalPoint &downstreamPoint);

    /**
     *  @brief  Get the primary network score
     *
     *  @return the primary network score
     */
    float GetPrimaryScore() const;

    /**
     *  @brief  Set the primary network score
     *
     *  @param  the primary network score
     */
    void SetPrimaryScore(const float primaryScore);

    /**
     *  @brief  Get the later tier network score
     *
     *  @return the later tier network score
     */
    float GetLaterTierScore() const;

    /**
     *  @brief  Set the later tier network score
     *
     *  @param  the later tier network score
     */
    void SetLaterTierScore(const float laterTierScore);

    /**
     *  @brief  Get whether the pfo has been assigned to a particle hierarchy
     *
     *  @return whether the pfo has been assigned to a particle hierarchy
     */
    bool GetIsInHierarchy() const;

    /**
     *  @brief  Set whether the pfo has been assigned to a particle hierarchy
     *
     *  @param  isInHierarchy whether the pfo has been assigned to a particle hierarchy
     */
    void SetIsInHierarchy(const bool isInHierarchy);

    /**
     *  @brief  HierarchyPfo == operator
     * 
     *  @param  rhs the HierarchyPfo to compare
     */
    bool operator==(const HierarchyPfo &rhs) const;

    /**
     *  @brief  HierarchyPfo == operator
     * 
     *  @param  rhs the pfo to compare
     */
    bool operator==(const pandora::ParticleFlowObject *const rhs) const;

private:
    const pandora::ParticleFlowObject *m_pPfo;                ///< a pointer to the corresponding pfo
    ThreeDSlidingFitResult m_slidingFitResult;                ///< the 3D sliding fit result of the pfo
    const pandora::ParticleFlowObject *m_pPredictedParentPfo; ///< a pointer to the best matched parent pfo
    const pandora::ParticleFlowObject *m_pParentPfo;          ///< a pointer to the assigned parent pfo
    pandora::PfoVector m_childPfoVector;                      ///< the vector of pointers to the assigned child pfos
    ExtremalPoint m_upstreamPoint;                            ///< The extremal point that lies closest to the neutrino vertex
    ExtremalPoint m_downstreamPoint;                          ///< The extremal point that lies furthest from the neutrino vertex
    float m_primaryScore;                                     ///< the primary network score
    float m_laterTierScore;                                   ///< the later tier network score
    bool m_isInHierarchy;                                     ///< whether the pfo has been assigned to a particle hierarchy
};

typedef std::vector<HierarchyPfo> HierarchyPfoVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ExtremalPoint::ExtremalPoint() :
    m_isSet(false),
    m_position(-999.f, -999.f, -999.f),
    m_direction(-999.f, -999.f, -999.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ExtremalPoint::GetPosition() const
{
    return m_position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ExtremalPoint::GetDirection() const
{
    return m_direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ExtremalPoint::Set(const pandora::CartesianVector &position, const pandora::CartesianVector &direction)
{
    m_isSet = true;
    m_position = position;
    m_direction = direction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ExtremalPoint::IsSet() const
{
    return m_isSet;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ExtremalPoint &ExtremalPoint::operator=(const ExtremalPoint &rhs)
{
    this->m_isSet = rhs.m_isSet;
    this->m_position = rhs.m_position;
    this->m_direction = rhs.m_direction;

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline HierarchyPfo::HierarchyPfo(const pandora::ParticleFlowObject *pPfo, const ThreeDSlidingFitResult &slidingFitResult,
    const ExtremalPoint &upstreamPoint, const ExtremalPoint &downstreamPoint) :
    m_pPfo(pPfo),
    m_slidingFitResult(slidingFitResult),
    m_pPredictedParentPfo(nullptr),
    m_pParentPfo(nullptr),
    m_childPfoVector(pandora::PfoVector()),
    m_upstreamPoint(upstreamPoint),
    m_downstreamPoint(downstreamPoint),
    m_primaryScore(-std::numeric_limits<float>::max()),
    m_laterTierScore(-std::numeric_limits<float>::max()),
    m_isInHierarchy(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject *HierarchyPfo::GetPfo() const
{
    return m_pPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetPfo(const pandora::ParticleFlowObject *pPfo)
{
    m_pPfo = pPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ThreeDSlidingFitResult &HierarchyPfo::GetSlidingFitResult() const
{
    return m_slidingFitResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject *HierarchyPfo::GetPredictedParentPfo() const
{
    return m_pPredictedParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetPredictedParentPfo(const pandora::ParticleFlowObject *pPredictedParentPfo)
{
    m_pPredictedParentPfo = pPredictedParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject *HierarchyPfo::GetParentPfo() const
{
    return m_pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetParentPfo(const pandora::ParticleFlowObject *pParentPfo)
{
    m_pParentPfo = pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::PfoVector &HierarchyPfo::GetChildPfoVector() const
{
    return m_childPfoVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::AddChildPfo(const pandora::ParticleFlowObject *const pChildPfo)
{
    m_childPfoVector.emplace_back(pChildPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ExtremalPoint &HierarchyPfo::GetUpstreamPoint() const
{
    return m_upstreamPoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetUpstreamPoint(const ExtremalPoint &upstreamPoint)
{
    m_upstreamPoint = upstreamPoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ExtremalPoint &HierarchyPfo::GetDownstreamPoint() const
{
    return m_downstreamPoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetDownstreamPoint(const ExtremalPoint &downstreamPoint)
{
    m_downstreamPoint = downstreamPoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float HierarchyPfo::GetPrimaryScore() const
{
    return m_primaryScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetPrimaryScore(const float primaryScore)
{
    m_primaryScore = primaryScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float HierarchyPfo::GetLaterTierScore() const
{
    return m_laterTierScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetLaterTierScore(const float laterTierScore)
{
    m_laterTierScore = laterTierScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool HierarchyPfo::GetIsInHierarchy() const
{
    return m_isInHierarchy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetIsInHierarchy(const bool isInHierarchy)
{
    m_isInHierarchy = isInHierarchy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool HierarchyPfo::operator==(const HierarchyPfo &rhs) const
{
    return this->GetPfo() == rhs.GetPfo();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool HierarchyPfo::operator==(const pandora::ParticleFlowObject *const rhs) const
{
    return this->GetPfo() == rhs;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content

#endif // #ifndef LAR_HIERARCHY_PFO_H
