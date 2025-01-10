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

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

using namespace lar_content;

namespace lar_dl_content
{

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
     *  @param  isTrack whether the input pfo is track-like
     *  @param  pPfo pointer to the input pfo
     *  @param  threeDSlidingFitResult 3D sliding fit of the input pfo
     *  @param  upstreamVertex the particle endpoint that is closest to the neutrino vertex
     *  @param  upstreamDirection the direction at the upstream vertex (pointing into the particle)
     *  @param  downstreamVertex the particle endpoint that is furthest from the neutrino vertex
     *  @param  downstreamDirection the direction at the downstream vertex (pointing into the particle)
     */
    HierarchyPfo(const bool isTrack, const pandora::ParticleFlowObject *pPfo, const ThreeDSlidingFitResult &threeDSlidingFitResult, 
        const pandora::CartesianVector &upstreamVertex, const pandora::CartesianVector &upstreamDirection, 
        const pandora::CartesianVector &downstreamVertex, const pandora::CartesianVector &downstreamDirection);

    /**
     *  @brief  HierarchyPfo == operator
     * 
     *  @param  otherHierarchyPfo the HierarchyPfo to compare
     */
    bool operator== (const HierarchyPfo &otherHierarchyPfo) const;

    /**
     *  @brief  Return whether the pfo is track-like
     *
     *  @return whether the pfo is track-like
     */
    bool GetIsTrack() const;

    /**
     *  @brief  Set whether the pfo is track-like
     *
     *  @param  isTrack whether the pfo is track-like
     */
    void SetIsTrack(const bool isTrack);

    /**
     *  @brief  Get the pfo
     *
     *  @return a pointer to the pfo
     */
    const pandora::ParticleFlowObject* GetPfo() const;

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
    const ThreeDSlidingFitResult& GetSlidingFitResult() const;

    /**
     *  @brief  Get the best matched parent pfo
     *
     *  @return a pointer to the best matched parent pfo
     */
    const pandora::ParticleFlowObject* GetPredictedParentPfo() const;

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
    const pandora::ParticleFlowObject* GetParentPfo() const;

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
    const pandora::PfoVector& GetChildPfoVector() const;

    /**
     *  @brief  Add a child pfo to the child pfo vector
     *
     *  @param  pChildPfo the pointer to the child pfo to add
     */
    void AddChildPfo(const pandora::ParticleFlowObject *const pChildPfo);

    /**
     *  @brief  Get the upstream vertex
     *
     *  @return the upstream vertex
     */
    const pandora::CartesianVector& GetUpstreamVertex() const;

    /**
     *  @brief  Set the upstream vertex
     *
     *  @param  upstreamVertex the upstream vertex
     */
    void SetUpstreamVertex(const pandora::CartesianVector &upstreamVertex);

    /**
     *  @brief  Get the upstream direction
     *
     *  @return the upstream direction
     */
    const pandora::CartesianVector& GetUpstreamDirection() const;

    /**
     *  @brief  Set the upstream direction
     *
     *  @param  upstreamDirection the upstream direction
     */
    void SetUpstreamDirection(const pandora::CartesianVector &upstreamDirection);

    /**
     *  @brief  Get the downstream vertex
     *
     *  @return the downstream vertex
     */
    const pandora::CartesianVector& GetDownstreamVertex() const;

    /**
     *  @brief  Set the downstream vertex
     *
     *  @param  downstreamVertex the downstream vertex
     */
    void SetDownstreamVertex(const pandora::CartesianVector &downstreamVertex);

    /**
     *  @brief  Get the downstream direction
     *
     *  @return the downstream direction
     */
    const pandora::CartesianVector& GetDownstreamDirection() const;

    /**
     *  @brief  Set the downstream direction
     *
     *  @param  downstreamDirection the downstream direction
     */
    void SetDownstreamDirection(const pandora::CartesianVector &downstreamDirection);

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

private:

    bool m_isTrack;                                           ///< whether the pfo is track-like
    const pandora::ParticleFlowObject *m_pPfo;                ///< a pointer to the corresponding pfo
    ThreeDSlidingFitResult m_slidingFitResult;                ///< the 3D sliding fit result of the pfo
    const pandora::ParticleFlowObject *m_pPredictedParentPfo; ///< a pointer to the best matched parent pfo
    const pandora::ParticleFlowObject *m_pParentPfo;          ///< a pointer to the assigned parent pfo
    pandora::PfoVector m_childPfoVector;                      ///< the vector of pointers to the assigned child pfos
    pandora::CartesianVector m_upstreamVertex;                ///< the particle endpoint that is closest to the neutrino vertex
    pandora::CartesianVector m_upstreamDirection;             ///< the direction at the upstream vertex (pointing into the particle)
    pandora::CartesianVector m_downstreamVertex;              ///< the particle endpoint that is furthest from the neutrino vertex
    pandora::CartesianVector m_downstreamDirection;           ///< the direction at the downstream vertex (pointing into the particle)
    float m_primaryScore;                                     ///< the primary network score 
    float m_laterTierScore;                                   ///< the later tier network score 
    bool m_isInHierarchy;                                     ///< whether the pfo has been assigned to a particle hierarchy
};

typedef std::map<const pandora::ParticleFlowObject*, HierarchyPfo> HierarchyPfoMap;

//------------------------------------------------------------------------------------------------------------------------------------------

inline HierarchyPfo::HierarchyPfo(const bool isTrack, const pandora::ParticleFlowObject *pPfo, const ThreeDSlidingFitResult &slidingFitResult,
    const pandora::CartesianVector &upstreamVertex, const pandora::CartesianVector &upstreamDirection, const pandora::CartesianVector &downstreamVertex, 
    const pandora::CartesianVector &downstreamDirection) :
        m_isTrack(isTrack),
        m_pPfo(pPfo),
        m_slidingFitResult(slidingFitResult),
        m_pPredictedParentPfo(nullptr),
        m_pParentPfo(nullptr),
        m_childPfoVector(pandora::PfoVector()),
        m_upstreamVertex(upstreamVertex),
        m_upstreamDirection(upstreamDirection),
        m_downstreamVertex(downstreamVertex),
        m_downstreamDirection(downstreamDirection),
        m_primaryScore(-std::numeric_limits<float>::max()),
        m_laterTierScore(-std::numeric_limits<float>::max()),
        m_isInHierarchy(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool HierarchyPfo::operator== (const HierarchyPfo &otherHierarchyPfo) const
{
    return this->GetPfo() == otherHierarchyPfo.GetPfo();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool HierarchyPfo::GetIsTrack() const
{
    return m_isTrack;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetIsTrack(const bool isTrack)
{
    m_isTrack = isTrack;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject* HierarchyPfo::GetPfo() const
{
    return m_pPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetPfo(const pandora::ParticleFlowObject* pPfo)
{
    m_pPfo = pPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ThreeDSlidingFitResult& HierarchyPfo::GetSlidingFitResult() const
{
    return m_slidingFitResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject* HierarchyPfo::GetPredictedParentPfo() const
{
    return m_pPredictedParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetPredictedParentPfo(const pandora::ParticleFlowObject* pPredictedParentPfo)
{
    m_pPredictedParentPfo = pPredictedParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject* HierarchyPfo::GetParentPfo() const
{
    return m_pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetParentPfo(const pandora::ParticleFlowObject* pParentPfo)
{
    m_pParentPfo = pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::PfoVector& HierarchyPfo::GetChildPfoVector() const
{
    return m_childPfoVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::AddChildPfo(const pandora::ParticleFlowObject *const pChildPfo)
{
    m_childPfoVector.push_back(pChildPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector& HierarchyPfo::GetUpstreamVertex() const
{
    return m_upstreamVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetUpstreamVertex(const pandora::CartesianVector &upstreamVertex)
{
    m_upstreamVertex = upstreamVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector& HierarchyPfo::GetUpstreamDirection() const
{
    return m_upstreamDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetUpstreamDirection(const pandora::CartesianVector &upstreamDirection)
{
    m_upstreamDirection = upstreamDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector& HierarchyPfo::GetDownstreamVertex() const
{
    return m_downstreamVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetDownstreamVertex(const pandora::CartesianVector &downstreamVertex)
{
    m_downstreamVertex = downstreamVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector& HierarchyPfo::GetDownstreamDirection() const
{
    return m_downstreamDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetDownstreamDirection(const pandora::CartesianVector &downstreamDirection)
{
    m_downstreamDirection = downstreamDirection;
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

} // namespace lar_dl_content

#endif // #ifndef LAR_HIERARCHY_PFO_H
