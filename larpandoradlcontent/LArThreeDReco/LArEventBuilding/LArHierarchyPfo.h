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
  *  @brief  HierarchyPfo class
  */
class HierarchyPfo
{
    public:
    /**
     *  @brief  Default constructor
     */
    HierarchyPfo();

    HierarchyPfo(const bool isTrack, const pandora::ParticleFlowObject *pPfo, const ThreeDSlidingFitResult &threeDSlidingFitResult, 
        const pandora::CartesianVector &upstreamVertex, const pandora::CartesianVector &upstreamDirection, 
        const pandora::CartesianVector &downstreamVertex, const pandora::CartesianVector &downstreamDirection);

    bool operator== (const HierarchyPfo &otherHierarchyPfo) const;

    bool GetIsTrack() const;
    void SetIsTrack(const bool isTrack);
    const pandora::ParticleFlowObject* GetPfo() const;
    void SetPfo(const pandora::ParticleFlowObject *pPfo);
    const pandora::ParticleFlowObject* GetPredictedParentPfo() const;
    const ThreeDSlidingFitResult& GetSlidingFitResult() const;
    void SetPredictedParentPfo(const pandora::ParticleFlowObject *pPredictedParentPfo);
    const pandora::ParticleFlowObject* GetParentPfo() const;
    void SetParentPfo(const pandora::ParticleFlowObject *pParentPfo);
    const pandora::PfoVector& GetSortedChildPfoVector();
    void AddChildPfo(const pandora::ParticleFlowObject *const pChildPfo);
    const pandora::CartesianVector& GetUpstreamVertex() const;
    void SetUpstreamVertex(const pandora::CartesianVector &upstreamVertex);
    const pandora::CartesianVector& GetUpstreamDirection() const;
    void SetUpstreamDirection(const pandora::CartesianVector &upstreamDirection);
    const pandora::CartesianVector& GetDownstreamVertex() const;
    void SetDownstreamVertex(const pandora::CartesianVector &downstreamVertex);
    const pandora::CartesianVector& GetDownstreamDirection() const;
    void SetDownstreamDirection(const pandora::CartesianVector &downstreamDirection);

    float GetPrimaryScore() const;
    void SetPrimaryScore(const float primaryScore);
    float GetLaterTierScore() const;
    void SetLaterTierScore(const float laterTierScore);
    int GetParentOrientation() const;
    void SetParentOrientation(const int parentOrientation);
    int GetChildOrientation() const;
    void SetChildOrientation(const int childOrientation);
    bool GetIsInHierarchy() const;
    void SetIsInHierarchy(const bool isInHierarchy);

    private:

    bool m_isTrack;
    const pandora::ParticleFlowObject *m_pPfo;
    ThreeDSlidingFitResult m_slidingFitResult;
    const pandora::ParticleFlowObject *m_pPredictedParentPfo;
    const pandora::ParticleFlowObject *m_pParentPfo;
    pandora::PfoVector m_childPfoVector;

    pandora::CartesianVector m_upstreamVertex; // extremal point closest to vertex
    pandora::CartesianVector m_upstreamDirection; // direction at extremal point closest to vertex
    pandora::CartesianVector m_downstreamVertex; // extremal point furthest from vertex
    pandora::CartesianVector m_downstreamDirection; // direction at extremal point furthest from vertex

    float m_primaryScore;
    float m_laterTierScore;
    int m_parentOrientation; // 0 - start is not closest to neutrino vertex
    int m_childOrientation; // 0 - start is not closest to neutrino vertex
    bool m_isInHierarchy;

};

typedef std::map<const pandora::ParticleFlowObject*, HierarchyPfo> HierarchyPfoMap;
typedef std::pair<const pandora::ParticleFlowObject*, HierarchyPfo> HierarchyPfoMapEntry;

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
        m_parentOrientation(-1),
        m_childOrientation(-1),
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

inline const pandora::PfoVector& HierarchyPfo::GetSortedChildPfoVector()
{
    std::sort(m_childPfoVector.begin(), m_childPfoVector.end(), LArPfoHelper::SortByNHits);

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

inline int HierarchyPfo::GetParentOrientation() const
{
    return m_parentOrientation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetParentOrientation(const int parentOrientation)
{
    m_parentOrientation = parentOrientation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int HierarchyPfo::GetChildOrientation() const
{
    return m_childOrientation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void HierarchyPfo::SetChildOrientation(const int childOrientation)
{
    m_childOrientation = childOrientation;
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
