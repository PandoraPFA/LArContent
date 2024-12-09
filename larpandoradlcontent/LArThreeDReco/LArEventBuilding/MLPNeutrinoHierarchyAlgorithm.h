/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/MLPNeutrinoHierarchyAlgorithm.h
 *
 *  @brief  Header file for the MLP neutrino hierarchy algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MLP_NEUTRINO_HIERARCHY_ALGORITHM_H
#define LAR_MLP_NEUTRINO_HIERARCHY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

using namespace lar_content;

namespace lar_dl_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MLPNeutrinoHierarchyAlgorithm class
 */
class MLPNeutrinoHierarchyAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MLPNeutrinoHierarchyAlgorithm();

private:
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

    HierarchyPfo(const pandora::ParticleFlowObject *pPfo, const pandora::CartesianVector &upstreamVertex,
        const pandora::CartesianVector &upstreamDirection, const pandora::CartesianVector &downstreamVertex, 
        const pandora::CartesianVector &downstreamDirection);

    bool operator== (const HierarchyPfo &otherHierarchyPfo) const;

    const pandora::ParticleFlowObject* GetPfo() const;
    void SetPfo(const pandora::ParticleFlowObject *pPfo);
    const pandora::ParticleFlowObject* GetPredictedParentPfo() const;
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

    const pandora::ParticleFlowObject *m_pPfo;
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

    pandora::StatusCode Run();

    void DetermineIsobelID();

    void PrintHierarchy();

    bool GetNeutrinoPfo();

    void FillTrackShowerVectors();

    bool GetExtremalVerticesAndDirections(const pandora::ParticleFlowObject *const pPfo, pandora::CartesianVector &upstreamVertex, 
        pandora::CartesianVector &upstreamDirection, pandora::CartesianVector &downstreamVertex, pandora::CartesianVector &downstreamDirection);

    bool GetShowerDirection(const pandora::ParticleFlowObject *const pPfp, const pandora::CartesianVector &vertex, const float searchRegion, 
        pandora::CartesianVector &direction);

    void SetPrimaryScores();

    void SetPrimaryScoreTrack(HierarchyPfo &trackPfo);

    void SetPrimaryScoreShower(HierarchyPfo &showerPfo);

    void BuildPrimaryTierPass1();

    void SetLaterTierScores();

    void BuildLaterTierPass1();

    float GetLaterTierScoreTrackToTrack(const HierarchyPfo &parentPfo, const HierarchyPfo &childPfo, 
        int &parentOrientation, int &childOrientation) const;
    float GetLaterTierScoreTrackToShower(const HierarchyPfo &parentPfo, const HierarchyPfo &childPfo,
        int &parentOrientation, int &childOrientation) const;

    void BuildPandoraHierarchy();

    void PrintPandoraHierarchy();

    void CheckForOrphans();

    float GetRandomNumber() const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    //typedef std::vector<PfoRelationTool *> PfoRelationToolVector;

    std::string m_neutrinoPfoListName;
    pandora::StringVector m_pfoListNames;

    std::string m_primaryTrackBranchModelName;
    std::string m_primaryTrackClassifierModelName;
    std::string m_primaryShowerClassifierModelName;

    LArDLHelper::TorchModel m_primaryTrackBranchModel;
    LArDLHelper::TorchModel m_primaryTrackClassifierModel;
    LArDLHelper::TorchModel m_primaryShowerClassifierModel;

    float m_primaryThresholdTrackPass1;
    float m_primaryThresholdShowerPass1;
    float m_laterTierThresholdTrackPass1;
    float m_laterTierThresholdShowerPass1;

    ///////////////////////////
    std::map<const pandora::ParticleFlowObject*, int> m_isobelID;
    ///////////////////////////

    const pandora::ParticleFlowObject* m_pNeutrinoPfo;
    std::map<const pandora::ParticleFlowObject*, HierarchyPfo> m_trackPfos;
    std::map<const pandora::ParticleFlowObject*, HierarchyPfo> m_showerPfos;
    std::vector<pandora::PfoVector> m_hierarchy;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject* MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetPfo() const
{
    return m_pPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetPfo(const pandora::ParticleFlowObject* pPfo)
{
    m_pPfo = pPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject* MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetPredictedParentPfo() const
{
    return m_pPredictedParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetPredictedParentPfo(const pandora::ParticleFlowObject* pPredictedParentPfo)
{
    m_pPredictedParentPfo = pPredictedParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject* MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetParentPfo() const
{
    return m_pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetParentPfo(const pandora::ParticleFlowObject* pParentPfo)
{
    m_pParentPfo = pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::PfoVector& MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetSortedChildPfoVector()
{
    std::sort(m_childPfoVector.begin(), m_childPfoVector.end(), LArPfoHelper::SortByNHits);

    return m_childPfoVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::AddChildPfo(const pandora::ParticleFlowObject *const pChildPfo)
{
    m_childPfoVector.push_back(pChildPfo);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector& MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetUpstreamVertex() const
{
    return m_upstreamVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetUpstreamVertex(const pandora::CartesianVector &upstreamVertex)
{
    m_upstreamVertex = upstreamVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector& MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetUpstreamDirection() const
{
    return m_upstreamDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetUpstreamDirection(const pandora::CartesianVector &upstreamDirection)
{
    m_upstreamDirection = upstreamDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector& MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetDownstreamVertex() const
{
    return m_downstreamVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetDownstreamVertex(const pandora::CartesianVector &downstreamVertex)
{
    m_downstreamVertex = downstreamVertex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector& MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetDownstreamDirection() const
{
    return m_downstreamDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetDownstreamDirection(const pandora::CartesianVector &downstreamDirection)
{
    m_downstreamDirection = downstreamDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetPrimaryScore() const
{
    return m_primaryScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetPrimaryScore(const float primaryScore)
{
    m_primaryScore = primaryScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetLaterTierScore() const
{
    return m_laterTierScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetLaterTierScore(const float laterTierScore)
{
    m_laterTierScore = laterTierScore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetParentOrientation() const
{
    return m_parentOrientation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetParentOrientation(const int parentOrientation)
{
    m_parentOrientation = parentOrientation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetChildOrientation() const
{
    return m_childOrientation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetChildOrientation(const int childOrientation)
{
    m_childOrientation = childOrientation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetIsInHierarchy() const
{
    return m_isInHierarchy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetIsInHierarchy(const bool isInHierarchy)
{
    m_isInHierarchy = isInHierarchy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_NEUTRINO_HIERARCHY_ALGORITHM_H
