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

    HierarchyPfo(const pandora::ParticleFlowObject *pPfo);

    bool operator== (const HierarchyPfo &otherHierarchyPfo) const;

    const pandora::ParticleFlowObject* GetPfo() const;
    void SetPfo(const pandora::ParticleFlowObject *pPfo);
    const pandora::ParticleFlowObject* GetParentPfo() const;
    void SetParentPfo(const pandora::ParticleFlowObject *pParentPfo);
    const pandora::ParticleFlowObject* GetPredictedParentPfo() const;
    void SetPredictedParentPfo(const pandora::ParticleFlowObject *pPredictedParentPfo);
    float GetPrimaryScore() const;
    void SetPrimaryScore(const float primaryScore);
    float GetLaterTierScore() const;
    void SetLaterTierScore(const float laterTierScore);
    int GetParentOrientation() const;
    void SetParentOrientation(const int parentOrientation);
    int GetChildOrientation() const;
    void SetChildOrientation(const int childOrientation);
    bool GetIsSet() const;
    void SetIsSet(const bool isSet);

    private:

    const pandora::ParticleFlowObject *m_pPfo;
    const pandora::ParticleFlowObject *m_pParentPfo;
    const pandora::ParticleFlowObject *m_pPredictedParentPfo;
    float m_primaryScore;
    float m_laterTierScore;
    int m_parentOrientation; // 0 - start is not closest to neutrino vertex
    int m_childOrientation; // 0 - start is not closest to neutrino vertex
    bool m_isSet;
    };

    pandora::StatusCode Run();

    void FillTrackShowerVectors();

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

    float GetRandomNumber() const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    //typedef std::vector<PfoRelationTool *> PfoRelationToolVector;

    pandora::StringVector m_pfoListNames;

    float m_primaryThresholdTrackPass1;
    float m_primaryThresholdShowerPass1;
    float m_laterTierThresholdTrackPass1;
    float m_laterTierThresholdShowerPass1;

    std::map<const pandora::ParticleFlowObject*, HierarchyPfo> m_trackPfos;
    std::map<const pandora::ParticleFlowObject*, HierarchyPfo> m_showerPfos;
    std::vector<std::vector<const pandora::ParticleFlowObject*>> m_hierarchy;
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

inline bool MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::GetIsSet() const
{
    return m_isSet;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void MLPNeutrinoHierarchyAlgorithm::HierarchyPfo::SetIsSet(const bool isSet)
{
    m_isSet = isSet;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_NEUTRINO_HIERARCHY_ALGORITHM_H
