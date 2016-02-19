/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/NeutrinoDaughterConsolidationAlgorithm.h
 *
 *  @brief  Header file for the neutrino daughter consolidation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_DAUGHTER_CONSOLIDATION_ALGORITHM_H
#define LAR_NEUTRINO_DAUGHTER_CONSOLIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  NeutrinoDaughterConsolidationAlgorithm class
 */
class NeutrinoDaughterConsolidationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    NeutrinoDaughterConsolidationAlgorithm();

private:
    /**
     *  @brief  PfoAssociation class
     */
    class PfoAssociation
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pParentPfo the address of the parent pfo
         *  @param  pDaughterPfo the address of the daughter pfo
         *  @param  score the association score
         */
        PfoAssociation(const pandora::Pfo *const pParentPfo, const pandora::Pfo *const pDaughterPfo, const float score);

        /**
         *  @brief  Get the address of the parent pfo
         *
         *  @return the address of the daughter pfo
         */
        const pandora::Pfo *GetParentPfo() const;

        /**
         *  @brief  Get the address of the daughter pfo
         *
         *  @return the address of the daughter pfo
         */
        const pandora::Pfo *GetDaughterPfo() const;

        /**
         *  @brief  Get the association score
         *
         *  @return the association score
         */
        float GetScore() const;

        /**
         *  @brief  operator<
         * 
         *  @param  rhs the pfo association object for comparison
         * 
         *  @return boolean
         */
        bool operator< (const PfoAssociation &rhs) const;

    private:
        const pandora::Pfo     *m_pParentPfo;               ///< The address of the parent pfo
        const pandora::Pfo     *m_pDaughterPfo;             ///< The address of the daughter pfo
        float                   m_score;                    ///< The association score
    };

    typedef std::vector<PfoAssociation> PfoAssociationList;

    pandora::StatusCode Run();

    /**
     *  @brief 
     * 
     *  @param  pNeutrinoPfo
     *  @param  pNeutrinoVertex
     */
    void GetNeutrinoProperties(const pandora::ParticleFlowObject *&pNeutrinoPfo, const pandora::Vertex *&pNeutrinoVertex) const;

    /**
     *  @brief  
     * 
     *  @param  pNeutrinoPfo
     *  @param  pNeutrinoVertex
     *  @param  parentCandidates
     *  @param  daughtersToMerge
     */
    void ClassifyPrimaryParticles(const pandora::ParticleFlowObject *const pNeutrinoPfo, const pandora::Vertex *const pNeutrinoVertex,
        pandora::PfoList &parentCandidates, pandora::PfoList &daughtersToMerge) const;

    /**
     *  @brief  
     * 
     *  @param  
     *  @param  
     *  @param  
     *  @param  pfoAssociationList
     */
    void GetPfoAssociations(const pandora::Vertex *const pNeutrinoVertex, const pandora::PfoList &parentCandidates,
        const pandora::PfoList &daughtersToMerge, PfoAssociationList &pfoAssociationList) const;

    /**
     *  @brief  
     * 
     *  @param  pfoAssociationList
     * 
     *  @return boolean
     */
    bool ProcessPfoAssociations(const PfoAssociationList &pfoAssociationList) const;

    /**
     *  @brief  Merge and delete a pair of pfos, with a specific set of conventions for cluster merging, vertex use, etc.
     * 
     *  @param  pPfoToEnlarge the address of the pfo to enlarge
     *  @param  pPfoToDelete the address of the pfo to delete (will become a dangling pointer)
     */
    void MergeAndDeletePfos(const pandora::ParticleFlowObject *const pPfoToEnlarge, const pandora::ParticleFlowObject *const pPfoToDelete) const;

    /**
     *  @brief  Select the parent cluster (same hit type and most hits) using a provided cluster list and hit type
     * 
     *  @param  clusterList the cluster list
     *  @param  hitType the hit type
     * 
     *  @return the address of the parent cluster
     */
    const pandora::Cluster *GetParentCluster(const pandora::ClusterList &clusterList, const pandora::HitType hitType) const;

    /**
     *  @brief  Find the name of the list hosting a specific object
     * 
     *  @param  pT the address of the object
     * 
     *  @return the name of the list
     */
    template <typename T>
    const std::string &GetListName(const T *const pT) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_neutrinoPfoListName;          ///< The neutrino pfo list name
    pandora::StringVector   m_daughterListNames;            ///< The list of potential daughter object list names
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *NeutrinoDaughterConsolidationAlgorithm::Factory::CreateAlgorithm() const
{
    return new NeutrinoDaughterConsolidationAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline NeutrinoDaughterConsolidationAlgorithm::PfoAssociation::PfoAssociation(const pandora::Pfo *const pParentPfo, const pandora::Pfo *const pDaughterPfo,
        const float score) :
    m_pParentPfo(pParentPfo),
    m_pDaughterPfo(pDaughterPfo),
    m_score(score)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Pfo *NeutrinoDaughterConsolidationAlgorithm::PfoAssociation::GetParentPfo() const
{
    return m_pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Pfo *NeutrinoDaughterConsolidationAlgorithm::PfoAssociation::GetDaughterPfo() const
{
    return m_pDaughterPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float NeutrinoDaughterConsolidationAlgorithm::PfoAssociation::GetScore() const
{
    return m_score;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool NeutrinoDaughterConsolidationAlgorithm::PfoAssociation::operator< (const PfoAssociation &rhs) const
{
    return (this->GetScore() > rhs.GetScore());
}

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_DAUGHTER_CONSOLIDATION_ALGORITHM_H
