/**
 *  @file   LArContent/LArThreeDReco/LArEventBuilding/NeutrinoHierarchyAlgorithm.h
 * 
 *  @brief  Header file for the neutrino hierarchy algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_HIERARCHY_ALGORITHM_H
#define LAR_NEUTRINO_HIERARCHY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include <unordered_map>

namespace lar_content
{

class PfoRelationTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  NeutrinoHierarchyAlgorithm class
 */
class NeutrinoHierarchyAlgorithm : public pandora::Algorithm
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
    NeutrinoHierarchyAlgorithm();

    /**
     *  @brief  PfoInfo class
     */
    class PfoInfo
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pPfo the address of the pfo
         *  @param  halfWindowLayers the number of layers to use for half-window of sliding fit
         *  @param  layerPitch the sliding fit z pitch, units cm
         */
        PfoInfo(const pandora::ParticleFlowObject *const pPfo, const unsigned int halfWindowLayers, const float layerPitch);

        /**
         *  @brief  Copy constructor
         * 
         *  @param  rhs the pfo info to copy
         */
        PfoInfo(const PfoInfo &rhs);

        /**
         *  @brief  Assignment operator
         * 
         *  @param  rhs the pfo info to assign
         */
        PfoInfo &operator=(const PfoInfo &rhs);

        /**
         *  @brief  Destructor
         */
        ~PfoInfo();

        /**
         *  @brief  Get the address of the pfo
         * 
         *  @return the address of the pfo
         */
        const pandora::ParticleFlowObject *GetThisPfo() const;

        /**
         *  @brief  Get the address of the three dimensional cluster
         * 
         *  @return the address of the three dimensional cluster
         */
        const pandora::Cluster *GetCluster3D() const;

        /**
         *  @brief  Get the address of the three dimensional sliding fit result
         * 
         *  @return the address of the three dimensional sliding fit result
         */
        const ThreeDSlidingFitResult *GetSlidingFitResult3D() const;

        /**
         *  @brief  Whether the pfo is associated with the neutrino vertex
         * 
         *  @return boolean
         */
        bool IsNeutrinoVertexAssociated() const;

        /**
         *  @brief  If associated, whether association to parent (vtx or pfo) is at sliding fit inner layer
         * 
         *  @return boolean
         */
        bool IsInnerLayerAssociated() const;

        /**
         *  @brief  Get the address of the parent pfo
         * 
         *  @return the address of the parent pfo
         */
        const pandora::ParticleFlowObject *GetParentPfo() const;

        /**
         *  @brief  Get the daughter pfo list
         * 
         *  @return the daughter pfo list
         */
        const pandora::PfoList &GetDaughterPfoList() const;

        /**
         *  @brief  Set the neutrino vertex association flag
         * 
         *  @param  isNeutrinoVertexAssociated the neutrino vertex association flag
         */
        void SetNeutrinoVertexAssociation(const bool isNeutrinoVertexAssociated);

        /**
         *  @brief  Set the inner layer association flag
         * 
         *  @param  isInnerLayerAssociated the inner layer association flag
         */
        void SetInnerLayerAssociation(const bool isInnerLayerAssociated);

        /**
         *  @brief  Set the parent pfo
         * 
         *  @param  pParentPfo the address of the parent pfo
         */
        void SetParentPfo(const pandora::ParticleFlowObject *const pParentPfo);

        /**
         *  @brief  Remove the parent pfo
         */
        void RemoveParentPfo();

        /**
         *  @brief  Add a daughter pfo
         * 
         *  @param  pDaughterPfo the address of the daughter pfo to add
         */
        void AddDaughterPfo(const pandora::ParticleFlowObject *const pDaughterPfo);

        /**
         *  @brief  Remove a daughter pfo
         * 
         *  @param  pDaughterPfo the address of the daughter pfo to remove
         */
        void RemoveDaughterPfo(const pandora::ParticleFlowObject *const pDaughterPfo);

    private:
        const pandora::ParticleFlowObject  *m_pThisPfo;                     ///< The address of the pfo
        const pandora::Cluster             *m_pCluster3D;                   ///< The address of the three dimensional cluster
        const pandora::Vertex              *m_pVertex3D;                    ///< The address of the three dimensional vertex
        ThreeDSlidingFitResult             *m_pSlidingFitResult3D;          ///< The three dimensional sliding fit result

        bool                                m_isNeutrinoVertexAssociated;   ///< Whether the pfo is associated with the neutrino vertex
        bool                                m_isInnerLayerAssociated;       ///< If associated, whether association to parent (vtx or pfo) is at sliding fit inner layer
        const pandora::ParticleFlowObject  *m_pParentPfo;                   ///< The address of the parent pfo
        pandora::PfoList                    m_daughterPfoList;              ///< The daughter pfo list
    };

    typedef std::unordered_map<const pandora::ParticleFlowObject*, PfoInfo*> PfoInfoMap;

    /**
     *  @brief  Query the pfo info map and separate/extract pfos currently either acting as parents or associated with the neutrino vertex
     * 
     *  @param  pfoInfoMap the pfo info map
     *  @param  assignedPfos to receive the sorted vector of assigned pfos
     *  @param  unassignedPfos to receive the sorted vector of unassigned pfos
     */
    void SeparatePfos(const NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap, pandora::PfoVector &assignedPfos, pandora::PfoVector &unassignedPfos) const;

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Get the address of the input neutrino pfo - enforces only one pfo present in input list; can return NULL if no neutrino exists
     * 
     *  @param  to receive the address of the input neutrino pfo
     */
    void GetNeutrinoPfo(const pandora::ParticleFlowObject *&pNeutrinoPfo) const;

    /**
     *  @brief  Get the list of candidate daughter pfos
     * 
     *  @param  candidateDaughterPfoList to receive the candidate daughter pfo list
     */
    void GetCandidateDaughterPfoList(pandora::PfoList &candidateDaughterPfoList) const;

    /**
     *  @brief  Process a provided pfo list and populate an initial pfo info map
     * 
     *  @param  pfoList the provided pfo list
     *  @param  pfoInfoMap to receive the initial pfo info map
     */
    void GetInitialPfoInfoMap(const pandora::PfoList &pfoList, PfoInfoMap &pfoInfoMap) const;

    /**
     *  @brief  Process the information in a pfo info map, creating pfo parent/daughter links
     * 
     *  @param  pNeutrinoPfo the address of the (original) parent neutrino pfo
     *  @param  candidateDaughterPfoList the list of candidate daughter pfos
     *  @param  pfoInfoMap the pfo info map
     */
    void ProcessPfoInfoMap(const pandora::ParticleFlowObject *const pNeutrinoPfo, const pandora::PfoList &candidateDaughterPfoList,
        const PfoInfoMap &pfoInfoMap) const;

    /**
     *  @brief  Display the information in a pfo info map, visualising pfo parent/daughter links
     * 
     *  @param  pNeutrinoPfo the address of the (original) parent neutrino pfo
     *  @param  pfoInfoMap the pfo info map
     */
    void DisplayPfoInfoMap(const pandora::ParticleFlowObject *const pNeutrinoPfo, const PfoInfoMap &pfoInfoMap) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<PfoRelationTool*> PfoRelationToolList;
    PfoRelationToolList             m_algorithmToolList;        ///< The algorithm tool list

    std::string                     m_neutrinoPfoListName;      ///< The neutrino pfo list name
    pandora::StringVector           m_daughterPfoListNames;     ///< The list of daughter pfo list names

    unsigned int                    m_halfWindowLayers;         ///< The number of layers to use for half-window of sliding fit
    bool                            m_displayPfoInfoMap;        ///< Whether to display the pfo info map (if monitoring is enabled)
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  PfoRelationTool class
 */
class PfoRelationTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pNeutrinoVertex the address of the three dimensional neutrino interaction vertex
     *  @param  pfoInfoMap mapping from pfos to three dimensional clusters, sliding fits, vertices, etc.
     */
    virtual void Run(NeutrinoHierarchyAlgorithm *const pAlgorithm, const pandora::Vertex *const pNeutrinoVertex, NeutrinoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *NeutrinoHierarchyAlgorithm::Factory::CreateAlgorithm() const
{
    return new NeutrinoHierarchyAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject *NeutrinoHierarchyAlgorithm::PfoInfo::GetThisPfo() const
{
    return m_pThisPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *NeutrinoHierarchyAlgorithm::PfoInfo::GetCluster3D() const
{
    return m_pCluster3D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ThreeDSlidingFitResult *NeutrinoHierarchyAlgorithm::PfoInfo::GetSlidingFitResult3D() const
{
    return m_pSlidingFitResult3D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool NeutrinoHierarchyAlgorithm::PfoInfo::IsNeutrinoVertexAssociated() const
{
    return m_isNeutrinoVertexAssociated;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool NeutrinoHierarchyAlgorithm::PfoInfo::IsInnerLayerAssociated() const
{
    return m_isInnerLayerAssociated;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject *NeutrinoHierarchyAlgorithm::PfoInfo::GetParentPfo() const
{
    return m_pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::PfoList &NeutrinoHierarchyAlgorithm::PfoInfo::GetDaughterPfoList() const
{
    return m_daughterPfoList;
}

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_HIERARCHY_ALGORITHM_H
