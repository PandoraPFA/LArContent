/**
 *  @file   LArContent/include/LArThreeDReco/LArEventBuilding/PfoHierarchyAlgorithm.h
 * 
 *  @brief  Header file for the pfo hierarchy algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_PFO_HIERARCHY_ALGORITHM_H
#define LAR_PFO_HIERARCHY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArObjects/LArThreeDSlidingFitResult.h"

#include <unordered_map>

namespace lar_content
{

class PfoRelationTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  PfoHierarchyAlgorithm class
 */
class PfoHierarchyAlgorithm : public pandora::Algorithm
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
    PfoHierarchyAlgorithm();

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
     *  @param  assignedPfos to receive the list of assigned pfos
     *  @param  unassignedPfos to receive the list of unassigned pfos
     */
    void SeparatePfos(const PfoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap, pandora::PfoList &assignedPfos, pandora::PfoList &unassignedPfos) const;

private:
    pandora::StatusCode Run();

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
     *  @param  pfoInfoMap the pfo info map
     */
    void ProcessPfoInfoMap(const pandora::ParticleFlowObject *const pNeutrinoPfo, const PfoInfoMap &pfoInfoMap) const;

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

    std::string                     m_neutrinoListName;         ///< The input list of pfo list names

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
    virtual void Run(PfoHierarchyAlgorithm *const pAlgorithm, const pandora::Vertex *const pNeutrinoVertex, PfoHierarchyAlgorithm::PfoInfoMap &pfoInfoMap) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *PfoHierarchyAlgorithm::Factory::CreateAlgorithm() const
{
    return new PfoHierarchyAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject *PfoHierarchyAlgorithm::PfoInfo::GetThisPfo() const
{
    return m_pThisPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *PfoHierarchyAlgorithm::PfoInfo::GetCluster3D() const
{
    return m_pCluster3D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ThreeDSlidingFitResult *PfoHierarchyAlgorithm::PfoInfo::GetSlidingFitResult3D() const
{
    return m_pSlidingFitResult3D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool PfoHierarchyAlgorithm::PfoInfo::IsNeutrinoVertexAssociated() const
{
    return m_isNeutrinoVertexAssociated;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool PfoHierarchyAlgorithm::PfoInfo::IsInnerLayerAssociated() const
{
    return m_isInnerLayerAssociated;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::ParticleFlowObject *PfoHierarchyAlgorithm::PfoInfo::GetParentPfo() const
{
    return m_pParentPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::PfoList &PfoHierarchyAlgorithm::PfoInfo::GetDaughterPfoList() const
{
    return m_daughterPfoList;
}

} // namespace lar_content

#endif // #ifndef LAR_PFO_HIERARCHY_ALGORITHM_H
