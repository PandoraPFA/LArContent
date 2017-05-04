/**
 *  @file   larpandoracontent/LArUtility/NeutrinoParentAlgorithm.h
 * 
 *  @brief  Header file for the neutrino parent algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_NEUTRINO_PARENT_ALGORITHM_H
#define LAR_NEUTRINO_PARENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

namespace lar_content
{

class SlicingTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  NeutrinoParentAlgorithm class
 */
class NeutrinoParentAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    NeutrinoParentAlgorithm();

    /**
     *  @brief  Slice class
     */
    class Slice
    {
    public:
        pandora::CaloHitList    m_caloHitListU;                     ///< The u calo hit list
        pandora::CaloHitList    m_caloHitListV;                     ///< The v calo hit list
        pandora::CaloHitList    m_caloHitListW;                     ///< The w calo hit list
    };

    typedef std::vector<Slice> SliceList;

    typedef std::vector<pandora::HitType> HitTypeList;
    typedef std::map<pandora::HitType, std::string> HitTypeToNameMap;

    /**
     *  @brief  Copy all the input hits in an event into a single slice
     *
     *  @param  sliceList the slice list to receive the single new slice
     */
    void CopyAllHitsToSingleSlice(SliceList &sliceList) const;

private:
    pandora::StatusCode Initialize();
    pandora::StatusCode Run();

    /**
     *  @brief  Use first-pass 3D event reconstruction to slice events into separate, distinct interactions for processing
     *
     *  @param  sliceList the slice list to receive the slice list
     */
    void PerformSlicing(SliceList &sliceList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    HitTypeList                 m_hitTypeList;                      ///< The hit type list
    HitTypeToNameMap            m_caloHitListNames;                 ///< The hit type to calo hit list name map
    HitTypeToNameMap            m_clusterListNames;                 ///< The hit type to cluster list name map

    std::string                 m_caloHitListNameU;                 ///< The name of the u input calo hit list
    std::string                 m_caloHitListNameV;                 ///< The name of the v input calo hit list
    std::string                 m_caloHitListNameW;                 ///< The name of the w input calo hit list
        
    std::string                 m_clusterListNameU;                 ///< The name of the u working cluster list
    std::string                 m_clusterListNameV;                 ///< The name of the v working cluster list
    std::string                 m_clusterListNameW;                 ///< The name of the w working cluster list

    bool                        m_shouldPerformSlicing;             ///< Whether to slice events into separate, distinct interactions for processing
    SlicingTool                *m_pSlicingTool;                     ///< The address of the slicing tool
    std::string                 m_listDeletionAlgorithm;            ///< The name of the list deletion algorithm
        
    std::string                 m_clusteringAlgorithm;              ///< The name of the two dimensional clustering algorithm
    std::string                 m_listMovingAlgorithm;              ///< The name of the list moving algorithm

    pandora::StringVector       m_twoDAlgorithms;                   ///< The names of the two dimensional reconstruction algorithms
    pandora::StringVector       m_threeDAlgorithms;                 ///< The names of the three dimensional reconstruction algorithms
    pandora::StringVector       m_threeDHitAlgorithms;              ///< The names of the three dimensional hit creation algorithms
    pandora::StringVector       m_vertexAlgorithms;                 ///< The names of the vertex reconstruction algorithms
    pandora::StringVector       m_twoDMopUpAlgorithms;              ///< The names of the two dimensional mop-up algorithms
    pandora::StringVector       m_threeDMopUpAlgorithms;            ///< The names of the three dimensional mop-up algorithms
    pandora::StringVector       m_neutrinoAlgorithms;               ///< The names of the neutrino building algorithms
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  SlicingTool class
 */
class SlicingTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  caloHitListNames the hit type to calo hit list name map
     *  @param  clusterListNames the hit type to cluster list name map
     *  @param  sliceList to receive the populated slice list
     */
    virtual void Slice(const NeutrinoParentAlgorithm *const pAlgorithm, const NeutrinoParentAlgorithm::HitTypeToNameMap &caloHitListNames,
        const NeutrinoParentAlgorithm::HitTypeToNameMap &clusterListNames, NeutrinoParentAlgorithm::SliceList &sliceList) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_NEUTRINO_PARENT_ALGORITHM_H
