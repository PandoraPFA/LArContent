/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerFragments/ThreeDRemnantsAlgorithm.h
 *
 *  @brief  Header file for the three dimensional remnants algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_REMNANTS_ALGORITHM_H
#define LAR_THREE_D_REMNANTS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArOverlapTensor.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/ThreeViewMatchingAlgorithm.h"

namespace lar_content
{

class RemnantTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDRemnantsAlgorithm class
 */
class ThreeDRemnantsAlgorithm : public ThreeViewMatchingAlgorithm<float>
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDRemnantsAlgorithm();

    void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;

private:
    void CalculateOverlapResult(const pandora::Cluster *const pClusterU, const pandora::Cluster *const pClusterV, const pandora::Cluster *const pClusterW);
    void ExamineTensor();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<RemnantTensorTool*> RemnantTensorToolVector;
    RemnantTensorToolVector     m_algorithmToolVector;      ///< The algorithm tool list

    unsigned int                m_nMaxTensorToolRepeats;    ///< The maximum number of repeat loops over tensor tools
    unsigned int                m_minClusterCaloHits;       ///< The selection cut on the number of cluster calo hits
    float                       m_xOverlapWindow;           ///< The sampling pitch in the x coordinate
    float                       m_pseudoChi2Cut;            ///< The selection cut on the matched chi2
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  RemnantTensorTool class
 */
class RemnantTensorTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeDRemnantsAlgorithm::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeDRemnantsAlgorithm *const pAlgorithm, TensorType &overlapTensor) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_REMNANTS_ALGORITHM_H
