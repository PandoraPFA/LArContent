/**
 *  @file   LArContent/include/LArThreeDReco/LArShowerFragments/ThreeDRemnantsAlgorithm.h
 *
 *  @brief  Header file for the three dimensional remnants algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_REMNANTS_ALGORITHM_H
#define LAR_THREE_D_REMNANTS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "LArObjects/LArOverlapTensor.h"

#include "LArThreeDReco/LArThreeDBase/ThreeDBaseAlgorithm.h"

namespace lar
{

class RemnantTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDRemnantsAlgorithm class
 */
class ThreeDRemnantsAlgorithm : public ThreeDBaseAlgorithm<float>
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

    virtual void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;

private:
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);
    void ExamineTensor();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_nMaxTensorToolRepeats;            ///< The maximum number of repeat loops over tensor tools
    float           m_minXOverlap;                      ///< The minimum X overlap required for clusters
    float           m_minXOverlapFraction;              ///< The minimum X overlap fraction required for clusters
    float           m_maxXDisplacement;                 ///< The maximum X extent for a group of associated hits
    float           m_maxMatchedChi2;                   ///< The maximum chi2 for a group of associated hits
    float           m_minMatchedFraction;               ///< The minimum fraction of matched hits required for a good match

    typedef std::vector<RemnantTensorTool*> RemnantTensorToolList;
    RemnantTensorToolList  m_algorithmToolList;         ///< The algorithm tool list
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
    virtual bool Run(ThreeDRemnantsAlgorithm *pAlgorithm, TensorType &overlapTensor) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDRemnantsAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDRemnantsAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_REMNANTS_ALGORITHM_H
