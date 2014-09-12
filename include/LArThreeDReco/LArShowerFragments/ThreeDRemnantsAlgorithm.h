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

namespace lar_content
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

    void SelectInputClusters(const pandora::ClusterList *const pInputClusterList, pandora::ClusterList &selectedClusterList) const;
    void SetPfoParameters(const ProtoParticle &protoParticle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const;

private:
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);
    void ExamineTensor();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_nMaxTensorToolRepeats;            ///< The maximum number of repeat loops over tensor tools

    unsigned int    m_minClusterCaloHits;               ///< The selection cut on the number of cluster calo hits
    float           m_xOverlapWindow;                   ///< The sampling pitch in the x coordinate
    float           m_pseudoChi2Cut;                    ///< The selection cut on the matched chi2

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

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_REMNANTS_ALGORITHM_H
