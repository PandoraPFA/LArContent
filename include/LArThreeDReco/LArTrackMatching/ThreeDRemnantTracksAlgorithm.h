/**
 *  @file   LArContent/include/LArThreeDReco/LArTrackMatching/ThreeDRemnantTracksAlgorithm.h
 *
 *  @brief  Header file for the three dimensional remnant tracks algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_REMNANT_TRACKS_ALGORITHM_H
#define LAR_THREE_D_REMNANT_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArObjects/LArOverlapTensor.h"
#include "LArObjects/LArTrackOverlapResult.h"

#include "LArThreeDReco/ThreeDBaseAlgorithm.h"

namespace lar
{

class RemnantTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDRemnantTracksAlgorithm class
 */
class ThreeDRemnantTracksAlgorithm : public ThreeDBaseAlgorithm<float> // TODO
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

private:
    void CalculateOverlapResult(pandora::Cluster *pClusterU, pandora::Cluster *pClusterV, pandora::Cluster *pClusterW);
    void ExamineTensor();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    unsigned int    m_nMaxTensorToolRepeats;            ///< The maximum number of repeat loops over tensor tools

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
    typedef ThreeDRemnantTracksAlgorithm::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeDRemnantTracksAlgorithm *pAlgorithm, TensorType &overlapTensor) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDRemnantTracksAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDRemnantTracksAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_REMNANT_TRACKS_ALGORITHM_H
