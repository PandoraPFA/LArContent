/**
 *  @file   LArContent/include/LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTrackFragmentsAlg.h
 *
 *  @brief  Header file for the three dimensional longitudinal track fragments algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_TRANSVERSE_TRACK_FRAGMENTS_ALGORITHM_H
#define LAR_THREE_D_TRANSVERSE_TRACK_FRAGMENTS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "LArThreeDReco/LArThreeDBase/ThreeDTracksBaseAlgorithm.h"

namespace lar
{

class TransverseFragmentTensorTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDTransverseTrackFragmentsAlg class
 */
class ThreeDTransverseTrackFragmentsAlg : public ThreeDTracksBaseAlgorithm<float> // TODO
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

    unsigned int        m_nMaxTensorToolRepeats;            ///< The maximum number of repeat loops over tensor tools

    typedef std::vector<TransverseFragmentTensorTool*> TensorToolList;
    TensorToolList      m_algorithmToolList;                ///< The algorithm tool list
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TransverseFragmentTensorTool class
 */
class TransverseFragmentTensorTool : public pandora::AlgorithmTool
{
public:
    typedef ThreeDTransverseTrackFragmentsAlg::TensorType TensorType;
    typedef std::vector<TensorType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapTensor the overlap tensor
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(ThreeDTransverseTrackFragmentsAlg *pAlgorithm, TensorType &overlapTensor) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDTransverseTrackFragmentsAlg::Factory::CreateAlgorithm() const
{
    return new ThreeDTransverseTrackFragmentsAlg();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_TRANSVERSE_TRACK_FRAGMENTS_ALGORITHM_H
