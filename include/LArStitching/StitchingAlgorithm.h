/**
 *  @file   LArContent/include/LArStitching/StitchingAlgorithm.h
 * 
 *  @brief  Header file for the Stitching algorithm class.
 * 
 *  $Log: $
 */
#ifndef PANDORA_STITCHING_ALGORITHM_H
#define PANDORA_STITCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "LArStitching/MultiPandora.h"

namespace lar_content
{

class StitchingTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  StitchingAlgorithm class
 */
class StitchingAlgorithm : public pandora::Algorithm
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
     *  @brief  StitchingInfo class
     */
    class StitchingInfo
    {
    public:
        // Placeholder
    };

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<StitchingTool*> StitchingToolList;
    StitchingToolList       m_algorithmToolList;        ///< The algorithm tool list
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  StitchingTool class
 */
class StitchingTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  multiPandora the multi pandora details
     *  @param  stitchingInfo the source for additional, local, stitching information
     */
    virtual void Run(const StitchingAlgorithm *const pAlgorithm, const MultiPandora &multiPandora, StitchingAlgorithm::StitchingInfo &stitchingInfo) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *StitchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new StitchingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef PANDORA_STITCHING_ALGORITHM_H
