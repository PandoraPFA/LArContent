/**
 *  @file   larpandoracontent/LArStitching/StitchingAlgorithm.h
 * 
 *  @brief  Header file for the Stitching algorithm class.
 * 
 *  $Log: $
 */
#ifndef PANDORA_STITCHING_ALGORITHM_H
#define PANDORA_STITCHING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

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

    typedef std::unordered_map<const pandora::ParticleFlowObject*, int> PfoToVolumeIdMap;

    /**
     *  @brief  StitchingInfo class
     */
    class StitchingInfo
    {
    public:
        // Placeholder
        PfoToVolumeIdMap    m_pfoToVolumeIdMap;         ///< The pfo to volume id map
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
     *  @param  stitchingInfo the source for additional, local, stitching information
     */
    virtual void Run(const StitchingAlgorithm *const pAlgorithm, StitchingAlgorithm::StitchingInfo &stitchingInfo) = 0;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *StitchingAlgorithm::Factory::CreateAlgorithm() const
{
    return new StitchingAlgorithm();
}

} // namespace lar_content

#endif // #ifndef PANDORA_STITCHING_ALGORITHM_H
