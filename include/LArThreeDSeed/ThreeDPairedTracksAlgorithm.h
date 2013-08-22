/**
 *  @file   LArContent/include/LArThreeDSeed/ThreeDPairedTracksAlgorithm.h
 * 
 *  @brief  Header file for the 3D paired tracks algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_THREE_D_PAIRED_TRACKS_ALGORITHM_H
#define LAR_THREE_D_PAIRED_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  ThreeDPairedTracksAlgorithm class
 */
class ThreeDPairedTracksAlgorithm : public pandora::Algorithm
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
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);



  
    std::string             m_inputClusterListNameU;           ///< Input cluster list name for U view
    std::string             m_inputClusterListNameV;           ///< Input cluster list name for V view
    std::string             m_inputClusterListNameW;           ///< Input cluster list name for W view
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDPairedTracksAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDPairedTracksAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_PAIRED_TRACKS_ALGORITHM_H
