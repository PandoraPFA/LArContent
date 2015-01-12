/**
 *  @file   LArContent/include/LArMonitoring/ParticleAnalysisAlgorithm.h
 *
 *  @brief  Header file for the particle analysis algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_PARTICLE_ANALYSIS_ALGORITHM_H
#define LAR_PARTICLE_ANALYSIS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ParticleAnalysisAlgorithm class
 */
class ParticleAnalysisAlgorithm: public pandora::Algorithm
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
    ParticleAnalysisAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~ParticleAnalysisAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Collect up all downstream pfos from an input pfo list
     *
     *  @param  pParentPfo the input Pfo list
     *  @param  pfoVector the output pfo vector
     */
    void GetConnectedPfos(const pandora::PfoList *const pPfoList, pandora::PfoVector &pfoVector) const;

    std::string     m_pfoListName;              ///< Name of input Pfo list
    std::string     m_fileName;                 ///< Name of output file
    std::string     m_treeName;                 ///< Name of output tree
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ParticleAnalysisAlgorithm::Factory::CreateAlgorithm() const
{
    return new ParticleAnalysisAlgorithm();
}

} // namespace lar_content

#endif // LAR_PARTICLE_ANALYSIS_ALGORITHM_H
