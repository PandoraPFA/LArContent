/**
 *  @file   LArContent/include/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h
 *
 *  @brief  Header file for the three dimensional hit creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_HIT_CREATION_ALGORITHM_H
#define LAR_THREE_D_HIT_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

namespace lar_content
{

class HitCreationBaseTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDHitCreationAlgorithm::Algorithm class
 */
class ThreeDHitCreationAlgorithm : public pandora::Algorithm
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
     *  @brief  Get the list of two dimensional calo hits (in a specified pfo) for which three dimensional hits have been created
     *
     *  @param  pPfo the address of the pfo
     *  @param  caloHitList to receive the list of two dimensional calo hits for which three dimensional hits have been created
     */
    void GetUsedTwoDHits(const pandora::ParticleFlowObject *const pPfo, pandora::CaloHitList &caloHitList) const;

    /**
     *  @brief  Get the list of two dimensional calo hits (in a specified pfo) for which no three dimensional hits have been created
     *
     *  @param  pPfo the address of the pfo
     *  @param  caloHitList to receive the list of two dimensional calo hits for which no three dimensional hits have been created
     */
    void GetRemainingTwoDHits(const pandora::ParticleFlowObject *const pPfo, pandora::CaloHitList &caloHitList) const;

    /**
     *  @brief  Get the subset of a provided calo hit list corresponding to a specified hit type
     *
     *  @param  inputCaloHitList the input calo hit list
     *  @param  hitType the hit type to filter upon
     *  @param  outputCaloHitList to receive the output calo hit list
     */
    void FilterCaloHitsByType(const pandora::CaloHitList &inputCaloHitList, const pandora::HitType hitType, pandora::CaloHitList &outputCaloHitList) const;

    /**
     *  @brief  Add a specified list of three dimensional hits to a cluster in a pfo, creating the new cluster if required
     *
     *  @param  pPfo the address of the pfo
     *  @param  caloHitList the list of three dimensional hits
     *  @param  pCluster3D to receive the address of the (new) three dimensional cluster
     */
    void AddThreeDHitsToPfo(pandora::ParticleFlowObject *const pPfo, pandora::CaloHitList &caloHitList, pandora::Cluster *&pCluster3D) const;

    /**
     *  @brief  Create a new three dimensional hit from a two dimensional hit
     *
     *  @param  pCaloHit2D the address of the two dimensional calo hit, for which a new three dimensional hit is to be created
     *  @param  position3D the position vector for the new three dimensional calo hit
     *  @param  pCaloHit3D to receive the address of the new three dimensional calo hit
     */
    void CreateThreeDHit(pandora::CaloHit *pCaloHit2D, const pandora::CartesianVector &position3D, pandora::CaloHit *&pCaloHit3D) const;

    /**
     *  @brief  Check that a new three dimensional position is not unphysical
     *
     *  @param  position3D the position vector for the new three dimensional calo hit
     *
     *  @param  boolean  
     */
    bool CheckThreeDHit(const pandora::CartesianVector &position3D) const;

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Get the address of the three dimensional cluster in a specified pfo
     *
     *  @param  pPfo the address of the pfo
     *
     *  @return the address of the three dimensional cluster
     */
    pandora::Cluster *GetThreeDCluster(pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Create a new three dimensional cluster, using a list of provided three dimensional hits, and add it to a specified pfo
     *
     *  @param  caloHitList the list of three dimensional hits
     *  @param  pCluster to receive the address of the new cluster
     */
    void CreateThreeDCluster(const pandora::CaloHitList &caloHitList, pandora::Cluster *&pCluster) const;

    /**
     *  @brief  Get the list of 2D calo hits in a pfo for which 3D hits have and have not been created
     *
     *  @param  pPfo the address of the pfo
     *  @param  usedHits to receive the list of two dimensional calo hits for which three dimensional hits have been created
     *  @param  remainingHits to receive the list of two dimensional calo hits for which three dimensional hits have not been created
     */
    void SeparateTwoDHits(const pandora::ParticleFlowObject *const pPfo, pandora::CaloHitList &usedHits,
        pandora::CaloHitList &remainingHits) const;

    typedef std::vector<HitCreationBaseTool*> HitCreationToolList;
    HitCreationToolList     m_algorithmToolList;        ///< The algorithm tool list

    std::string             m_inputPfoListName;         ///< The name of the input pfo list
    std::string             m_outputCaloHitListName;    ///< The name of the output calo hit list
    std::string             m_outputClusterListName;    ///< The name of the output cluster list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDHitCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDHitCreationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_HIT_CREATION_ALGORITHM_H
