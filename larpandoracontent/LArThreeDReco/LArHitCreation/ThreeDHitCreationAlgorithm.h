/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h
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
     *  @brief  Get the subset of a provided calo hit vector corresponding to a specified hit type
     *
     *  @param  inputCaloHitVector the input calo hit vector
     *  @param  hitType the hit type to filter upon
     *  @param  outputCaloHitVector to receive the output calo hit vector
     */
    void FilterCaloHitsByType(const pandora::CaloHitVector &inputCaloHitVector, const pandora::HitType hitType, pandora::CaloHitVector &outputCaloHitVector) const;

    /**
     *  @brief  Add a specified vector of three dimensional hits to a cluster in a pfo, creating the new cluster if required
     *
     *  @param  pPfo the address of the pfo
     *  @param  caloHitVector the vector of three dimensional hits
     *  @param  pCluster3D to receive the address of the (new) three dimensional cluster
     */
    void AddThreeDHitsToPfo(const pandora::ParticleFlowObject *const pPfo, const pandora::CaloHitVector &caloHitVector, const pandora::Cluster *&pCluster3D) const;

    /**
     *  @brief  Create a new three dimensional hit from a two dimensional hit
     *
     *  @param  pCaloHit2D the address of the two dimensional calo hit, for which a new three dimensional hit is to be created
     *  @param  position3D the position vector for the new three dimensional calo hit
     *  @param  pCaloHit3D to receive the address of the new three dimensional calo hit
     */
    void CreateThreeDHit(const pandora::CaloHit *const pCaloHit2D, const pandora::CartesianVector &position3D, const pandora::CaloHit *&pCaloHit3D) const;

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

    /**
     *  @brief  Create a new three dimensional cluster, using a vector of provided three dimensional hits, and add it to a specified pfo
     *
     *  @param  caloHitVector the vector of three dimensional hits
     *  @param  pCluster to receive the address of the new cluster
     */
    void CreateThreeDCluster(const pandora::CaloHitVector &caloHitVector, const pandora::Cluster *&pCluster) const;

    /**
     *  @brief  Get the list of 2D calo hits in a pfo for which 3D hits have and have not been created
     *
     *  @param  pPfo the address of the pfo
     *  @param  usedHitVector to receive the vector of two dimensional calo hits for which three dimensional hits have been created
     *  @param  remainingHitVector to receive the vector of two dimensional calo hits for which three dimensional hits have not been created
     */
    void SeparateTwoDHits(const pandora::ParticleFlowObject *const pPfo, pandora::CaloHitVector &usedHitVector,
        pandora::CaloHitVector &remainingHitVector) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<HitCreationBaseTool*> HitCreationToolVector;
    HitCreationToolVector   m_algorithmToolVector;      ///< The algorithm tool vector

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
