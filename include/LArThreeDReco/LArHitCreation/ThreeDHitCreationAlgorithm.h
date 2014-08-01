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

namespace lar
{

class HitCreationTool;

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
     *  @brief  Create a new three dimensional cluster, using a list of provided three dimensional hits
     *
     *  @param  pCaloHit2D the address of the two dimensional calo hit, for which a new three dimensional hit is to be created
     *  @param  position3D the position vector for the new three dimensional calo hit
     *  @param  pCaloHit3D to receive the address of the new three dimensional calo hit
     */
    void CreateThreeDHit(pandora::CaloHit *pCaloHit2D, const pandora::CartesianVector &position3D, pandora::CaloHit *&pCaloHit3D) const;

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

    std::string             m_inputPfoListName;         ///< The name of the input pfo list
    std::string             m_outputCaloHitListName;    ///< The name of the output calo hit list
    std::string             m_outputClusterListName;    ///< The name of the output cluster list

    typedef std::vector<HitCreationTool*> HitCreationToolList;
    HitCreationToolList     m_algorithmToolList;        ///< The algorithm tool list
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  HitCreationTool class
 */
class HitCreationTool : public pandora::AlgorithmTool
{
public:

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  pPfo the address of the pfo
     *  @param  inputTwoDHits the list of input two dimensional hits
     *  @param  newThreeDHits to receive the new three dimensional hits
     */
    virtual void Run(ThreeDHitCreationAlgorithm *pAlgorithm, const pandora::ParticleFlowObject *const pPfo, const pandora::CaloHitList &inputTwoDHits,
        pandora::CaloHitList &newThreeDHits) = 0;

protected:
    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and candidate fit positions from the other two views
     *
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  hitType1 the hit type identifying the first view
     *  @param  hitType2 the hit type identifying the second view
     *  @param  fitPositionList1 the candidate sliding fit position in the first view
     *  @param  fitPositionList2 the candidate sliding fit position in the second view
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    void GetBestPosition3D(const pandora::CaloHit *const pCaloHit2D, const pandora::HitType hitType1, const pandora::HitType hitType2,
        const pandora::CartesianPointList &fitPositionList1, const pandora::CartesianPointList &fitPositionList2, pandora::CartesianVector &position3D, float &chiSquared) const;

    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and candidate fit positions from the other two views
     *
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  hitType1 the hit type identifying the first view
     *  @param  hitType2 the hit type identifying the second view
     *  @param  fitPosition1 the candidate sliding fit position in the first view
     *  @param  fitPosition2 the candidate sliding fit position in the second view
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    void GetPosition3D(const pandora::CaloHit *const pCaloHit2D, const pandora::HitType hitType1, const pandora::HitType hitType2,
        const pandora::CartesianVector &fitPosition1, const pandora::CartesianVector &fitPosition2, pandora::CartesianVector &position3D, float &chiSquared) const;

    /**
     *  @brief  Get the three dimensional position using a provided two dimensional calo hit and candidate fit positions from the other two views
     *
     *  @param  pCaloHit2D address of the two dimensional calo hit
     *  @param  hitType the hit type identifying the other view
     *  @param  fitPosition the candidate sliding fit position in the other view
     *  @param  position3D to receive the three dimensional position
     *  @param  chiSquared to receive the chi squared value
     */
    void GetPosition3D(const pandora::CaloHit *const pCaloHit2D, const pandora::HitType hitType, const pandora::CartesianVector &fitPosition,
        pandora::CartesianVector &position3D, float &chiSquared) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool    m_useChiSquaredApproach;    ///< Whether to obtain y, z positions via chi2 approach, or projected position approach
    bool    m_useDeltaXCorrection;      ///< Whether to add a term to chi squared accounting for hit combination delta x values
    float   m_sigmaX;                   ///< Resolution in x dimension, used for delta x correction to chi squared
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDHitCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDHitCreationAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_THREE_D_HIT_CREATION_ALGORITHM_H
