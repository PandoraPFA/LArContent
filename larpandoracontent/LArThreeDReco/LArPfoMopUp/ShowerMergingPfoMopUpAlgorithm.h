/**
*  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/ShowerMergingPfoMopUpAlgorithm.h
*
*  @brief  Header file for  the shower merging pfo mop up algorithm class.
*  
*  $Log: $
*/
#ifndef LAR_SHOWER_MERGING_ALGORITHM_H
#define LAR_SHOWER_MERGING_ALGORITHM_H 1

#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief ShowerMergingPfoMopUpAlgorithm class 
 */
class ShowerMergingPfoMopUpAlgorithm : public PfoMopUpBaseAlgorithm
{
public:
    /**
     *  @brief Default constructor
     */
    ShowerMergingPfoMopUpAlgorithm();

private:
    pandora::StatusCode Run();
    
    /**
     *  @brief Get the sorted list of input particle flow objects (PFOs)
     *
     *  @param sortedPfos to receive the sorted list of input pfos 
     */
    void GetSortedPfos(pandora::PfoVector &sortedPfos) const;
    
     /**
      *  @brief Check if the directions of pfo1 and pfo2 are aligned 
      *
      *  @param pPfo1 Pointer to the first particle flow object
      *  @param pPfo2 Pointer to the second particle flow object
      *  return boolean 
      */
    bool AreDirectionsAligned(const pandora::ParticleFlowObject *const pPfo1, const pandora::ParticleFlowObject *const pPfo2) const;
   
    /**
     *  @brief Check if two pfos can be associated as a upstream stub and a downstream shower pair 
     *
     *  @param pPfo1 Pointer to the first particle flow object
     *  @param pPfo2 Pointer to the second particle flow object
     *  return boolean 
     */
    bool IsPfoVertexAssociated(const pandora::ParticleFlowObject *const pPfo1, const pandora::ParticleFlowObject *const pPfo2) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_inputPfoListNames;   ///< The input pfo list names 
    float m_alignmentAngle;                      ///< The maximum angular separation between two pfos
    float m_maxVtxPfosSeparation;                ///< The maximum distance between the neutrino vertex and start of the pfo 
    float m_stubShowerSeparationSquared;         ///< The square of the maximum distance between the upstream stub and the downstream shower 


};

}

#endif
