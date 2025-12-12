#ifndef LAR_SHOWER_MERGING_ALGORITHM_H
#define LAR_SHOWER_MERGING_ALGORITHM_H 1

#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"


namespace lar_content
{

class ShowerMergingPfoMopUpAlgorithm : public PfoMopUpBaseAlgorithm
{

private:

pandora::StatusCode Run();

void GetSortedPfos(pandora::PfoVector &sortedPfos) const;

bool AreDirectionsAligned(const pandora::ParticleFlowObject *const pPfo1, const pandora::ParticleFlowObject *const pPfo2) const;

bool HaveSameVertex(const pandora::ParticleFlowObject *const pPfo1, const pandora::ParticleFlowObject *const pPfo2, bool &invert) const;

pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

pandora::StringVector m_inputPfoListNames;
float m_alignmentAngle;
float m_stubShowerSeperation;


};

}

#endif
