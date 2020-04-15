/**
 *  @file   larpandoracontent/LArTest/TestAlgorithm.cc
 *
 *  @brief  Implementation of the test algorithm merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTest/TestAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"


using namespace pandora;

namespace lar_content
{


TestAlgorithm::TestAlgorithm() :
  m_clusterListName(),
  m_caloHitListName()
{
}

StatusCode  TestAlgorithm::Run()
{


    const ClusterList *pClusterList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);


    unsigned int maxHits(0);
    const Cluster *clusterToEnlarge(nullptr);

    for(const Cluster *const pCluster : *pClusterList)
    {
        if(pCluster->GetNCaloHits() > maxHits)
	{
   	    maxHits = pCluster->GetNCaloHits();
	    clusterToEnlarge = pCluster;
	}
    }

    ClusterList initialCluster;
    initialCluster.push_back(clusterToEnlarge);

    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &initialCluster, "INITIAL CLUSTER", BLACK, 2);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    for(const CaloHit *const pCaloHit : *pCaloHitList)
    {

      if(pCaloHit->GetParentAddress() == clusterToEnlarge)
	    continue;



      if(PandoraContentApi::IsAvailable(*this, pCaloHit))
      {
            std::cout << "CALO HIT IS AVAILABLE" << std::endl;
	    PandoraContentApi::AddToCluster(*this, clusterToEnlarge, pCaloHit);
      }
    }

    ClusterList finalCluster;
    finalCluster.push_back(clusterToEnlarge);

    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &finalCluster, "FINAL CLUSTER", RED, 2);

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    /*

    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);

    typedef std::vector<const BoxGap*> BoxGapVector;
    typedef std::vector<const LineGap*> LineGapVector;
    BoxGapVector boxGapVector;
    LineGapVector lineGapVector;

    for (const DetectorGap *const pDetectorGap : this->GetPandora().GetGeometry()->GetDetectorGapList())
    {
        const BoxGap *pBoxGap(nullptr);
	pBoxGap = dynamic_cast<const BoxGap*>(pDetectorGap);

	if(pBoxGap)
	    boxGapVector.push_back(pBoxGap);


        const LineGap *pLineGap(nullptr);
	pLineGap = dynamic_cast<const LineGap*>(pDetectorGap);

	if(pLineGap)
	    lineGapVector.push_back(pLineGap);

    }


    std::cout << "BOX GAP VECTOR SIZE: " << boxGapVector.size() << std::endl;
    std::cout << "LINE GAP VECTOR SIZE: " << lineGapVector.size() << std::endl;


    for (const LineGap *const pLineGap : lineGapVector)
    {

        float startX(pLineGap->GetLineStartX());
        float endX(pLineGap->GetLineEndX());

	CartesianVector startA(startX, 0, 0);
	CartesianVector startB(startX, 0, 1390);

	CartesianVector endA(endX, 0, 0);
	CartesianVector endB(endX, 0, 1390);
	
	PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &startA, &startB, "Start X", BLUE, 2, 1);
	PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &endA, &endB, "END X", BLUE, 2, 1);


	float m_maxXDistanceFromGap(3);

	CartesianVector startC(startX-m_maxXDistanceFromGap, 0, 0);
	CartesianVector startD(startX-m_maxXDistanceFromGap, 0, 1390);

	CartesianVector endC(endX+m_maxXDistanceFromGap, 0, 0);
	CartesianVector endD(endX+m_maxXDistanceFromGap, 0, 1390);

	PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &startC, &startD, "Start X", RED, 2, 1);
	PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &endC, &endD, "END X", RED, 2, 1);


    }


    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    return STATUS_CODE_SUCCESS;

}

StatusCode TestAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    
    return STATUS_CODE_SUCCESS;
}

}