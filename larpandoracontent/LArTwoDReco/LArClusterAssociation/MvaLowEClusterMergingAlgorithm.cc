/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/MvaLowEClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of the mva lowe cluster merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Helpers/MCParticleHelper.h"
#include "Objects/MCParticle.h"
#include "Objects/CartesianVector.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArTrackShowerId/ShowerGrowingAlgorithm.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/MvaLowEClusterMergingAlgorithm.h"

#include <set>
#include <cmath>
#include <numeric>

using namespace pandora;

namespace lar_content
{
template <typename T>
MvaLowEClusterMergingAlgorithm<T>::MvaLowEClusterMergingAlgorithm() :
    m_trainingSetMode{false},
    m_enableProbability{true},
    m_minProbabilityCut{0.5f},
    m_event{-1},
    m_maxClusterFraction{0.9f},
    m_minNCaloHits{1},
    m_writeTree{true},
    m_upperHitThreshold{100},
    m_filePathEnvironmentVariable{"FW_SEARCH_PATH"},
    m_countHitsThreshold{0},
    m_vertexListName{"NeutrinoVertices3D"},
    m_contactThreshold{2},
    m_proximityThreshold{5},
    m_divisions{8},
    m_sectorTolerance{0.03}

{
}

//------------------------------------------------------------------------------------------------------------------------------------------
template <typename T>
MvaLowEClusterMergingAlgorithm<T>::~MvaLowEClusterMergingAlgorithm()
{
    if (m_writeTree)
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "RECREATE"));
         }

}

//------------------------------------------------------------------------------------------------------------------------------------------
template <typename T>
StatusCode MvaLowEClusterMergingAlgorithm<T>::Run()
{
    ++m_event;
    std::map<const Cluster*, ClusterVector> clustersToMerge;

    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        try
        {
            const ClusterList *pClusterList = nullptr;
            PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

            if (!pClusterList || pClusterList->empty())
            {
                if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo()) 
                    std::cout << "ClusterMergingAlgorithm: unable to find cluster list " << clusterListName << std::endl;

                continue;
            }

            PandoraContentApi::ParticleFlowObject::Metadata lowEClustersMetadata;
            this->EdgeHitComparer(pClusterList, clusterListName);  
	}
        
        catch (StatusCodeException &statusCodeException)
        {
            throw statusCodeException;
        }
        
	if (m_writeTree)
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "RECREATE"))
         }

    }

    return STATUS_CODE_SUCCESS;
}

//--------------------------------------------------------------------------------------------//

template <typename T>
const MCParticle* MvaLowEClusterMergingAlgorithm<T>::GetMCForCluster(const Cluster *const cluster, std::map<const Cluster*,
    const MCParticle*> &clusterToMCMap) const
{
    const MCParticle* clusterMC{nullptr};

    if (clusterToMCMap.count(cluster) > 0)
    {
        clusterMC = clusterToMCMap.at(cluster);
    }
    else
    {
        try
        {
            clusterMC = MCParticleHelper::GetMainMCParticle(cluster);
            clusterToMCMap[cluster] = clusterMC;
        }
        catch (StatusCodeException e)
        {
            std::cout << "Failed to get MC particle for cluster of " << cluster->GetOrderedCaloHitList().size()
                      << " : " << e.ToString() << std::endl;
        }
    }

    return clusterMC;
}

//---------------------------------------------------------------------------------------------//

template <typename T>
bool MvaLowEClusterMergingAlgorithm<T>::IsValidToUse(const Cluster *const cluster, std::map<const Cluster*, bool> &clusterIsUsed) const
{
    
    if (clusterIsUsed.count(cluster) > 0)
        return false;

    if (cluster->GetNCaloHits() < m_minNCaloHits)
        return false;

    if (cluster->GetNCaloHits() > m_upperHitThreshold)
        return false;

    return true;
}

//--------------------------------------------------------------------------------------------------------------------
template <typename T>
StatusCode MvaLowEClusterMergingAlgorithm<T>::EdgeHitComparer(const pandora::ClusterList *const pClusterList, const std::string &listName ) const
{
    std::map<const Cluster*, const MCParticle*> clusterToMCParticleMap;
    std::map<const Cluster*, bool> clusterIsUsed;
    std::map<const Cluster*, ClusterVector> clustersToMerge;
    int sigBac{0};
    int countHits{0};

    for (const Cluster *const pCluster : *pClusterList)
    {
	const Cluster *const cluster(*iter);
	countHits = countHits + (cluster->GetNCaloHits());
    }

    for (auto iter = pClusterList->begin(); iter != pClusterList->end(); ++iter)
    {
       const Cluster *const cluster(*iter);
       const MCParticle *clusterMC(this->GetMCForCluster(cluster, clusterToMCParticleMap));
       const int nClusters(pClusterList->size());
       
       if (!this->IsValidToUse(cluster, clusterIsUsed))
            continue;

	CaloHitList clusterEdgeHits;
	CartesianVector centroid11(0.f, 0.f, 0.f);
        CartesianVector centroid12(0.f, 0.f, 0.f);
        try {
       	centroid11 = cluster->GetCentroid(cluster->GetInnerPseudoLayer());
        centroid12 = cluster->GetCentroid(cluster->GetOuterPseudoLayer());
        }
        catch (...) {
        CaloHitList hits;
        cluster->GetOrderedCaloHitList().FillCaloHitList(hits);
        centroid11 = hits.front()->GetPositionVector();
        centroid12 = hits.front()->GetPositionVector();
        }

        const CartesianVector centroid1((centroid11 + centroid12) * 0.5f);
	
	this->EdgeHitFinder(cluster, clusterEdgeHits); 

        const unsigned int clusterNHits{cluster->GetNCaloHits()};
        const unsigned long int clusterNEdgeHits{clusterEdgeHits.size()};

	for (auto iter2 = std::next(iter) ; iter2 != pClusterList->end() ; ++iter2)
	{
	    const Cluster *const otherCluster(*iter2);
	    const MCParticle *otherClusterMC(this->GetMCForCluster(otherCluster, clusterToMCParticleMap));
            sigBac = 0;

            if (clusterMC == otherClusterMC)
            {
                   sigBac = 1;
            }

	    CaloHitList otherClusterEdgeHits;
	    CartesianVector centroid21(0.f, 0.f, 0.f);
	    CartesianVector centroid22(0.f, 0.f, 0.f);
	    try {
	    centroid21 = otherCluster->GetCentroid(otherCluster->GetInnerPseudoLayer());
            centroid22 = otherCluster->GetCentroid(otherCluster->GetOuterPseudoLayer());
	    }
	    catch (...) {
	    CaloHitList hits;
            cluster->GetOrderedCaloHitList().FillCaloHitList(hits);
            centroid21 = hits.front()->GetPositionVector();
            centroid22 = hits.front()->GetPositionVector();
	    }
	    const CartesianVector centroid2((centroid21 + centroid22) * 0.5f);
	    const float centroidSeparation{std::sqrt(centroid1.GetDistanceSquared(centroid2))};
	    const CartesianVector centroidVector(centroid2.GetX() - centroid1.GetX(), 0.f, centroid2.GetZ() - centroid1.GetZ());
	    
	    const VertexList *pVertexList{nullptr};
            float vtxClusterAngle{-1};
	    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,                       vertexListName, pVertexList));
	    
	    if (pVertexList && !pVertexList->empty())
            {
                const Vertex* pVertex{nullptr}; 
   	        pVertex = pVertexList->front();
	        const HitType hitType(LArClusterHelper::GetClusterHitType(cluster));
                CartesianVector vertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition() , hitType));

		const CartesianVector vtxToCluster(centroid1.GetX() - vertexPosition.GetX(), 0.f, centroid1.GetZ() - vertexPosition.GetZ());
                const CartesianVector vtxToOtherCluster(centroid2.GetX() - vertexPosition.GetX(), 0.f, centroid2.GetZ() - vertexPosition.GetZ());
		const float D1(vtxToCluster.GetMagnitude());
		const float D2(vtxToOtherCluster.GetMagnitude());
	        if ( D1*D2 != 0)
		{
		    vtxClusterAngle = (vtxToCluster.GetDotProduct(vtxToOtherCluster) / (D1 * D2));
		}
	    }

	    this->EdgeHitFinder(otherCluster, otherClusterEdgeHits);

	    const unsigned int otherClusterNHits{otherCluster->GetNCaloHits()};
      	    const unsigned long int otherClusterNEdgeHits{otherClusterEdgeHits.size()};
	    
            CaloHitList largestList, smallestList;
	    CartesianVector vectorTool(0.f, 0.f, 0.f);
	    if ( clusterEdgeHits.size() > otherClusterEdgeHits.size() )
            {
		    largestList = clusterEdgeHits;
		    smallestList = otherClusterEdgeHits;
		    vectorTool = centroid2;
            }
	    else 
            {
                largestList = otherClusterEdgeHits;
                smallestList = clusterEdgeHits;
		vectorTool = centroid1;
	    }

            std::vector<float> distanceDistribution;

	    float finalDistance{0}, maxEdgeHitSeparation{-std::numeric_limits<float>::max()}, minEdgeHitSeparation{std::numeric_limits<float>::max()}, adc1{0}, adc2{0};
            int contact{0}, proximity{0};;

	    for ( const CaloHit *const caloHit1 : *largestList)
            {
                  CartesianVector hit1Position(caloHit1->GetPositionVector());
                  float nearestHit{-1};
		  adc1 += caloHit1->GetInputEnergy();

		  for (const CaloHit *const caloHit2 : *smallestList)
                  {
		      CartesianVector hit2Position(caloHit2->GetPositionVector());
                      float distance{std::sqrt(hit1Position.GetDistanceSquared(hit2Position))};
 		      CartesianVector vector1(hit2Position.GetX() - hit1Position.GetX(), 0.f, hit2Position.GetZ() - hit1Position.GetZ());
                      double smallestAngle{1000};
                      CartesianVector tangent(centroidVector.GetZ(), 0.f, -centroidVector.GetX());
		      const float c(tangent.GetX() * vectorTool.GetX() + tangent.GetZ() * vectorTool.GetZ());
		      const float intercept(tangent.GetX() * hit2Position.GetX() + tangent.GetZ() * hit2Position.GetZ());
                      adc2 += caloHit2->GetInputEnergy();

		      if ( intercept < c)
                          continue;

		      if ( distance > maxEdgeHitSeparation )
		      {
		          maxEdgeHitSeparation = distance;
		      }

		      if ( distance < minEdgeHitSeparation)
	              {
                          minEdgeHitSeparation = distance;
	              }

		      if ( distance < nearestHit )
		      {
		          nearestHit = distance;

		      }
		      
		      if (std::abs( vector1.GetDotProduct(centroidVector) ) > smallestAngle)
		         continue;
                    
                      if (caloHit1->GetPseudoLayer() == caloHit2->GetPseudoLayer())
	              {
			      if ( distance < m_proximityThreshold)
			      { 
				      ++proximity;
				      if ( distance < m_contactThreshold)
				      {
					      ++contact;
					      --proximity;
				      }
		               }

	              }
		      smallestAngle = std::abs( vector1.GetDotProduct(centroidVector));
		      finalDistance = distance;
                  }
		  distanceDistribution.push_back(finalDistance);
	    }

            float combinedNClusterHits(clusterNHits + otherClusterNHits);
	    float combinedNClusterEdgeHits(clusterNEdgeHits + otherClusterNEdgeHits);
 
            const double avgDistance{!distanceDistribution.empty() ? std::accumulate(distanceDistribution.begin(), distanceDistribution.end(), 0.0) / distanceDistribution.size() : -1.0};

	    std::vector<std::string> featureOrder;
	    const bool mcMatch(clusterMC == otherClusterMC);

	    featureOrder.emplace_back("VertexClusterAngle");
            featureOrder.emplace_back("MinEdgeHitSeparation");
            featureOrder.emplace_back("Cluster1NHits");
            featureOrder.emplace_back("Cluster2NHits");
            featureOrder.emplace_back("NHitsInContact");
            featureOrder.emplace_back("NHitsInProximity");
            featureOrder.emplace_back("CentroidSeparation");
            featureOrder.emplace_back("AvgDistance");
	    featureOrder.emplace_back("CentroidVectorX");
	    featureOrder.emplace_back("CentroidVectorZ");
	    featureOrder.emplace_back("Cluster1ADC");
	    featureOrder.emplace_back("Cluster2ADC");

	    if (m_trainingSetMode)
            {
                            
                LArMvaHelper::MvaFeatureVector featureVector;
                
		featureVector.emplace_back(static_cast<double>(vtxClusterAngle));
		featureVector.emplace_back(static_cast<double>(minEdgeHitSeparation));
                featureVector.emplace_back(static_cast<double>(clusterNHits));
                featureVector.emplace_back(static_cast<double>(otherClusterNHits));
                featureVector.emplace_back(static_cast<double>(contact));
		featureVector.emplace_back(static_cast<double>(proximity));
 		featureVector.emplace_back(static_cast<double>(centroidSeparation));
		featureVector.emplace_back(static_cast<double>(avgDistance));
		featureVector.emplace_back(static_cast<double>(centroidVector.GetX()));
		featureVector.emplace_back(static_cast<double>(centroidVector.GetZ()));
		featureVector.emplace_back(static_cast<double>(adc1));
		featureVector.emplace_back(static_cast<double>(adc2));

                LArMvaHelper::ProduceTrainingExample(m_trainingOutputFile, mcMatch, featureVector);
 
            }
	   
	    if (! m_trainingSetMode)
	    {	    
	        LArMvaHelper::MvaFeatureMap featureMap;
	    	    
                featureMap["VertexClusterAngle"] = vtxClusterAngle;
	        featureMap["MinEdgeHitSeparation"] = minEdgeHitSeparation;
                featureMap["Cluster1NHits"] = clusterNHits;
	        featureMap["Cluster2NHits"] = otherClusterNHits;
	        featureMap["NHitsInContact"] = contact;
                featureMap["NHitsInProximity"] = proximity;
	        featureMap["CentroidSeparation"] = centroidSeparation;
	        featureMap["AvgDistance"] = avgDistance;
		featureMap["CentroidVectorX"] = centroidVector.GetX();
		featureMap["CentroidVectorZ"] = centroidVector.GetZ();
		featureMap["Cluster1ADC"] = adc1;
		featureMap["Cluster2ADC"] = adc2;
           
	        const bool areClustersToMerge(this->ClusterTool(featureOrder, featureMap));

                if (areClustersToMerge && countHits > m_countHitsThreshold)
                {
                   clusterIsUsed[cluster] = true;
                   clusterIsUsed[otherCluster] = true;
                   clustersToMerge[cluster].push_back(otherCluster);
                }
	    }
            
	    if (m_writeTree)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "EventNumber", m_event));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SignalvBackground", sigBac));
                if (pVertexList &&  ! pVertexList->empty()){
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "VertexClusterAngle", vtxClusterAngle));}
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "maxEdgeHitSeparation", maxEdgeHitSeparation));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NHitsInContact", contact));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NHitsInProximity", proximity));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "minEdgeHitSeparation", minEdgeHitSeparation));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "CentroidSeparation", centroidSeparation));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NClusters", nClusters));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "CentroidVectorX", centroidVector.GetX()));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "CentroidVectorZ", centroidVector.GetZ()));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TotalNClusterHits", (int) combinedNClusterHits));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TotalNClusterEdgeHits",(int) combinedNClusterEdgeHits));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
                
	    }
	     
        }
    }

    for (auto clusterToMergePair : clustersToMerge)
    {
        const Cluster *currentCluster = clusterToMergePair.first;
        const auto clusters = clusterToMergePair.second;

        for (auto clusterToMerge : clusters)
        {
            if (! clusterToMerge->IsAvailable())
                continue;

            try
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, currentCluster, clusterToMerge, listName, listName));
            } catch (StatusCodeException) {}
        }
    }


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool MvaLowEClusterMergingAlgorithm<T>::ClusterTool(std::vector<std::string> featureOrder, LArMvaHelper::MvaFeatureMap featureMap) const
{

    if (!m_enableProbability)
    {
        return LArMvaHelper::Classify(m_mva, featureOrder, featureMap);
    }
    else
    {
	return (LArMvaHelper::CalculateProbability(m_mva, featureOrder, featureMap) > m_minProbabilityCut);
    }
}

//------------------------------------------------------------------------------------------------------------------

template <typename T>
const CaloHitList MvaLowEClusterMergingAlgorithm<T>::EdgeHitFinder(const pandora::Cluster *const cluster, pandora::CaloHitList &clusterEdgeHits) const
{
        CartesianVector vector1(0.f,0.f,0.f);
        CartesianVector vector2(0.f,0.f,0.f);
        CaloHitList clusterCaloHits;
        CartesianVector clusterCentroid(0.f,0.f,0.f);
        std::map<const Cluster*, bool> clusterIsUsed;
        cluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHits);
        const CartesianVector centroid1(cluster->GetCentroid(cluster->GetInnerPseudoLayer()));
        const CartesianVector centroid2(cluster->GetCentroid(cluster->GetOuterPseudoLayer()));
        const CartesianVector centroid((centroid1 + centroid2) * 0.5f);
        std::map<const CaloHit*, bool> hitIsUsed;


        for (int i = 0 ; i < m_divisions ; i++)
        {
                float pi(3.14159265);
                float phi((2 * pi * i) / m_divisions);
                CartesianVector vec(std::cos(phi), std::sin(phi),0.f);
                CaloHitList sectorHits;
                sectorHits.clear();
                float maxMag{0};

                for (const CaloHit *pCaloHit: clusterCaloHits)
                {
                                            if (hitIsUsed.count(pCaloHit) == 0)
                        {
                            const float x(pCaloHit->GetPositionVector().GetX());
                            const float z(pCaloHit->GetPositionVector().GetZ());
                            CartesianVector centroidToHitVec(x - centroid.GetX(), 0.f, z - centroid.GetZ());
                            float mag(centroidToHitVec.GetMagnitude());
                            CartesianVector normCentroidVec((x - centroid.GetX()) / mag, 0.f , (z - centroid.GetZ()) / mag );
                            float dotProduct(vec.GetDotProduct(normCentroidVec));

                            if (std::abs(dotProduct) < m_sectorTolerance )
                            {
                                if (maxMag < mag)
                                {
                                        maxMag = mag;
                                        sectorHits.clear();
                                        hitIsUsed[pCaloHit] = true;
                                        sectorHits.push_back(pCaloHit);
                                }
                            }
                         }



                }
                for (const CaloHit *pCaloHit: sectorHits)
                {
                       clusterEdgeHits.push_back(pCaloHit);
                }
                
		if (clusterCaloHits.size() < clusterEdgeHits.size())
                {
                        std::cout << "Error in Edge Hit Algorithm: Incorrect edge hit classification !!!!" << std::endl;
                }
        }
        return clusterEdgeHits;
}

//-------------------------------------------------------------------------------------------------

template <typename T>
double MvaLowEClusterMergingAlgorithm<T>::Angle(const CartesianVector vector1, const CartesianVector vector2) const
{
    const double dx{vector1.GetX() - vector2.GetX()};
    const double dz{vector1.GetZ() - vector2.GetZ()};
    double angle{tan(dx / dz)};

    return angle;
}

//------------------------------------------------------------------------------------------------------------------------------------------
template <typename T>
StatusCode MvaLowEClusterMergingAlgorithm<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingSetMode", m_trainingSetMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaFileName", m_mvaFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MvaName", m_mvaName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EnableProbability", m_enableProbability));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinProbabilityCut", m_minProbabilityCut));
 
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));


    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UpperHitThreshold", m_upperHitThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinNCaloHits", m_minNCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CountHitsThreshold", m_countHitsThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ContactThreshold", m_contactThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ProximityThreshold", m_proxmityThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "EdgeHitToolDivisions", m_divisions));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SectorTolerance", m_sectorTolerance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    if (m_writeTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
    }   
    
    
    
    if (m_trainingSetMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
    }
    else
    {
        if (m_mvaFileName.empty() || m_mvaName.empty())
        {
            std::cout << "MvaLowEClusterMergingAlgorithm: MvaFileName and MvaName must be set if in classification mode " << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }

        const std::string fullMvaFileName(LArFileHelper::FindFileInPath(m_mvaFileName, m_filePathEnvironmentVariable));
        m_mva.Initialize(fullMvaFileName, m_mvaName);

    }


    return STATUS_CODE_SUCCESS;
}
//--------------------------------------------------------------------------------------------------------------------------------------

template class MvaLowEClusterMergingAlgorithm<AdaBoostDecisionTree>;


} // namespace lar_content
