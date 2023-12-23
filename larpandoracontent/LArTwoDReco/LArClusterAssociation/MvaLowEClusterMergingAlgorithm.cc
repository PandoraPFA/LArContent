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
    m_minNCaloHits{1},
    m_writeTree{true},
    m_upperHitThreshold{100},
    m_filePathEnvironmentVariable{"FW_SEARCH_PATH"},
    m_countHitsThreshold{0}

{
}

//------------------------------------------------------------------------------------------------------------------------------------------
template <typename T>
MvaLowEClusterMergingAlgorithm<T>::~MvaLowEClusterMergingAlgorithm()
{
    if (m_writeTree)
        {
           // PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "RECREATE"));
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName2.c_str(), m_fileName.c_str(), "RECREATE"));
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
        
	if (!clustersToMerge.empty())
        {
	       	//PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, m_clustersToMerge, clustersToMerge));
		std::cout << std::endl;
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
    const MCParticle* clusterMC = nullptr;

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

   // if (!cluster->IsAvailable())
   //     return false;
    
    if (cluster->GetNCaloHits() < m_minNCaloHits)
        return false;

    if (cluster->GetNCaloHits() > m_upperHitThreshold)
        return false;

    if (clusterIsUsed.count(cluster) > 0)
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
    int sigBac;
    std::cout << listName << std::endl;
    int countHits{0};

    for (auto iter = pClusterList->begin(); iter != pClusterList->end(); ++iter)
    {
	const Cluster *const cluster(*iter);
	countHits = countHits + (cluster->GetNCaloHits());
    }
    std::cout << " ************************ HITS **************** " << countHits << std::endl;

    for (auto iter = pClusterList->begin(); iter != pClusterList->end(); ++iter)
    {
       const Cluster *const cluster(*iter);
       const MCParticle *clusterMC(this->GetMCForCluster(cluster, clusterToMCParticleMap));
       //const MCParticle *const clusterMC(MCParticleHelper::GetMainMCParticle(cluster));
       const int nClusters(pClusterList->size());
       std::cout << "Number of clusters to check: " << nClusters << std::endl;	
       
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

        //const CartesianVector centroid11(cluster->GetCentroid(cluster->GetInnerPseudoLayer()));
        //const CartesianVector centroid12(cluster->GetCentroid(cluster->GetOuterPseudoLayer()));
        const CartesianVector centroid1((centroid11 + centroid12) * 0.5f);
	
	this->EdgeHitFinder(cluster, clusterEdgeHits); 

        const unsigned int clusterNHits{cluster->GetNCaloHits()};
        const unsigned long int clusterNEdgeHits{clusterEdgeHits.size()};

	std::cout << "Running 1st Loop" << std::endl;
	for (auto iter2 = std::next(iter) ; iter2 != pClusterList->end() ; ++iter2)
	{
	    const Cluster *const otherCluster(*iter2);
	    const MCParticle *otherClusterMC(this->GetMCForCluster(otherCluster, clusterToMCParticleMap));
            sigBac = 0;

            if (clusterMC == otherClusterMC)
            {
                   sigBac = 1;

            }

	    std::cout << "********* CLUSTER MC ********" << std::endl;
            std::cout << "Cluster 1: " << clusterMC <<std::endl;
            std::cout << "Cluster 2: " << otherClusterMC << std::endl;
            std::cout << "Bool: " << sigBac << std::endl;    
	    std::cout << "*****************************" << std::endl;

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
	    const float centroidSeparation{Distance(centroid1, centroid2)};
	    const CartesianVector centroidVector(centroid2.GetX() - centroid1.GetX(), 0.f, centroid2.GetZ() - centroid1.GetZ());
	    
	    std::string vertexListName("NeutrinoVertices3D");
	    const VertexList *pVertexList = nullptr;
            float vtxClusterAngle{-1};
	    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this,                       vertexListName, pVertexList));
	    
	    if (pVertexList && !pVertexList->empty())
            {
                const Vertex* pVertex{nullptr}; 
   	        pVertex = pVertexList->front();
	        const HitType hitType(LArClusterHelper::GetClusterHitType(cluster));
                CartesianVector vertexPosition(0.f, 0.f, 0.f);

		vertexPosition = LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition() , hitType);
		const CartesianVector vtxToCluster(centroid1.GetX() - vertexPosition.GetX(), 0.f, centroid1.GetZ() - vertexPosition.GetZ());
                const CartesianVector vtxToOtherCluster(centroid2.GetX() - vertexPosition.GetX(), 0.f, centroid2.GetZ() - vertexPosition.GetZ());
		const float D1(std::sqrt((vtxToCluster.GetX() * vtxToCluster.GetX()) + (vtxToCluster.GetZ() * vtxToCluster.GetZ())));
		const float D2(std::sqrt((vtxToOtherCluster.GetX() * vtxToOtherCluster.GetX()) + (vtxToOtherCluster.GetZ() * vtxToOtherCluster.GetZ())));
	        vtxClusterAngle = (vtxToCluster.GetDotProduct(vtxToOtherCluster) / (D1 * D2));
	     
	     	std::cout << "Vertex Cluster Angle: " << vtxClusterAngle << std::endl;
	    }
	    std::cout << "Vertex List: " << pVertexList << std::endl;
	    std::cout << "Centroid Separation is: " << centroidSeparation << std::endl;
            std::cout << "Centroid Vector is: " << centroidVector << std::endl;
            std::cout << "Cluster1 pos: " << centroid1 << std::endl;
	    std::cout << "Cluster2 pos: " << centroid2 << std::endl;

	    this->EdgeHitFinder(otherCluster, otherClusterEdgeHits);
            // cluster->GetOrderedCaloHitList().FillCaloHitList(clusterEdgeHits);
	    // otherCluster->GetOrderedCaloHitList().FillCaloHitList(otherClusterEdgeHits);

	    const unsigned int otherClusterNHits{otherCluster->GetNCaloHits()};
      	    const unsigned long int otherClusterNEdgeHits{otherClusterEdgeHits.size()};
	    
	   // std::cout << "Comparing Cluster 1 and 2, Cluster 1 has " << clusterNEdgeHits << " edge hits whilst Cluster 2 has " <<                         otherClusterNEdgeHits << std::endl;
	    std::cout << "Running 2nd Loop" << std::endl;
            CaloHitList largestList;
	    CaloHitList smallestList;
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
            std::vector<double>  angleDistribution;
	    std::vector<float> nearestHitDistribution;
	    std::vector<float> nearestHitAngleDistribution;

	    float finalDistance{0};
            double finalAngle{0};
            float maxEdgeHitSeparation{0};
	    float minEdgeHitSeparation{10000};
            int contact{0};
	    int proximity{0};
	    float contactThreshold{2};
	    float proximityThreshold{5}; 
	    float adc1(0);
	    float adc2(0);

	    for (auto iter3 = largestList.begin() ; iter3 != largestList.end() ; ++iter3)
            {
 		  std::cout << "Running 3rd Loop" << std::endl;
		  std::cout << largestList.size() << " hits" << std::endl;
		  const CaloHit *const caloHit1(*iter3);
                  CartesianVector hit1Position(caloHit1->GetPositionVector());
                  float nearestHit{-1};
		  double nearestHitAngle{-1};
		  float hitAdc1(caloHit1->GetInputEnergy());
		  adc1 += hitAdc1;
		  

		  for (auto iter4 = smallestList.begin() ; iter4 != smallestList.end() ; ++iter4)
                  {
		      std::cout << "Running 4th Loop" << std::endl;
		      std::cout << smallestList.size() << " hits" << std::endl;
                      const CaloHit *const caloHit2(*iter4);
		      CartesianVector hit2Position(caloHit2->GetPositionVector());
                      float distance1{Distance(hit1Position, hit2Position)};
 		      CartesianVector vector1(hit2Position.GetX() - hit1Position.GetX(), 0.f, hit2Position.GetZ() - hit1Position.GetZ());
                      double smallestAngle{1000};
                      CartesianVector tangent(centroidVector.GetZ(), 0.f, -centroidVector.GetX());
		      const float c(tangent.GetX() * vectorTool.GetX() + tangent.GetZ() * vectorTool.GetZ());
		      const float intercept(tangent.GetX() * hit2Position.GetX() + tangent.GetZ() * hit2Position.GetZ());
                      std::cout << c << intercept << std::endl;
		      std::cout << vector1.GetDotProduct(centroidVector) << std::endl;
		      float hitAdc2(caloHit2->GetInputEnergy());
                      adc2 += hitAdc2;
                      
                      nearestHitDistribution.push_back(nearestHit);
                      nearestHitAngleDistribution.push_back(nearestHitAngle);

		      if ( intercept < c)
                          continue;


		      if ( distance1 > maxEdgeHitSeparation )
		      {
		          maxEdgeHitSeparation = distance1;

		      }
		      

		      if ( distance1 < minEdgeHitSeparation)
	              {
                          minEdgeHitSeparation = distance1;
	              }

		      if ( distance1 < nearestHit )
		      {
		          nearestHit = distance1;
			  nearestHitAngle = Angle(vector1, CartesianVector (1., 0., 0.));

		      }
		      
		      if (std::abs( vector1.GetDotProduct(centroidVector) ) > smallestAngle)
		         continue;
                    
                      if (caloHit1->GetPseudoLayer() == caloHit2->GetPseudoLayer())
	              {
			      if ( distance1 < proximityThreshold)
			      { 
				      ++proximity;
				      if ( distance1 < contactThreshold)
				      {
					      ++contact;
					      --proximity;
				      }
		               }

	              }
		      smallestAngle = std::abs( vector1.GetDotProduct(centroidVector));
		      finalDistance = distance1;
		      finalAngle = smallestAngle;           
                      std::cout << "final parameters: UPDATE ***" << std::endl;

                  }

		  distanceDistribution.push_back(finalDistance);
		  angleDistribution.push_back(finalAngle);

            //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "DistanceBetweenHits", finalDistance));
            //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "AngleBetweenHits", finalAngle));

	    }

	     long unsigned int i;
	     float varDistance{0};
	     float varAngle{0};
	     float varNearestHitDistance{0};
	     float varNearestHitAngle{0};
	     const double avgDistance{!distanceDistribution.empty() ? std::accumulate(distanceDistribution.begin(), distanceDistribution.end(), 0.0) / distanceDistribution.size() : -1.0};
             const double avgAngle{!angleDistribution.empty() ? std::accumulate(angleDistribution.begin(), angleDistribution.end(), 0.0) / angleDistribution.size() : -1.0};
	     const float avgNearestHitDistance(std::accumulate(nearestHitDistribution.begin(), nearestHitDistribution.end(), 0.0) / nearestHitDistribution.size());
             const float avgNearestHitAngle(std::accumulate(nearestHitAngleDistribution.begin(), nearestHitAngleDistribution.end(), 0.0) / nearestHitAngleDistribution.size());

	     for (i = 0 ; i < distanceDistribution.size() ; ++i)
             {
                 varDistance += ((avgDistance - distanceDistribution[i]) * (avgDistance - distanceDistribution[i]));
	     }
	     for (i = 0 ; i < angleDistribution.size() ; ++i)
             {
                 varAngle += ((avgAngle - angleDistribution[i]) * (avgAngle - angleDistribution[i]));
             }
             for (i = 0 ; i < nearestHitDistribution.size() ; ++i)
             {
                 varNearestHitDistance += ((avgNearestHitDistance - nearestHitDistribution[i]) * (avgNearestHitDistance - nearestHitDistribution[i]));
             }
             for (i = 0 ; i < nearestHitAngleDistribution.size() ; ++i)
             {
                 varNearestHitAngle += ((avgNearestHitAngle - nearestHitAngleDistribution[i]) * (avgNearestHitAngle - nearestHitAngleDistribution[i]));
             }

             varAngle = std::sqrt(varAngle / angleDistribution.size());
             varDistance = std::sqrt(varDistance / distanceDistribution.size());
             varNearestHitDistance = std::sqrt(varNearestHitDistance / nearestHitDistribution.size());
	     varNearestHitAngle = std::sqrt(varNearestHitAngle / nearestHitAngleDistribution.size());

            float combinedNClusterHits(clusterNHits + otherClusterNHits);
	    float combinedNClusterEdgeHits(clusterNEdgeHits + otherClusterNEdgeHits);
  
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
                   // lowEClustersMetadata = m_merge;
                   clusterIsUsed[cluster] = true;
                   clusterIsUsed[otherCluster] = true;
                   clustersToMerge[cluster].push_back(otherCluster);
                }
	    }
            
	    if (m_writeTree)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "EventNumber", m_event));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "SignalvBackground", sigBac));
                if (pVertexList &&  ! pVertexList->empty()){
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "VertexClusterAngle", vtxClusterAngle));}
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "maxEdgeHitSeparation", maxEdgeHitSeparation));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "NHitsInContact", contact));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "NHitsInProximity", proximity));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "minEdgeHitSeparation", minEdgeHitSeparation));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "AverageDistance", avgDistance));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "AverageAngle", avgAngle));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "AverageNearestHitDistance", avgNearestHitDistance));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "AverageNearestHitAngle", avgNearestHitAngle));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "DistanceVariance", varDistance));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "AngleVariance", varAngle));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "NearestHitDistanceVariance", varNearestHitDistance));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "NearestHitAngleVariance", varNearestHitAngle));

                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "CentroidSeparation", centroidSeparation));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "NClusters", nClusters));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "CentroidVectorX", centroidVector.GetX()));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "CentroidVectorZ", centroidVector.GetZ()));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "TotalNClusterHits", (int) combinedNClusterHits));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName2.c_str(), "TotalNClusterEdgeHits",(int) combinedNClusterEdgeHits));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName2.c_str()));
                
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
	std::cout << "CUT VALUE USED: " << m_minProbabilityCut << std::endl;
	std::cout << "PROBABILITY: " << LArMvaHelper::CalculateProbability(m_mva, featureOrder, featureMap) << std::endl;
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
        float divisions(8);
        std::map<const CaloHit*, bool> hitIsUsed;


        for (int i(0) ; i < divisions ; i++)
        {
                float pi(3.14159265);
                float phi((2 * pi * i) / divisions);
                CartesianVector vec(std::cos(phi), std::sin(phi),0.f);
                //std::cout << "Sweep Angle: " << vec << std::endl;
                CaloHitList sectorHits;
                sectorHits.clear();
                float maxMag{0};
                //std::cout << "CLUSTER CALO HIT SIZE: " << clusterCaloHits.size() << std::endl;

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
                            //std::cout << " Centroid1: " << centroid1.GetX() << " || " << centroid1.GetZ() <<  " Centroid2: " << centroid2.                              GetX() << " || " << centroid2.GetZ() << pMCParticle << listName << std::endl;
                            //  std::cout << "CENTROID: " << centroid.GetX() << " || " << centroid.GetZ() << std::endl;

                            if (std::abs(dotProduct) < 0.03 )
                            {
                                //std::cout << "Max Mag: " << maxMag << std::endl;
                                if (maxMag < mag)
                                {
                                        maxMag = mag;
                                        sectorHits.clear();
                                        hitIsUsed[pCaloHit] = true;
                                        sectorHits.push_back(pCaloHit);
                                        //std::cout << "DOT PRODUCT: " << dotProduct << "Magnitude: " << mag << std::endl;
                                }
                            }
                         }



                }
                //std::cout << "Last Hit to add: " << sectorHits.size() << std::endl;
                for (const CaloHit *pCaloHit: sectorHits)
                {
                       clusterEdgeHits.push_back(pCaloHit);
                }
                //std::cout << "Number of hits found: " << clusterEdgeHits.size() << std::endl;
                if (clusterCaloHits.size() < clusterEdgeHits.size())
                {
                        std::cout << "Error in Edge Hit Algorithm: Incorrect edge hit classification !!!!" << std::endl;

                }
        }




        return clusterEdgeHits;
}

//----------------------------------------------------------------------------------------------------

template <typename T>
float MvaLowEClusterMergingAlgorithm<T>::Distance(const CartesianVector vector1, const CartesianVector vector2) const
{
    const float dx{std::abs(vector1.GetX() - vector2.GetX())};
    const float dz{vector1.GetZ() - vector2.GetZ()};
    float distance{std::sqrt(dx * dx + dz * dz)};

    return distance;
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));

    if (m_writeTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));
        std::cout << "Input Tree Name for Cluster Merging: " << m_treeName << std::endl;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName2", m_treeName2));
        std::cout << "Input Tree Name for Edge Hit Merging: " << m_treeName2 << std::endl;

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
