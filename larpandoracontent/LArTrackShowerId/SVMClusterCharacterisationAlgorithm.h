/**
 *  @file   larpandoracontent/LArTrackShowerId/SVMClusterCharacterisationAlgorithm.h
 *
 *  @brief  Header file for the cluster characterisation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SVM_CLUSTER_CHARACTERISATION_ALGORITHM_H
#define LAR_SVM_CLUSTER_CHARACTERISATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "larpandoracontent/LArHelpers/LArSVMHelper.h"

namespace lar_content
{

/**
 *  @brief  SVMClusterCharacterisationAlgorithm class
 */
class SVMClusterCharacterisationAlgorithm : public pandora::Algorithm
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
    SVMClusterCharacterisationAlgorithm();


    typedef SVMFeatureTool<const SVMClusterCharacterisationAlgorithm *const, const pandora::Cluster * const>  ClusterCharacterisationFeatureTool; 
    
    //--------------------------------------------------------------------------------------------------------------------------------------
    
    /**
     *  @brief Event feature info class
     */
    class ClusterFeatureInfo
    {
    public:
      /**
       *  @brief  Constructor
       * 
       *  @param  showerFitWidth 
       *  @param  showerFitWidth 
       */
      ClusterFeatureInfo(/*const float hitType, const float nHits, */const float nGoodHits, /*const float straightLineLength, */const float straightLineLengthLarge,  
        const float showerFitWidth,  const float showerFitGapLength, const float diffWithStraigthLine, 
        const float widthDirectionX, const float nNearbyClusters, const float mipEnergy);
    
 float m_hitType;  
 float           m_nHits;             ///<  
      float           m_nGoodHits;             ///<      
      float           m_straightLineLength;                   ///
    float m_straightLineLengthLarge;                  
      float           m_showerFitWidth;        ///< 
      float           m_showerFitGapLength;         ///< 
    float m_diffWithStraigthLine;
    float m_widthDirectionX;
    float m_nNearbyClusters;
    float m_mipEnergy;
    
    };
    
     typedef std::map<const pandora::Cluster * const, ClusterFeatureInfo> ClusterFeatureInfoMap;
  

private:
    pandora::StatusCode Run();


  void PopulateClusterFeatureInfoMap(const pandora::Cluster *const pCluster, ClusterFeatureInfoMap &clusterFeatureInfoMap) const;
  SupportVectorMachine::DoubleVector GenerateFeatureList(const pandora::Cluster *const pCluster, ClusterFeatureInfoMap &clusterFeatureInfoMap);
  
    ClusterCharacterisationFeatureTool::FeatureToolVector m_featureToolVector; ///< The feature tool map
    SupportVectorMachine                 m_svMachine;         ///< The support vector machine
  
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool                    m_isPfoLevel;                   ///< Whether to use configuration for pfo vs cluster characterisation
    bool                    m_postBranchAddition;           ///< Whether to use configuration for shower clusters post branch addition

    pandora::StringVector  m_inputPfoListNames;
    std::string            m_trackPfoListName;
    std::string            m_showerPfoListName;
    pandora::StringVector   m_inputClusterListNames;        ///< The names of the input cluster lists
    
    std::string             m_trainingOutputFile;
    bool                    m_overwriteExistingId;          ///< Whether to consider any clusters that already have an assigned particle id
    bool                    m_useUnavailableClusters;       ///< Whether to consider clusters that are already constituents of a pfo
    unsigned int            m_minCaloHitsCut;               ///< The minimum number of calo hits to qualify as a track
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SVMClusterCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new SVMClusterCharacterisationAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline SVMClusterCharacterisationAlgorithm::ClusterFeatureInfo::ClusterFeatureInfo(/*const float hitType, const float nHits, */const float nGoodHits, /*const float straightLineLength,*/ 
    const float straightLineLengthLarge,const float showerFitWidth, const float showerFitGapLength,  const float diffWithStraigthLine, 
            const float widthDirectionX, const float nNearbyClusters, const float mipEnergy) :
   // m_hitType(hitType),
   // m_nHits(nHits),
    m_nGoodHits(nGoodHits), 
     //   m_straightLineLength(straightLineLength),
    m_straightLineLengthLarge(straightLineLengthLarge),
   m_showerFitWidth(showerFitWidth),
    m_showerFitGapLength(showerFitGapLength),
    m_diffWithStraigthLine(diffWithStraigthLine),
    m_widthDirectionX(widthDirectionX),
    m_nNearbyClusters(nNearbyClusters),
    m_mipEnergy(mipEnergy)
{
}

} // namespace lar_content

#endif // #ifndef LAR_SVM_CLUSTER_CHARACTERISATION_ALGORITHM_H
