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
      ClusterFeatureInfo(const float straightLineLengthLarge,  
        const float showerFitWidthRatio,  const float showerFitGapLengthRatio, const float diffWithStraigthLineRatio, 
        const float widthDirectionXRatio, const float nNearbyClustersRatio);
    
 float m_hitType;  
 float           m_nHits;             ///<  
      float           m_nGoodHits;             ///<      
      float           m_straightLineLength;                   ///
    float m_straightLineLengthLarge;                  
      float           m_showerFitWidthRatio;        ///< 
      float           m_showerFitGapLengthRatio;         ///< 
    float m_diffWithStraigthLineRatio;
    float m_widthDirectionXRatio;
    float m_nNearbyClustersRatio;
    float m_mipEnergy;
    
    };
    
     typedef std::map<const pandora::Cluster *const, ClusterFeatureInfo> ClusterFeatureInfoMap;
  

private:
    pandora::StatusCode Run();


  void PopulateClusterFeatureInfoMap(const pandora::Cluster *const pCluster, ClusterFeatureInfoMap &clusterFeatureInfoMap) const;
  SupportVectorMachine::DoubleVector GenerateFeatureList(const pandora::Cluster *const pCluster, ClusterFeatureInfoMap &clusterFeatureInfoMap);
  
    ClusterCharacterisationFeatureTool::FeatureToolVector m_featureToolVector; ///< The feature tool map
    SupportVectorMachine                 m_svMachine;         ///< The support vector machine
  
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool                    m_produceSamplesMode;           ///< Whether to run for training samples production
    bool                    m_isPfoLevel;                   ///< Whether to use configuration for pfo vs cluster characterisation
    bool                    m_postBranchAddition;           ///< Whether to use configuration for shower clusters post branch addition
    
    std::string             m_trainingOutputFile;           ///< The training output file
    std::string             m_parameterInputFile;           ///< The parameter input file
    std::string             m_svmName;                      ///< The name of the SVM to find

    pandora::StringVector   m_inputPfoListNames;            ///< The names of the input pfos lists
    std::string             m_trackPfoListName;             ///< The names of the input track pfos lists
    std::string             m_showerPfoListName;            ///< The names of the input shower pfos lists
    pandora::StringVector   m_inputClusterListNames;        ///< The names of the input cluster lists
    
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

inline SVMClusterCharacterisationAlgorithm::ClusterFeatureInfo::ClusterFeatureInfo(const float straightLineLengthLarge, const float showerFitWidthRatio, const float showerFitGapLengthRatio, 
            const float diffWithStraigthLineRatio, const float widthDirectionXRatio, const float nNearbyClustersRatio) :
   // m_hitType(hitType),
   // m_nHits(nHits),
   // m_nGoodHits(nGoodHits), 
     //   m_straightLineLength(straightLineLength),
    m_straightLineLengthLarge(straightLineLengthLarge),
    m_showerFitWidthRatio(showerFitWidthRatio),
    m_showerFitGapLengthRatio(showerFitGapLengthRatio),
    m_diffWithStraigthLineRatio(diffWithStraigthLineRatio),
    m_widthDirectionXRatio(widthDirectionXRatio),
    m_nNearbyClustersRatio(nNearbyClustersRatio)
{
}



} // namespace lar_content

#endif // #ifndef LAR_SVM_CLUSTER_CHARACTERISATION_ALGORITHM_H
