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
    

private:
    pandora::StatusCode Run();


    SupportVectorMachine::DoubleVector GenerateFeatureVector(const pandora::Cluster *const pCluster) const;
    SupportVectorMachine::DoubleVector ProcessFeatureVector(const SupportVectorMachine::DoubleVector &featureVector) const; 
    bool IsClearTrack(const SupportVectorMachine::DoubleVector &featureVector) const;   
  
    ClusterCharacterisationFeatureTool::FeatureToolVector m_featureToolVector; ///< The feature tool map
    SupportVectorMachine                 m_svMachine;         ///< The support vector machine
  
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool                    m_produceSamplesMode;           ///< Whether to run for training samples production
    bool                    m_isPfoLevel;                   ///< Whether to use configuration for pfo vs cluster characterisation
    bool                    m_postBranchAddition;           ///< Whether to use configuration for shower clusters post branch addition
    bool                    m_ratioVariables;               ///< Whether to use variables as they are or their ratio over the straight line length
    
    std::string             m_trainingOutputFile;           ///< The training output file
    std::string             m_svmFileName;                  ///< The SVM input file
    std::string             m_svmName;                      ///< The name of the SVM to find

    pandora::StringVector   m_inputPfoListNames;            ///< The names of the input pfos lists
    std::string             m_trackPfoListName;             ///< The names of the input track pfos lists
    std::string             m_showerPfoListName;            ///< The names of the input shower pfos lists
    pandora::StringVector   m_inputClusterListNames;        ///< The names of the input cluster lists
    
    bool                    m_overwriteExistingId;          ///< Whether to consider any clusters that already have an assigned particle id
    bool                    m_updateClusterIds;             ///< Whether to update daughter cluster particle id labels to match pfo id
    bool                    m_useUnavailableClusters;       ///< Whether to consider clusters that are already constituents of a pfo
    unsigned int            m_minTrackLikeViews;            ///< The minimum number of clusters (views) to be track-like to label the pfo as track
    unsigned int            m_minCaloHitsCut;               ///< The minimum number of calo hits to qualify as a track
    
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *SVMClusterCharacterisationAlgorithm::Factory::CreateAlgorithm() const
{
    return new SVMClusterCharacterisationAlgorithm();
}


} // namespace lar_content

#endif // #ifndef LAR_SVM_CLUSTER_CHARACTERISATION_ALGORITHM_H
