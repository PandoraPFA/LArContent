/**
 *  @file   larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.h
 * 
 *  @brief  Header file for the track shower id feature tools 
 * 
 *  $Log: $
 */
#ifndef LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H
#define LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H 1

#include "larpandoracontent/LArTrackShowerId/SVMClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/BranchGrowingAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h" 

namespace lar_content
{

/**
*  @brief  TrackShowerIdFeatureTool class
*/
class TrackShowerIdFeatureTool 
{
    public:

//#ifndef LAR_SHOWER_FIT_WIDTH_FEATURE_TOOL_H
//#define LAR_SHOWER_FIT_WIDTH_FEATURE_TOOL_H 1

/**
*   @brief  Feature tool to calculate shower fit width
*/
class ShowerFitWidthFeatureTool : public SVMClusterCharacterisationAlgorithm::ClusterCharacterisationFeatureTool
{
public:
   /**
    *  @brief  Factory class for instantiating algorithm tool
    */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

   /**
    *  @brief  Default constructor
    */
    ShowerFitWidthFeatureTool();

   /**
    *  @brief  Run the tool
    * 
    *  @param  pAlgorithm address of the calling algorithm
    *  @param  pCluster, the cluster we are characterizing
    *  @param  pClusterList, list of clusters in the relevant view (for some of the variables)
    * 
    *  @return the global asymmetry feature
    */
    void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float CalculateShowerFitWidth(const pandora::Cluster *const pCluster) const;

    float   m_slidingShowerFitWindow;       ///< 
};

//#endif
//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

class NHitsFeatureTool : public SVMClusterCharacterisationAlgorithm::ClusterCharacterisationFeatureTool
{
public:

    /**
    *  @brief  Factory class for instantiating algorithm tool
    */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
    pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };
    
    /**
    *  @brief  Default constructor
    */
    NHitsFeatureTool();

    /**
    *  @brief  Run the tool
    * 
    *   @param  pAlgorithm address of the calling algorithm
    *  @param  pCluster, the cluster we are characterizing
    *  @param  pClusterList, list of clusters in the relevant view (for some of the variables)
    * 
    *  @return the global asymmetry feature
    */
    void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster *const pCluster);

    private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    int CalculateNHits(const pandora::Cluster * const pCluster) const;
    int CalculateNGoodHits(const pandora::Cluster * const pCluster,const SVMClusterCharacterisationAlgorithm * const pAlgorithm) const;
    std::string             m_mcParticleListName;           ///< Name of input MC particle list 
};

//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
class ShowerFitGapLengthFeatureTool : public SVMClusterCharacterisationAlgorithm::ClusterCharacterisationFeatureTool
{
public:

    /**
    *  @brief  Factory class for instantiating algorithm tool
    */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
    pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };
//--------------------------------------------------------------------------------------------------------------------------------------

    /**
    *  @brief  Default constructor
*/
ShowerFitGapLengthFeatureTool();

/**
*  @brief  Run the tool
* 
*  @param  pAlgorithm address of the calling algorithm
*  @param  pCluster, the cluster we are characterizing
*  @param  pClusterList, list of clusters in the relevant view (for some of the variables)
* 
*  @return the global asymmetry feature
*/
void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster);

private:
pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
float CalculateShowerFitGapLength(const pandora::Cluster * const pCluster) const;
float     m_slidingShowerFitWindow; 
};

//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

class StraightLineLengthFeatureTool : public SVMClusterCharacterisationAlgorithm::ClusterCharacterisationFeatureTool
{
public:

    /**
    *  @brief  Factory class for instantiating algorithm tool
    */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
    pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };
    
    /**
    *  @brief  Default constructor
    */
    StraightLineLengthFeatureTool();
    
/**
*  @brief  Run the tool
* 
*  @param  pAlgorithm address of the calling algorithm
*  @param  pCluster, the cluster we are characterizing
*  @param  pClusterList, list of clusters in the relevant view (for some of the variables)
* 
*  @return the global asymmetry feature
*/
void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster);

private:
pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
void CalculateVariablesSlidingFit(const pandora::Cluster * const pCluster, float &straightLineLength, float &straightLineLengthLarge, float &diffWithStraigthLine, float &widthDirectionX) const;

float     m_slidingLinearFitWindow; 
bool      m_addDiffWithStraightLine;
bool      m_addWidthDirectionX;
};



//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

class PointsOfContactFeatureTool : public SVMClusterCharacterisationAlgorithm::ClusterCharacterisationFeatureTool
{
public:

    /**
    *  @brief  Factory class for instantiating algorithm tool
    */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
    pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };
    
    /**
    *  @brief  Default constructor
    */
    PointsOfContactFeatureTool();
    
/**
*  @brief  Run the tool
* 
*  @param  pAlgorithm address of the calling algorithm
*  @param  pCluster, the cluster we are characterizing
*  @param  pClusterList, list of clusters in the relevant view (for some of the variables)
* 
*  @return the global asymmetry feature
*/
void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster);
private:
pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
int CalculatePointsOfContact(const pandora::Cluster * const pCluster, const SVMClusterCharacterisationAlgorithm *const pAlgorithm) const;
void FindAssociatedClusters(const pandora::Cluster *const pParticleSeed, pandora::ClusterVector &candidateClusters, 
     BranchGrowingAlgorithm::ClusterUsageMap &forwardUsageMap, BranchGrowingAlgorithm::ClusterUsageMap &backwardUsageMap, const SVMClusterCharacterisationAlgorithm *const pAlgorithm) const;
void IdentifyClusterMerges(const pandora::ClusterVector &particleSeedVector, const BranchGrowingAlgorithm::ClusterUsageMap &backwardUsageMap, 
     BranchGrowingAlgorithm::SeedAssociationList &seedAssociationList) const;
BranchGrowingAlgorithm::AssociationType AreClustersAssociated(const pandora::Cluster *const pClusterSeed, 
     const pandora::Cluster *const pCluster, const SVMClusterCharacterisationAlgorithm *const pAlgorithm) const; 
     
  typedef std::unordered_map<const pandora::Cluster*, LArVertexHelper::ClusterDirection> ClusterDirectionMap;                                                
  mutable ClusterDirectionMap m_clusterDirectionMap;          ///< The cluster direction map 

pandora::StringVector             m_clusterListNames;           ///< Name of input MC particle list 
   float                       m_nearbyClusterDistance;        ///< The nearby cluster distance, used for determining cluster associations                
    float                       m_remoteClusterDistance;        ///< The remote cluster distance, used for determining cluster associations           
    float                       m_directionTanAngle;            ///< Direction determination, look for vertex inside triangle with apex shifted along the clus
    float                       m_directionApexShift;  


};


//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

class NNearbyClustersFeatureTool : public SVMClusterCharacterisationAlgorithm::ClusterCharacterisationFeatureTool
{
public:

    /**
    *  @brief  Factory class for instantiating algorithm tool
    */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
    pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };
    
    /**
    *  @brief  Default constructor
    */
    NNearbyClustersFeatureTool();
    
/**
*  @brief  Run the tool
* 
*  @param  pAlgorithm address of the calling algorithm
*  @param  pCluster, the cluster we are characterizing
*  @param  pClusterList, list of clusters in the relevant view (for some of the variables)
* 
*  @return the global asymmetry feature
*/
void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster);
private:
pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
int CalculatePointsOfContact(const pandora::Cluster * const pCluster, const SVMClusterCharacterisationAlgorithm *const pAlgorithm) const;
bool IsClusterNearby(const pandora::Cluster *const pCluster,const pandora::Cluster *const pCandidateCluster) const;    

pandora::StringVector             m_clusterListNames;           ///< Name of input MC particle list 
float m_nearbyClusterDistance;
};



//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

class MipEnergyFeatureTool : public SVMClusterCharacterisationAlgorithm::ClusterCharacterisationFeatureTool
{
public:

    /**
    *  @brief  Factory class for instantiating algorithm tool
    */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
    pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };
    
    /**
    *  @brief  Default constructor
    */
    MipEnergyFeatureTool();
    
    /**
*  @brief  Run the tool
* 
*  @param  pAlgorithm address of the calling algorithm
*  @param  pCluster, the cluster we are characterizing
*  @param  pClusterList, list of clusters in the relevant view (for some of the variables)
* 
*  @return the global asymmetry feature
*/
void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster); 
private:
pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
float CalculateMipEnergy(const pandora::Cluster * const pCluster) const;
float m_mipCorrectionPerHit;

};

};

/*inline pandora::Algorithm *TrackShowerIdFeatureTool::Factory::CreateAlgorithm() const
{
    return new TrackShowerIdFeatureTool();
}*/

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::ShowerFitWidthFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::ShowerFitWidthFeatureTool();
}

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::NHitsFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::NHitsFeatureTool();
}

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::ShowerFitGapLengthFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::ShowerFitGapLengthFeatureTool();
}

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::StraightLineLengthFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::StraightLineLengthFeatureTool();
}

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::PointsOfContactFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::PointsOfContactFeatureTool();
}

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::NNearbyClustersFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::NNearbyClustersFeatureTool();
}

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::MipEnergyFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::MipEnergyFeatureTool();
}

} // namespace lar_content

#endif // #ifndef LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H
