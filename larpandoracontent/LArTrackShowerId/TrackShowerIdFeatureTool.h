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
*   @brief  
*/
class ShowerFitFeatureTool : public SVMClusterCharacterisationAlgorithm::ClusterCharacterisationFeatureTool
{
public:
   /**
    *  @brief  Factory class for instantiating algorithm tool to calculate variables related to sliding shower fit 
    */
    class Factory : public pandora::AlgorithmToolFactory
    {
    public:
        pandora::AlgorithmTool *CreateAlgorithmTool() const;
    };

   /**
    *  @brief  Default constructor
    */
    ShowerFitFeatureTool();

   /**
    *  @brief  Run the tool, add the feature(s) calculated through sliding shower fit to the featureVector
    * 
    *  @param  pAlgorithm address of the calling algorithm
    *  @param  pCluster, the cluster we are characterizing
    * 
    */
    void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


   /**
    *  @brief  Calculation of the shower fit width variable
    * 
    *  @param  pAlgorithm, address of the calling algorithm
    *  @param  pCluster, the cluster we are characterizing
    * 
    *  @return shower fit width
    */
    float CalculateShowerFitWidth(const SVMClusterCharacterisationAlgorithm *const pAlgorithm, const pandora::Cluster *const pCluster) const;

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
class VertexDistanceFeatureTool : public SVMClusterCharacterisationAlgorithm::ClusterCharacterisationFeatureTool
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
VertexDistanceFeatureTool();

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
float CalculateVertexDistance(const SVMClusterCharacterisationAlgorithm *const pAlgorithm, const pandora::Cluster * const pCluster) const;
bool m_addVertexDistance;
};

//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

class LinearFitFeatureTool : public SVMClusterCharacterisationAlgorithm::ClusterCharacterisationFeatureTool
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
    LinearFitFeatureTool();
    
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
void CalculateVariablesSlidingLinearFit(const pandora::Cluster * const pCluster, float &straightLineLengthLarge, float &diffWithStraigthLine, float &dTdLWidth, float &maxFitGapLength) const;

float     m_slidingLinearFitWindow; 
bool      m_addDiffWithStraightLine;
bool      m_adddTdLWidth;
bool      m_addMaxFitGapLength;
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


//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------


/*inline pandora::Algorithm *TrackShowerIdFeatureTool::Factory::CreateAlgorithm() const
{
    return new TrackShowerIdFeatureTool();
}*/

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::ShowerFitFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::ShowerFitFeatureTool();
}

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::NHitsFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::NHitsFeatureTool();
}

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::ShowerFitGapLengthFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::ShowerFitGapLengthFeatureTool();
}

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::LinearFitFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::LinearFitFeatureTool();
}

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::NNearbyClustersFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::NNearbyClustersFeatureTool();
}

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::MipEnergyFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::MipEnergyFeatureTool();
}

inline pandora::AlgorithmTool *TrackShowerIdFeatureTool::VertexDistanceFeatureTool::Factory::CreateAlgorithmTool() const
{
    return new TrackShowerIdFeatureTool::VertexDistanceFeatureTool();
}

} // namespace lar_content

#endif // #ifndef LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H
