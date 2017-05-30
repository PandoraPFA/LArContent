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
*  @brief  TrackShowerIdFeatureTool class, encompassing different tools calculating cunning variables 
*/
class TrackShowerIdFeatureTool 
{
    public:

//#ifndef LAR_SHOWER_FIT_WIDTH_FEATURE_TOOL_H
//#define LAR_SHOWER_FIT_WIDTH_FEATURE_TOOL_H 1

/**
*   @brief  ShowerFitFeatureTool to calculate variables related to sliding shower fit
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

    float 			m_slidingShowerFitWindow;		///<  The sliding shower fit window
};

//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
/**
*   @brief  NHitsFeatureTool class for the calculation of number of hits
*/
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
    *  @param  pAlgorithm address of the calling algorithm
    *  @param  pCluster, the cluster we are characterizing
    * 
    */
    void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster *const pCluster);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
	
	/**
    *  @brief  Calculation of number of hits
    * 
    *  @param  pCluster, the cluster we are characterizing
    * 
    *  @return number of hits
    */
    int CalculateNHits(const pandora::Cluster * const pCluster) const;
	
	/**
    *  @brief  Calculation of number of good hits: ATTN - this will go, in MCHelper now
    * 
    *  @param  pCluster, the cluster we are characterizing
    * 
    *  @return number of good hits
    */
    int CalculateNGoodHits(const pandora::Cluster * const pCluster,const SVMClusterCharacterisationAlgorithm * const pAlgorithm) const;
	
    std::string			m_mcParticleListName;		///< Name of input MC particle list 
};

//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
/**
*   @brief  ShowerFitGapLengthFeatureTool class for the shower fit gap length
*/
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
    /**
    *  @brief  Default constructor
	*/
	ShowerFitGapLengthFeatureTool();

	/**
	*  @brief  Run the tool
	* 
	*  @param  pAlgorithm address of the calling algorithm
	*  @param  pCluster, the cluster we are characterizing
	*/
	void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster);

private:
	pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

	/**
    *  @brief  Calculation of shower fit gap length 
    * 
    *  @param  pCluster, the cluster we are characterizing
    * 
    *  @return shower fit gap length 
    */	
	float CalculateShowerFitGapLength(const pandora::Cluster * const pCluster) const;
	
	float 			m_slidingShowerFitWindow; 		///<  The sliding shower fit window
};

//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
/**
*   @brief  VertexDistanceFeatureTool class for the calculation of distance to neutrino vertex
*/
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

    /**
    *  @brief  Default constructor
	*/
	VertexDistanceFeatureTool();

	/**
	*  @brief  Run the tool
	* 
	*  @param  pAlgorithm address of the calling algorithm
	*  @param  pCluster, the cluster we are characterizing
	*/
	void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster);

private: 
	pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
	
	/**
    *  @brief  Calculation of vertex distance
    * 
    *  @param  pCluster, the cluster we are characterizing
    * 
    *  @return distance to neutrino vertex 
    */	
	float CalculateVertexDistance(const SVMClusterCharacterisationAlgorithm *const pAlgorithm, const pandora::Cluster * const pCluster) const;
	
	bool 			m_addVertexDistance; 		///<  decide whether to add the vertex distance variable to the feature list 
};

//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
/**
*   @brief  LinearFitFeatureTool class for the calculation of variables related to sliding linear fit
*/
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
	*/
	void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster);

private:
	pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
			
	/**
    *  @brief  Calculation of several variables related to sliding linear fit
    * 
    *  @param  pCluster, the cluster we are characterizing
    */	
	void CalculateVariablesSlidingLinearFit(const pandora::Cluster * const pCluster, float &straightLineLengthLarge, float &diffWithStraightLine,float &diffWithStraightLineMean, float &diffWithStraightLineSigma, float &dTdLWidth, float &maxFitGapLength, float &rmsSlidingLinearFit) const;

	float			m_slidingLinearFitWindow;     	///<  The sliding linear fit window
	bool			m_addDiffWithStraightLine; 		///<  decide whether to add the difference with straight line variable to the feature list
	bool			m_adddTdLWidth; 				///<  decide whether to add the dTdL variable to the feature list
	bool			m_addMaxFitGapLength;			///<  decide whether to add the max fit gap length variable to the feature list
	bool                    m_addRMSLinearFit;                      ///<  decide whether to add the RMS from the linear fit
};

//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
/**
*   @brief  NNearbyClustersFeatureTool class for the calculation of number of clusters nearby
*/
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
	*/
	void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster);
private:
	pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

	/**
    *  @brief  Calculation of number of points of contact, i.e. number of clusters nearby
    * 
    *  @param  pCluster, the cluster we are characterizing
	*  @param  pAlgorithm address of the calling algorithm
	*     
    *  @return number of clusters nearby pCluster
    */	
	int CalculatePointsOfContact(const pandora::Cluster * const pCluster, const SVMClusterCharacterisationAlgorithm *const pAlgorithm) const;

	/**
    *  @brief  Address whether pCandidateCluster can be considered nearby pCluster
    * 
    *  @param  pCluster, the cluster we are characterizing
	*  @param  pCandidateCluster, the cluster we are asking whether it is nearby pCluster 
	*     
    *  @return yes or no
    */
	bool IsClusterNearby(const pandora::Cluster *const pCluster,const pandora::Cluster *const pCandidateCluster) const;    

	pandora::StringVector             m_clusterListNames;           ///< Name of input MC particle list 
	float 							  m_nearbyClusterDistance;      ///< Distance to decide whether a cluster is considered nearby
};


//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
/**
*   @brief  MipEnergyFeatureTool class for the calculation of mip energy
*/
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
	*/
	void Run(SupportVectorMachine::DoubleVector &featureVector, const SVMClusterCharacterisationAlgorithm * const pAlgorithm, const pandora::Cluster * const pCluster); 

private:
	pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
	
	/**
    *  @brief  Calculate the mip energy of the cluster
    * 
    *  @param  pCluster, the cluster we are characterizing
	*     
    *  @return the mip energy of the cluster 
    */
	float CalculateMipEnergy(const pandora::Cluster * const pCluster) const;

	float		m_mipCorrectionPerHit;		///< Mip correction per hit

};

};


//--------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

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
