/**
 *  @file   larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.h
 *
 *  @brief  Header file for the track shower id feature tools
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H
#define LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H 1

//#include "larpandoracontent/LArObjects/LArSupportVectorMachine.h"
#include "larpandoracontent/LArObjects/LArAdaBoostDecisionTree.h"
#include <Eigen/Dense>

namespace lar_content
{

typedef MvaFeatureTool<const pandora::Algorithm *const, const pandora::Cluster *const>  ClusterCharacterisationFeatureTool;
typedef MvaFeatureTool<const pandora::Algorithm *const, const pandora::ParticleFlowObject *const>  PfoCharacterisationFeatureTool;

/**
 *   @brief  ShowerFitFeatureTool to calculate variables related to sliding shower fit
 */
class TwoDShowerFitFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoDShowerFitFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

   /**
    *  @brief  Calculation of the shower fit width variable
    *
    *  @param  pAlgorithm address of the calling algorithm
    *  @param  pCluster the cluster we are characterizing
    *
    *  @return shower fit width
    */
    float CalculateShowerFitWidth(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster) const;

    unsigned int    m_slidingShowerFitWindow;        ///< The sliding shower fit window
    unsigned int    m_slidingLinearFitWindow;        ///< The sliding linear fit window
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  LinearFitFeatureTool class for the calculation of variables related to sliding linear fit
 */
class TwoDLinearFitFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoDLinearFitFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Calculation of several variables related to sliding linear fit
     *
     *  @param  pCluster the cluster we are characterizing
     *  @param  straightLineLengthLarge to receive to length reported by the straight line fit
     *  @param  diffWithStraigthLineMean to receive the difference with straight line mean variable
     *  @param  diffWithStraightLineSigma to receive the difference with straight line sigma variable
     *  @param  dTdLWidth to receive the dTdL width variable
     *  @param  maxFitGapLength to receive the max fit gap length variable
     *  @param  rmsSlidingLinearFit to receive the RMS from the linear fit
     */
    void CalculateVariablesSlidingLinearFit(const pandora::Cluster *const pCluster, float &straightLineLengthLarge, float &diffWithStraigthLineMean,
        float &diffWithStraightLineSigma, float &dTdLWidth, float &maxFitGapLength, float &rmsSlidingLinearFit) const;

    unsigned int    m_slidingLinearFitWindow;       ///< The sliding linear fit window
    unsigned int    m_slidingLinearFitWindowLarge;  ///< The sliding linear fit window - should be large, providing a simple linear fit
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  VertexDistanceFeatureTool class for the calculation of distance to neutrino vertex
 */
class TwoDVertexDistanceFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoDVertexDistanceFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Calculation of vertex distance
     *
     *  @param  pCluster the cluster we are characterizing
     *
     *  @return distance to neutrino vertex
     */
    float CalculateVertexDistance(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster) const;

    unsigned int    m_slidingLinearFitWindow;       ///< The sliding linear fit window
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
/**
 *   @brief  class for the calculation of curvature/"wiggliness" of pfos
 */
class PfoHierarchyFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    PfoHierarchyFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

 /**
   *  @brief  Calculation of how much a pfo wiggles
   *
   *  @param  pAlgorithm                   address of the calling algorithm
   *  @param  pInputPfo                    PFO that we are characterising      
   */

};
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------


/**
 *   @brief  LinearFitFeatureTool class for the calculation of variables related to sliding linear fit
 */
class ThreeDLinearFitFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDLinearFitFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Calculation of several variables related to sliding linear fit
     *
     *  @param  pCluster the cluster we are characterizing
     *  @param  straightLineLengthLarge to receive to length reported by the straight line fit
     *  @param  diffWithStraigthLineMean to receive the difference with straight line mean variable
     *  @param  diffWithStraightLineSigma to receive the difference with straight line sigma variable
     *  @param  dTdLWidth to receive the dTdL width variable
     *  @param  maxFitGapLength to receive the max fit gap length variable
     *  @param  rmsSlidingLinearFit to receive the RMS from the linear fit
     */
    void CalculateVariablesSlidingLinearFit(const pandora::Cluster *const pCluster, float &straightLineLengthLarge, float &diffWithStraigthLineMean,
        float &maxFitGapLength, float &rmsSlidingLinearFit) const;

    unsigned int    m_slidingLinearFitWindow;       ///< The sliding linear fit window
    unsigned int    m_slidingLinearFitWindowLarge;  ///< The sliding linear fit window - should be large, providing a simple linear fit
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  VertexDistanceFeatureTool class for the calculation of distance to neutrino vertex
 */
class ThreeDVertexDistanceFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDVertexDistanceFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  VertexDistanceFeatureTool class for the calculation of distance to neutrino vertex
 */
class ThreeDOpeningAngleFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDOpeningAngleFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void Divide3DCaloHitList(const pandora::Algorithm *const pAlgorithm, pandora::CaloHitList &threeDCaloHitList,
        pandora::CartesianPointVector &pointVectorStart, pandora::CartesianPointVector &pointVectorEnd);

    float OpeningAngle(const pandora::CartesianVector &principal, const pandora::CartesianVector &secondary, const pandora::CartesianVector &eigenValues) const;

    float m_hitFraction;           ///< fraction of hits in start and end of pfo
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  PCA class for the calculation of PCA-related variables
 */
class ThreeDPCAFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDPCAFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  ChargeFeatureTool class for the calculation of concentration
 */
class ThreeDChargeFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDChargeFeatureTool();

    /**
     *  @brief  VertexComparator class for comparison of two points wrt neutrino vertex position
     */
    class VertexComparator
    {
    public:
        /**
         *  @brief  Constructor
         */
        VertexComparator(const pandora::CartesianVector vertexPosition2D);

        /**
         *  @brief  operator <
         *
         *  @param  rhs object for comparison
         *
         *  @return boolean
         */
        bool operator()(const pandora::CaloHit *const left, const pandora::CaloHit *const right) const;

        pandora::CartesianVector   m_neutrinoVertex;    //The neutrino vertex used to sort

    };

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    /**
     *  @brief  Calculation of the charge variables
     *
     *  @param  pAlgorithm, the algorithm
     *  @param  pCluster the cluster we are characterizing
     *  @param  totalCharge, to receive the total charge
     *  @param  chargeSigma, to receive the charge sigma
     *  @param  chargeMean, to receive the charge mean
     *  @param  startCharge, to receive the charge in the initial 10% hits
     *  @param  endCharge, to receive the charge in the last 10% hits
     */
    void CalculateChargeVariables(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, float &totalCharge, float &chargeSigma,
        float &chargeMean, float &endCharge);

    /**
     *  @brief  Function to order the calo hit list by distance to neutrino vertex
     *
     *  @param  pAlgorithm, the algorithm
     *  @param  pCluster the cluster we are characterizing
     *  @param  caloHitList to receive the ordered calo hit list
     *
     */
    void OrderCaloHitsByDistanceToVertex(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, pandora::CaloHitList &caloHitList);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float           m_endChargeFraction;           ///< Fraction of hits that will be considered to calculate end charge (default 10%)
};

} // namespace lar_content

#endif // #ifndef LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H
