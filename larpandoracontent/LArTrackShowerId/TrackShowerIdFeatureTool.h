/**
 *  @file   larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.h
 *
 *  @brief  Header file for the track shower id feature tools
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H
#define LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H 1

#include "larpandoracontent/LArObjects/LArSupportVectorMachine.h"
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
class TwoDCurvatureFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoDCurvatureFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

 /**
   *  @brief  Calculation of how much a pfo wiggles
   *
   *  @param  pAlgorithm                   address of the calling algorithm
   *  @param  pInputPfo                    PFO that we are characterising      
   */

    unsigned int    m_slidingLinearFitWindow;    ///< The sliding linear fit window to calculate the direction
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
 *   @brief  Class for the calculation of variables using PCA coordinates
 */
 class ThreeDPCAVariablesFeatureTool : public PfoCharacterisationFeatureTool
{
 public:
  /**
   *  @brief  Default constructor
   */
  ThreeDPCAVariablesFeatureTool();

  void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

 private:
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

  /**
   *  @brief  Calculation of concentration and conicalness
   *
   *  @param  threeDCaloHitList            CaloHitList
   *  @param  eVecsTranspose               Transpose of matrix of eigenvectors from principal compoenent analysis (PCA) 
   *  @param  PcaCentroid                  Centroid from PCA
   *  @param  paStart                      Start of pfo on principal axis
   *  @param  pfoLengthOnPA                Length of pfo on prinicpal axis
   *  @param  concentration                To receive valus of concentration
   *  @param  concentration2               To receive valus of concentration2 (alternative version of concentration)
   *  @param  conicalness                  To receive value of conicalness
   */
  void CalculateConcentrationConicalness(pandora::CaloHitList &threeDCaloHitList, Eigen::Matrix3f &eVecsTranspose, Eigen::Vector3f &PcaCentroid,
                                         const float paStart, const float pfoengthOnPA, float &concentration, float &concentration2, float &conicalness);

  /**
   *  @brief  Calculation of inverse hit density
   *
   *  @param  paStart                      Start of pfo on principal axis
   *  @param  pfoLengthOnPA                Length of pfo on prinicpal axis 
   *  @param  pcaPositions                 Map with PCA postions of hits ordered by position on principal axis
   *  @param  invHitDensityY               To receive inverse hit density in the y PCA coordinate
   *  @param  invHitDensityZ               To receive inverse hit density in the z PCA coordinate
   *  @param  invHitDensityYRatio          To receive ratio of inverse hit density between ebd and start in the y PCA coordinate
   *  @param  invHitDensityZRatio          To receive ratio of inverse hit density between ebd and start in the z PCA coordinate 
   */
  void CalculateDensity(const float paStart, const float pfoLengthOnPA, std::map<float,Eigen::Vector3f> &pcaPositions,
                        float &invHitDensityY, float &invHitDensityZ, float &invHitDensityYRatio, float &invHitDensityZRatio);

  /**
   *  @brief  Calculation of concentration and conicalness
   *
   *  @param  paStart                      Start of pfo on principal axis
   *  @param  pcaPositions                 Map with PCA postions of hits ordered by position on principal axis
   *  @param  nSegmentsDoubleHits          To receive number of segments with double hits transverse to principal axis
   *  @param  nSegmentsDoubleHitsRatio     To receive ratio of number of segments with double hits transverse to principal axis to number of segments
   */
  void SegmentsWithDoubleHits(const float paStart, std::map<float,Eigen::Vector3f> &pcaPositions, int &nSegmentsDoubleHits, float &nSegmentsDoubleHitsRatio);

  float        m_pfoFraction;              ///< Fraction of length of pfo used to define start and end
  unsigned int m_minHits;                  ///< Minimum number of hits for concentration, conicalness and hit densities
  float        m_segmentWidth;             ///< Width of segments of pfo used for inverse hit density and double hits transverse to principal axis  
  float        m_MoliereRadius;            ///< Moliere radius of liquid argon (10.1 cm)
  float        m_MoliereFraction;          ///< Fraction of Moliere radius used in definition of double hits transverse to principal axis
  int          m_slidingLinearFitWindow;   ///< Number of windows to use for slidinglinearfit
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
