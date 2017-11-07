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

namespace lar_content
{

typedef SvmFeatureTool<const pandora::Algorithm *const, const pandora::Cluster *const>  ClusterCharacterisationFeatureTool;
typedef SvmFeatureTool<const pandora::Algorithm *const, const pandora::ParticleFlowObject *const>  PfoCharacterisationFeatureTool;

/**
 *   @brief  ShowerFitFeatureTool to calculate variables related to sliding shower fit
 */
class ShowerFitFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ShowerFitFeatureTool();

    void Run(SupportVectorMachine::DoubleVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

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

    bool            m_ratioVariables;                ///< Whether to divide the vertex distance by the straight line length
    unsigned int    m_slidingShowerFitWindow;        ///< The sliding shower fit window
    unsigned int    m_slidingLinearFitWindow;        ///< The sliding linear fit window
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  NHitsFeatureTool class for the calculation of number of hits
 */
class NHitsFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    void Run(SupportVectorMachine::DoubleVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  LinearFitFeatureTool class for the calculation of variables related to sliding linear fit
 */
class LinearFitFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    LinearFitFeatureTool();

    void Run(SupportVectorMachine::DoubleVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

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
    bool            m_ratioVariables;               ///< Whether to divide all variables by the straight line length
    bool            m_addStraightLineLength;        ///< Decide whether to add the straight line length to the feature list
    bool            m_addDiffWithStraightLineMean;  ///< Decide whether to add the difference with straight line mean variable to the feature list
    bool            m_addDiffWithStraightLineSigma; ///< Decide whether to add the difference with straight line sigma variable to the feature list
    bool            m_addDTDLWidth;                 ///< Decide whether to add the dTdL width variable to the feature list
    bool            m_addMaxFitGapLength;           ///< Decide whether to add the max fit gap length variable to the feature list
    bool            m_addRMSLinearFit;              ///< Decide whether to add the RMS from the linear fit
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  NNearbyClustersFeatureTool class for the calculation of number of clusters nearby
 */
class NNearbyClustersFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    /**
    *  @brief  Default constructor
    */
    NNearbyClustersFeatureTool();

    void Run(SupportVectorMachine::DoubleVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Calculation of number of points of contact, i.e. number of clusters nearby
     *
     *  @param  pCluster the cluster we are characterizing
     *  @param  pAlgorithm address of the calling algorithm
     *
     *  @return number of clusters nearby pCluster
     */
    int CalculatePointsOfContact(const pandora::Cluster *const pCluster, const pandora::Algorithm *const pAlgorithm) const;

    /**
     *  @brief  Address whether pCandidateCluster can be considered nearby pCluster
     *
     *  @param  pCluster the cluster we are characterizing
     *  @param  pCandidateCluster the cluster we are asking whether it is nearby pCluster
     *
     *  @return boolean
     */
    bool IsClusterNearby(const pandora::Cluster *const pCluster,const pandora::Cluster *const pCandidateCluster) const;

    pandora::StringVector   m_clusterListNames;             ///< Name of input MC particle list
    unsigned int            m_minClusterCaloHits;           ///< Minimum number of hits per cluster considered
    float                   m_nearbyClusterDistance;        ///< Distance to decide whether a cluster is considered nearby
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  MipEnergyFeatureTool class for the calculation of mip energy
 */
class MipEnergyFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    MipEnergyFeatureTool();

    void Run(SupportVectorMachine::DoubleVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Calculate the mip energy of the cluster
     *
     *  @param  pCluster the cluster we are characterizing
     *
     *  @return the mip energy of the cluster
     */
    float CalculateMipEnergy(const pandora::Cluster *const pCluster) const;

    float   m_mipCorrectionPerHit;      ///< Mip correction per hit
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  VertexDistanceFeatureTool class for the calculation of distance to neutrino vertex
 */
class VertexDistanceFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexDistanceFeatureTool();

    void Run(SupportVectorMachine::DoubleVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

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

    bool            m_ratioVariables;               ///< Whether to divide the vertex distance by the straight line length
    unsigned int    m_slidingLinearFitWindow;       ///< The sliding linear fit window
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

    void Run(SupportVectorMachine::DoubleVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

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
    bool            m_ratioVariables;               ///< Whether to divide all variables by the straight line length
    bool            m_addStraightLineLength;        ///< Decide whether to add the straight line length to the feature list
    bool            m_addDiffWithStraightLineMean;  ///< Decide whether to add the difference with straight line mean variable to the feature list
    bool            m_addDiffWithStraightLineSigma; ///< Decide whether to add the difference with straight line sigma variable to the feature list
    bool            m_addDTDLWidth;                 ///< Decide whether to add the dTdL width variable to the feature list
    bool            m_addMaxFitGapLength;           ///< Decide whether to add the max fit gap length variable to the feature list
    bool            m_addRMSLinearFit;              ///< Decide whether to add the RMS from the linear fit
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

    void Run(SupportVectorMachine::DoubleVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool            m_ratioVariables;               ///< Whether to divide the vertex distance by the straight line length
    unsigned int    m_slidingLinearFitWindow;       ///< The sliding linear fit window
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  VertexDistanceFeatureTool class for the calculation of distance to neutrino vertex
 */
class OpeningAngleFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    OpeningAngleFeatureTool();

    void Run(SupportVectorMachine::DoubleVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    void Divide3DCaloHitList(const pandora::Algorithm *const pAlgorithm, pandora::CaloHitList threeDCaloHitList,
        pandora::CartesianPointVector &pointVectorStart, pandora::CartesianPointVector &pointVectorEnd);

    //static bool SortByDistanceToVertex(const pandora::CaloHit *const left, const pandora::CaloHit *const right, const pandora::CartesianVector &nuVertex);
    float OpeningAngle(const pandora::CartesianVector &principal, const pandora::CartesianVector &secondary, const pandora::CartesianVector &eigenValues) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  PCA class for the calculation of PCA-related variables
 */
class PCAFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    PCAFeatureTool();

    void Run(SupportVectorMachine::DoubleVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool            m_ratioVariables;               ///< Whether to use ratios of secondary to principal and tertiary to principal
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  ChargeFeatureTool class for the calculation of charge related variables
 */
class ChargeFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ChargeFeatureTool();

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
    static void CalculateChargeVariables(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, float &totalCharge, float &chargeSigma,
        float &chargeMean, float &startCharge, float &endCharge, float &incrementCharge, float endChargeFraction);

    /**
     *  @brief  Function to order the calo hit list by distance to neutrino vertex
     *
     *  @param  pAlgorithm, the algorithm
     *  @param  pCluster the cluster we are characterizing
     *
     */
    static pandora::CaloHitList OrderCaloHitsByDistanceToVertex(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

    /**
     *  @brief  Static function deciding the sorting
     *
     *  @param  left, the calohit on the left
     *  @param  right, the calohit on the right
     *  @param  nuVertex2D, the 2D projection of the position of the neutrino vertex
     *
     */
    static bool SortByDistanceToVertex(const pandora::CaloHit *const left, const pandora::CaloHit *const right, const pandora::CartesianVector &nuVertex2D);

    void Run(SupportVectorMachine::DoubleVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool            m_ratioVariables;              ///< Whether to divide the charge variables by the total or mean charge
    bool            m_addTotalCharge;              ///< Decide whether to add total charge to the feature vector
    bool            m_addChargeSigma;              ///< Decide whether to add charge sigma to the feature vector
    bool            m_addChargeIncrements;         ///< Decide whether to add charge increments to the feature vector
    bool            m_addStartCharge;              ///< Decide whether to add charge in the first 10% hits to the feature vector
    bool            m_addEndCharge;                ///< Decide whether to add charge in the last 10% hits to the feature vector
    float           m_endChargeFraction;           ///< Fraction of hits that will be considered to calculate end charge (default 10%)
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

    void Run(SupportVectorMachine::DoubleVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool            m_ratioVariables;              ///< Whether to divide the charge variables by the total or mean charge
    bool            m_addTotalCharge;              ///< Decide whether to add total charge to the feature vector
    bool            m_addChargeSigma;              ///< Decide whether to add charge sigma to the feature vector
    bool            m_addChargeIncrements;         ///< Decide whether to add charge increments to the feature vector
    bool            m_addStartCharge;              ///< Decide whether to add charge in the first 10% hits to the feature vector
    bool            m_addEndCharge;                ///< Decide whether to add charge in the last 10% hits to the feature vector
    float           m_endChargeFraction;           ///< Fraction of hits that will be considered to calculate end charge (default 10%)
};

} // namespace lar_content

#endif // #ifndef LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H
