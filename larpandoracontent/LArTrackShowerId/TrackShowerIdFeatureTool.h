/**
 *  @file   larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.h
 *
 *  @brief  Header file for the track shower id feature tools
 *
 *  $Log: $
 */
#ifndef LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H
#define LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H 1

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "Pandora/PandoraInternal.h"

namespace lar_content
{

typedef MvaFeatureTool<const pandora::Algorithm *const, const pandora::Cluster *const> ClusterCharacterisationFeatureTool;
typedef MvaFeatureTool<const pandora::Algorithm *const, const pandora::ParticleFlowObject *const> PfoCharacterisationFeatureTool;

// Helper function to get the position of a point with respect to the cathode in ICARUS
enum PositionInCryostat {AboveCathode, BelowCathode, WithinCathode};
PositionInCryostat LocatePointInCryostat_ICARUS(const float point_X); 

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  TwoDShowerFitFeatureTool to calculate variables related to sliding shower fit
 */
class TwoDShowerFitFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoDShowerFitFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);
    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

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

    unsigned int m_slidingShowerFitWindow; ///< The sliding shower fit window
    unsigned int m_slidingLinearFitWindow; ///< The sliding linear fit window
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  TwoDLinearFitFeatureTool class for the calculation of variables related to 2d sliding linear fit
 */
class TwoDLinearFitFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoDLinearFitFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);
    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

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

    unsigned int m_slidingLinearFitWindow;      ///< The sliding linear fit window
    unsigned int m_slidingLinearFitWindowLarge; ///< The sliding linear fit window - should be large, providing a simple linear fit
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  TwoDVertexDistanceFeatureTool class for the calculation of 2d distance to neutrino vertex
 */
class TwoDVertexDistanceFeatureTool : public ClusterCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoDVertexDistanceFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);
    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster);

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

    unsigned int m_slidingLinearFitWindow; ///< The sliding linear fit window
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  PfoHierarchyFeatureTool for calculation of features relating to reconstructed particle hierarchy
 */
class PfoHierarchyFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    PfoHierarchyFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);
    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  ThreeDLinearFitFeatureTool class for the calculation of variables related to 3d sliding linear fit
 */
class ThreeDLinearFitFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDLinearFitFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);
    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

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
    void CalculateVariablesSlidingLinearFit(const pandora::Cluster *const pCluster, float &straightLineLengthLarge,
        float &diffWithStraigthLineMean, float &maxFitGapLength, float &rmsSlidingLinearFit) const;

    unsigned int m_slidingLinearFitWindow;      ///< The sliding linear fit window
    unsigned int m_slidingLinearFitWindowLarge; ///< The sliding linear fit window - should be large, providing a simple linear fit
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  ThreeDVertexDistanceFeatureTool class for the calculation of 3d distance to neutrino vertex
 */
class ThreeDVertexDistanceFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDVertexDistanceFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);
    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  ConeChargeFeatureTool class for the calculation of charge distribution and conicalness
 */
class ConeChargeFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ConeChargeFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);
    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Configurable parameters to calculate cone charge variables
     *
     *  @param conMinHits: minimum hit requirement at start and end of pfo to calculate conicalness
     *  @param minCharge: minimum charge requirement at start and end of pfo to calculate conicalness
     *  @param conFracRange: conincal fractional range to determine start/end of pfo
     *  @param MoliereRadius: 10.1 cm to determine halo/core of pfo
     *  @param MoliereRadiusFrac: fraction of Moliere radius, default = 0.2
     */
    unsigned int m_conMinHits;
    float m_minCharge;
    float m_conFracRange;
    float m_MoliereRadius;
    float m_MoliereRadiusFrac;

    /**
     *  @brief Calculate charge distribution in relation to the Moeliere radius
     *
     *  @param caloHitList: the calo hit list of plane w
     *  @param pfoStart: start position of the pfo
     *  @param pfoDir: direction of pfo from the principle vector of pca
     *  @param chargeCore: to receive sum of charge within Moeliete radius * fraction
     *  @param chargeHalo: to receive sum of charge outside of Moeliere radius * fraction
     *  @param chargeCon: to receive weighted sum of total charge
     */
    void CalculateChargeDistribution(const pandora::CaloHitList &caloHitList, const pandora::CartesianVector &pfoStart,
        const pandora::CartesianVector &pfoDir, float &chargeCore, float &chargeHalo, float &chargeCon);

    /**
     *  @brief Calculate conicalness as a ratio of charge distribution at the end and start of pfo
     *
     *  @param caloHitList: the calo hit list of plane w
     *  @param pfoStart: start position of the pfo
     *  @param pfoDir: direction of pfo from the principle vector of pca
     *  @param pfoLength: length of the whole pfo
     *  return conicalness
     */
    float CalculateConicalness(const pandora::CaloHitList &caloHitList, const pandora::CartesianVector &pfoStart,
        const pandora::CartesianVector &pfoDir, const float pfoLength);
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  ConeChargeFeatureTool_ICARUS class for the calculation of charge distribution and conicalness, in ICARUS
 */
class ConeChargeFeatureTool_ICARUS : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ConeChargeFeatureTool_ICARUS();

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

        pandora::CartesianVector m_neutrinoVertex; //The neutrino vertex used to sort
    };

    /**                                                                                                                                                                                                    
     *  @brief  DistanceToVertexComparator class for comparison of two points with different hit type (left and right) wrt neutrino vertex position                                                        
     */
    class DistanceToVertexComparator
    {
    public:
        /**                                                                                                                                                                                                  
         * @brief Constructor                                                                                                                                                                                
         */
        DistanceToVertexComparator(const pandora::CartesianVector vertexPosition2D_A, const pandora::CartesianVector vertexPosition2D_B, 
        const pandora::HitType hitType_A, const pandora::HitType hitType_B);

        /**                                                                                                                                                                                                  
         @brief operator <                                                                                                                                                                                 
            *                                                                                                                                                                                                 
            * @param  rhs object for comparison                                                                                                                                                               
            *                                                                                                                                                                                                 
            * @return boolean                                                                                                                                                                                 
            */
        bool operator()(const pandora::CaloHit *const left, const pandora::CaloHit *const right) const;

        pandora::CartesianVector m_neutrinoVertex_A; // The neutrino vertex used to sort for hit type A, i.e. position is projected onto the correct view based on the hit type                              
        pandora::CartesianVector m_neutrinoVertex_B; // same for hit type B                                                                                                                                  
        pandora::HitType m_hitType_A;
        pandora::HitType m_hitType_B;
    };
    
    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);
    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_useICARUSCollectionPlane; ///< Whether to change the view logic in order to match the ICARUS Collection plane

    /**
     *  @brief Configurable parameters to calculate cone charge variables
     *
     *  @param conMinHits: minimum hit requirement at start and end of pfo to calculate conicalness
     *  @param minCharge: minimum charge requirement at start and end of pfo to calculate conicalness
     *  @param conFracRange: conincal fractional range to determine start/end of pfo
     *  @param MoliereRadius: 10.1 cm to determine halo/core of pfo
     *  @param MoliereRadiusFrac: fraction of Moliere radius, default = 0.2 
     */
    unsigned int m_conMinHits;
    float m_minCharge;
    float m_conFracRange;
    float m_MoliereRadius;
    float m_MoliereRadiusFrac;

    /**
     *  @brief  Functions to order the calo hit list by distance to neutrino vertex
     *
     *  @param  pAlgorithm, the algorithm
     *  @param  pCluster the cluster we are characterizing
     *  @param  caloHitList to receive the ordered calo hit list
     */
    void OrderCaloHitsByDistanceToVertex(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, 
        pandora::CaloHitList &caloHitList);
    void CombineCaloHitListsToHaveCollection(const pandora::Algorithm *const pAlgorithm,
        const pandora::CaloHitList &orderedCaloHitList1, const pandora::CaloHitList &orderedCaloHitList2, pandora::CaloHitList &mergedCaloHitList);

    /**
     *  @brief Calculate charge distribution in relation to the Moeliere radius
     *
     *  @param caloHitList: the calo hit list of plane w
     *  @param pfoStart: start position of the pfo
     *  @param pfoDir: direction of pfo from the principle vector of pca
     *  @param chargeCore: to receive sum of charge within Moeliete radius * fraction
     *  @param chargeHalo: to receive sum of charge outside of Moeliere radius * fraction
     *  @param chargeCon: to receive weighted sum of total charge
     */
    void CalculateChargeDistribution(const pandora::CaloHitList &caloHitList, const pandora::CartesianVector &pfoStart,
        const pandora::CartesianVector &pfoDir, float &chargeCore, float &chargeHalo, float &chargeCon);

    /**
     *  @brief Calculate conicalness as a ratio of charge distribution at the end and start of pfo
     *
     *  @param caloHitList: the calo hit list of plane w
     *  @param pfoStart: start position of the pfo
     *  @param pfoDir: direction of pfo from the principle vector of pca
     *  @param pfoLength: length of the whole pfo
     *  return conicalness
     */
    float CalculateConicalness(const pandora::CaloHitList &caloHitList, const pandora::CartesianVector &pfoStart,
        const pandora::CartesianVector &pfoDir, const float pfoLength);
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  ThreeDOpeningAngleFeatureTool class for the calculation of distance to neutrino vertex
 */
class ThreeDOpeningAngleFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDOpeningAngleFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);
    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief  Obtain positions at the vertex and non-vertex end of a list of three dimensional calo hits
     *
     *  @param  threeDCaloHitList the list of three dimensional calo hits
     *  @param  pointVectorStart to receive the positions at the start/vertex region
     *  @param  pointVectorEnd to receive the positions at the end region (opposite end to vertex)
     */
    void Divide3DCaloHitList(const pandora::Algorithm *const pAlgorithm, const pandora::CaloHitList &threeDCaloHitList,
        pandora::CartesianPointVector &pointVectorStart, pandora::CartesianPointVector &pointVectorEnd);

    /**
     *  @brief  Use the results of principal component analysis to calculate an opening angle
     *
     *  @param  principal the principal axis
     *  @param  secondary the secondary axis
     *  @param  eigenValues the eigenvalues
     *
     *  @return the opening angle
     */
    float OpeningAngle(const pandora::CartesianVector &principal, const pandora::CartesianVector &secondary,
        const pandora::CartesianVector &eigenValues) const;

    float m_hitFraction;  ///< Fraction of hits in start and end of pfo
    float m_defaultValue; ///< Default value to return, in case calculation not feasible
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  ThreeDPCAFeatureTool class for the calculation of PCA-related variables
 */
class ThreeDPCAFeatureTool : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDPCAFeatureTool();

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);
    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  ThreeDChargeFeatureTool class for the calculation of charge-related features
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

        pandora::CartesianVector m_neutrinoVertex; //The neutrino vertex used to sort
    };

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);
    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

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
    void CalculateChargeVariables(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, float &totalCharge,
        float &chargeSigma, float &chargeMean, float &endCharge);

    /**
     *  @brief  Function to order the calo hit list by distance to neutrino vertex
     *
     *  @param  pAlgorithm, the algorithm
     *  @param  pCluster the cluster we are characterizing
     *  @param  caloHitList to receive the ordered calo hit list
     *
     */
    void OrderCaloHitsByDistanceToVertex(
        const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, pandora::CaloHitList &caloHitList);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_endChargeFraction; ///< Fraction of hits that will be considered to calculate end charge (default 10%)
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *   @brief  ThreeDChargeFeatureTool_ICARUS class for the calculation of charge-related features, in ICARUS
 */
class ThreeDChargeFeatureTool_ICARUS : public PfoCharacterisationFeatureTool
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDChargeFeatureTool_ICARUS();

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

        pandora::CartesianVector m_neutrinoVertex; //The neutrino vertex used to sort
    };

    /**                                                                                                                                                                                                    
     *  @brief  DistanceToVertexComparator class for comparison of two points with different hit type (left and right) wrt neutrino vertex position                                                        
     */
    class DistanceToVertexComparator
    {
    public:
        /**                                                                                                                                                                                                  
         * @brief Constructor                                                                                                                                                                                
         */
        DistanceToVertexComparator(const pandora::CartesianVector vertexPosition2D_A, const pandora::CartesianVector vertexPosition2D_B, 
        const pandora::HitType hitType_A, const pandora::HitType hitType_B);

        /**                                                                                                                                                                                                  
         @brief operator <                                                                                                                                                                                 
            *                                                                                                                                                                                                 
            * @param  rhs object for comparison                                                                                                                                                               
            *                                                                                                                                                                                                 
            * @return boolean                                                                                                                                                                                 
            */
        bool operator()(const pandora::CaloHit *const left, const pandora::CaloHit *const right) const;

        pandora::CartesianVector m_neutrinoVertex_A; // The neutrino vertex used to sort for hit type A, i.e. position is projected onto the correct view based on the hit type                              
        pandora::CartesianVector m_neutrinoVertex_B; // same for hit type B                                                                                                                                  
        pandora::HitType m_hitType_A;
        pandora::HitType m_hitType_B;
    };

    void Run(LArMvaHelper::MvaFeatureVector &featureVector, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);
    void Run(LArMvaHelper::MvaFeatureMap &featureMap, pandora::StringVector &featureOrder, const std::string &featureToolName,
        const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo);

private:
    /**
     *  @brief  Calculation of the charge variables
     *
     *  @param  orderedCaloHitList, that contains the list of collection hits ordered by distance to vertex 
     *  @param  totalCharge, to receive the total charge
     *  @param  chargeSigma, to receive the charge sigma
     *  @param  chargeMean, to receive the charge mean
     *  @param  endCharge, to receive the charge in the last 10% hits
     */
    void CalculateChargeVariables(const pandora::CaloHitList &orderedCaloHitList, float &totalCharge, 
        float &chargeSigma, float &chargeMean, float &endCharge);

    /**
     *  @brief  Functions to order the calo hit list by distance to neutrino vertex
     *
     *  @param  pAlgorithm, the algorithm
     *  @param  pCluster the cluster we are characterizing
     *  @param  caloHitList to receive the ordered calo hit list
     */
    void OrderCaloHitsByDistanceToVertex(const pandora::Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, 
        pandora::CaloHitList &caloHitList);
    void CombineCaloHitListsToHaveCollection(const pandora::Algorithm *const pAlgorithm,
        const pandora::CaloHitList &orderedCaloHitList1, const pandora::CaloHitList &orderedCaloHitList2, pandora::CaloHitList &mergedCaloHitList);
    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_useICARUSCollectionPlane; ///< Whether to change the view logic in order to match the ICARUS Collection plane
    float m_endChargeFraction;       ///< Fraction of hits that will be considered to calculate end charge (default 10%)
};

} // namespace lar_content

#endif // #ifndef LAR_TRACK_SHOWER_ID_FEATURE_TOOLS_H
