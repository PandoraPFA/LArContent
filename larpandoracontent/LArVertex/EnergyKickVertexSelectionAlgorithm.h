/**
 *  @file   larpandoracontent/LArVertex/EnergyKickVertexSelectionAlgorithm.h
 * 
 *  @brief  Header file for the energy kick vertex selection algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_ENERGY_KICK_VERTEX_SELECTION_ALGORITHM_H
#define LAR_ENERGY_KICK_VERTEX_SELECTION_ALGORITHM_H 1

#include "larpandoracontent/LArVertex/VertexSelectionBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  EnergyKickVertexSelectionAlgorithm class
 */
class EnergyKickVertexSelectionAlgorithm : public VertexSelectionBaseAlgorithm
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
    EnergyKickVertexSelectionAlgorithm();

private:
     /**
     *  @brief Sliding fit data class.
     */
    class SlidingFitData
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pCluster pointer to the cluster
         *  @param  slidingFitWindow the sliding fit window
         *  @param  slidingFitPitch the sliding fit pitch
         */
        SlidingFitData(const pandora::Cluster *const pCluster, const int slidingFitWindow, const float slidingFitPitch);
        
        /**
         *  @brief  Get the min layer direction
         * 
         *  @return the min layer direction
         */
        const pandora::CartesianVector &GetMinLayerDirection() const;

        /**
         *  @brief  Get the max layer direction
         * 
         *  @return the max layer direction
         */
        const pandora::CartesianVector &GetMaxLayerDirection() const;

        /**
         *  @brief  Get the min layer position
         * 
         *  @return the min layer position
         */
        const pandora::CartesianVector &GetMinLayerPosition() const;

        /**
         *  @brief  Get the max layer position
         * 
         *  @return the max layer position
         */
        const pandora::CartesianVector &GetMaxLayerPosition() const;
        
        /**
         *  @brief  Get a pointer to the corresponding cluster
         * 
         *  @return pointer to the corresponding cluster
         */
        const pandora::Cluster *GetCluster() const;

    private:
        pandora::CartesianVector    m_minLayerDirection;    ///< The direction of the fit at the min layer
        pandora::CartesianVector    m_maxLayerDirection;    ///< The direction of the fit at the min layer
        pandora::CartesianVector    m_minLayerPosition;     ///< The position of the fit at the max layer
        pandora::CartesianVector    m_maxLayerPosition;     ///< The position of the fit at the max layer
        const pandora::Cluster     *m_pCluster;             ///< Pointer to the corresponding cluster
    };

    typedef std::vector<SlidingFitData> SlidingFitDataList;

    void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const;

    /**
     *  @brief  Calculate the sliding fits data objects for the clusters in a given view
     * 
     *  @param  hitType the view hit type
     *  @param  slidingFitPitch the sliding fit pitch
     *  @param  slidingFitDataList to receive the list of sliding fit data objects
     */
    void CalculateClusterSlidingFits(const pandora::HitType hitType, const float slidingFitPitch, SlidingFitDataList &slidingFitDataList) const;

    /**
     *  @brief  Get the energy score for a given vertex
     * 
     *  @param  pVertex pointer to the vertex
     *  @param  beamConstants the beam constants
     *  @param  slidingFitDataListU the list of sliding fit data objects in the U view
     *  @param  slidingFitDataListV the list of sliding fit data objects in the V view
     *  @param  slidingFitDataListW the list of sliding fit data objects in the W view
     *  
     *  @return the energy score
     */
    float GetEnergyScore(const pandora::Vertex *const pVertex, const BeamConstants &beamConstants, const SlidingFitDataList &slidingFitDataListU,
        const SlidingFitDataList &slidingFitDataListV, const SlidingFitDataList &slidingFitDataListW) const;
    
    /**
     *  @brief  Increment the energy kick and energy asymmetry for a given vertex in a given view
     * 
     *  @param  pVertex pointer to the vertex
     *  @param  energyKick the energy kick to increment
     *  @param  energyAsymmetry the energy asymmetry to increment
     *  @param  hitType the view hit type
     *  @param  slidingFitDataList the list of sliding fit data objects in this view
     *  
     *  @return the energy score
     */
    void IncrementEnergyScoresForView(const pandora::Vertex *const pVertex, float &energyKick, float &energyAsymmetry, const pandora::HitType hitType,
        const SlidingFitDataList &slidingFitDataList) const;
                                                      
    /**
     *  @brief  Increment the parameters used to calculate the energy kick for a given cluster and a given vertex
     * 
     *  @param  distanceToStart distance from the vertex to the min layer of the cluster's sliding fit
     *  @param  distanceToEnd distance from the vertex to the max layer of the cluster's sliding fit
     *  @param  pointToFitStartVector vector from vertex position to min layer of cluster's sliding fit
     *  @param  pointToFitEndVector vector from vertex position to max layer of cluster's sliding fit
     *  @param  axisDirection the cluster axis direction
     *  @param  useEnergyMetrics whether to use the energy metrics, or to revert to hit-based metrics
     *  @param  totEnergyKick the total energy kick to increment
     *  @param  totEnergy the total energy to increment
     *  @param  totHitKick the total hit kick to increment
     *  @param  totHits the total number of hits to increment
     *  @param  pCluster pointer to the vertex
     */
    void IncrementEnergyKickParameters(const float distanceToStart, const float distanceToEnd, const pandora::CartesianVector &pointToFitStartVector, 
        const pandora::CartesianVector &pointToFitEndVector, const pandora::CartesianVector &axisDirection, bool &useEnergyMetrics, float &totEnergyKick,
        float &totEnergy, float &totHitKick, unsigned int &totHits, const pandora::Cluster *const pCluster) const;
                                                 
    /**
     *  @brief  Increment the parameters used to calculate the energy asymmetry for a given cluster and a given vertex
     * 
     *  @param  vertexPosition2D the projection of the vertex's position into this view
     *  @param  axisDirection the cluster axis direction
     *  @param  useEnergyMetrics whether to use the energy metrics, or to revert to hit-based metrics
     *  @param  localEvtAxisDirEnergy current local event axis using energy weighting
     *  @param  localEvtAxisDirHits the current local event axis using hit weighting
     *  @param  pCluster pointer to the vertex
     *  @param  isViable whether a set of clusters is producing a viable energy asymmetry score (i.e. at no more than 5 deg angle)
     *  @param  asymmetryConsideredClusters the clusters considered in the asymmetry calculation (either one or two)
     */
    void IncrementEnergyAsymmetryParameters(const pandora::CartesianVector &vertexPosition2D, const pandora::CartesianVector &axisDirection, 
        bool &useEnergyMetrics, pandora::CartesianVector &localEvtAxisDirEnergy, pandora::CartesianVector &localEvtAxisDirHits,
        const pandora::Cluster *const pCluster, bool &isViable, pandora::ClusterList &asymmetryConsideredClusters) const;
    
    /**
     *  @brief  Calculate the energy asymmetry for a vertex in a given view using the calculated parameters
     * 
     *  @param  consideredClusters the list of clusters considered for the energy asymmetry calculation
     *  @param  vertexPosition2D the projection of the vertex's position into this view
     *  @param  useEnergyMetrics whether to use the energy metrics, or to revert to hit-based metrics
     *  @param  localEvtAxisDirEnergy current local event axis using energy weighting
     *  @param  localEvtAxisDirHits the current local event axis using hit weighting
     *  @param  isViable whether a set of clusters is producing a viable energy asymmetry score (i.e. at no more than 5 deg angle)
     */
    float CalculateEnergyAsymmetry(const pandora::ClusterList &consideredClusters, const pandora::CartesianVector &vertexPosition2D, 
        const bool useEnergyMetrics, const pandora::CartesianVector &localEvtAxisDirEnergy, const pandora::CartesianVector &localEvtAxisDirHits,
        bool isViable) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    
    unsigned int    m_slidingFitWindow;             ///< The layer window for the sliding linear fits
    float           m_rOffset;                      ///< The r offset parameter in the energy score
    float           m_xOffset;                      ///< The x offset parameter in the energy score
    float           m_epsilon;                      ///< The epsilon parameter in the energy score
    float           m_asymmetryConstant;            ///< The asymmetry constant parameter in the energy score
    unsigned int    m_minNHits;                     ///< The min number of hits parameter in the energy score
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *EnergyKickVertexSelectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new EnergyKickVertexSelectionAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &EnergyKickVertexSelectionAlgorithm::SlidingFitData::GetMinLayerDirection() const
{
    return m_minLayerDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &EnergyKickVertexSelectionAlgorithm::SlidingFitData::GetMaxLayerDirection() const
{
    return m_maxLayerDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &EnergyKickVertexSelectionAlgorithm::SlidingFitData::GetMinLayerPosition() const
{
    return m_minLayerPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &EnergyKickVertexSelectionAlgorithm::SlidingFitData::GetMaxLayerPosition() const
{
    return m_maxLayerPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster *EnergyKickVertexSelectionAlgorithm::SlidingFitData::GetCluster() const
{
    return m_pCluster;
}

} // namespace lar_content

#endif // #ifndef LAR_ENERGY_KICK_VERTEX_SELECTION_ALGORITHM_H