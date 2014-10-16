/**
 *  @file   LArContent/include/LArUtility/VertexSelectionAlgorithm.h
 * 
 *  @brief  Header file for the vertex selection algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_SELECTION_ALGORITHM_H
#define LAR_VERTEX_SELECTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  VertexSelectionAlgorithm::Algorithm class
 */
class VertexSelectionAlgorithm : public pandora::Algorithm
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
    VertexSelectionAlgorithm();

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Use hits in clusters (in the provided named list) to fill a provided histogram with hit-vertex relationship information
     * 
     *  @param  pVertex the address of the vertex
     *  @param  hitType the relevant hit type
     *  @param  clusterListName the cluster list name
     *  @param  histogram to receive the populated histogram
     */
    void FillHistogram(const pandora::Vertex *const pVertex, const pandora::HitType hitType, const std::string &clusterListName,
        pandora::Histogram &histogram) const;

    /**
     *  @brief  Use a provided vertex position and cluster to fill a provided histogram with hit-vertex relationship information
     * 
     *  @param  vertexPosition2D the projected vertex position
     *  @param  pCluster the address of the cluster
     *  @param  histogram to receive the populated histogram
     */
    void FillHistogram(const pandora::CartesianVector &vertexPosition2D, const pandora::Cluster *const pCluster,
        pandora::Histogram &histogram) const;

    /**
     *  @brief  Get the figure of merit for a trio of histograms
     * 
     *  @param  histogramU the histogram for the u view
     *  @param  histogramV the histogram for the v view
     *  @param  histogramW the histogram for the w view
     * 
     *  @return the figure of merit
     */
    float GetFigureOfMerit(const pandora::Histogram &histogramU, const pandora::Histogram &histogramV, const pandora::Histogram &histogramW) const;

    /**
     *  @brief  Get the figure of merit contribution for a single histogram
     * 
     *  @param  histogram the histogram
     * 
     *  @return the figure of merit
     */
    float GetFigureOfMerit(const pandora::Histogram &histogram) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_inputClusterListNameU;        ///< The name of the view U cluster list
    std::string             m_inputClusterListNameV;        ///< The name of the view V cluster list
    std::string             m_inputClusterListNameW;        ///< The name of the view W cluster list
    std::string             m_outputVertexListName;         ///< The name under which to save the output vertex list
    bool                    m_replaceCurrentVertexList;     ///< Whether to replace the current vertex list with the output list

    unsigned int            m_histogramNPhiBins;            ///< The number of histogram bins in phi
    float                   m_histogramPhiMin;              ///< The histogram lower phi bound
    float                   m_histogramPhiMax;              ///< The histogram upper phi bound

    float                   m_maxHitVertexDisplacement;     ///< Max hit-vertex displacement for contribution to histograms
    float                   m_hitDeweightingPower;          ///< The hit power used for distance-weighting hit contributions to histograms
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *VertexSelectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new VertexSelectionAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_VERTEX_SELECTION_ALGORITHM_H
