/**
 *  @file   larpandoracontent/LArVertex/MvaVertexSelectionAlgorithm.h
 *
 *  @brief  Header file for the mva vertex selection algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MVA_VERTEX_SELECTION_ALGORITHM_H
#define LAR_MVA_VERTEX_SELECTION_ALGORITHM_H 1

#include "Api/PandoraContentApi.h"

#include "larpandoracontent/LArObjects/LArAdaBoostDecisionTree.h"
#include "larpandoracontent/LArObjects/LArSupportVectorMachine.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArVertex/TrainedVertexSelectionAlgorithm.h"

#include <random>

namespace lar_content
{

template <typename, unsigned int>
class KDTreeLinkerAlgo;
template <typename, unsigned int>
class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MvaVertexSelectionAlgorithm class
 */
template <typename T>
class MvaVertexSelectionAlgorithm : public TrainedVertexSelectionAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MvaVertexSelectionAlgorithm();

protected:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

private:
    /**
     *  @brief  Get the vertex score list
     *
     *  @param  vertexVector the vector of vertices
     *  @param  beamConstants the beam constants
     *  @param  kdTreeU the hit kd tree for the U view
     *  @param  kdTreeV the hit kd tree for the V view
     *  @param  kdTreeW the hit kd tree for the W view
     *  @param  vertexScoreList the vertex score list to fill
     */
    void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const;

    /**
     *  @brief  Used a binary classifier to compare a set of vertices and pick the best one
     *
     *  @param  vertexVector the vector of vertices
     *  @param  vertexFeatureInfoMap the vertex feature info map
     *  @param  eventFeatureList the event feature list
     *  @param  t the mva
     *  @param  useRPhi whether to include the r/phi feature
     *
     *  @return address of the best vertex
     */
    const pandora::Vertex *CompareVertices(const pandora::VertexVector &vertexVector, const VertexFeatureInfoMap &vertexFeatureInfoMap,
        const LArMvaHelper::MvaFeatureVector &eventFeatureList, const T &t, const bool useRPhi) const;

    std::string m_filePathEnvironmentVariable; ///< The environment variable providing a list of paths to mva files
    std::string m_mvaFileName;                 ///< The mva file name
    std::string m_regionMvaName;               ///< The name of the region mva to find
    std::string m_vertexMvaName;               ///< The name of the vertex mva to find
    T m_mvaRegion;                             ///< The region mva
    T m_mvaVertex;                             ///< The vertex mva
};

typedef MvaVertexSelectionAlgorithm<AdaBoostDecisionTree> BdtVertexSelectionAlgorithm;
typedef MvaVertexSelectionAlgorithm<SupportVectorMachine> SvmVertexSelectionAlgorithm;

} // namespace lar_content

#endif // #ifndef LAR_MVA_VERTEX_SELECTION_ALGORITHM_H
