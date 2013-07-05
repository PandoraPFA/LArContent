/**
 *  @file   LArContent/include/LArVertex/VertexSplittingAlgorithm.h
 * 
 *  @brief  Header file for the vertex splitting algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_VERTEX_SPLITTING_ALGORITHM_H
#define LAR_VERTEX_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"



namespace lar
{

/**
 *  @brief  VertexSplittingAlgorithm class
 */
class VertexSplittingAlgorithm : public pandora::Algorithm
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

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    void GetListOfCleanClusters(const pandora::ClusterList* const pClusterList, pandora::ClusterVector &clusterVector);  

    std::string         m_inputClusterListName;        ///< The name of the cluster list
      
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *VertexSplittingAlgorithm::Factory::CreateAlgorithm() const
{
    return new VertexSplittingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_VERTEX_SPLITTING_ALGORITHM_H
