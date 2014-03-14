/**
 *  @file   LArContent/include/LArMonitoring/NtupleWritingAlgorithm.h
 * 
 *  @brief  Header file for the ntuple writing algorithm
 * 
 *  $Log: $
 */
#ifndef LAR_NTUPLE_WRITING_ALGORITHM_H
#define LAR_NTUPLE_WRITING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar
{

/**
 *  @brief  NtupleWritingAlgorithm class
 */
class NtupleWritingAlgorithm : public pandora::Algorithm
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


    typedef std::vector<int> IntVector;
    typedef std::vector<float> FloatVector;

    class NtupleWriter {

    public:
        NtupleWriter();
        NtupleWriter( std::string m_fileName, std::string m_treeName );
        ~NtupleWriter();

        void NewEntry( int eventNumber, int viewType );

        void SetVertex(const pandora::CartesianVector& theVertex);

        void AddCluster(pandora::Cluster* pCluster, bool IsSeed);

        void WriteEntry();

    private:

        std::string m_fileName;

        std::string m_treeNameU;
        std::string m_treeNameV;
        std::string m_treeNameW;

        bool m_treeExistsU;
        bool m_treeExistsV;
        bool m_treeExistsW;

        int m_event;
        int m_view;

        int m_nHits;
        int m_nClusters;
        int m_nSeedClusters;

        int m_foundVertex;
        float m_vertexPosX;
        float m_vertexPosY;
        float m_vertexPosZ;

        IntVector m_hitID;
        IntVector m_hitClusterID;
        FloatVector m_hitPosX;
        FloatVector m_hitPosY;
        FloatVector m_hitPosZ;
        FloatVector m_hitEnergy;

        IntVector m_clusterID;
        IntVector m_clusterIsSeed;
        IntVector m_clusterDirectionInZ;
        IntVector m_clusterLayers;
        FloatVector m_clusterStartPosX;
        FloatVector m_clusterStartPosY;
        FloatVector m_clusterStartPosZ;
        FloatVector m_clusterStartDirX;
        FloatVector m_clusterStartDirY;
        FloatVector m_clusterStartDirZ;
        FloatVector m_clusterEndPosX;
        FloatVector m_clusterEndPosY;
        FloatVector m_clusterEndPosZ;
        FloatVector m_clusterEndDirX;
        FloatVector m_clusterEndDirY;
        FloatVector m_clusterEndDirZ;
        FloatVector m_clusterEnergy;
        FloatVector m_clusterLength;
    };



    std::string        m_seedClusterListName;
    std::string        m_nonSeedClusterListName;
    std::string        m_vertexName;

    std::string        m_outputFileName;
    std::string        m_outputTreeName;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *NtupleWritingAlgorithm::Factory::CreateAlgorithm() const
{
    return new NtupleWritingAlgorithm();
}

} // namespace lar

#endif // #ifndef LAR_NTUPLE_WRITING_ALGORITHM_H
