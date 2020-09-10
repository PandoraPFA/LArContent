/**
 *  @file   ExampleContent/include/ExampleAlgorithms/DirectionAnalysisAlgorithm.h
 * 
 *  @brief  Header file for the access lists algorithm class.
 * 
 *  $Log: $
 */
#ifndef DIRECTION_ANALYSIS_ALGORITHM_H
#define DIRECTION_ANALYSIS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "larpandoracontent/LArControlFlow/TrackDirectionTool.h"

namespace lar_content
{

/**
 *  @brief  DirectionAnalysisAlgorithm class
 */
class DirectionAnalysisAlgorithm : public pandora::Algorithm
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

    DirectionAnalysisAlgorithm();

    ~DirectionAnalysisAlgorithm();

private:
    pandora::StatusCode     Run();

    void                    WriteMCInformation(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList);
    void                    CheckEventType(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, pandora::PfoVector &pfoVector, int &targetNumberMuons, int &targetNumberProtons, int &nOthers, int &targetNumberPfos);
    void                    WriteVertexInformation(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, const pandora::Vertex* const pVertex, pandora::PfoVector &pfoVector);
    float                   GetVertexDR(const pandora::MCParticleList *pMCParticleList, const pandora::Vertex* const pVertex, bool enableSpaceChargeCorrection);
    void                    WritePfoInformation(pandora::PfoVector &pfoVector);
    bool                    IsGoodPfo(const pandora::ParticleFlowObject* pPfo);
    void                    WriteClusterAndHitInformation(pandora::ClusterVector &clusterVector);

    const pandora::Cluster* GetTargetClusterFromPFO(const pandora::ParticleFlowObject* pPfo);
    void                    IsClusterTwoParticles(const pandora::Cluster *const pCluster, TrackDirectionTool::HitChargeVector forwardsFitCharges, TrackDirectionTool::HitChargeVector backwardsFitCharges, bool &isTwoParticles);
    bool                    IsParticleContained(const pandora::MCParticle* pMCParticle);
    bool                    IsRecoParticleContained(const pandora::ParticleFlowObject* pPfo);
    bool                    IsRecoParticleContained(const pandora::Cluster* pCluster);
    bool                    IntersectsYFace(const pandora::MCParticle* pMCParticle);
    bool                    RecoIntersectsYFace(TrackDirectionTool::DirectionFitObject &fitResult);

    void                    WriteToTree(const pandora::Cluster* pCluster, TrackDirectionTool::DirectionFitObject &fitResult, std::string &treeName);
    void                    WriteHitToTree(TrackDirectionTool::HitCharge &hitCharge, std::string &treeName, int mcDirection);

    pandora::StatusCode     ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_mcParticleListName;
    std::string             m_inputHitListName;
    std::string             m_clusterListName;
    std::string             m_vertexListName;
    std::string             m_pfoListName;

    int                     m_targetParticlePDG;

    bool                    m_particleContained;
    bool                    m_cosmic;
    bool                    m_data;
    
    bool                    m_drawFit;

    bool                    m_writeToTree;
    std::string             m_treeName;
    std::string             m_fileName;

    int                     m_fileIdentifier;
    int                     m_eventNumber;

    TrackDirectionTool      *m_pTrackDirectionTool;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *DirectionAnalysisAlgorithm::Factory::CreateAlgorithm() const
{
    return new DirectionAnalysisAlgorithm();
}

} // namespace lar_content

#endif // #ifndef DIRECTION_ANALYSIS_ALGORITHM_H

