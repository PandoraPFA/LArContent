/**
 *  @file   larpandoracontent/LArMonitoring/TwoViewTransverseTracksValidationTool.h
 *
 *  @brief  Header file for the two view transverse tracks validation tool
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_TRANSVERSE_TRACKS_VALIDATION_TOOL_H
#define TWO_VIEW_TRANSVERSE_TRACKS_VALIDATION_TOOL_H 1

#include <memory>

#include "PandoraMonitoringApi.h"

#include "Objects/Cluster.h"

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArDiscreteProbabilityVector.h"
#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"
namespace lar_content
{

/**
 *  @brief  TwoViewTransverseTracksValidationTool class
 */
class TwoViewTransverseTracksValidationTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewTransverseTracksValidationTool();

    ~TwoViewTransverseTracksValidationTool();

    bool Run(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        const DiscreteProbabilityVector &discreteProbabilityVector1,
        const DiscreteProbabilityVector &discreteProbabilityVector2,
        const TwoViewTransverseOverlapResult &overlapResult);

private:
    class AbstractDatum
    {
//        public:
//            virtual void RegisterDatumWithTree(const pandora::Pandora *const p_Pandora, const std::string treeName, 
//                const std::string branchName) = 0;
//
    };

    template <typename T> 
    class TreeDatum : public AbstractDatum
    {
        public:
            TreeDatum(T t, const TwoViewTransverseTracksValidationTool *const pTool, const std::string datumName);

        private:
//            void RegisterDatumWithTree(const pandora::Pandora *const p_Pandora, const std::string treeName, 
//                const std::string branchName) /*override*/;

            T m_Datum;
    };

    typedef std::map<std::string, std::unique_ptr<AbstractDatum> > TreeDataBox;

    template <typename T>
    void StoreAndRegisterDatum(const T &t, const std::string datumName, TreeDataBox &treeDataBox);

    template <typename T>
    void RegisterTreeDatum(const T &t, const std::string datumName) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_treeName;            ///< name of output tree
    std::string m_outputFileName;            ///< name of output file
};

template <typename T>
TwoViewTransverseTracksValidationTool::TreeDatum<T>::TreeDatum(T t, const TwoViewTransverseTracksValidationTool *const pTool, 
    const std::string datumName) :
    m_Datum(t)
{
    
    pTool->RegisterTreeDatum(t, datumName);
    //RegisterDatumWithTree(p_Pandora, treeName, branchName, m_Datum);
    //PANDORA_MONITORING_API(SetTreeVariable(pandoraInstance, treeName, branchName, m_Datum));
}

template <typename T>
void TwoViewTransverseTracksValidationTool::StoreAndRegisterDatum(const T &t, const std::string datumName, TreeDataBox &treeDataBox)
{
    treeDataBox.try_emplace(datumName, (new TreeDatum<T>(t, this, datumName)));
    return;
}


template <typename T>
void TwoViewTransverseTracksValidationTool::RegisterTreeDatum(const T &t, const std::string datumName) const
{
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), datumName.c_str(), t));
}

//void TwoViewTransverseTracksValidationTool::TreeDatum::RegisterDatumWithTree(const pandora::Pandora *const p_Pandora, 
//    const std::string treeName, const std::string branchName)
//{
//    PANDORA_MONITORING_API(SetTreeVariable(p_Pandora, treeName, branchName, m_Datum));
//}

} // namespace lar_content

#endif // #ifndef TWO_VIEW_TRANSVERSE_TRACKS_VALIDATION_TOOL_H
