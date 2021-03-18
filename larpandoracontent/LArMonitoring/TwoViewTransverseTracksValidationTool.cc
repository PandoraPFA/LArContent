/**
 *  @file   larpandoracontent/LArMonitoring/TwoViewTransverseTracksValidationTool.cc
 *
 *  @brief  Implementation of the two view transverse tracks validation tool.
 *
 *  $Log: $
 */

//#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/TwoViewTransverseTracksValidationTool.h"

using namespace pandora;

namespace lar_content
{

TwoViewTransverseTracksValidationTool::TwoViewTransverseTracksValidationTool() :
    m_treeName("mytree"),
    m_outputFileName("output.root")
{
    //PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "matchtree", "ksuv", var));
}

TwoViewTransverseTracksValidationTool::~TwoViewTransverseTracksValidationTool()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_outputFileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewTransverseTracksValidationTool::Run(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
    const DiscreteProbabilityVector &discreteProbabilityVector1, const DiscreteProbabilityVector &discreteProbabilityVector2,
    const TwoViewTransverseOverlapResult &overlapResult)
{

    TreeDataBox treeDataBox;
    //RegisterTreeDatum(1.f, treeDataBox);

    StoreAndRegisterDatum(static_cast<int>(discreteProbabilityVector1.GetSize()), "distSize1", treeDataBox);

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
   
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewTransverseTracksValidationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
