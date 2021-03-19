/**
 *  @file   larpandoracontent/LArMonitoring/TwoViewTransverseTracksValidationTool.cc
 *
 *  @brief  Implementation of the two view transverse tracks validation tool.
 *
 *  $Log: $
 */
#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArMonitoring/TwoViewTransverseTracksValidationTool.h"

#include <limits>

using namespace pandora;

namespace lar_content
{

TwoViewTransverseTracksValidationTool::TwoViewTransverseTracksValidationTool() :
    m_treeName("mytree"), 
    m_outputFileName("output.root")
{
}

TwoViewTransverseTracksValidationTool::~TwoViewTransverseTracksValidationTool()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_outputFileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewTransverseTracksValidationTool::Run(const Cluster *const pCluster1, const Cluster *const pCluster2,
    const DiscreteProbabilityVector &discreteProbabilityVector1, const DiscreteProbabilityVector &discreteProbabilityVector2,
    const TwoViewTransverseOverlapResult &overlapResult)
{

    TreeDataBox treeDataBox;

    this->CollectTruthInformation(pCluster1, treeDataBox, "1");
    this->CollectTruthInformation(pCluster2, treeDataBox, "2");

    StoreAndRegisterDatum(static_cast<int>(this->IsSameMCParticle(pCluster1, pCluster2)), "isSameMCParticle", treeDataBox);

    StoreAndRegisterDatum(static_cast<int>(discreteProbabilityVector1.GetSize()), "discProbVecSize1", treeDataBox);
    StoreAndRegisterDatum(static_cast<int>(discreteProbabilityVector2.GetSize()), "discProbVecSize2", treeDataBox);
    StoreAndRegisterDatum(static_cast<int>(overlapResult.GetNSamplingPoints()), "nSamplingPts", treeDataBox);
    StoreAndRegisterDatum(static_cast<int>(overlapResult.GetNMatchedSamplingPoints()), "nMatchedSamplingPts", treeDataBox);
    StoreAndRegisterDatum(static_cast<int>(overlapResult.GetNReUpsampledSamplingPoints()), "nReUpsampSamplingPts", treeDataBox);
    StoreAndRegisterDatum(static_cast<int>(overlapResult.GetNMatchedReUpsampledSamplingPoints()), "nMatchedReUpsampSamplingPts", treeDataBox);
    StoreAndRegisterDatum(overlapResult.GetCorrelationCoefficient(), "correlationCoef", treeDataBox);
    StoreAndRegisterDatum(overlapResult.GetLocallyMatchedFraction(), "locallyMatchedFrac", treeDataBox);
    StoreAndRegisterDatum(overlapResult.GetMatchingScore(), "matchingScore", treeDataBox);
    StoreAndRegisterDatum(overlapResult.GetTwoViewXOverlap().GetTwoViewXOverlapSpan(), "xOverlapSpan", treeDataBox);
    StoreAndRegisterDatum(overlapResult.GetTwoViewXOverlap().GetXOverlapFraction0(), "xOverlapFrac1", treeDataBox);
    StoreAndRegisterDatum(overlapResult.GetTwoViewXOverlap().GetXOverlapFraction1(), "xOverlapFrac2", treeDataBox);

    CollectChargeProfileInformation(discreteProbabilityVector1, treeDataBox, "1");
    CollectChargeProfileInformation(discreteProbabilityVector2, treeDataBox, "2");

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewTransverseTracksValidationTool::IsSameMCParticle(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    const MCParticle *const pMCParticle1(MCParticleHelper::GetMainMCParticle(pCluster1));
    const MCParticle *const pMCParticle2(MCParticleHelper::GetMainMCParticle(pCluster2));

    if (pMCParticle1 && pMCParticle2)
    {
        return pMCParticle1->GetUid() == pMCParticle2->GetUid();
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewTransverseTracksValidationTool::CollectTruthInformation(const Cluster *const pCluster, TreeDataBox &treeDataBox, std::string clusterName)
{
    constexpr float floatDefault(std::numeric_limits<int>::lowest());

    int foundMCParticle(0);
    int pdgMCParticle(std::numeric_limits<int>::lowest());
    int isPrimaryMCParticle(0);
    float energyMCParticle(floatDefault);
    float momXMCParticle(floatDefault), momYMCParticle(floatDefault), momZMCParticle(floatDefault);
    float posXMCParticle(floatDefault), posYMCParticle(floatDefault), posZMCParticle(floatDefault);

    const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));
    if (pMCParticle)
    {
        foundMCParticle = 1;
        pdgMCParticle = pMCParticle->GetParticleId();
        isPrimaryMCParticle = static_cast<int>(pMCParticle->IsRootParticle());
        energyMCParticle = pMCParticle->GetEnergy();
        momXMCParticle = pMCParticle->GetMomentum().GetX();
        momYMCParticle = pMCParticle->GetMomentum().GetY();
        momZMCParticle = pMCParticle->GetMomentum().GetZ();
        posXMCParticle = pMCParticle->GetVertex().GetX();
        posYMCParticle = pMCParticle->GetVertex().GetY();
        posZMCParticle = pMCParticle->GetVertex().GetZ();
    }

    StoreAndRegisterDatum(foundMCParticle, ("foundMCParticle" + clusterName), treeDataBox);
    StoreAndRegisterDatum(pdgMCParticle, ("pdgMCParticle" + clusterName), treeDataBox);
    StoreAndRegisterDatum(isPrimaryMCParticle, ("isPrimaryMCParticle" + clusterName), treeDataBox);
    StoreAndRegisterDatum(energyMCParticle, ("energyMCParticle" + clusterName), treeDataBox);
    StoreAndRegisterDatum(momXMCParticle, ("momXMCParticle" + clusterName), treeDataBox);
    StoreAndRegisterDatum(momYMCParticle, ("momYMCParticle" + clusterName), treeDataBox);
    StoreAndRegisterDatum(momZMCParticle, ("momZMCParticle" + clusterName), treeDataBox);
    StoreAndRegisterDatum(posXMCParticle, ("posXMCParticle" + clusterName), treeDataBox);
    StoreAndRegisterDatum(posYMCParticle, ("posYMCParticle" + clusterName), treeDataBox);
    StoreAndRegisterDatum(posZMCParticle, ("posZMCParticle" + clusterName), treeDataBox);

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewTransverseTracksValidationTool::CollectChargeProfileInformation(
    const DiscreteProbabilityVector &discreteProbabilityVector, TreeDataBox &treeDataBox, std::string clusterName)
{
    std::string xName("xValues" + clusterName);
    std::string probabilityDensityName("probDensityValues" + clusterName);
    std::string probabilityName("probValues" + clusterName);
    std::string cumulativeProbabilityName("cumulProbValues" + clusterName);
    std::string widthName("widthValues" + clusterName);

    StoreAndRegisterDatum(std::vector<float>(), xName, treeDataBox);
    StoreAndRegisterDatum(std::vector<float>(), probabilityDensityName, treeDataBox);
    StoreAndRegisterDatum(std::vector<float>(), probabilityName, treeDataBox);
    StoreAndRegisterDatum(std::vector<float>(), cumulativeProbabilityName, treeDataBox);
    StoreAndRegisterDatum(std::vector<float>(), widthName, treeDataBox);

    float x, probabilityDensity, cumulativeProbability, width;
    for (unsigned int iBin = 0; iBin < discreteProbabilityVector.GetSize(); iBin++)
    {
        discreteProbabilityVector.GetAllAtIndex(iBin, x, probabilityDensity, cumulativeProbability, width);
        treeDataBox[xName]->emplace_back(x);
        treeDataBox[probabilityDensityName]->emplace_back(probabilityDensity);
        treeDataBox[probabilityName]->emplace_back(discreteProbabilityVector.GetProbability(iBin));
        treeDataBox[cumulativeProbabilityName]->emplace_back(cumulativeProbability);
        treeDataBox[widthName]->emplace_back(width);
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewTransverseTracksValidationTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OutputFileName", m_outputFileName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
