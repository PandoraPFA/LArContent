/**
 *  @file   larpandoracontent/LArMonitoring/SliceIdMonitoringTool.cc
 *
 *  @brief  Implementation of the Slice Id monitoring tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"

#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArMonitoring/SliceIdMonitoringTool.h"

#include <math.h>

using namespace pandora;

namespace lar_content
{

SliceIdMonitoringTool::SliceIdMonitoringTool() : m_minPurity(0.95), m_minSignificance(0.1)
{
}

SliceIdMonitoringTool::~SliceIdMonitoringTool()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ttreed", "outputd.root", "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << " Unable to write tree  to file " << std::endl;
    }
}

void SliceIdMonitoringTool::SelectOutputPfos(const Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
    const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos, const PfoToFloatMap &pfoToProbabilityMap, const SliceVector &sliceVector)
{

    if (this->GetPandora().GetSettings()->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << " with probability map size "
                  << pfoToProbabilityMap.size() << std::endl;

    const CaloHitList *pAllReconstructedCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pAllReconstructedCaloHitList));

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));

    // Obtain map: [mc particle -> primary mc particle]
    LArMCParticleHelper::MCRelationMap mcToPrimaryMCMap;
    LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToPrimaryMCMap);

    // Remove non-reconstructable hits, e.g. those downstream of a neutron
    CaloHitList reconstructableCaloHitList;
    LArMCParticleHelper::PrimaryParameters parameters;
    LArMCParticleHelper::SelectCaloHits(pAllReconstructedCaloHitList, mcToPrimaryMCMap, reconstructableCaloHitList,
        parameters.m_selectInputHits, parameters.m_maxPhotonPropagation);

    const CaloHitList nuNHitsTotal(this->CountNeutrinoInducedHits(reconstructableCaloHitList));
    std::cout << "nu hits in mc " << nuNHitsTotal.size() << std::endl;

    std::list<float> compList;
    std::vector<CaloHitList> caloListList;
    std::vector<bool> isSelectedBestList;
    float completenessTotal{0.f};
    float totalSlicesSize{0};

    const unsigned int nSlices(sliceVector.size());
    std::cout << "nSlices " << nSlices << std::endl;
    for (unsigned int sliceIndex = 0; sliceIndex < nSlices; ++sliceIndex)
    {
        const auto nuHypothesis = nuSliceHypotheses.at(sliceIndex);
        const auto crHypothesis = crSliceHypotheses.at(sliceIndex); // Check that something was reconstructed

        totalSlicesSize = nuHypothesis.size() + crHypothesis.size();

        if (nuHypothesis.empty() && crHypothesis.empty())
            continue;

        // Work out if the slice was identified as neutrino or CR
        const auto isSelectedAsNu =
            !nuHypothesis.empty() && std::find(selectedPfos.begin(), selectedPfos.end(), nuHypothesis.front()) == selectedPfos.begin();

        //get all the hits from one slice
        const CaloHitList caloHitList = sliceVector[sliceIndex];

        //gets the parent hits, the ones from the same worker instence as the mc particles
        CaloHitList parentCaloHitList;
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const CaloHit *const pParentHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
            parentCaloHitList.push_back(pParentHit);
        }

        const CaloHitList nuNHitsHere(this->CountNeutrinoInducedHits(parentCaloHitList));
        std::cout << "total hits in slice " << parentCaloHitList.size() << std::endl;
        std::cout << "nu hits in slice " << nuNHitsHere.size() << std::endl;

        int sharedHits{0};

        for (const CaloHit *const pCaloHit : nuNHitsHere)
        {
            if (std::find(nuNHitsTotal.begin(), nuNHitsTotal.end(), pCaloHit) != nuNHitsTotal.end())
                ++sharedHits;
        }
        std::cout << "shared hits: " << sharedHits << std::endl;

        float completeness = static_cast<float>(sharedHits) / nuNHitsTotal.size();
        float purity = static_cast<float>(sharedHits) / nuNHitsHere.size();

        std::cout << "completeness = " << completeness << std::endl;
        std::cout << "purity = " << purity << std::endl;
        std::cout << "isSelectedAsNu = " << isSelectedAsNu << std::endl;
        std::cout << "-----------------------------------" << std::endl;

        compList.push_back(completeness);
        caloListList.push_back(caloHitList);
        isSelectedBestList.push_back(isSelectedAsNu);
        completenessTotal += completeness;
    }

    //find max completeness
    int isBestSelected{0};
    int isa0nan{0};
    if (totalSlicesSize != 0)
    {
        std::list<float>::iterator maxcomp = std::max_element(compList.begin(), compList.end());
        auto index = std::distance(compList.begin(), maxcomp);

        int indexvalue = index;
        CaloHitList bestCaloHitList = caloListList[indexvalue];
        isBestSelected = isSelectedBestList[indexvalue];

        std::cout << "sum comp = " << completenessTotal << std::endl;
        if (completenessTotal == 0.0 || isnan(completenessTotal) == true)
        {
            isBestSelected = 0;
            isa0nan = 1;
            std::cout << "This is a 0/nan: " << isa0nan << std::endl;
        }
    }

    std::cout << "true = " << isBestSelected << std::endl;

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "wasright", isBestSelected));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ttreed", "isa0nan", isa0nan));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ttreed"));
}

CaloHitList SliceIdMonitoringTool::CountNeutrinoInducedHits(const CaloHitList &caloHitList) const
{
    CaloHitList nNuHits(0);
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        try
        {
            if (LArMCParticleHelper::IsNeutrino(LArMCParticleHelper::GetParentMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit))))
                nNuHits.push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }

    return nNuHits;
}

StatusCode SliceIdMonitoringTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinPurity", m_minPurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSignificance", m_minSignificance));
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
