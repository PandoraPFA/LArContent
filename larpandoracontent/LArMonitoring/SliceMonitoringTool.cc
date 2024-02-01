/**
 *  @file   larpandoracontent/LArMonitoring/SliceMonitoringTool.cc
 *
 *  @brief  Implementation of the slice monitoring tool.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/SliceMonitoringTool.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include <algorithm>
#include <iterator>

using namespace pandora;

namespace lar_content
{

SliceMonitoringTool::SliceMonitoringTool() {}

//------------------------------------------------------------------------------------------------------------------------------------------

SliceMonitoringTool::~SliceMonitoringTool()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceMonitoringTool::ProcessSlices(const pandora::Algorithm *const pAlgorithm, SlicingAlgorithm::SliceList &inputSliceList)
{

    const MCParticleList *pMCParticleList{nullptr};
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));

    const CaloHitList *pCompleteCaloHitList{nullptr};
    PANDORA_THROW_RESULT_IF(
        STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::GetList(*pAlgorithm, "FullHitList", pCompleteCaloHitList));

    // Populate the complete calo hit list, based on every hit in every slice.
    CaloHitList fullSliceCaloHitList{};
    for (const auto &slice : inputSliceList)
    {
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListU)
            fullSliceCaloHitList.push_back(pSliceCaloHit);
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListV)
            fullSliceCaloHitList.push_back(pSliceCaloHit);
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListW)
            fullSliceCaloHitList.push_back(pSliceCaloHit);
    }

    const MCParticle *pTrueNeutrino{nullptr};
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;

    // Since the slice hits, and the full, unfiltered hits may not match up, first
    // perform a quick match between these two sets of hits.
    std::map<const std::tuple<float, float, float>, const CaloHit*> caloHitListMatchMap;
    auto getHitKey = [](const CaloHit* pCaloHit) -> std::tuple<float, float, float> {
        const auto pos = pCaloHit->GetPositionVector();
        return {pos.GetX(), pos.GetZ(), pCaloHit->GetHadronicEnergy()};
    };
    CaloHitList matchedCaloHitList;

    for (const CaloHit* pCaloHit : fullSliceCaloHitList)
        caloHitListMatchMap[getHitKey(pCaloHit)] = pCaloHit;

    for (const CaloHit *pCaloHit : *pCompleteCaloHitList) {

        // If there is a match for this hit...we want to point the slice-based hit,
        // to the complete calo hit list hit instead.
        // Otherwise, we can't really make connections between the two.
        if (caloHitListMatchMap.count(getHitKey(pCaloHit)) > 0)
            caloHitListMatchMap[getHitKey(pCaloHit)] = pCaloHit;

        // Populate MC info.
        LArCaloHit *pLArCaloHit{
            const_cast<LArCaloHit *>(dynamic_cast<const LArCaloHit *>(pCaloHit))};
        const auto mcWeights = pLArCaloHit->GetMCParticleWeightMap();

        const MCParticle *largestContributor{nullptr};
        float weight = -1;

        if (mcWeights.empty())
            continue;

        for (const auto &mcWeight : mcWeights) {
            const MCParticle *mc{mcWeight.first};
            const auto parent(LArMCParticleHelper::GetParentMCParticle(mc));

            if (mcWeight.second > weight)
                largestContributor = mc;

            if (LArMCParticleHelper::IsNeutrino(parent)) {
                pTrueNeutrino = parent;
                largestContributor = parent;
            }
        }

        if (largestContributor == nullptr)
            continue;

        mcToTrueHitListMap[largestContributor].push_back(pCaloHit);
    }

    if (m_trainingMode)
        WriteOutHits(inputSliceList, mcToTrueHitListMap);

    if (pTrueNeutrino) {
        const float trueNuEnergy{pTrueNeutrino->GetEnergy()};
        const int success{1};

        CaloHitList mcHits;
        if (mcToTrueHitListMap.count(pTrueNeutrino) > 0)
            mcHits = mcToTrueHitListMap.at(pTrueNeutrino);
        mcHits.sort();

        // First, assess the slices. We want to know which is the most-neutrino-filled, largest, etc.
        std::pair<unsigned int, unsigned int> bestSlice({0, 0});
        std::map<unsigned int, CaloHitList> matchedSliceHits;
        for (unsigned int sliceNumber = 0; sliceNumber < inputSliceList.size(); ++sliceNumber)
        {
            const auto slice(inputSliceList[sliceNumber]);
            CaloHitList sliceCaloHits({});

            for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListU)
                sliceCaloHits.push_back(caloHitListMatchMap[getHitKey(pSliceCaloHit)]);
            for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListV)
                sliceCaloHits.push_back(caloHitListMatchMap[getHitKey(pSliceCaloHit)]);
            for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListW)
                sliceCaloHits.push_back(caloHitListMatchMap[getHitKey(pSliceCaloHit)]);

            CaloHitList sliceNuHits;
            sliceCaloHits.sort();
            std::set_intersection(sliceCaloHits.begin(), sliceCaloHits.end(),
                                  mcHits.begin(), mcHits.end(),
                                  std::inserter(sliceNuHits, sliceNuHits.end()));

            if (sliceNuHits.size() > bestSlice.first)
                bestSlice = {sliceNuHits.size(), sliceNumber};

            matchedSliceHits[sliceNumber] = sliceCaloHits;
        }

        // Lets calculate some slice properties.
        // Perform this for every given input slice.
        for (unsigned int sliceNumber = 0; sliceNumber < inputSliceList.size(); ++sliceNumber)
        {
            const auto slice(inputSliceList[sliceNumber]);
            CaloHitList sliceCaloHits(matchedSliceHits[sliceNumber]);
            std::cout << "There is " << sliceCaloHits.size() << " hits to pick from..." << std::endl;

            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "sliceNumber", (int) sliceNumber));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "mostNuHitSlice", (int) bestSlice.second));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "isBestSlice", (int) (sliceNumber == bestSlice.second)));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "trueNuEnergy", trueNuEnergy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(),
                                                   m_treename.c_str(), "trueNuHits",
                                                   (float)mcHits.size()));

            const std::map<HitType, std::string> allViews(
                {{TPC_VIEW_U, "_U"}, {TPC_VIEW_V, "_V"}, {TPC_VIEW_W, "_W"}});

            for (const auto &viewNamePair : allViews) {
                HitType view(viewNamePair.first);
                auto viewName(viewNamePair.second);

                CaloHitList totalNuHitsInView;
                std::copy_if(mcHits.begin(), mcHits.end(),
                             std::back_inserter(totalNuHitsInView),
                             [&](const pandora::CaloHit *hit) {
                             return hit->GetHitType() == view;
                             });
                totalNuHitsInView.sort();

                CaloHitList allCaloHitsInView;
                std::copy_if(sliceCaloHits.begin(), sliceCaloHits.end(),
                             std::back_inserter(allCaloHitsInView),
                             [&](const pandora::CaloHit *hit) {
                             return hit->GetHitType() == view;
                             });
                allCaloHitsInView.sort();

                // Now we have all the hits + all the MC-based (i.e. neutrino hits), get
                // the overlap / missing / extra. First, the hits that match in the slice.
                CaloHitList sliceNuHits;
                std::set_intersection(totalNuHitsInView.begin(), totalNuHitsInView.end(),
                                      allCaloHitsInView.begin(), allCaloHitsInView.end(),
                                      std::inserter(sliceNuHits, sliceNuHits.end()));

                // Now, the MC hits that are missing (i.e. the neutrino hits that are
                // missing).
                CaloHitList missingNuHits;
                std::set_difference(totalNuHitsInView.begin(), totalNuHitsInView.end(),
                                    allCaloHitsInView.begin(), allCaloHitsInView.end(),
                                    std::inserter(missingNuHits, missingNuHits.end()));

                // Finally, the inverse, the hits that are in the slice that aren't
                // associated with the neutrino (cosmic contamination).
                CaloHitList sliceCRHits;
                std::set_difference(allCaloHitsInView.begin(), allCaloHitsInView.end(),
                                    totalNuHitsInView.begin(), totalNuHitsInView.end(),
                                    std::inserter(sliceCRHits, sliceCRHits.end()));

                // Calculate slice completeness (how much of the neutrino is here), and
                // purity (how much CR contamination is there).
                const float containsNeutrinoHits = sliceNuHits.size() > 0;
                float nuComp(0.f);
                float nuPurity(1.f);

                if (containsNeutrinoHits) {
                    nuComp = sliceNuHits.size() / (float)totalNuHitsInView.size();
                    nuPurity = sliceNuHits.size() / (float)allCaloHitsInView.size();
                }

                std::cout << sliceNumber << ": " << nuComp << " / " << nuPurity <<
                        "(" << (sliceNumber == bestSlice.second) <<
                        " / " << (sliceNuHits.size()) <<
                        " / " << (sliceCRHits.size()) <<
                        " / " << (allCaloHitsInView.size()) <<
                        ")" << std::endl;

                PANDORA_MONITORING_API(SetTreeVariable(
                    this->GetPandora(), m_treename.c_str(),
                    "totalNuHitsInView" + viewName, (float)totalNuHitsInView.size()));
                PANDORA_MONITORING_API(SetTreeVariable(
                    this->GetPandora(), m_treename.c_str(), "sliceHits" + viewName,
                    (float)allCaloHitsInView.size()));
                PANDORA_MONITORING_API(
                    SetTreeVariable(this->GetPandora(), m_treename.c_str(),
                                    "sliceNuHits" + viewName, (float)sliceNuHits.size()));
                PANDORA_MONITORING_API(
                    SetTreeVariable(this->GetPandora(), m_treename.c_str(),
                                    "sliceCRHits" + viewName, (float)sliceCRHits.size()));
                PANDORA_MONITORING_API(SetTreeVariable(
                    this->GetPandora(), m_treename.c_str(), "missingNuHits" + viewName,
                    (float)missingNuHits.size()));
                PANDORA_MONITORING_API(
                    SetTreeVariable(this->GetPandora(), m_treename.c_str(),
                                    "containsNeutrino" + viewName, containsNeutrinoHits));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(),
                                                       m_treename.c_str(),
                                                       "nuSliceComp" + viewName, nuComp));
                PANDORA_MONITORING_API(
                    SetTreeVariable(this->GetPandora(), m_treename.c_str(),
                                    "nuSlicePur" + viewName, nuPurity));
            }

            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    } else {
        // Lets calculate some slice properties.
        for (unsigned int sliceNumber = 0; sliceNumber < inputSliceList.size(); ++sliceNumber)
        {
            const int success{0};
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "success", success));
            PANDORA_MONITORING_API(SetTreeVariable(
                this->GetPandora(), m_treename.c_str(), "sliceNumber", (int) sliceNumber));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
        }
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SliceMonitoringTool::WriteOutHits(SlicingAlgorithm::SliceList &inputSliceList, const LArMCParticleHelper::MCContributionMap &mcToTrueHitListMap)
{

    // Build up a full list of all the slicing input hits, whilst retaining the
    // slice that Pandora currently chose.
    std::map<HitType, CaloHitList> perViewSliceHits({});
    std::map<const CaloHit*, int> hitToSliceNumber;

    for (unsigned int sliceNumber = 0; sliceNumber < inputSliceList.size(); ++sliceNumber)
    {
        const auto slice(inputSliceList[sliceNumber]);

        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListU) {
            perViewSliceHits[TPC_VIEW_U].push_back(pSliceCaloHit);
            hitToSliceNumber.insert({pSliceCaloHit, sliceNumber});
        }
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListV) {
            perViewSliceHits[TPC_VIEW_V].push_back(pSliceCaloHit);
            hitToSliceNumber.insert({pSliceCaloHit, sliceNumber});
        }
        for (const CaloHit *const pSliceCaloHit : slice.m_caloHitListW) {
            perViewSliceHits[TPC_VIEW_W].push_back(pSliceCaloHit);
            hitToSliceNumber.insert({pSliceCaloHit, sliceNumber});
        }
    }

    std::map<const CaloHit*, int> hitToPdgCode;
    for (const auto &mcHitListPair : mcToTrueHitListMap)
        for (const auto &hit : mcHitListPair.second)
            hitToPdgCode.insert({hit, mcHitListPair.first->GetParticleId()});

    const std::map<HitType, std::string> allViews({{TPC_VIEW_U, "_U_View"}, {TPC_VIEW_V, "_V_View"}, {TPC_VIEW_W, "_W_View"}});


    for (const auto &viewNamePair : allViews)
    {
        HitType view(viewNamePair.first);
        auto viewName(viewNamePair.second);
        const auto sliceCaloHits(perViewSliceHits[view]);

        LArMvaHelper::MvaFeatureVector featureVector;
        featureVector.emplace_back(static_cast<double>(inputSliceList.size()));
        featureVector.emplace_back(static_cast<double>(sliceCaloHits.size()));

        for (const CaloHit *pCaloHit : sliceCaloHits)
        {
            const float x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()}, adc{pCaloHit->GetMipEquivalentEnergy()};
            const float particlePdg(hitToPdgCode.count(pCaloHit) != 0 ? hitToPdgCode[pCaloHit] : 0);

            featureVector.emplace_back(static_cast<double>(x));
            featureVector.emplace_back(static_cast<double>(z));
            featureVector.emplace_back(static_cast<double>(adc));
            featureVector.emplace_back(static_cast<double>(particlePdg));
            featureVector.emplace_back(static_cast<double>(hitToSliceNumber[pCaloHit]));
        }

        const std::string trainingFileName{m_trainingOutputFile + viewName + ".csv"};
        LArMvaHelper::ProduceTrainingExample(trainingFileName, true, featureVector);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SliceMonitoringTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Filename", m_filename));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "Treename", m_treename));

    if (m_trainingMode)
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
