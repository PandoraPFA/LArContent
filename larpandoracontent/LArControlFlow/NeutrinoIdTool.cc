/**
 *  @file   larpandoracontent/LArControlFlow/NeutrinoIdTool.cc
 *
 *  @brief  Implementation of the neutrino id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/NeutrinoIdTool.h"

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArSvmHelper.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

NeutrinoIdTool::NeutrinoIdTool() :
    m_useTrainingMode(false),
    m_selectNuanceCode(false),
    m_minPurity(0.9f),
    m_minCompleteness(0.9f),
    m_minProbability(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SelectOutputPfos(const MasterAlgorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    bool passesTrainingSelection(true);
    unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());

    if (m_useTrainingMode)
    {
        // TODO put all of this into a separate function
        // Find the best slice and check to see if it passes the selection cuts
        float purity(-std::numeric_limits<float>::max());
        float completeness(-std::numeric_limits<float>::max());
        bestSliceIndex = this->GetBestSliceIndex(nuSliceHypotheses, crSliceHypotheses, purity, completeness);
        
        if (purity < m_minPurity || completeness < m_minCompleteness)
            passesTrainingSelection = false;

        if (m_selectNuanceCode && (this->GetNuanceCode(pAlgorithm) != m_nuance))
            passesTrainingSelection = false;
    }

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)                                
    {
        SliceFeatures features(nuSliceHypotheses.at(sliceIndex), crSliceHypotheses.at(sliceIndex), this);

        if (features.IsFeatureVectorAvailable())
        {
            if (m_useTrainingMode && passesTrainingSelection)
                LArSvmHelper::ProduceTrainingExample(m_trainingOutputFile, sliceIndex == bestSliceIndex, features.GetFeatureVector());
        }
        else
        {
            std::cout << "Can't use slice " << sliceIndex << ", one or more features are incalculable" << std::endl;
        }
        // TODO Use trained SVM to calculate neutrino probability
    }

    
    // TODO Fill the selected Pfos based on the SVM probabilities
    (void) selectedPfos;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int NeutrinoIdTool::GetBestSliceIndex(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, float &purity, float &completeness) const
{
    unsigned int bestSliceIndex(0);
    unsigned int nNuHitsInBestSlice(0);
    unsigned int nHitsInBestSlice(0);
    unsigned int nuNHitsTotal(0);
   
    // ATTN for events with no neutrino induced hits, default neutrino purity and completeness to zero
    purity = 0;
    completeness = 0;

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        CaloHitList reconstructedHits;
        this->Collect2DHits(nuSliceHypotheses.at(sliceIndex), reconstructedHits);
        this->Collect2DHits(crSliceHypotheses.at(sliceIndex), reconstructedHits);

        const unsigned int nNuHits(this->CountNeutrinoInducedHits(reconstructedHits));
        nuNHitsTotal += nNuHits;
        
        // TODO may need to be careful with reproducibility here?
        if (nNuHits > nNuHitsInBestSlice)
        {
            nNuHitsInBestSlice = nNuHits;
            nHitsInBestSlice = reconstructedHits.size();
            bestSliceIndex = sliceIndex; 
        }
    }

    if (nHitsInBestSlice > 0)
        purity = static_cast<float>(nNuHitsInBestSlice) / static_cast<float>(nHitsInBestSlice);

    if (nuNHitsTotal > 0)
        completeness = static_cast<float>(nNuHitsInBestSlice) / static_cast<float>(nuNHitsTotal);

    return bestSliceIndex;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::Collect2DHits(const PfoList &pfos, CaloHitList &hitList) const
{
    CaloHitList collectedHits;
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_U, collectedHits);
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_V, collectedHits);
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_W, collectedHits);

    for (const CaloHit *const pCaloHit : collectedHits)
    {
        // ATTN hits collected from Pfos are copies of hits passed from master instance, we need to access their parent to use MC info
        const CaloHit *pParentCaloHit(static_cast<const CaloHit *>(pCaloHit->GetParentAddress()));

        // Ensure no hits have been double counted
        if (std::find(hitList.begin(), hitList.end(), pParentCaloHit) == hitList.end())
            hitList.push_back(pParentCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int NeutrinoIdTool::CountNeutrinoInducedHits(const CaloHitList &hitList) const
{
    unsigned int nNuHits(0);
    for (const CaloHit *const pCaloHit : hitList)
    {
        try
        {
            if (LArMCParticleHelper::IsNeutrino(LArMCParticleHelper::GetParentMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit))))
                nNuHits++;
        }
        catch (StatusCodeException &)
        {
        }
    }

    return nNuHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int NeutrinoIdTool::GetNuanceCode(const MasterAlgorithm *const pAlgorithm) const
{
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));

    MCParticleVector trueNeutrinos;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, trueNeutrinos);

    if (trueNeutrinos.size() != 1)
    {
        std::cout << "NeutrinoIdTool::GetNuanceCode - Error: number of true neutrinos in event must be exactly one" << std::endl;    
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
    }

    return LArMCParticleHelper::GetNuanceCode(trueNeutrinos.front());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

// TODO think about how to make this function cleaner when features are more established
NeutrinoIdTool::SliceFeatures::SliceFeatures(const PfoList &nuPfos, const PfoList &crPfos, const NeutrinoIdTool *const pTool) : 
    m_isAvailable(false), m_pTool(pTool)
{
    const ParticleFlowObject *const pNeutrino(this->GetNeutrino(nuPfos));
    const CartesianVector nuVertex(LArPfoHelper::GetVertex(pNeutrino)->GetPosition());

    // Neutrino features
    const PfoList &nuFinalStates(pNeutrino->GetDaughterPfoList());
    CartesianVector nuWeightedDirTotal(0, 0, 0);
    unsigned int nuNHitsTotal(0);
    for (const ParticleFlowObject *const pPfo : nuFinalStates)
    {
        CartesianPointVector spacePoints;
        this->GetSpacePoints(pPfo, spacePoints);

        // TODO make this number configurable
        if (spacePoints.size() < 5) continue;

        nuNHitsTotal += spacePoints.size();

        CartesianVector dir(this->GetDirectionFromVertex(spacePoints, nuVertex));
        nuWeightedDirTotal += dir * static_cast<float>(spacePoints.size());
    }

    if (nuNHitsTotal == 0) return;
    CartesianVector nuWeightedDir(nuWeightedDirTotal * (1.f / static_cast<float>(nuNHitsTotal)));

    float nuNFinalStatePfos(static_cast<float>(nuFinalStates.size()));
    float nuVertexY(nuVertex.GetY());
    float nuWeightedDirZ(nuWeightedDir.GetZ());
    
    // Cosmic-ray features
    unsigned int nCRHitsMax(0);
    unsigned int nCRHitsTotal(0);
    float crLongestTrackDirY(std::numeric_limits<float>::max());
    CartesianVector ceiling(0, std::numeric_limits<float>::max(), 0);
    for (const ParticleFlowObject *const pPfo : crPfos)
    {
        CartesianPointVector spacePoints;                                                                                                    
        this->GetSpacePoints(pPfo, spacePoints);
        
        nCRHitsTotal += spacePoints.size();
        
        // TODO make this number configurable
        if (spacePoints.size() < 5) continue;
        
        if (spacePoints.size() > nCRHitsMax)
        {
            nCRHitsMax = spacePoints.size();
            crLongestTrackDirY = this->GetDirectionFromVertex(spacePoints, ceiling).GetY();
        }
    }

    if (nCRHitsTotal == 0) return;
    float crFracHitsInLongestTrack = static_cast<float>(nCRHitsMax)/static_cast<float>(nCRHitsTotal);

    // Push the features to the feature vector
    m_featureVector.push_back(nuNFinalStatePfos);
    m_featureVector.push_back(nuNHitsTotal);
    m_featureVector.push_back(nuVertexY);
    m_featureVector.push_back(nuWeightedDirZ);
    m_featureVector.push_back(crLongestTrackDirY);
    m_featureVector.push_back(crFracHitsInLongestTrack);

    m_isAvailable = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoIdTool::SliceFeatures::IsFeatureVectorAvailable() const
{
    return m_isAvailable;
}

//------------------------------------------------------------------------------------------------------------------------------------------
        
SupportVectorMachine::DoubleVector NeutrinoIdTool::SliceFeatures::GetFeatureVector() const
{
    if (!m_isAvailable)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return m_featureVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ParticleFlowObject *NeutrinoIdTool::SliceFeatures::GetNeutrino(const PfoList &nuPfos)
{
    // ATTN we should only ever have one neutrino reconstructed per slice
    if (nuPfos.size() != 1)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
        
    return nuPfos.front();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SliceFeatures::GetSpacePoints(const ParticleFlowObject *const pPfo, CartesianPointVector &spacePoints) const
{
    ClusterList clusters3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);
    
    if (clusters3D.size() > 1)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
                                                                                                                          
    if (clusters3D.size() == 0)
        return;

    CaloHitList caloHits;
    clusters3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHits);

    for (const CaloHit *const pCaloHit : caloHits)
        spacePoints.push_back(pCaloHit->GetPositionVector());
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector NeutrinoIdTool::SliceFeatures::GetDirectionFromVertex(const CartesianPointVector &spacePoints, const CartesianVector &vertex) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pFirstLArTPC(m_pTool->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float layerPitch(pFirstLArTPC->GetWirePitchW());

    const ThreeDSlidingFitResult fit(&spacePoints, 5, layerPitch);                                                                   
    const CartesianVector endMin(fit.GetGlobalMinLayerPosition());                                                                         
    const CartesianVector endMax(fit.GetGlobalMaxLayerPosition());                                                                         
    const CartesianVector dirMin(fit.GetGlobalMinLayerDirection());                                                                        
    const CartesianVector dirMax(fit.GetGlobalMaxLayerDirection());     

    const bool isMinStart((endMin - vertex).GetMagnitude() < (endMax - vertex).GetMagnitude());
    const CartesianVector startPoint(isMinStart ? endMin : endMax);
    const CartesianVector endPoint(isMinStart ? endMax : endMin);
    const CartesianVector startDir(isMinStart ? dirMin : dirMax);

    const bool shouldFlip((endPoint - startPoint).GetUnitVector().GetDotProduct(startDir) < 0);
    return (shouldFlip ? startDir*(-1) : startDir);
}

/*

//------------------------------------------------------------------------------------------------------------------------------------------

NeutrinoIdTool::~NeutrinoIdTool()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "BestSlices", "Slices.root", "UPDATE"));
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "AllSlices", "Slices.root", "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::SelectOutputPfos(const MasterAlgorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    m_eventId++;

    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        const PfoList &sliceOutput((m_selectAllNeutrinos || (m_selectOnlyFirstSliceNeutrinos && (0 == sliceIndex))) ?
            nuSliceHypotheses.at(sliceIndex) : crSliceHypotheses.at(sliceIndex));

        selectedPfos.insert(selectedPfos.end(), sliceOutput.begin(), sliceOutput.end());
    }
        
    // -------------------------------------------------------------------------------------------------------------------------------------

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*pAlgorithm, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*pAlgorithm, "CaloHitList2D", pCaloHitList));
    
    // Identify reconstructable MCParticles, and get mappings to their good hits
    LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, LArMCParticleHelper::PrimaryParameters(), LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);
    
    int nuanceCode(-1);
    if (!nuMCParticlesToGoodHitsMap.empty())
        nuanceCode = LArMCParticleHelper::GetNuanceCode(LArMCParticleHelper::GetParentNeutrino(nuMCParticlesToGoodHitsMap.begin()->first));

    unsigned int nMuon(this->CountParticlesByPdg(nuMCParticlesToGoodHitsMap, MU_MINUS));
    unsigned int nProton(this->CountParticlesByPdg(nuMCParticlesToGoodHitsMap, PROTON));
    unsigned int nPiPlus(this->CountParticlesByPdg(nuMCParticlesToGoodHitsMap, PI_PLUS));
    unsigned int nPiMinus(this->CountParticlesByPdg(nuMCParticlesToGoodHitsMap, PI_MINUS));
    unsigned int nPhoton(this->CountParticlesByPdg(nuMCParticlesToGoodHitsMap, PHOTON));

    // -------------------------------------------------------------------------------------------------------------------------------------
    
    float totalNeutrinoWeight(0);
    unsigned int bestSliceIndex(0);
    float bestPurity(0); 
    float bestNeutrinoWeight(-std::numeric_limits<float>::max());

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        CaloHitList reconstructedHitsInSlice;
        this->GetReconstructedHitsInSlice(nuSliceHypotheses.at(sliceIndex), crSliceHypotheses.at(sliceIndex), reconstructedHitsInSlice);

        float neutrinoWeight, totalWeight;
        LArMCParticleHelper::GetNeutrinoWeight(&reconstructedHitsInSlice, neutrinoWeight, totalWeight);

        if (totalWeight <= std::numeric_limits<float>::epsilon())
            continue;

        if (neutrinoWeight > bestNeutrinoWeight)
        {
            bestNeutrinoWeight = neutrinoWeight;
            bestPurity = neutrinoWeight / totalWeight;
            bestSliceIndex = sliceIndex;
        }

        totalNeutrinoWeight += neutrinoWeight;
    }

    float bestCompleteness(-1);
    if (totalNeutrinoWeight > std::numeric_limits<float>::epsilon())
    {        
        bestCompleteness = bestNeutrinoWeight / totalNeutrinoWeight;
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BestSlices", "mc_nuance", nuanceCode));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BestSlices", "mc_nMuon", static_cast<int>(nMuon)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BestSlices", "mc_nProton", static_cast<int>(nProton)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BestSlices", "mc_nPiPlus", static_cast<int>(nPiPlus)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BestSlices", "mc_nPiMinus", static_cast<int>(nPiMinus)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BestSlices", "mc_nPhoton", static_cast<int>(nPhoton)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BestSlices", "purity", bestPurity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "BestSlices", "completeness", bestCompleteness));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), "BestSlices"));
   
    if (bestPurity < 0.9 || bestCompleteness < 0.9)
        std::cout << "CAUTION - THIS EVENT DOESN'T HAVE A CLEAN SLICE!" << std::endl;
    
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pFirstLArTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float layerPitch(pFirstLArTPC->GetWirePitchW());

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "mc_nuance", nuanceCode));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "mc_nMuon", static_cast<int>(nMuon)));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "mc_nProton", static_cast<int>(nProton)));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "mc_nPiPlus", static_cast<int>(nPiPlus)));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "mc_nPiMinus", static_cast<int>(nPiMinus)));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "mc_nPhoton", static_cast<int>(nPhoton)));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "bestSlicePurity", bestPurity));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "bestSliceCompleteness", bestCompleteness));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "isBestSlice", (sliceIndex == bestSliceIndex) ? 1 : 0));

        if (sliceIndex == bestSliceIndex)
            std::cout << "NEUTRINO SLICE!" << std::endl;

        const PfoList &nuPfos(nuSliceHypotheses.at(sliceIndex));
        const PfoList &crPfos(crSliceHypotheses.at(sliceIndex));

        // ATTN we should only ever have one neutrino reconstructed per slice
        if (nuPfos.size() != 1)
            throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
        
        const ParticleFlowObject *const pNeutrino(nuPfos.front());

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_nPfos", static_cast<int>(pNeutrino->GetNDaughterPfos())));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "cr_nPfos", static_cast<int>(crPfos.size())));

        CartesianVector vertex(LArPfoHelper::GetVertex(pNeutrino)->GetPosition());

        // Neutrino primaries loop
        CartesianVector nu_weightedDir(0, 0, 0);
        CartesianVector nu_unweightedDir(0, 0, 0);
        unsigned int nHitsTotal(0);
        unsigned int nPfos(0);
        float maxStartPointDistance(-std::numeric_limits<float>::max());
        for (const ParticleFlowObject *const pPfo : pNeutrino->GetDaughterPfoList())
        {
            ClusterList clusters3D;
            LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

            // ATTN Only uses first cluster with hits of type TPC_3D
            if (clusters3D.empty())
                continue;
            
            unsigned int n3DHits(clusters3D.front()->GetNCaloHits());
            if (n3DHits < 15)
                continue;

            ThreeDSlidingFitResult fit(clusters3D.front(), 5, layerPitch);
            CartesianVector endMin(fit.GetGlobalMinLayerPosition());
            CartesianVector endMax(fit.GetGlobalMaxLayerPosition());
            CartesianVector dirMin(fit.GetGlobalMinLayerDirection());
            CartesianVector dirMax(fit.GetGlobalMaxLayerDirection());

            bool isMinStart((endMin - vertex).GetMagnitude() < (endMax - vertex).GetMagnitude());
            CartesianVector startPoint(isMinStart ? endMin : endMax);
            CartesianVector endPoint(isMinStart ? endMax : endMin);
            CartesianVector startDir(isMinStart ? dirMin : dirMax);
            CartesianVector endDir(isMinStart ? dirMax : dirMin);

            float startPointDistance((startPoint - vertex).GetMagnitude());
            if (startPointDistance > maxStartPointDistance)
            {
                maxStartPointDistance = startPointDistance;
            }

            if (startPointDistance < 10)
            {
                nu_weightedDir += startDir * static_cast<float>(n3DHits);
                nu_unweightedDir += startDir;
                nHitsTotal += n3DHits;
                nPfos++;
            }
        }

        unsigned int nGenerations(this->GetNGenerations(pNeutrino, 0));
        std::cout << "nGenerations = " << nGenerations << std::endl;

        CartesianVector nu_meanDir(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        float anisotropy(std::numeric_limits<float>::max());
        float linearity(std::numeric_limits<float>::max());

        if (nHitsTotal > 0)
        {
            nu_meanDir = nu_weightedDir.GetUnitVector();
            anisotropy = nu_weightedDir.GetMagnitude() / static_cast<float>(nHitsTotal);
        }

        if (nPfos > 0)
        {
            linearity = nu_unweightedDir.GetMagnitude() / static_cast<float>(nPfos);
        }

        // ---------------------------------------------------------------------------------------------------------------------------------

        unsigned int maxHits(0);
        float cr_cosThetaY(std::numeric_limits<float>::max());
        for (const ParticleFlowObject *const pPfo : crPfos)
        {
            ClusterList clusters3D;
            LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

            if (clusters3D.size() > 1)
                throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
            
            if (clusters3D.size() == 0) continue;

            CaloHitList threeDHits;
            clusters3D.front()->GetOrderedCaloHitList().FillCaloHitList(threeDHits);
            
            ThreeDSlidingFitResult fit(clusters3D.front(), 15, layerPitch);
            CartesianVector endMin(fit.GetGlobalMinLayerPosition());
            CartesianVector endMax(fit.GetGlobalMaxLayerPosition());
            CartesianVector dirMin(fit.GetGlobalMinLayerDirection());
            CartesianVector dirMax(fit.GetGlobalMaxLayerDirection());

            bool isMinUpper(endMin.GetY() > endMax.GetY());
            CartesianVector dir(isMinUpper ? dirMin : dirMax);

            if (threeDHits.size() > maxHits)
            {
                maxHits = threeDHits.size();
                cr_cosThetaY = dir.GetY();
            }
        }

        // ---------------------------------------------------------------------------------------------------------------------------------
        

            
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_nHitsTotal", static_cast<int>(nHitsTotal)));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_cosWeightedDirZ", nu_meanDir.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_anisotropy", anisotropy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_linearity", linearity));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_maxStartPointDistance", maxStartPointDistance));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_nGenerations", static_cast<int>(nGenerations)));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_vertexX", vertex.GetX()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_vertexY", vertex.GetY()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_vertexZ", vertex.GetZ()));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "cr_cosThetaY", cr_cosThetaY));

    
        // Calculate the hit density
        unsigned int hitsInSphere_10(0);
        unsigned int hitsInSphere_20(0);
        unsigned int hitsInSphere_30(0);
        unsigned int hitsInSphere_40(0);
        unsigned int hitsInSphere_50(0);
        for (const ParticleFlowObject *const pPfo : pNeutrino->GetDaughterPfoList())
        {
            ClusterList clusters3D;
            LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

            if (clusters3D.size() > 1)
                throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);
            
            if (clusters3D.size() == 0) continue;

            CaloHitList threeDHits;
            clusters3D.front()->GetOrderedCaloHitList().FillCaloHitList(threeDHits);
            for (const CaloHit *const pCaloHit : threeDHits)
            {
                const float distToVertex((pCaloHit->GetPositionVector() - vertex).GetMagnitude());
                if (distToVertex < 10.) hitsInSphere_10++;
                if (distToVertex < 20.) hitsInSphere_20++;
                if (distToVertex < 30.) hitsInSphere_30++;
                if (distToVertex < 40.) hitsInSphere_40++;
                if (distToVertex < 50.) hitsInSphere_50++;
            }
        }

        float hitDensity_10(static_cast<float>(hitsInSphere_10) / (10*10*10));
        float hitDensity_20(static_cast<float>(hitsInSphere_20) / (20*20*20));
        float hitDensity_30(static_cast<float>(hitsInSphere_30) / (30*30*30));
        float hitDensity_40(static_cast<float>(hitsInSphere_40) / (40*40*40));
        float hitDensity_50(static_cast<float>(hitsInSphere_50) / (50*50*50));
        
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_hitDensity_10", hitDensity_10));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_hitDensity_20", hitDensity_20));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_hitDensity_30", hitDensity_30));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_hitDensity_40", hitDensity_40));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "AllSlices", "nu_hitDensity_50", hitDensity_50));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), "AllSlices"));


        // BEGIN SCREENSHOTS
        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &nuPfos, "Nu Slice " + sliceIndex, AUTO));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        std::string fileNameNu("screenshot__event-" + std::to_string(m_eventId) + "__slice-" + std::to_string(sliceIndex)  + "__truth-" + ((sliceIndex == bestSliceIndex) ? "best" : "other") + "__reco-nu__passCuts-" + ((bestPurity < 0.9 || bestCompleteness < 0.9) ? "n" : "y") + ".png");
        system(("~asmith/neutrinoId/pandora/test/screenshots/takeScreenshot.sh " + fileNameNu).c_str());

        PANDORA_MONITORING_API(VisualizeParticleFlowObjects(this->GetPandora(), &crPfos, "CR Slice " + sliceIndex, AUTO));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        std::string fileNameCR("screenshot__event-" + std::to_string(m_eventId) + "__slice-" + std::to_string(sliceIndex)  + "__truth-" + ((sliceIndex == bestSliceIndex) ? "best" : "other") + "__reco-cr__passCuts-" + ((bestPurity < 0.9 || bestCompleteness < 0.9) ? "n" : "y") + ".png");
        system(("~asmith/neutrinoId/pandora/test/screenshots/takeScreenshot.sh " + fileNameCR).c_str());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
unsigned int NeutrinoIdTool::CountParticlesByPdg(const LArMCParticleHelper::MCContributionMap &nuMCParticlesToGoodHitsMap, const ParticleType &type) const
{
    unsigned int count(0);
    for (const auto &entry : nuMCParticlesToGoodHitsMap)
        if (entry.first->GetParticleId() == type) count++;

    return count;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int NeutrinoIdTool::GetNGenerations(const ParticleFlowObject *const pPfo, unsigned int generation)
{
    generation++;
    unsigned int maxGeneration(generation);
    for (const ParticleFlowObject *const pDaughter : pPfo->GetDaughterPfoList())
    {
        unsigned int daughterGeneration(this->GetNGenerations(pDaughter, generation));
        maxGeneration = std::max(maxGeneration, daughterGeneration);
    }

    return maxGeneration;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void NeutrinoIdTool::GetReconstructedHitsInSlice(const PfoList &nuPfos, const PfoList &crPfos, CaloHitList &reconstructedHitsInSlice) const
{
    CaloHitList collectedHits;
    LArPfoHelper::GetCaloHits(nuPfos, TPC_VIEW_U, collectedHits);
    LArPfoHelper::GetCaloHits(nuPfos, TPC_VIEW_V, collectedHits);
    LArPfoHelper::GetCaloHits(nuPfos, TPC_VIEW_W, collectedHits);
    LArPfoHelper::GetCaloHits(crPfos, TPC_VIEW_U, collectedHits);
    LArPfoHelper::GetCaloHits(crPfos, TPC_VIEW_V, collectedHits);
    LArPfoHelper::GetCaloHits(crPfos, TPC_VIEW_W, collectedHits);

    // Ensure no hits have been double counted
    for (const CaloHit *const pCaloHit : collectedHits)
    {
        if (std::find(reconstructedHitsInSlice.begin(), reconstructedHitsInSlice.end(), pCaloHit) == reconstructedHitsInSlice.end())
        {
            // ATTN Hits collected from Pfos are copies of hit passed from master instance, we need to access their parent to use MC information
            const CaloHit *pParentCaloHit(static_cast<const CaloHit *>(pCaloHit->GetParentAddress()));

            reconstructedHitsInSlice.push_back(pParentCaloHit);
        }
    }
}
*/
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode NeutrinoIdTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseTrainingMode", m_useTrainingMode));

    if (m_useTrainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "TrainingOutputFileName", m_trainingOutputFile));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumPurity", m_minPurity));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumCompleteness", m_minCompleteness));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectNuanceCode", m_selectNuanceCode));

    if (m_selectNuanceCode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "NuanceCode", m_nuance));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumNeutrinoProbability", m_minProbability));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
