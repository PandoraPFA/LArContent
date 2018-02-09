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

    if (nuSliceHypotheses.size() == 0) return;

    // If using training mode, find the slice with the most neutrino induced hits and see if the event passes the training selection cuts
    // ATTN if not using training mode then these variables are unused.
    unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());
    const bool passesTrainingSelection(m_useTrainingMode ? this->GetBestSliceIndex(pAlgorithm, nuSliceHypotheses, crSliceHypotheses, bestSliceIndex) : false);

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)                                
    {
        SliceFeatures features(nuSliceHypotheses.at(sliceIndex), crSliceHypotheses.at(sliceIndex), this);
        if (!features.IsFeatureVectorAvailable()) continue;

        if (m_useTrainingMode && passesTrainingSelection)
            LArSvmHelper::ProduceTrainingExample(m_trainingOutputFile, sliceIndex == bestSliceIndex, features.GetFeatureVector());

        // TODO Use trained SVM to calculate neutrino probability
    }

    
    // TODO Fill the selected Pfos based on the SVM probabilities
    (void) selectedPfos;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoIdTool::GetBestSliceIndex(const MasterAlgorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, unsigned int &bestSliceIndex) const
{
    unsigned int nHitsInBestSlice(0);
    unsigned int nNuHitsInBestSlice(0);
    unsigned int nuNHitsTotal(0);
   
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

    // ATTN for events with no neutrino induced hits, default neutrino purity and completeness to zero
    float purity(0);
    float completeness(0);

    if (nHitsInBestSlice > 0)
        purity = static_cast<float>(nNuHitsInBestSlice) / static_cast<float>(nHitsInBestSlice);

    if (nuNHitsTotal > 0)
        completeness = static_cast<float>(nNuHitsInBestSlice) / static_cast<float>(nuNHitsTotal);

    return this->PassesQualityCuts(pAlgorithm, purity, completeness);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool NeutrinoIdTool::PassesQualityCuts(const MasterAlgorithm *const pAlgorithm, const float purity, const float completeness) const
{
    if (purity < m_minPurity || completeness < m_minCompleteness)
        return false;

    if (m_selectNuanceCode && (this->GetNuanceCode(pAlgorithm) != m_nuance))
        return false;

    return true;
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
    unsigned int nuNHitsUsedTotal(0);
    unsigned int nuNHitsTotal(0);
    CartesianPointVector nuAllSpacePoints;
    for (const ParticleFlowObject *const pPfo : nuFinalStates)
    {
        CartesianPointVector spacePoints;
        this->GetSpacePoints(pPfo, spacePoints);

        nuAllSpacePoints.insert(nuAllSpacePoints.end(), spacePoints.begin(), spacePoints.end());
        nuNHitsTotal += spacePoints.size();

        // TODO make this number configurable
        if (spacePoints.size() < 5) continue;

        CartesianVector dir(this->GetDirectionFromVertex(spacePoints, nuVertex));
        nuWeightedDirTotal += dir * static_cast<float>(spacePoints.size());
        nuNHitsUsedTotal += spacePoints.size();
    }

    if (nuNHitsUsedTotal == 0) return;
    CartesianVector nuWeightedDir(nuWeightedDirTotal * (1.f / static_cast<float>(nuNHitsUsedTotal)));

    float nuNFinalStatePfos(static_cast<float>(nuFinalStates.size()));
    float nuVertexY(nuVertex.GetY());
    float nuWeightedDirZ(nuWeightedDir.GetZ());
    float nuNSpacePointsInSphere(static_cast<float>(this->GetNPointsInSphere(nuAllSpacePoints, nuVertex, 10))); // TODO Make this configurable
    
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

    if (nCRHitsMax == 0) return;
    if (nCRHitsTotal == 0) return;

    float crFracHitsInLongestTrack = static_cast<float>(nCRHitsMax)/static_cast<float>(nCRHitsTotal);

    // Push the features to the feature vector
    m_featureVector.push_back(nuNFinalStatePfos);
    m_featureVector.push_back(nuNHitsTotal);
    m_featureVector.push_back(nuVertexY);
    m_featureVector.push_back(nuWeightedDirZ);
    m_featureVector.push_back(nuNSpacePointsInSphere);
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

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int NeutrinoIdTool::SliceFeatures::GetNPointsInSphere(const CartesianPointVector &spacePoints, const CartesianVector &vertex, float radius) const
{
    unsigned int nPointsInSphere(0);
    for (const CartesianVector &point : spacePoints)
        if ((point - vertex).GetMagnitudeSquared() <= radius*radius) nPointsInSphere++;

    return nPointsInSphere;
}

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
