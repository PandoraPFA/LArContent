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

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

template <typename T>
NeutrinoIdTool<T>::NeutrinoIdTool() :
    m_useTrainingMode(false),
    m_selectNuanceCode(false),
    m_nuance(-std::numeric_limits<int>::max()),
    m_minPurity(0.9f),
    m_minCompleteness(0.9f),
    m_minProbability(0.0f),
    m_maxNeutrinos(1),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NeutrinoIdTool<T>::SelectOutputPfos(const Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
    const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const unsigned int nSlices(nuSliceHypotheses.size());
    if (nSlices == 0)
        return;

    SliceFeaturesVector sliceFeaturesVector;
    this->GetSliceFeatures(this, nuSliceHypotheses, crSliceHypotheses, sliceFeaturesVector);

    if (m_useTrainingMode)
    {
        // ATTN in training mode, just return everything as a cosmic-ray
        this->SelectAllPfos(pAlgorithm, crSliceHypotheses, selectedPfos);

        unsigned int bestSliceIndex(std::numeric_limits<unsigned int>::max());
        if (!this->GetBestMCSliceIndex(pAlgorithm, nuSliceHypotheses, crSliceHypotheses, bestSliceIndex))
            return;

        for (unsigned int sliceIndex = 0; sliceIndex < nSlices; ++sliceIndex)
        {
            const SliceFeatures &features(sliceFeaturesVector.at(sliceIndex));
            if (!features.IsFeatureVectorAvailable())
                continue;

            LArMvaHelper::MvaFeatureVector featureVector;
            features.GetFeatureVector(featureVector);
            LArMvaHelper::ProduceTrainingExample(m_trainingOutputFile, sliceIndex == bestSliceIndex, featureVector);
        }

        return;
    }

    this->SelectPfosByProbability(pAlgorithm, nuSliceHypotheses, crSliceHypotheses, sliceFeaturesVector, selectedPfos);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NeutrinoIdTool<T>::GetSliceFeatures(const NeutrinoIdTool<T> *const pTool, const SliceHypotheses &nuSliceHypotheses,
    const SliceHypotheses &crSliceHypotheses, SliceFeaturesVector &sliceFeaturesVector) const
{
    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
        sliceFeaturesVector.push_back(SliceFeatures(nuSliceHypotheses.at(sliceIndex), crSliceHypotheses.at(sliceIndex), pTool));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool NeutrinoIdTool<T>::GetBestMCSliceIndex(const Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
    const SliceHypotheses &crSliceHypotheses, unsigned int &bestSliceIndex) const
{
    unsigned int nHitsInBestSlice(0), nNuHitsInBestSlice(0);

    // Get all hits in all slices to find true number of mc hits
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

    const int nuNHitsTotal(this->CountNeutrinoInducedHits(reconstructableCaloHitList));
    const CaloHitSet reconstructableCaloHitSet(reconstructableCaloHitList.begin(), reconstructableCaloHitList.end());

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        CaloHitList reconstructedCaloHitList;
        this->Collect2DHits(crSliceHypotheses.at(sliceIndex), reconstructedCaloHitList, reconstructableCaloHitSet);

        for (const ParticleFlowObject *const pNeutrino : nuSliceHypotheses.at(sliceIndex))
        {
            const PfoList &nuFinalStates(pNeutrino->GetDaughterPfoList());
            this->Collect2DHits(nuFinalStates, reconstructedCaloHitList, reconstructableCaloHitSet);
        }

        const unsigned int nNuHits(this->CountNeutrinoInducedHits(reconstructedCaloHitList));

        if (nNuHits > nNuHitsInBestSlice)
        {
            nNuHitsInBestSlice = nNuHits;
            nHitsInBestSlice = reconstructedCaloHitList.size();
            bestSliceIndex = sliceIndex;
        }
    }

    // ATTN for events with no neutrino induced hits, default neutrino purity and completeness to zero
    const float purity(nHitsInBestSlice > 0 ? static_cast<float>(nNuHitsInBestSlice) / static_cast<float>(nHitsInBestSlice) : 0.f);
    const float completeness(nuNHitsTotal > 0 ? static_cast<float>(nNuHitsInBestSlice) / static_cast<float>(nuNHitsTotal) : 0.f);
    return this->PassesQualityCuts(pAlgorithm, purity, completeness);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool NeutrinoIdTool<T>::PassesQualityCuts(const Algorithm *const pAlgorithm, const float purity, const float completeness) const
{
    if (purity < m_minPurity || completeness < m_minCompleteness)
        return false;
    if (m_selectNuanceCode && (this->GetNuanceCode(pAlgorithm) != m_nuance))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NeutrinoIdTool<T>::Collect2DHits(const PfoList &pfos, CaloHitList &reconstructedCaloHitList, const CaloHitSet &reconstructableCaloHitSet) const
{
    CaloHitList collectedHits;
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_U, collectedHits);
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_V, collectedHits);
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_W, collectedHits);

    for (const CaloHit *const pCaloHit : collectedHits)
    {
        const CaloHit *const pParentHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
        if (!reconstructableCaloHitSet.count(pParentHit))
            continue;

        // Ensure no hits have been double counted
        if (std::find(reconstructedCaloHitList.begin(), reconstructedCaloHitList.end(), pParentHit) == reconstructedCaloHitList.end())
            reconstructedCaloHitList.push_back(pParentHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
unsigned int NeutrinoIdTool<T>::CountNeutrinoInducedHits(const CaloHitList &caloHitList) const
{
    unsigned int nNuHits(0);
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        try
        {
            if (LArMCParticleHelper::IsNeutrino(LArMCParticleHelper::GetParentMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit))))
                nNuHits++;
        }
        catch (const StatusCodeException &)
        {
        }
    }

    return nNuHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
int NeutrinoIdTool<T>::GetNuanceCode(const Algorithm *const pAlgorithm) const
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

template <typename T>
void NeutrinoIdTool<T>::SelectAllPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &hypotheses, PfoList &selectedPfos) const
{
    for (const PfoList &pfos : hypotheses)
    {
        for (const ParticleFlowObject *const pPfo : pfos)
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["NuScore"] = -1.f;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
        }

        this->SelectPfos(pfos, selectedPfos);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NeutrinoIdTool<T>::SelectPfosByProbability(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses,
    const SliceHypotheses &crSliceHypotheses, const SliceFeaturesVector &sliceFeaturesVector, PfoList &selectedPfos) const
{
    // Calculate the probability of each slice that passes the minimum probability cut
    std::vector<UintFloatPair> sliceIndexProbabilityPairs;
    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        const float nuProbability(sliceFeaturesVector.at(sliceIndex).GetNeutrinoProbability(m_mva));

        for (const ParticleFlowObject *const pPfo : crSliceHypotheses.at(sliceIndex))
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["NuScore"] = nuProbability;
	    
	    std::map<std::string, double> featureMap;
	    sliceFeaturesVector.at(sliceIndex).GetFeatureMap(featureMap);

	    for(auto const& [ name, value ] : featureMap)
	      metadata.m_propertiesToAdd[name] = value;
	    
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
        }

        for (const ParticleFlowObject *const pPfo : nuSliceHypotheses.at(sliceIndex))
        {
            object_creation::ParticleFlowObject::Metadata metadata;
            metadata.m_propertiesToAdd["NuScore"] = nuProbability;

	    std::map<std::string, double> featureMap;
	    sliceFeaturesVector.at(sliceIndex).GetFeatureMap(featureMap);

	    for(auto const& [ name, value ] : featureMap)
	      metadata.m_propertiesToAdd[name] = value;

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*pAlgorithm, pPfo, metadata));
        }

        if (nuProbability < m_minProbability)
        {
            this->SelectPfos(crSliceHypotheses.at(sliceIndex), selectedPfos);
            continue;
        }

        sliceIndexProbabilityPairs.push_back(UintFloatPair(sliceIndex, nuProbability));
    }

    // Sort the slices by probability
    std::sort(sliceIndexProbabilityPairs.begin(), sliceIndexProbabilityPairs.end(),
        [](const UintFloatPair &a, const UintFloatPair &b) { return (a.second > b.second); });

    // Select the first m_maxNeutrinos as neutrinos, and the rest as cosmic
    unsigned int nNuSlices(0);
    for (const UintFloatPair &slice : sliceIndexProbabilityPairs)
    {
        if (nNuSlices < m_maxNeutrinos)
        {
            this->SelectPfos(nuSliceHypotheses.at(slice.first), selectedPfos);
            nNuSlices++;
            continue;
        }

        this->SelectPfos(crSliceHypotheses.at(slice.first), selectedPfos);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NeutrinoIdTool<T>::SelectPfos(const PfoList &pfos, PfoList &selectedPfos) const
{
    selectedPfos.insert(selectedPfos.end(), pfos.begin(), pfos.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

// TODO think about how to make this function cleaner when features are more established
template <typename T>
NeutrinoIdTool<T>::SliceFeatures::SliceFeatures(const PfoList &nuPfos, const PfoList &crPfos, const NeutrinoIdTool<T> *const pTool) :
    m_isAvailable(false),
    m_pTool(pTool)
{
    try
    {
        const ParticleFlowObject *const pNeutrino(this->GetNeutrino(nuPfos));
        const CartesianVector &nuVertex(LArPfoHelper::GetVertex(pNeutrino)->GetPosition());
        const PfoList &nuFinalStates(pNeutrino->GetDaughterPfoList());

        // Neutrino features
        CartesianVector nuWeightedDirTotal(0.f, 0.f, 0.f);
        unsigned int nuNHitsUsedTotal(0);
        unsigned int nuNHitsTotal(0);
        CartesianPointVector nuAllSpacePoints;
        for (const ParticleFlowObject *const pPfo : nuFinalStates)
        {
            CartesianPointVector spacePoints;
            this->GetSpacePoints(pPfo, spacePoints);

            nuAllSpacePoints.insert(nuAllSpacePoints.end(), spacePoints.begin(), spacePoints.end());
            nuNHitsTotal += spacePoints.size();

            if (spacePoints.size() < 5)
                continue;

            const CartesianVector dir(this->GetDirectionFromVertex(spacePoints, nuVertex));
            nuWeightedDirTotal += dir * static_cast<float>(spacePoints.size());
            nuNHitsUsedTotal += spacePoints.size();
        }

        if (nuNHitsUsedTotal == 0)
            return;
        const CartesianVector nuWeightedDir(nuWeightedDirTotal * (1.f / static_cast<float>(nuNHitsUsedTotal)));

        CartesianPointVector pointsInSphere;
        this->GetPointsInSphere(nuAllSpacePoints, nuVertex, 10, pointsInSphere);

        CartesianVector centroid(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        LArPcaHelper::EigenValues eigenValues(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
        LArPcaHelper::EigenVectors eigenVectors;
        LArPcaHelper::RunPca(pointsInSphere, centroid, eigenValues, eigenVectors);

        const float nuNFinalStatePfos(static_cast<float>(nuFinalStates.size()));
        const float nuVertexY(nuVertex.GetY());
        const float nuWeightedDirZ(nuWeightedDir.GetZ());
        const float nuNSpacePointsInSphere(static_cast<float>(pointsInSphere.size()));

        if (eigenValues.GetX() <= std::numeric_limits<float>::epsilon())
            return;
        const float nuEigenRatioInSphere(eigenValues.GetY() / eigenValues.GetX());

        // Cosmic-ray features
        unsigned int nCRHitsMax(0);
        unsigned int nCRHitsTotal(0);
        float crLongestTrackDirY(std::numeric_limits<float>::max());
        float crLongestTrackDeflection(-std::numeric_limits<float>::max());

        for (const ParticleFlowObject *const pPfo : crPfos)
        {
            CartesianPointVector spacePoints;
            this->GetSpacePoints(pPfo, spacePoints);

            nCRHitsTotal += spacePoints.size();

            if (spacePoints.size() < 5)
                continue;

            if (spacePoints.size() > nCRHitsMax)
            {
                nCRHitsMax = spacePoints.size();
                const CartesianVector upperDir(this->GetUpperDirection(spacePoints));
                const CartesianVector lowerDir(this->GetLowerDirection(spacePoints));

                crLongestTrackDirY = upperDir.GetY();
                crLongestTrackDeflection = 1.f - upperDir.GetDotProduct(lowerDir * (-1.f));
            }
        }

        if (nCRHitsMax == 0)
            return;
        if (nCRHitsTotal == 0)
            return;

        const float crFracHitsInLongestTrack = static_cast<float>(nCRHitsMax) / static_cast<float>(nCRHitsTotal);

        // Push the features to the feature vector
        m_featureVector.push_back(nuNFinalStatePfos);
        m_featureVector.push_back(nuNHitsTotal);
        m_featureVector.push_back(nuVertexY);
        m_featureVector.push_back(nuWeightedDirZ);
        m_featureVector.push_back(nuNSpacePointsInSphere);
        m_featureVector.push_back(nuEigenRatioInSphere);
        m_featureVector.push_back(crLongestTrackDirY);
        m_featureVector.push_back(crLongestTrackDeflection);
        m_featureVector.push_back(crFracHitsInLongestTrack);
        m_featureVector.push_back(nCRHitsMax);

        m_featureMap["NuNFinalStatePfos"] = nuNFinalStatePfos;
        m_featureMap["NuNHitsTotal"] = nuNHitsTotal;
        m_featureMap["NuVertexY"] = nuVertexY;
        m_featureMap["NuWeightedDirZ"] = nuWeightedDirZ;
        m_featureMap["NuNSpacePointsInSphere"] = nuNSpacePointsInSphere;
        m_featureMap["NuEigenRatioInSphere"] = nuEigenRatioInSphere;
        m_featureMap["CRLongestTrackDirY"] = crLongestTrackDirY;
        m_featureMap["CRLongestTrackDeflection"] = crLongestTrackDeflection;
        m_featureMap["CRFracHitsInLongestTrack"] = crFracHitsInLongestTrack;
        m_featureMap["CRNHitsMax"] = nCRHitsMax;

        m_isAvailable = true;
    }
    catch (StatusCodeException &)
    {
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
bool NeutrinoIdTool<T>::SliceFeatures::IsFeatureVectorAvailable() const
{
    return m_isAvailable;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NeutrinoIdTool<T>::SliceFeatures::GetFeatureVector(LArMvaHelper::MvaFeatureVector &featureVector) const
{
    if (!m_isAvailable)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    featureVector.insert(featureVector.end(), m_featureVector.begin(), m_featureVector.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NeutrinoIdTool<T>::SliceFeatures::GetFeatureMap(std::map<std::string, double> &featureMap) const
{
    if (!m_isAvailable)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    featureMap.insert(m_featureMap.begin(), m_featureMap.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
float NeutrinoIdTool<T>::SliceFeatures::GetNeutrinoProbability(const T &t) const
{
    // ATTN if one or more of the features can not be calculated, then default to calling the slice a cosmic ray
    if (!this->IsFeatureVectorAvailable())
        return 0.f;

    LArMvaHelper::MvaFeatureVector featureVector;
    this->GetFeatureVector(featureVector);
    return LArMvaHelper::CalculateProbability(t, featureVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
const ParticleFlowObject *NeutrinoIdTool<T>::SliceFeatures::GetNeutrino(const PfoList &nuPfos) const
{
    // ATTN we should only ever have one neutrino reconstructed per slice
    if (nuPfos.size() != 1)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    return nuPfos.front();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NeutrinoIdTool<T>::SliceFeatures::GetSpacePoints(const ParticleFlowObject *const pPfo, CartesianPointVector &spacePoints) const
{
    ClusterList clusters3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

    if (clusters3D.size() > 1)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    if (clusters3D.empty())
        return;

    CaloHitList caloHits;
    clusters3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHits);

    for (const CaloHit *const pCaloHit : caloHits)
        spacePoints.push_back(pCaloHit->GetPositionVector());
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
CartesianVector NeutrinoIdTool<T>::SliceFeatures::GetDirectionFromVertex(const CartesianPointVector &spacePoints, const CartesianVector &vertex) const
{
    return this->GetDirection(spacePoints, [&](const CartesianVector &pointA, const CartesianVector &pointB) {
        return ((pointA - vertex).GetMagnitude() < (pointB - vertex).GetMagnitude());
    });
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
CartesianVector NeutrinoIdTool<T>::SliceFeatures::GetUpperDirection(const CartesianPointVector &spacePoints) const
{
    return this->GetDirection(
        spacePoints, [&](const CartesianVector &pointA, const CartesianVector &pointB) { return (pointA.GetY() > pointB.GetY()); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
CartesianVector NeutrinoIdTool<T>::SliceFeatures::GetLowerDirection(const CartesianPointVector &spacePoints) const
{
    return this->GetDirection(
        spacePoints, [&](const CartesianVector &pointA, const CartesianVector &pointB) { return (pointA.GetY() < pointB.GetY()); });
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
CartesianVector NeutrinoIdTool<T>::SliceFeatures::GetDirection(const CartesianPointVector &spacePoints,
    std::function<bool(const CartesianVector &pointA, const CartesianVector &pointB)> fShouldChooseA) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pFirstLArTPC(m_pTool->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float layerPitch(pFirstLArTPC->GetWirePitchW());

    const ThreeDSlidingFitResult fit(&spacePoints, 5, layerPitch);
    const CartesianVector endMin(fit.GetGlobalMinLayerPosition());
    const CartesianVector endMax(fit.GetGlobalMaxLayerPosition());
    const CartesianVector dirMin(fit.GetGlobalMinLayerDirection());
    const CartesianVector dirMax(fit.GetGlobalMaxLayerDirection());

    const bool isMinStart(fShouldChooseA(endMin, endMax));
    const CartesianVector startPoint(isMinStart ? endMin : endMax);
    const CartesianVector endPoint(isMinStart ? endMax : endMin);
    const CartesianVector startDir(isMinStart ? dirMin : dirMax);

    const bool shouldFlip((endPoint - startPoint).GetUnitVector().GetDotProduct(startDir) < 0.f);
    return (shouldFlip ? startDir * (-1.f) : startDir);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void NeutrinoIdTool<T>::SliceFeatures::GetPointsInSphere(const CartesianPointVector &spacePoints, const CartesianVector &vertex,
    const float radius, CartesianPointVector &spacePointsInSphere) const
{
    for (const CartesianVector &point : spacePoints)
    {
        if ((point - vertex).GetMagnitudeSquared() <= radius * radius)
            spacePointsInSphere.push_back(point);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
StatusCode NeutrinoIdTool<T>::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseTrainingMode", m_useTrainingMode));

    if (m_useTrainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinimumPurity", m_minPurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinimumCompleteness", m_minCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SelectNuanceCode", m_selectNuanceCode));

    if (m_selectNuanceCode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NuanceCode", m_nuance));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinimumNeutrinoProbability", m_minProbability));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaximumNeutrinos", m_maxNeutrinos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    if (!m_useTrainingMode)
    {
        std::string mvaName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MvaName", mvaName));

        std::string mvaFileName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MvaFileName", mvaFileName));

        const std::string fullMvaFileName(LArFileHelper::FindFileInPath(mvaFileName, m_filePathEnvironmentVariable));
        m_mva.Initialize(fullMvaFileName, mvaName);
    }

    return STATUS_CODE_SUCCESS;
}

template class NeutrinoIdTool<AdaBoostDecisionTree>;
template class NeutrinoIdTool<SupportVectorMachine>;

} // namespace lar_content
