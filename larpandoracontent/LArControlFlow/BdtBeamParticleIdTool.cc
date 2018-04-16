/**
 *  @file   larpandoracontent/LArControlFlow/BdtBeamParticleIdTool.cc
 *
 *  @brief  Implementation of the beam particle id tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/BdtBeamParticleIdTool.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

using namespace pandora;

namespace lar_content
{

BdtBeamParticleIdTool::BdtBeamParticleIdTool() :
    m_useTrainingMode(false),
    m_trainingOutputFile(""),
    m_minPurity(0.8f),
    m_minCompleteness(0.8f),
    m_selectInputHits(true),
    m_maxPhotonPropagation(2.5f),
    m_adaBoostDecisionTree(AdaBoostDecisionTree()),
    m_filePathEnvironmentVariable("FW_SEARCH_PATH"),
    m_maxNeutrinos(std::numeric_limits<int>::max()),
    m_minAdaBDTScore(0.f)
{
    m_pSliceFeatureParamters = new SliceFeatureParamters();
}

//------------------------------------------------------------------------------------------------------------------------------------------

BdtBeamParticleIdTool::BdtBeamParticleIdTool(const BdtBeamParticleIdTool &rhs) : 
    m_useTrainingMode(rhs.m_useTrainingMode),
    m_trainingOutputFile(rhs.m_trainingOutputFile),
    m_minPurity(rhs.m_minPurity),
    m_minCompleteness(rhs.m_minCompleteness),
    m_selectInputHits(rhs.m_selectInputHits),
    m_maxPhotonPropagation(rhs.m_maxPhotonPropagation),
    m_adaBoostDecisionTree(rhs.m_adaBoostDecisionTree),
    m_filePathEnvironmentVariable(rhs.m_filePathEnvironmentVariable),
    m_maxNeutrinos(rhs.m_maxNeutrinos),
    m_minAdaBDTScore(rhs.m_minAdaBDTScore)
{
    m_pSliceFeatureParamters = new SliceFeatureParamters(*(rhs.m_pSliceFeatureParamters));
}

//------------------------------------------------------------------------------------------------------------------------------------------

BdtBeamParticleIdTool &BdtBeamParticleIdTool::operator=(const BdtBeamParticleIdTool &rhs)
{
    if (this != &rhs)
    {
        m_useTrainingMode = rhs.m_useTrainingMode;
        m_trainingOutputFile = rhs.m_trainingOutputFile;
        m_minPurity = rhs.m_minPurity;
        m_minCompleteness = rhs.m_minCompleteness;
        m_selectInputHits = rhs.m_selectInputHits;
        m_maxPhotonPropagation = rhs.m_maxPhotonPropagation;
        m_adaBoostDecisionTree = rhs.m_adaBoostDecisionTree;
        m_filePathEnvironmentVariable = rhs.m_filePathEnvironmentVariable;
        m_maxNeutrinos = rhs.m_maxNeutrinos;
        m_minAdaBDTScore = rhs.m_minAdaBDTScore;

        m_pSliceFeatureParamters = new SliceFeatureParamters(*(rhs.m_pSliceFeatureParamters));
    }

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

BdtBeamParticleIdTool::~BdtBeamParticleIdTool()
{
    delete m_pSliceFeatureParamters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BdtBeamParticleIdTool::Initialize()
{   
    // Get global TPC geometry information
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);
    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    m_pSliceFeatureParamters->SetTPCGeometryInformation(parentMinX, parentMaxX, parentMinY, parentMaxY, parentMinZ, parentMaxZ);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SelectOutputPfos(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, PfoList &selectedPfos)
{
    if (nuSliceHypotheses.size() != crSliceHypotheses.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const unsigned int nSlices(nuSliceHypotheses.size());
    if (nSlices == 0) return;

    SliceFeaturesVector sliceFeaturesVector;
    this->GetSliceFeatures(this, nuSliceHypotheses, crSliceHypotheses, sliceFeaturesVector);

    if (m_useTrainingMode)
    {
        // ATTN in training mode, just return everything as a cosmic-ray
        this->SelectAllPfos(crSliceHypotheses, selectedPfos);

        pandora::IntVector bestSliceIndicies;
        this->GetBestMCSliceIndicies(pAlgorithm, nuSliceHypotheses, crSliceHypotheses, bestSliceIndicies);

        for (unsigned int sliceIndex = 0; sliceIndex < nSlices; ++sliceIndex)
        {
            const SliceFeatures &features(sliceFeaturesVector.at(sliceIndex));
            if (!features.IsFeatureVectorAvailable()) continue;

            LArMvaHelper::MvaFeatureVector featureVector;
            features.GetFeatureVector(featureVector);

            bool isGoodTrainingSlice(false);
            if (std::find(bestSliceIndicies.begin(), bestSliceIndicies.end(), sliceIndex) != bestSliceIndicies.end())
                isGoodTrainingSlice = true;

            LArMvaHelper::ProduceTrainingExample(m_trainingOutputFile, isGoodTrainingSlice, featureVector);
        }

        return;
    }

    this->SelectPfosByAdaBDTScore(nuSliceHypotheses, crSliceHypotheses, sliceFeaturesVector, selectedPfos);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::GetSliceFeatures(const BdtBeamParticleIdTool *const pTool, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, SliceFeaturesVector &sliceFeaturesVector) const
{
    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
        sliceFeaturesVector.push_back(SliceFeatures(nuSliceHypotheses.at(sliceIndex), crSliceHypotheses.at(sliceIndex), pTool, m_pSliceFeatureParamters));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SelectAllPfos(const SliceHypotheses &hypotheses, PfoList &selectedPfos) const
{   
    for (const PfoList &pfos : hypotheses)
        this->SelectPfos(pfos, selectedPfos);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SelectPfos(const PfoList &pfos, PfoList &selectedPfos) const
{
    selectedPfos.insert(selectedPfos.end(), pfos.begin(), pfos.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::GetBestMCSliceIndicies(const pandora::Algorithm *const pAlgorithm, const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, pandora::IntVector &bestSliceIndicies) const
{
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
    LArMCParticleHelper::SelectCaloHits(pAllReconstructedCaloHitList, mcToPrimaryMCMap, reconstructableCaloHitList, m_selectInputHits, m_maxPhotonPropagation);

    MCParticleToIntMap mcParticleToReconstructableHitsMap;
    this->PopulateMCParticleToHitsMap(mcParticleToReconstructableHitsMap, reconstructableCaloHitList);

    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        // All hits in a slice - No double counting
        CaloHitList reconstructedHits;
        this->Collect2DHits(crSliceHypotheses.at(sliceIndex), reconstructedHits, reconstructableCaloHitList);

        if (nuSliceHypotheses.at(sliceIndex).size() == 1)
        {
            const PfoList &nuFinalStates(nuSliceHypotheses.at(sliceIndex).front()->GetDaughterPfoList());
            this->Collect2DHits(nuFinalStates, reconstructedHits, reconstructableCaloHitList);
        }

        int nRecoHits(reconstructedHits.size());

        // MCParticle to hits in slice map 
        MCParticleToIntMap mcParticleToHitsInSliceMap;
        this->PopulateMCParticleToHitsMap(mcParticleToHitsInSliceMap, reconstructedHits);

        if (mcParticleToHitsInSliceMap.size() == 0)
            continue;

        // Get best mc particle for slice
        const MCParticle *pBestMCParticle(nullptr);
        int nSharedHits(0);

        for (const auto &iter : mcParticleToHitsInSliceMap)
        {
            if (iter.second > nSharedHits)
            {
                pBestMCParticle = iter.first;
                nSharedHits = iter.second;
            }
        }

        // Only consider if target beam particles
        if (2001 != LArMCParticleHelper::GetNuanceCode(pBestMCParticle))
            continue;

        int nMCHits(mcParticleToReconstructableHitsMap.at(pBestMCParticle));
        float purity(nRecoHits > 0 ? static_cast<float>(nSharedHits) / static_cast<float>(nRecoHits) : 0.f);
        float completeness(nMCHits > 0 ? static_cast<float>(nSharedHits) / static_cast<float>(nMCHits) : 0.f);

        if (this->PassesQualityCuts(purity, completeness))
            bestSliceIndicies.push_back(sliceIndex);
    }
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::PopulateMCParticleToHitsMap(MCParticleToIntMap &mcParticleToIntMap, const pandora::CaloHitList &caloHitList) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        try
        {
            const MCParticle *pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit)));

            if (mcParticleToIntMap.find(pParentMCParticle) != mcParticleToIntMap.end())
            {
                mcParticleToIntMap.at(pParentMCParticle) = mcParticleToIntMap.at(pParentMCParticle) + 1;
            }
            else
            {
                mcParticleToIntMap.insert(MCParticleToIntMap::value_type(pParentMCParticle, 1));
            }
        }
        catch (StatusCodeException &statusCodeException)
        {
            if (STATUS_CODE_NOT_INITIALIZED != statusCodeException.GetStatusCode())
                throw statusCodeException;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::Collect2DHits(const PfoList &pfos, CaloHitList &hitList, const CaloHitList &reconstructableCaloHitList) const
{
    CaloHitList collectedHits;
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_U, collectedHits);
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_V, collectedHits);
    LArPfoHelper::GetCaloHits(pfos, TPC_VIEW_W, collectedHits);

    for (const CaloHit *const pCaloHit : collectedHits)
    {
        // ATTN hits collected from Pfos are copies of hits passed from master instance, we need to access their parent to use MC info
        const CaloHit *pParentCaloHit(static_cast<const CaloHit *>(pCaloHit->GetParentAddress()));

        if (std::find(reconstructableCaloHitList.begin(), reconstructableCaloHitList.end(), pParentCaloHit) == reconstructableCaloHitList.end())
            continue;

        // Ensure no hits have been double counted
        if (std::find(hitList.begin(), hitList.end(), pParentCaloHit) == hitList.end())
            hitList.push_back(pParentCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BdtBeamParticleIdTool::PassesQualityCuts(const float purity, const float completeness) const
{
    if (purity < m_minPurity || completeness < m_minCompleteness)
    {
        return false;
    }
    else
    {
        return true;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SelectPfosByAdaBDTScore(const SliceHypotheses &nuSliceHypotheses, const SliceHypotheses &crSliceHypotheses, const SliceFeaturesVector &sliceFeaturesVector, PfoList &selectedPfos) const
{
    // Calculate the probability of each slice that passes the minimum probability cut
    std::vector<UintFloatPair> sliceIndexAdaBDTScorePairs;
    for (unsigned int sliceIndex = 0, nSlices = nuSliceHypotheses.size(); sliceIndex < nSlices; ++sliceIndex)
    {
        const float nuAdaBDTScore(sliceFeaturesVector.at(sliceIndex).GetAdaBoostDecisionTreeScore(m_adaBoostDecisionTree));

        if (nuAdaBDTScore < m_minAdaBDTScore)
        {
            this->SelectPfos(crSliceHypotheses.at(sliceIndex), selectedPfos);
            continue;
        }

        sliceIndexAdaBDTScorePairs.push_back(UintFloatPair(sliceIndex, nuAdaBDTScore));
    }

    // Sort the slices by probability
    std::sort(sliceIndexAdaBDTScorePairs.begin(), sliceIndexAdaBDTScorePairs.end(), [] (const UintFloatPair &a, const UintFloatPair &b) 
    {
        return (a.second > b.second);
    });

    // Select the first m_maxNeutrinos as neutrinos, and the rest as cosmic
    unsigned int nNuSlices(0);
    for (const UintFloatPair &slice : sliceIndexAdaBDTScorePairs)
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
//------------------------------------------------------------------------------------------------------------------------------------------

BdtBeamParticleIdTool::Plane::Plane(const CartesianVector &normal, const CartesianVector &point) : 
    m_unitNormal(normal.GetUnitVector()),
    m_point(point),
    m_d(-1.f * (normal.GetDotProduct(point)))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector BdtBeamParticleIdTool::Plane::GetLineIntersection(const CartesianVector &a0, const CartesianVector &a) const
{
    if (std::fabs(a.GetDotProduct(m_unitNormal)) < std::numeric_limits<float>::min())
        return CartesianVector(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());

    const float denominator(a.GetDotProduct(m_unitNormal));

    if (std::fabs(denominator) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    const float t(-1.f * (a0.GetDotProduct(m_unitNormal) + m_d) / denominator);
    return (a0 + (a * t));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

BdtBeamParticleIdTool::SliceFeatureParamters::SliceFeatureParamters() : 
    m_beamTPCIntersection(0.f, 0.f, 0.f),
    m_beamDirection(0.f, 0.f, 0.f),
    m_selectedFraction(10.f),
    m_nSelectedHits(100)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SliceFeatureParamters::SetTPCGeometryInformation(float &tpcMinX, float &tpcMaxX, float &tpcMinY, float &tpcMaxY, float &tpcMinZ, float &tpcMaxZ) 
{
    m_tpcMinX = tpcMinX;
    m_tpcMaxX = tpcMaxX;
    m_tpcMinY = tpcMinY;
    m_tpcMaxY = tpcMaxY;
    m_tpcMinZ = tpcMinZ;
    m_tpcMaxZ = tpcMaxZ;

    const CartesianVector normalTop(0.f, 0.f, 1.f), pointTop(0.f, 0.f, m_tpcMaxZ);
    const CartesianVector normalBottom(0.f, 0.f, -1.f), pointBottom(0.f, 0.f, m_tpcMinZ);
    const CartesianVector normalRight(1.f, 0.f, 0.f), pointRight(m_tpcMaxX, 0.f, 0.f);
    const CartesianVector normalLeft(-1.f, 0.f, 0.f), pointLeft(m_tpcMinX, 0.f, 0.f);
    const CartesianVector normalBack(0.f, 1.f, 0.f), pointBack(0.f, m_tpcMaxY, 0.f);
    const CartesianVector normalFront(0.f, -1.f, 0.f), pointFront(0.f, m_tpcMinY, 0.f);

    m_tpcPlanes.emplace_back(normalTop, pointTop);
    m_tpcPlanes.emplace_back(normalBottom, pointBottom);
    m_tpcPlanes.emplace_back(normalRight, pointRight);
    m_tpcPlanes.emplace_back(normalLeft, pointLeft);
    m_tpcPlanes.emplace_back(normalBack, pointBack);
    m_tpcPlanes.emplace_back(normalFront, pointFront);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

BdtBeamParticleIdTool::SliceFeatures::SliceFeatures(const PfoList &pfosNu, const PfoList &pfosCr, const BdtBeamParticleIdTool *const pTool, const SliceFeatureParamters *pSliceFeatureParamters) :
    m_isAvailable(false),
    m_pSliceFeatureParamters(new SliceFeatureParamters(*pSliceFeatureParamters)),
    m_pTool(pTool)
{
    try
    {
        PfoList allConnectedPfoListNu, allConnectedPfoListCr;
        LArPfoHelper::GetAllConnectedPfos(pfosNu, allConnectedPfoListNu);
        LArPfoHelper::GetAllConnectedPfos(pfosCr, allConnectedPfoListCr);

        CaloHitList caloHitList3DNu, caloHitList3DCr;
        LArPfoHelper::GetCaloHits(allConnectedPfoListNu, TPC_3D, caloHitList3DNu);
        LArPfoHelper::GetCaloHits(allConnectedPfoListCr, TPC_3D, caloHitList3DCr);

        CaloHitList selectedCaloHitListNu, selectedCaloHitListCr;

        double closestDistanceNu(std::numeric_limits<double>::max()), closestDistanceCr(std::numeric_limits<double>::max());

        this->GetSelectedCaloHits(caloHitList3DNu, selectedCaloHitListNu, closestDistanceNu);
        this->GetSelectedCaloHits(caloHitList3DCr, selectedCaloHitListCr, closestDistanceCr);

        double supplementaryAngleToBeamNu(std::numeric_limits<double>::max()), supplementaryAngleToBeamCr(std::numeric_limits<double>::max());
        double separationNu(std::numeric_limits<double>::max()), separationCr(std::numeric_limits<double>::max());

        LArPcaHelper::EigenValues eigenValuesNu(0.f, 0.f, 0.f);
        LArPcaHelper::EigenValues eigenValuesCr(0.f, 0.f, 0.f);

        if (!selectedCaloHitListNu.empty() && !selectedCaloHitListCr.empty())
        {
            // Beam
            CartesianVector centroidNu(0.f, 0.f, 0.f);
            LArPcaHelper::EigenVectors eigenVecsNu;
            LArPcaHelper::RunPca(selectedCaloHitListNu, centroidNu, eigenValuesNu, eigenVecsNu);

            const CartesianVector &majorAxisNu(eigenVecsNu.front());
            supplementaryAngleToBeamNu = majorAxisNu.GetOpeningAngle(m_pSliceFeatureParamters->GetBeamDirection());

            CartesianVector interceptOneNu(0.f, 0.f, 0.f), interceptTwoNu(0.f, 0.f, 0.f);
            this->GetTPCIntercepts(centroidNu, majorAxisNu, interceptOneNu, interceptTwoNu);

            const double separationOneNu((interceptOneNu - m_pSliceFeatureParamters->GetBeamTPCIntersection()).GetMagnitude());
            const double separationTwoNu((interceptTwoNu - m_pSliceFeatureParamters->GetBeamTPCIntersection()).GetMagnitude());

            separationNu = std::min(separationOneNu, separationTwoNu);

            // Cosmic
            CartesianVector centroidCr(0.f, 0.f, 0.f);
            LArPcaHelper::EigenVectors eigenVecsCr;
            LArPcaHelper::RunPca(selectedCaloHitListCr, centroidCr, eigenValuesCr, eigenVecsCr);

            const CartesianVector &majorAxisCr(eigenVecsCr.front());
            supplementaryAngleToBeamCr = majorAxisCr.GetOpeningAngle(m_pSliceFeatureParamters->GetBeamDirection());

            CartesianVector interceptOneCr(0.f, 0.f, 0.f), interceptTwoCr(0.f, 0.f, 0.f);
            this->GetTPCIntercepts(centroidCr, majorAxisCr, interceptOneCr, interceptTwoCr);

            const double separationOneCr((interceptOneCr - m_pSliceFeatureParamters->GetBeamTPCIntersection()).GetMagnitude());
            const double separationTwoCr((interceptTwoCr - m_pSliceFeatureParamters->GetBeamTPCIntersection()).GetMagnitude());

            separationCr = std::min(separationOneCr, separationTwoCr);

            float maxYNu(-std::numeric_limits<float>::max()), maxYCr(-std::numeric_limits<float>::max());

            for (const CaloHit *pCaloHit : caloHitList3DNu)
            {
                if (maxYNu < pCaloHit->GetPositionVector().GetY())
                    maxYNu = pCaloHit->GetPositionVector().GetY();
            }

            for (const CaloHit *pCaloHit : caloHitList3DCr)
            {
                if (maxYCr < pCaloHit->GetPositionVector().GetY())
                    maxYCr = pCaloHit->GetPositionVector().GetY();
            }

            m_featureVector.push_back(closestDistanceNu);
            m_featureVector.push_back(supplementaryAngleToBeamNu);
            m_featureVector.push_back(separationNu);
            m_featureVector.push_back(eigenValuesNu.GetX());
            m_featureVector.push_back(eigenValuesNu.GetY());
            m_featureVector.push_back(eigenValuesNu.GetZ());
            m_featureVector.push_back(maxYNu);
            m_featureVector.push_back(allConnectedPfoListNu.size());

            m_featureVector.push_back(closestDistanceCr);
            m_featureVector.push_back(supplementaryAngleToBeamCr);
            m_featureVector.push_back(separationCr);
            m_featureVector.push_back(eigenValuesCr.GetX());
            m_featureVector.push_back(eigenValuesCr.GetY());
            m_featureVector.push_back(eigenValuesCr.GetZ());
            m_featureVector.push_back(maxYCr);
            m_featureVector.push_back(allConnectedPfoListCr.size());

            m_isAvailable = true;
        }
    }
    catch (StatusCodeException &)
    {
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

BdtBeamParticleIdTool::SliceFeatures::SliceFeatures(const SliceFeatures &rhs) : 
    m_isAvailable(rhs.m_isAvailable),
    m_featureVector(rhs.m_featureVector),
    m_pTool(rhs.m_pTool)
{
    m_pSliceFeatureParamters = new SliceFeatureParamters(*(rhs.m_pSliceFeatureParamters));
}

//------------------------------------------------------------------------------------------------------------------------------------------

BdtBeamParticleIdTool::SliceFeatures &BdtBeamParticleIdTool::SliceFeatures::operator=(const SliceFeatures &rhs)
{
    if (this != &rhs)
    {
        m_isAvailable = rhs.m_isAvailable;
        m_featureVector = rhs.m_featureVector;
        m_pTool = rhs.m_pTool;
        m_pSliceFeatureParamters = new SliceFeatureParamters(*(rhs.m_pSliceFeatureParamters));        
    }

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

BdtBeamParticleIdTool::SliceFeatures::~SliceFeatures()
{
    delete m_pSliceFeatureParamters;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SliceFeatures::GetSelectedCaloHits(const CaloHitList &inputCaloHitList, CaloHitList &outputCaloHitList,
    double &closestHitToFaceDistance) const
{
    if (inputCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    typedef std::pair<const CaloHit*, float> HitDistancePair;
    typedef std::vector<HitDistancePair> HitDistanceVector;
    HitDistanceVector hitDistanceVector;

    for (const CaloHit *const pCaloHit : inputCaloHitList)
        hitDistanceVector.emplace_back(pCaloHit, (pCaloHit->GetPositionVector() - m_pSliceFeatureParamters->GetBeamTPCIntersection()).GetMagnitudeSquared());

    std::sort(hitDistanceVector.begin(), hitDistanceVector.end(), [](const HitDistancePair &lhs, const HitDistancePair &rhs) -> bool {return (lhs.second < rhs.second);});
    closestHitToFaceDistance = std::sqrt(hitDistanceVector.front().second);

    const unsigned int nInputHits(inputCaloHitList.size());
    const unsigned int nSelectedCaloHits(nInputHits < m_pSliceFeatureParamters->GetNSelectedHits() ? nInputHits :
        static_cast<unsigned int>(std::round(static_cast<float>(nInputHits) * m_pSliceFeatureParamters->GetSelectedFraction() / 100.f + 0.5f)));

    for (const HitDistancePair &hitDistancePair : hitDistanceVector)
    {
        outputCaloHitList.push_back(hitDistancePair.first);

        if (outputCaloHitList.size() >= nSelectedCaloHits)
            break;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void BdtBeamParticleIdTool::SliceFeatures::GetTPCIntercepts(const CartesianVector &a0, const CartesianVector &lineDirection,
    CartesianVector &interceptOne, CartesianVector &interceptTwo) const
{
    CartesianPointVector intercepts;
    const CartesianVector lineUnitVector(lineDirection.GetUnitVector());

    for (const Plane &plane : m_pSliceFeatureParamters->GetPlanes())
    {
        const CartesianVector intercept(plane.GetLineIntersection(a0, lineUnitVector));

        if (this->IsContained(intercept))
            intercepts.push_back(intercept);
    }

    if (intercepts.size() == 2)
    {
        interceptOne = intercepts.at(0);
        interceptTwo = intercepts.at(1);
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BdtBeamParticleIdTool::SliceFeatures::IsContained(const CartesianVector &spacePoint) const
{
    if ((m_pSliceFeatureParamters->GetTPCMinX() - spacePoint.GetX() > std::numeric_limits<float>::epsilon()) || (spacePoint.GetX() - m_pSliceFeatureParamters->GetTPCMaxX() > std::numeric_limits<float>::epsilon()) ||
        (m_pSliceFeatureParamters->GetTPCMinY() - spacePoint.GetY() > std::numeric_limits<float>::epsilon()) || (spacePoint.GetY() - m_pSliceFeatureParamters->GetTPCMaxY() > std::numeric_limits<float>::epsilon()) ||
        (m_pSliceFeatureParamters->GetTPCMinZ() - spacePoint.GetZ() > std::numeric_limits<float>::epsilon()) || (spacePoint.GetZ() - m_pSliceFeatureParamters->GetTPCMaxZ() > std::numeric_limits<float>::epsilon()))
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool BdtBeamParticleIdTool::SliceFeatures::IsFeatureVectorAvailable() const
{
    return m_isAvailable;
}

//------------------------------------------------------------------------------------------------------------------------------------------
        
void BdtBeamParticleIdTool::SliceFeatures::GetFeatureVector(LArMvaHelper::MvaFeatureVector &featureVector) const
{
    if (!m_isAvailable)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    featureVector.insert(featureVector.end(), m_featureVector.begin(), m_featureVector.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------
        
float BdtBeamParticleIdTool::SliceFeatures::GetAdaBoostDecisionTreeScore(const AdaBoostDecisionTree &adaBoostDecisionTree) const
{
    // ATTN if one or more of the features can not be calculated, then default to calling the slice a cosmic ray.  -1.f if the minimum score 
    // possible for a weighted bdt.
    if (!this->IsFeatureVectorAvailable()) return -1.f;

    LArMvaHelper::MvaFeatureVector featureVector;
    this->GetFeatureVector(featureVector);
    return LArMvaHelper::CalculateClassificationScore(adaBoostDecisionTree, featureVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode BdtBeamParticleIdTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    // BDT Settings
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseTrainingMode", m_useTrainingMode));

    if (m_useTrainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "TrainingOutputFileName", m_trainingOutputFile));
    }
    else
    {
        std::string bdtName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "BdtName", bdtName));

        std::string bdtFileName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "BdtFileName", bdtFileName));

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
            "MinAdaBDTScore", m_minAdaBDTScore));

        const std::string fullBdtFileName(LArFileHelper::FindFileInPath(bdtFileName, m_filePathEnvironmentVariable));
        m_adaBoostDecisionTree.Initialize(fullBdtFileName, bdtName);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumPurity", m_minPurity));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinimumCompleteness", m_minCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectInputHits", m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPhotonPropagation", m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaximumNeutrinos", m_maxNeutrinos));

    // Geometry Information for training
    FloatVector beamTPCIntersection;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "BeamTPCIntersection", beamTPCIntersection));
 
    if (3 == beamTPCIntersection.size())
    {
        pandora::CartesianVector beamTPCIntersectionCartesianVector(beamTPCIntersection.at(0), beamTPCIntersection.at(1), beamTPCIntersection.at(2));
        m_pSliceFeatureParamters->SetBeamTPCIntersection(beamTPCIntersectionCartesianVector);
    }
    else if (!beamTPCIntersection.empty())
    {
        std::cout << "BdtBeamParticleIdTool::ReadSettings - invalid BeamTPCIntersection specified " << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    else
    {
        // Default for protoDUNE.
        pandora::CartesianVector beamTPCIntersectionCartesianVector(-33.051, 461.06, 0);
        m_pSliceFeatureParamters->SetBeamTPCIntersection(beamTPCIntersectionCartesianVector);
    }

    FloatVector beamDirection;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "BeamDirection", beamDirection));

    if (3 == beamDirection.size())
    {
        CartesianVector beamDirectionCartesianVector(beamDirection.at(0), beamDirection.at(1), beamDirection.at(2));
        m_pSliceFeatureParamters->SetBeamDirection(beamDirectionCartesianVector);
    }
    else if (!beamDirection.empty())
    {
        std::cout << "BdtBeamParticleIdTool::ReadSettings - invalid BeamDirection specified " << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }
    else
    {
        // Default for protoDUNE.
        const float thetaXZ0(-11.844f * M_PI / 180.f);
        CartesianVector beamDirectionCartesianVector(std::sin(thetaXZ0), 0, std::cos(thetaXZ0));
        m_pSliceFeatureParamters->SetBeamDirection(beamDirectionCartesianVector);
    }

    float selectedFraction(0.f);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectedFraction", selectedFraction));

    if (selectedFraction > std::numeric_limits<float>::epsilon())
        m_pSliceFeatureParamters->SetSelectedFraction(selectedFraction);
    
    unsigned int nSelectedHits(0);

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NSelectedHits", nSelectedHits));

    if (nSelectedHits > 0)
        m_pSliceFeatureParamters->SetNSelectedHits(nSelectedHits);

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
