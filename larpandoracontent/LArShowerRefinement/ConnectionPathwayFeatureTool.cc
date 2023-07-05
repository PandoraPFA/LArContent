/**
 *  @file   larpandoracontent/LArShowerRefinement/ConnectionPathwayFeatureTool.cc
 *
 *  @brief  Implementation of the connection pathway feature tool
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArShowerRefinement/ConnectionPathwayFeatureTool.h"
#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

using namespace pandora;

namespace lar_content
{

InitialRegionFeatureTool::InitialRegionFeatureTool() :
    m_defaultFloat(-10.f),
    m_nHitsToConsider(10),
    m_maxInitialGapSizeLimit(4.f),
    m_minLargestGapSizeLimit(2.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void InitialRegionFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
    const ParticleFlowObject *const /*pShowerPfo*/, const CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch,
    const CartesianPointVector & /*showerStarts3D*/)
{
    float initialGapSizeU(m_defaultFloat), initialGapSizeV(m_defaultFloat), initialGapSizeW(m_defaultFloat);
    float largestGapSizeU(m_defaultFloat), largestGapSizeV(m_defaultFloat), largestGapSizeW(m_defaultFloat);
    this->GetViewInitialRegionVariables(pAlgorithm, nuVertex3D, protoShowerMatch, TPC_VIEW_U, initialGapSizeU, largestGapSizeU);
    this->GetViewInitialRegionVariables(pAlgorithm, nuVertex3D, protoShowerMatch, TPC_VIEW_V, initialGapSizeV, largestGapSizeV);
    this->GetViewInitialRegionVariables(pAlgorithm, nuVertex3D, protoShowerMatch, TPC_VIEW_W, initialGapSizeW, largestGapSizeW);

    float minInitialGapSize(m_defaultFloat), middleInitialGapSize(m_defaultFloat), maxInitialGapSize(m_defaultFloat);
    float minLargestGapSize(m_defaultFloat), middleLargestGapSize(m_defaultFloat), maxLargestGapSize(m_defaultFloat);
    LArConnectionPathwayHelper::GetMinMiddleMax(initialGapSizeU, initialGapSizeV, initialGapSizeW, minInitialGapSize, middleInitialGapSize, maxInitialGapSize);
    LArConnectionPathwayHelper::GetMinMiddleMax(largestGapSizeU, largestGapSizeV, largestGapSizeW, minLargestGapSize, middleLargestGapSize, maxLargestGapSize);

    featureVector.push_back(std::min(maxInitialGapSize, m_maxInitialGapSizeLimit));
    featureVector.push_back(std::min(minLargestGapSize, m_minLargestGapSizeLimit));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void InitialRegionFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder, const std::string &featureToolName,
    const Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertex3D,
    const ProtoShowerMatch &protoShowerMatch, const CartesianPointVector &showerStarts3D)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pShowerPfo, nuVertex3D, protoShowerMatch, showerStarts3D);

    if (featureMap.find(featureToolName + "_initialGapSize") != featureMap.end())
    {
        std::cout << "Already wrote initialGapSize feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_initialGapSize");
    featureMap[featureToolName + "_initialGapSize"] = toolFeatureVec[0];

    if (featureMap.find(featureToolName + "_largestGapSize") != featureMap.end())
    {
        std::cout << "Already wrote largestGapSize feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_largestGapSize");
    featureMap[featureToolName + "_largestGapSize"] = toolFeatureVec[1];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void InitialRegionFeatureTool::GetViewInitialRegionVariables(const Algorithm *const pAlgorithm, const CartesianVector &nuVertex3D,
    const ProtoShowerMatch &protoShowerMatch, const HitType hitType, float &initialGapSize, float &largestGapSize) const
{
    const ProtoShower &protoShower(hitType == TPC_VIEW_U
                                       ? protoShowerMatch.GetProtoShowerU()
                                       : (hitType == TPC_VIEW_V ? protoShowerMatch.GetProtoShowerV() : protoShowerMatch.GetProtoShowerW()));
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertex3D, hitType));
    const CartesianVector &startDirection(protoShower.GetConnectionPathway().GetStartDirection());

    // Initial gap size
    CaloHitVector spineHitVector(protoShower.GetSpineHitList().begin(), protoShower.GetSpineHitList().end());

    std::sort(spineHitVector.begin(), spineHitVector.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(nuVertex2D));
    initialGapSize = (nuVertex2D - spineHitVector.front()->GetPositionVector()).GetMagnitude();

    // Largest Gap Size
    largestGapSize = m_defaultFloat;

    std::vector<float> longitudinalProjections;

    for (const CaloHit *const pCaloHit : protoShower.GetSpineHitList())
        longitudinalProjections.push_back(startDirection.GetDotProduct(pCaloHit->GetPositionVector() - nuVertex2D));

    std::sort(longitudinalProjections.begin(), longitudinalProjections.end());

    const unsigned int nIterations(std::min(longitudinalProjections.size(), static_cast<long unsigned int>(m_nHitsToConsider)) - 1);

    for (unsigned int i = 0; i < nIterations; ++i)
        largestGapSize = std::max(std::fabs(longitudinalProjections[i] - longitudinalProjections[i + 1]), largestGapSize);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode InitialRegionFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "nHitsToConsider", m_nHitsToConsider));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DefaultFloat", m_defaultFloat));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxInitialGapSizeLimit", m_maxInitialGapSizeLimit));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinLargestGapSizeLimit", m_minLargestGapSizeLimit));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ConnectionRegionFeatureTool::ConnectionRegionFeatureTool() :
    m_defaultFloat(-10.f),
    m_spineFitWindow(10),
    m_nLayersHalfWindow(5),
    m_pathwayLengthLimit(30.f),
    m_pathwayScatteringAngle2DLimit(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConnectionRegionFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
    const ParticleFlowObject *const /*pShowerPfo*/, const CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch,
    const CartesianPointVector &showerStarts3D)
{
    const float pathwayLength = (nuVertex3D - showerStarts3D.front()).GetMagnitude();
    const float pathwayScatteringAngle2D = this->Get2DKink(pAlgorithm, protoShowerMatch, showerStarts3D.back());

    featureVector.push_back(std::min(pathwayLength, m_pathwayLengthLimit));
    featureVector.push_back(std::min(pathwayScatteringAngle2D, m_pathwayScatteringAngle2DLimit));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConnectionRegionFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder,
    const std::string &featureToolName, const Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch, const CartesianPointVector &showerStarts3D)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pShowerPfo, nuVertex3D, protoShowerMatch, showerStarts3D);

    if (featureMap.find(featureToolName + "_pathwayLength") != featureMap.end())
    {
        std::cout << "Already wrote pathwayLength feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_pathwayLength");
    featureMap[featureToolName + "_pathwayLength"] = toolFeatureVec[0];

    if (featureMap.find(featureToolName + "_pathwayScatteringAngle2D") != featureMap.end())
    {
        std::cout << "Already wrote pathwayScatteringAngle2D feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_pathwayScatteringAngle2D");
    featureMap[featureToolName + "_pathwayScatteringAngle2D"] = toolFeatureVec[1];
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ConnectionRegionFeatureTool::Get2DKink(
    const Algorithm *const pAlgorithm, const ProtoShowerMatch &protoShowerMatch, const CartesianVector &showerStart3D) const
{
    try
    {
        CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

        for (HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            const ProtoShower &protoShower(
                hitType == TPC_VIEW_U ? protoShowerMatch.GetProtoShowerU()
                                      : (hitType == TPC_VIEW_V ? protoShowerMatch.GetProtoShowerV() : protoShowerMatch.GetProtoShowerW()));
            CartesianPointVector &spinePositions(hitType == TPC_VIEW_U ? spinePositionsU : (hitType == TPC_VIEW_V ? spinePositionsV : spinePositionsW));

            for (const CaloHit *const pCaloHit : protoShower.GetSpineHitList())
                spinePositions.push_back(pCaloHit->GetPositionVector());
        }

        const TwoDSlidingFitResult spineFitU(&spinePositionsU, m_spineFitWindow, LArGeometryHelper::GetWirePitch(pAlgorithm->GetPandora(), TPC_VIEW_U));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, m_spineFitWindow, LArGeometryHelper::GetWirePitch(pAlgorithm->GetPandora(), TPC_VIEW_V));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, m_spineFitWindow, LArGeometryHelper::GetWirePitch(pAlgorithm->GetPandora(), TPC_VIEW_W));

        const float scatterAngleU(this->GetLargest2DKinkFromView(pAlgorithm, spineFitU, TPC_VIEW_U, showerStart3D));
        const float scatterAngleV(this->GetLargest2DKinkFromView(pAlgorithm, spineFitV, TPC_VIEW_V, showerStart3D));
        const float scatterAngleW(this->GetLargest2DKinkFromView(pAlgorithm, spineFitW, TPC_VIEW_W, showerStart3D));

        float minScatterAngle(m_defaultFloat), middleScatterAngle(m_defaultFloat), maxScatterAngle(m_defaultFloat);
        LArConnectionPathwayHelper::GetMinMiddleMax(scatterAngleU, scatterAngleV, scatterAngleW, minScatterAngle, middleScatterAngle, maxScatterAngle);

        return middleScatterAngle;
    }
    catch (...)
    {
    }

    return m_defaultFloat;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ConnectionRegionFeatureTool::GetLargest2DKinkFromView(const Algorithm *const pAlgorithm, const TwoDSlidingFitResult &spineFit,
    const HitType hitType, const CartesianVector &showerStart3D) const
{
    const LayerFitResultMap &layerFitResultMap(spineFit.GetLayerFitResultMap());
    const int minLayer(layerFitResultMap.begin()->first), maxLayer(layerFitResultMap.rbegin()->first);
    const int nLayersSpanned(1 + maxLayer - minLayer);

    if (nLayersSpanned <= 2 * m_nLayersHalfWindow)
        return m_defaultFloat;

    const CartesianVector projectedShowerStart(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), showerStart3D, hitType));

    float showerStartL(0.f), showerStartT(0.f);
    spineFit.GetLocalPosition(projectedShowerStart, showerStartL, showerStartT);

    float maxCentralLayer(spineFit.GetLayer(showerStartL) - m_nLayersHalfWindow);
    float minCentralLayer(layerFitResultMap.begin()->first + m_nLayersHalfWindow + 1);

    float highestOpeningAngle(m_defaultFloat);
    CartesianVector kinkPosition(0.f, 0.f, 0.f);
    CartesianVector highestKinkPosition(0.f, 0.f, 0.f);

    for (LayerFitResultMap::const_iterator iter = layerFitResultMap.begin(), iterEnd = layerFitResultMap.end(); iter != iterEnd; ++iter)
    {
        const int layer(iter->first);

        if (layer < minCentralLayer)
            continue;

        if (layer > maxCentralLayer)
            continue;

        bool found(false);
        float openingAngle2D(std::numeric_limits<float>::max());

        float thisHighestOpeningAngle(m_defaultFloat);
        CartesianVector thisHighestKinkPosition(0.f, 0.f, 0.f);

        for (int i = 0; i < m_nLayersHalfWindow; ++i)
        {
            const int testLayer(layer + i);
            const float rL(spineFit.GetL(testLayer));
            const float rL1(spineFit.GetL(testLayer - m_nLayersHalfWindow));
            const float rL2(spineFit.GetL(testLayer + m_nLayersHalfWindow));

            CartesianVector firstPosition(0.f, 0.f, 0.f), centralPosition(0.f, 0.f, 0.f), secondPosition(0.f, 0.f, 0.f);

            if ((STATUS_CODE_SUCCESS != spineFit.GetGlobalFitPosition(rL1, firstPosition)) ||
                (STATUS_CODE_SUCCESS != spineFit.GetGlobalFitPosition(rL, centralPosition)) ||
                (STATUS_CODE_SUCCESS != spineFit.GetGlobalFitPosition(rL2, secondPosition)))
            {
                continue;
            }

            const CartesianVector firstDirection((centralPosition - firstPosition).GetUnitVector());
            const CartesianVector secondDirection((secondPosition - centralPosition).GetUnitVector());

            const float testOpeningAngle2D(secondDirection.GetOpeningAngle(firstDirection) * 180.f / M_PI);

            if (testOpeningAngle2D < openingAngle2D)
            {
                openingAngle2D = testOpeningAngle2D;
                found = true;
            }

            if (testOpeningAngle2D > thisHighestOpeningAngle)
            {
                thisHighestOpeningAngle = openingAngle2D;
                thisHighestKinkPosition = centralPosition;
            }
        }

        if (found)
        {
            if (openingAngle2D > highestOpeningAngle)
            {
                highestOpeningAngle = openingAngle2D;
                highestKinkPosition = thisHighestKinkPosition;
            }
        }
    }

    return highestOpeningAngle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConnectionRegionFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DefaultFloat", m_defaultFloat));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SpineFitWindow", m_spineFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NLayersHalfWindow", m_nLayersHalfWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PathwayLengthLimit", m_pathwayLengthLimit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "PathwayScatteringAngle2DLimit", m_pathwayScatteringAngle2DLimit));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ShowerRegionFeatureTool::ShowerRegionFeatureTool() :
    m_defaultFloat(-10.f),
    m_defaultRatio(-0.5f),
    m_spineFitWindow(20),
    m_showerRadius(14.f),
    m_showerFitWindow(1000),
    m_edgeStep(2.f),
    m_moliereFraction(0.9f),
    m_maxNHitsLimit(2000.f),
    m_maxFoundHitRatioLimit(1.5f),
    m_maxScatterAngleLimit(40.f),
    m_maxOpeningAngleLimit(20.f),
    m_maxNuVertexEnergyWeightedMeanRadialDistanceLimit(20.f),
    m_minShowerStartMoliereRadiusLimit(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRegionFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
    const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch,
    const CartesianPointVector &showerStarts3D)
{
    float nHitsU(m_defaultFloat), foundHitRatioU(m_defaultRatio), scatterAngleU(m_defaultFloat), openingAngleU(m_defaultFloat),
        nuVertexEnergyAsymmetryU(m_defaultRatio), nuVertexEnergyWeightedMeanRadialDistanceU(m_defaultFloat),
        showerStartEnergyAsymmetryU(m_defaultRatio), showerStartMoliereRadiusU(m_defaultFloat);

    this->GetViewShowerRegionVariables(pAlgorithm, pShowerPfo, nuVertex3D, protoShowerMatch, TPC_VIEW_U, showerStarts3D.at(0), nHitsU,
        foundHitRatioU, scatterAngleU, openingAngleU, nuVertexEnergyAsymmetryU, nuVertexEnergyWeightedMeanRadialDistanceU,
        showerStartEnergyAsymmetryU, showerStartMoliereRadiusU);

    float nHitsV(m_defaultFloat), foundHitRatioV(m_defaultRatio), scatterAngleV(m_defaultFloat), openingAngleV(m_defaultFloat),
        nuVertexEnergyAsymmetryV(m_defaultRatio), nuVertexEnergyWeightedMeanRadialDistanceV(m_defaultFloat),
        showerStartEnergyAsymmetryV(m_defaultRatio), showerStartMoliereRadiusV(m_defaultFloat);

    this->GetViewShowerRegionVariables(pAlgorithm, pShowerPfo, nuVertex3D, protoShowerMatch, TPC_VIEW_V, showerStarts3D.at(0), nHitsV,
        foundHitRatioV, scatterAngleV, openingAngleV, nuVertexEnergyAsymmetryV, nuVertexEnergyWeightedMeanRadialDistanceV,
        showerStartEnergyAsymmetryV, showerStartMoliereRadiusV);

    float nHitsW(m_defaultFloat), foundHitRatioW(m_defaultRatio), scatterAngleW(m_defaultFloat), openingAngleW(m_defaultFloat),
        nuVertexEnergyAsymmetryW(m_defaultRatio), nuVertexEnergyWeightedMeanRadialDistanceW(m_defaultFloat),
        showerStartEnergyAsymmetryW(m_defaultRatio), showerStartMoliereRadiusW(m_defaultFloat);

    this->GetViewShowerRegionVariables(pAlgorithm, pShowerPfo, nuVertex3D, protoShowerMatch, TPC_VIEW_W, showerStarts3D.at(0), nHitsW,
        foundHitRatioW, scatterAngleW, openingAngleW, nuVertexEnergyAsymmetryW, nuVertexEnergyWeightedMeanRadialDistanceW,
        showerStartEnergyAsymmetryW, showerStartMoliereRadiusW);

    float minNHits(m_defaultFloat), minFoundHitRatio(m_defaultRatio), minScatterAngle(m_defaultFloat), minOpeningAngle(m_defaultFloat),
        minNuVertexEnergyAsymmetry(m_defaultRatio), minNuVertexEnergyWeightedMeanRadialDistance(m_defaultFloat),
        minShowerStartEnergyAsymmetry(m_defaultRatio), minShowerStartMoliereRadius(m_defaultFloat);

    float middleNHits(m_defaultFloat), middleFoundHitRatio(m_defaultRatio), middleScatterAngle(m_defaultFloat), middleOpeningAngle(m_defaultFloat),
        middleNuVertexEnergyAsymmetry(m_defaultRatio), middleNuVertexEnergyWeightedMeanRadialDistance(m_defaultFloat),
        middleShowerStartEnergyAsymmetry(m_defaultRatio), middleShowerStartMoliereRadius(m_defaultFloat);

    float maxNHits(m_defaultFloat), maxFoundHitRatio(m_defaultRatio), maxScatterAngle(m_defaultFloat), maxOpeningAngle(m_defaultFloat),
        maxNuVertexEnergyAsymmetry(m_defaultRatio), maxNuVertexEnergyWeightedMeanRadialDistance(m_defaultFloat),
        maxShowerStartEnergyAsymmetry(m_defaultRatio), maxShowerStartMoliereRadius(m_defaultFloat);

    LArConnectionPathwayHelper::GetMinMiddleMax(nHitsU, nHitsV, nHitsW, minNHits, middleNHits, maxNHits);
    LArConnectionPathwayHelper::GetMinMiddleMax(foundHitRatioU, foundHitRatioV, foundHitRatioW, minFoundHitRatio, middleFoundHitRatio, maxFoundHitRatio);
    LArConnectionPathwayHelper::GetMinMiddleMax(scatterAngleU, scatterAngleV, scatterAngleW, minScatterAngle, middleScatterAngle, maxScatterAngle);
    LArConnectionPathwayHelper::GetMinMiddleMax(openingAngleU, openingAngleV, openingAngleW, minOpeningAngle, middleOpeningAngle, maxOpeningAngle);
    LArConnectionPathwayHelper::GetMinMiddleMax(nuVertexEnergyAsymmetryU, nuVertexEnergyAsymmetryV, nuVertexEnergyAsymmetryW,
        minNuVertexEnergyAsymmetry, middleNuVertexEnergyAsymmetry, maxNuVertexEnergyAsymmetry);
    LArConnectionPathwayHelper::GetMinMiddleMax(nuVertexEnergyWeightedMeanRadialDistanceU, nuVertexEnergyWeightedMeanRadialDistanceV,
        nuVertexEnergyWeightedMeanRadialDistanceW, minNuVertexEnergyWeightedMeanRadialDistance,
        middleNuVertexEnergyWeightedMeanRadialDistance, maxNuVertexEnergyWeightedMeanRadialDistance);
    LArConnectionPathwayHelper::GetMinMiddleMax(showerStartEnergyAsymmetryU, showerStartEnergyAsymmetryV, showerStartEnergyAsymmetryW,
        minShowerStartEnergyAsymmetry, middleShowerStartEnergyAsymmetry, maxShowerStartEnergyAsymmetry);
    LArConnectionPathwayHelper::GetMinMiddleMax(showerStartMoliereRadiusU, showerStartMoliereRadiusV, showerStartMoliereRadiusW,
        minShowerStartMoliereRadius, middleShowerStartMoliereRadius, maxShowerStartMoliereRadius);

    featureVector.push_back(std::min(maxNHits, m_maxNHitsLimit));
    featureVector.push_back(std::min(maxFoundHitRatio, m_maxFoundHitRatioLimit));
    featureVector.push_back(std::min(maxScatterAngle, m_maxScatterAngleLimit));
    featureVector.push_back(std::min(maxOpeningAngle, m_maxOpeningAngleLimit));
    featureVector.push_back(maxNuVertexEnergyAsymmetry);
    featureVector.push_back(std::min(maxNuVertexEnergyWeightedMeanRadialDistance, m_maxNuVertexEnergyWeightedMeanRadialDistanceLimit));
    featureVector.push_back(maxShowerStartEnergyAsymmetry);
    featureVector.push_back(std::min(minShowerStartMoliereRadius, m_minShowerStartMoliereRadiusLimit));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRegionFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder, const std::string &featureToolName,
    const Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertex3D,
    const ProtoShowerMatch &protoShowerMatch, const CartesianPointVector &showerStarts3D)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pShowerPfo, nuVertex3D, protoShowerMatch, showerStarts3D);

    if (featureMap.find(featureToolName + "_nShowerHits") != featureMap.end())
    {
        std::cout << "Already wrote nShowerHits feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_nShowerHits");
    featureMap[featureToolName + "_nShowerHits"] = toolFeatureVec[0];

    if (featureMap.find(featureToolName + "_foundHitRatio") != featureMap.end())
    {
        std::cout << "Already wrote foundHitRatio feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_foundHitRatio");
    featureMap[featureToolName + "_foundHitRatio"] = toolFeatureVec[1];

    if (featureMap.find(featureToolName + "_scatterAngle") != featureMap.end())
    {
        std::cout << "Already wrote scatterAngle feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_scatterAngle");
    featureMap[featureToolName + "_scatterAngle"] = toolFeatureVec[2];

    if (featureMap.find(featureToolName + "_openingAngle") != featureMap.end())
    {
        std::cout << "Already wrote openingAngle feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_openingAngle");
    featureMap[featureToolName + "_openingAngle"] = toolFeatureVec[3];

    if (featureMap.find(featureToolName + "_nuVertexEnergyAsymmetry") != featureMap.end())
    {
        std::cout << "Already wrote nuVertexEnergyAsymmetry feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_nuVertexEnergyAsymmetry");
    featureMap[featureToolName + "_nuVertexEnergyAsymmetry"] = toolFeatureVec[4];

    if (featureMap.find(featureToolName + "_nuVertexEnergyWeightedMeanRadialDistance") != featureMap.end())
    {
        std::cout << "Already wrote nuVertexEnergyWeightedMeanRadialDistance feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_nuVertexEnergyWeightedMeanRadialDistance");
    featureMap[featureToolName + "_nuVertexEnergyWeightedMeanRadialDistance"] = toolFeatureVec[5];

    if (featureMap.find(featureToolName + "_showerStartEnergyAsymmetry") != featureMap.end())
    {
        std::cout << "Already wrote showerStartEnergyAsymmetry feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_showerStartEnergyAsymmetry");
    featureMap[featureToolName + "_showerStartEnergyAsymmetry"] = toolFeatureVec[6];

    if (featureMap.find(featureToolName + "_showerStartMoliereRadius") != featureMap.end())
    {
        std::cout << "Already wrote showerStartMoliereRadius feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_showerStartMoliereRadius");
    featureMap[featureToolName + "_showerStartMoliereRadius"] = toolFeatureVec[7];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRegionFeatureTool::GetViewShowerRegionVariables(const Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch, const HitType hitType, const CartesianVector &showerStart3D,
    float &nHits, float &foundHitRatio, float &scatterAngle, float &openingAngle, float &nuVertexEnergyAsymmetry,
    float &nuVertexEnergyWeightedMeanRadialDistance, float &showerStartEnergyAsymmetry, float &showerStartMoliereRadius)
{
    CaloHitList viewHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewHitList);

    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertex3D, hitType));
    const CartesianVector showerStart2D(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), showerStart3D, hitType));
    const bool isDownstream(showerStart2D.GetZ() > nuVertex2D.GetZ());

    try
    {
        // Fit the spine to get the initial shower direction
        const CaloHitList &spineHitList(hitType == TPC_VIEW_U ? protoShowerMatch.GetProtoShowerU().GetSpineHitList()
                                                              : (hitType == TPC_VIEW_V ? protoShowerMatch.GetProtoShowerV().GetSpineHitList()
                                                                                       : protoShowerMatch.GetProtoShowerW().GetSpineHitList()));

        CartesianPointVector spinePositions;
        for (const CaloHit *const pCaloHit : spineHitList)
            spinePositions.push_back(pCaloHit->GetPositionVector());

        const TwoDSlidingFitResult spineFitResult(&spinePositions, m_spineFitWindow, LArGeometryHelper::GetWirePitch(pAlgorithm->GetPandora(), hitType));

        // Now can build the shower
        CaloHitList postShowerHitList;
        CartesianPointVector postShowerPositions;

        this->BuildViewShower(pShowerPfo, spineFitResult, hitType, showerStart2D, nuVertex2D, postShowerHitList, postShowerPositions);

        this->GetShowerHitVariables(spineHitList, postShowerHitList, pShowerPfo, hitType, nHits, foundHitRatio);

        // Fit the shower
        TwoDSlidingFitResult showerFitResult(
            &postShowerPositions, m_showerFitWindow, LArGeometryHelper::GetWirePitch(pAlgorithm->GetPandora(), hitType));

        const bool isShowerDownstream((showerStart2D - showerFitResult.GetGlobalMinLayerPosition()).GetMagnitude() <
                                      (showerStart2D - showerFitResult.GetGlobalMaxLayerPosition()).GetMagnitude());

        // Collect more hits
        for (const CaloHit *const pCaloHit : viewHitList)
        {
            if (std::find(postShowerHitList.begin(), postShowerHitList.end(), pCaloHit) != postShowerHitList.end())
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const CartesianVector displacement(hitPosition - showerStart2D);
            const float l(showerFitResult.GetGlobalMinLayerDirection().GetDotProduct(displacement));
            const float t(showerFitResult.GetGlobalMinLayerDirection().GetCrossProduct(displacement).GetMagnitude());

            if (((isDownstream && (l > 0.f)) || (!isDownstream && (l < 0.f))) && (t < m_showerRadius))
            {
                postShowerHitList.push_back(pCaloHit);
                postShowerPositions.push_back(pCaloHit->GetPositionVector());
            }
        }

        this->GetShowerHitVariables(spineHitList, postShowerHitList, pShowerPfo, hitType, nHits, foundHitRatio);

        this->CalculateViewScatterAngle(nuVertex2D, spineFitResult, showerStart2D, showerFitResult, scatterAngle);

        this->CalculateViewOpeningAngle(showerFitResult, postShowerHitList, showerStart2D, openingAngle);

        this->CalculateViewNuVertexConsistencyVariables(
            spineFitResult, postShowerHitList, isDownstream, nuVertex2D, nuVertexEnergyAsymmetry, nuVertexEnergyWeightedMeanRadialDistance);

        this->CalculateViewShowerStartConsistencyVariables(
            showerFitResult, postShowerHitList, isShowerDownstream, showerStartEnergyAsymmetry, showerStartMoliereRadius);
    }
    catch (...)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRegionFeatureTool::BuildViewShower(const ParticleFlowObject *const pShowerPfo, const TwoDSlidingFitResult &spineFit, const HitType hitType,
    const CartesianVector &showerStart2D, const CartesianVector &nuVertex2D, CaloHitList &postShowerHitList, CartesianPointVector &postShowerPositions)
{
    // Find initial shower direction from spine fit
    float lShowerStart(0.f), tShowerStart(0.f);
    spineFit.GetLocalPosition(showerStart2D, lShowerStart, tShowerStart);

    const LayerFitResultMap &layerFitResultMap(spineFit.GetLayerFitResultMap());
    const int showerStartLayer(spineFit.GetLayer(lShowerStart));
    int closestLayer(std::numeric_limits<int>::max());

    for (const auto &entry : layerFitResultMap)
    {
        if (std::fabs(entry.first - showerStartLayer) < std::fabs(entry.first - closestLayer))
            closestLayer = entry.first;
    }

    const float gradient(layerFitResultMap.at(closestLayer).GetGradient());
    CartesianVector showerDirection(0.f, 0.f, 0.f);

    spineFit.GetGlobalDirection(gradient, showerDirection);

    // Now find the shower
    const bool isDownstream(showerStart2D.GetZ() > nuVertex2D.GetZ());
    CaloHitList caloHitList;

    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, caloHitList);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector displacement(hitPosition - showerStart2D);
        const float l(showerDirection.GetDotProduct(displacement));
        const float t(showerDirection.GetCrossProduct(displacement).GetMagnitude());

        if (((isDownstream && (l > 0.f)) || (!isDownstream && (l < 0.f))) && (t < m_showerRadius))
        {
            postShowerHitList.push_back(pCaloHit);
            postShowerPositions.push_back(pCaloHit->GetPositionVector());
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRegionFeatureTool::GetShowerHitVariables(const CaloHitList &spineHitList, const CaloHitList &postShowerHitList,
    const ParticleFlowObject *const pShowerPfo, const HitType hitType, float &nHits, float &foundHitRatio)
{
    int foundHits(spineHitList.size());

    for (const CaloHit *const pCaloHit : postShowerHitList)
    {
        if (std::find(spineHitList.begin(), spineHitList.end(), pCaloHit) == spineHitList.end())
            ++foundHits;
    }

    CaloHitList viewHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewHitList);

    foundHitRatio = static_cast<float>(foundHits) / static_cast<float>(viewHitList.size());
    nHits = postShowerHitList.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRegionFeatureTool::CalculateViewScatterAngle(const CartesianVector &nuVertex2D, const TwoDSlidingFitResult &spineFitResult,
    const CartesianVector &showerStart2D, const TwoDSlidingFitResult &showerFitResult, float &scatterAngle)
{
    const bool isDownstream(showerStart2D.GetZ() > nuVertex2D.GetZ());
    const CartesianVector streamCorrectedConnectionPathwayDirection(
        isDownstream ? spineFitResult.GetGlobalMinLayerDirection() : spineFitResult.GetGlobalMaxLayerDirection() * (-1.f));

    const bool isShowerDownstream((showerStart2D - showerFitResult.GetGlobalMinLayerPosition()).GetMagnitude() <
                                  (showerStart2D - showerFitResult.GetGlobalMaxLayerPosition()).GetMagnitude());
    const CartesianVector streamCorrectedShowerDirection(
        isShowerDownstream ? showerFitResult.GetGlobalMinLayerDirection() : showerFitResult.GetGlobalMaxLayerDirection() * (-1.f));

    scatterAngle = streamCorrectedConnectionPathwayDirection.GetOpeningAngle(streamCorrectedShowerDirection) * 180.f / M_PI;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRegionFeatureTool::CalculateViewOpeningAngle(const TwoDSlidingFitResult &showerFitResult, const CaloHitList &postShowerHitList,
    const CartesianVector &showerStart2D, float &openingAngle)
{
    try
    {
        const CartesianVector directionAxis = showerFitResult.GetGlobalMinLayerDirection();
        const CartesianVector orthoAxis = directionAxis.GetCrossProduct(CartesianVector(0.f, 1.f, 0.f));

        std::map<int, float> positiveEdges, negativeEdges;

        for (const CaloHit *const pCaloHit : postShowerHitList)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - showerStart2D);
            const float thisT(directionAxis.GetCrossProduct(position).GetMagnitude());
            const float thisL(directionAxis.GetDotProduct(position));
            const float orthoL(orthoAxis.GetDotProduct(position));

            std::map<int, float> &edgeMap(orthoL > 0.f ? positiveEdges : negativeEdges);

            const int lIndex(std::floor(thisL / m_edgeStep));

            edgeMap[lIndex] = (edgeMap.find(lIndex) == edgeMap.end() ? thisT : std::max(edgeMap[lIndex], thisT));
        }

        CartesianPointVector positiveEdgePositions, negativeEdgePositions;

        for (auto &entry : positiveEdges)
            positiveEdgePositions.push_back(CartesianVector(entry.second, 0.f, entry.first));

        for (auto &entry : negativeEdges)
            negativeEdgePositions.push_back(CartesianVector(entry.second, 0.f, entry.first));

        const TwoDSlidingFitResult positiveEdgeFit(&positiveEdgePositions, m_showerFitWindow, showerFitResult.GetLayerPitch());
        const TwoDSlidingFitResult negativeEdgeFit(&negativeEdgePositions, m_showerFitWindow, showerFitResult.GetLayerPitch());

        const CartesianVector positiveMinLayer(positiveEdgeFit.GetGlobalMinLayerPosition());
        const CartesianVector positiveMaxLayer(positiveEdgeFit.GetGlobalMaxLayerPosition());
        const CartesianVector negativeMinLayer(negativeEdgeFit.GetGlobalMinLayerPosition());
        const CartesianVector negativeMaxLayer(negativeEdgeFit.GetGlobalMaxLayerPosition());

        const CartesianVector globalPositiveMinLayer(
            showerStart2D + (directionAxis * positiveMinLayer.GetZ()) + (orthoAxis * positiveMinLayer.GetX()));
        const CartesianVector globalPositiveMaxLayer(
            showerStart2D + (directionAxis * positiveMaxLayer.GetZ()) + (orthoAxis * positiveMaxLayer.GetX()));
        const CartesianVector globalNegativeMinLayer(
            showerStart2D + (directionAxis * negativeMinLayer.GetZ()) - (orthoAxis * negativeMinLayer.GetX()));
        const CartesianVector globalNegativeMaxLayer(
            showerStart2D + (directionAxis * negativeMaxLayer.GetZ()) - (orthoAxis * negativeMaxLayer.GetX()));

        const CartesianVector positiveEdgeVector((globalPositiveMaxLayer - globalPositiveMinLayer).GetUnitVector());
        const CartesianVector negativeEdgeVector((globalNegativeMaxLayer - globalNegativeMinLayer).GetUnitVector());

        const float positiveOpeningAngle = directionAxis.GetOpeningAngle(positiveEdgeVector) * 180.f / M_PI;
        const float negativeOpeningAngle = directionAxis.GetOpeningAngle(negativeEdgeVector) * 180.f / M_PI;

        openingAngle = std::max(positiveOpeningAngle, negativeOpeningAngle);
    }
    catch (...)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRegionFeatureTool::CalculateViewNuVertexConsistencyVariables(const TwoDSlidingFitResult &spineFitResult, const CaloHitList &postShowerHitList,
    const bool isDownstream, const CartesianVector &nuVertex2D, float &nuVertexEnergyAsymmetry, float &nuVertexEnergyWeightedMeanRadialDistance)
{
    const CartesianVector &directionAxis(isDownstream ? spineFitResult.GetGlobalMinLayerDirection() : spineFitResult.GetGlobalMaxLayerDirection());
    const CartesianVector orthoAxis = directionAxis.GetCrossProduct(CartesianVector(0.f, 1.f, 0.f));

    float totalEnergy(0.f);

    // nuVertexEnergyAsymmetry
    nuVertexEnergyAsymmetry = 0.f;

    for (const CaloHit *const pCaloHit : postShowerHitList)
    {
        const float hitEnergy(std::fabs(pCaloHit->GetElectromagneticEnergy()));

        totalEnergy += hitEnergy;

        const CartesianVector position(pCaloHit->GetPositionVector() - nuVertex2D);
        const float thisL(orthoAxis.GetDotProduct(position));

        nuVertexEnergyAsymmetry += (thisL < 0.f) ? (-1.f * hitEnergy) : hitEnergy;
    }

    nuVertexEnergyAsymmetry = (totalEnergy < std::numeric_limits<float>::epsilon()) ? m_defaultRatio : (nuVertexEnergyAsymmetry / totalEnergy);
    nuVertexEnergyAsymmetry = std::fabs(nuVertexEnergyAsymmetry);

    // nuVertexEnergyWeightedMeanRadialDistance
    nuVertexEnergyWeightedMeanRadialDistance = 0.f;

    for (const CaloHit *const pCaloHit : postShowerHitList)
    {
        const CartesianVector position(pCaloHit->GetPositionVector() - nuVertex2D);
        const float hitEnergy(std::fabs(pCaloHit->GetElectromagneticEnergy()));

        nuVertexEnergyWeightedMeanRadialDistance += (directionAxis.GetCrossProduct(position).GetMagnitude() * hitEnergy);
    }

    nuVertexEnergyWeightedMeanRadialDistance =
        (totalEnergy < std::numeric_limits<float>::epsilon()) ? m_defaultFloat : nuVertexEnergyWeightedMeanRadialDistance / totalEnergy;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerRegionFeatureTool::CalculateViewShowerStartConsistencyVariables(const TwoDSlidingFitResult &showerFitResult,
    const CaloHitList &postShowerHitList, const bool isShowerDownstream, float &showerStartEnergyAsymmetry, float &showerStartMoliereRadius)
{
    const CartesianVector fitShowerStart(isShowerDownstream ? showerFitResult.GetGlobalMinLayerPosition() : showerFitResult.GetGlobalMaxLayerPosition());
    const CartesianVector directionAxis(
        isShowerDownstream ? showerFitResult.GetGlobalMinLayerDirection() : showerFitResult.GetGlobalMaxLayerDirection());
    const CartesianVector orthoAxis = directionAxis.GetCrossProduct(CartesianVector(0.f, 1.f, 0.f));

    float totalEnergy(0.f);

    // showerStartEnergyAsymmetry
    showerStartEnergyAsymmetry = 0.f;

    for (const CaloHit *const pCaloHit : postShowerHitList)
    {
        const float hitEnergy(std::fabs(pCaloHit->GetElectromagneticEnergy()));

        totalEnergy += hitEnergy;

        const CartesianVector position(pCaloHit->GetPositionVector() - fitShowerStart);
        const float thisL(orthoAxis.GetDotProduct(position));

        showerStartEnergyAsymmetry += (thisL < 0.f) ? (-1.f * hitEnergy) : hitEnergy;
    }

    showerStartEnergyAsymmetry = (totalEnergy < std::numeric_limits<float>::epsilon()) ? m_defaultRatio : (showerStartEnergyAsymmetry / totalEnergy);
    showerStartEnergyAsymmetry = std::fabs(showerStartEnergyAsymmetry);

    // showerStartMoliereRadius
    showerStartMoliereRadius = m_defaultFloat;

    CaloHitVector showerStartPostShowerHitVector(postShowerHitList.begin(), postShowerHitList.end());

    std::sort(showerStartPostShowerHitVector.begin(), showerStartPostShowerHitVector.end(),
        [&fitShowerStart, &directionAxis](const CaloHit *const pCaloHitA, const CaloHit *const pCaloHitB) -> bool {
            const CartesianVector positionA(pCaloHitA->GetPositionVector() - fitShowerStart);
            const CartesianVector positionB(pCaloHitB->GetPositionVector() - fitShowerStart);

            const float tA(directionAxis.GetCrossProduct(positionA).GetMagnitude());
            const float tB(directionAxis.GetCrossProduct(positionB).GetMagnitude());

            return tA < tB;
        });

    float showerStartRunningEnergySum(0.f);

    for (const CaloHit *const pCaloHit : showerStartPostShowerHitVector)
    {
        const float hitEnergy(std::fabs(pCaloHit->GetElectromagneticEnergy()));
        showerStartRunningEnergySum += hitEnergy;

        if ((showerStartRunningEnergySum / totalEnergy) > m_moliereFraction)
        {
            const CartesianVector position(pCaloHit->GetPositionVector() - fitShowerStart);
            showerStartMoliereRadius = directionAxis.GetCrossProduct(position).GetMagnitude();
            break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerRegionFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DefaultFloat", m_defaultFloat));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DefaultRatio", m_defaultRatio));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SpineFitWindow", m_spineFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowerRadius", m_showerRadius));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowerFitWindow", m_showerFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EdgeStep", m_edgeStep));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MoliereFraction", m_moliereFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxNHitsLimit", m_maxNHitsLimit));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxFoundHitRatioLimit", m_maxFoundHitRatioLimit));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxScatterAngleLimit", m_maxScatterAngleLimit));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxOpeningAngleLimit", m_maxOpeningAngleLimit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxNuVertexEnergyWeightedMeanRadialDistanceLimit", m_maxNuVertexEnergyWeightedMeanRadialDistanceLimit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinShowerStartMoliereRadiusLimit", m_minShowerStartMoliereRadiusLimit));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

AmbiguousRegionFeatureTool::AmbiguousRegionFeatureTool() :
    m_defaultFloat(-10.f),
    m_caloHitListNameU("CaloHitListU"),
    m_caloHitListNameV("CaloHitListV"),
    m_caloHitListNameW("CaloHitListW"),
    m_maxTransverseDistance(0.75f),
    m_maxSampleHits(3),
    m_maxHitSeparation(1.f),
    m_maxTrackFraction(0.8f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AmbiguousRegionFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm,
    const ParticleFlowObject *const /*pShowerPfo*/, const CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch,
    const CartesianPointVector & /*showerStarts3D*/)
{
    float nAmbiguousViews(0.f);
    this->CalculateNAmbiguousViews(protoShowerMatch, nAmbiguousViews);

    float maxUnaccountedEnergy(m_defaultFloat);
    float unaccountedHitEnergyU(m_defaultFloat), unaccountedHitEnergyV(m_defaultFloat), unaccountedHitEnergyW(m_defaultFloat);

    if (this->GetViewAmbiguousHitVariables(pAlgorithm, protoShowerMatch, TPC_VIEW_U, nuVertex3D, unaccountedHitEnergyU))
        maxUnaccountedEnergy = std::max(maxUnaccountedEnergy, unaccountedHitEnergyU);

    if (this->GetViewAmbiguousHitVariables(pAlgorithm, protoShowerMatch, TPC_VIEW_V, nuVertex3D, unaccountedHitEnergyV))
        maxUnaccountedEnergy = std::max(maxUnaccountedEnergy, unaccountedHitEnergyV);

    if (this->GetViewAmbiguousHitVariables(pAlgorithm, protoShowerMatch, TPC_VIEW_W, nuVertex3D, unaccountedHitEnergyW))
        maxUnaccountedEnergy = std::max(maxUnaccountedEnergy, unaccountedHitEnergyW);

    featureVector.push_back(nAmbiguousViews);
    featureVector.push_back(maxUnaccountedEnergy);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AmbiguousRegionFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder,
    const std::string &featureToolName, const Algorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertex3D, const ProtoShowerMatch &protoShowerMatch, const CartesianPointVector &showerStarts3D)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pShowerPfo, nuVertex3D, protoShowerMatch, showerStarts3D);

    if (featureMap.find(featureToolName + "_nAmbiguousViews") != featureMap.end())
    {
        std::cout << "Already wrote nAmbiguousViews feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_nAmbiguousViews");
    featureMap[featureToolName + "_nAmbiguousViews"] = toolFeatureVec[0];

    if (featureMap.find(featureToolName + "_maxUnaccountedEnergy") != featureMap.end())
    {
        std::cout << "Already wrote maxUnaccountedEnergy feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_maxUnaccountedEnergy");
    featureMap[featureToolName + "_maxUnaccountedEnergy"] = toolFeatureVec[1];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AmbiguousRegionFeatureTool::CalculateNAmbiguousViews(const ProtoShowerMatch &protoShowerMatch, float &nAmbiguousViews)
{
    nAmbiguousViews = 0.f;

    const int nAmbiguousHitsU(protoShowerMatch.GetProtoShowerU().GetAmbiguousHitList().size());
    nAmbiguousViews += (nAmbiguousHitsU == 0) ? 0.f : 1.f;

    const int nAmbiguousHitsV(protoShowerMatch.GetProtoShowerV().GetAmbiguousHitList().size());
    nAmbiguousViews += (nAmbiguousHitsV == 0) ? 0.f : 1.f;

    const int nAmbiguousHitsW(protoShowerMatch.GetProtoShowerW().GetAmbiguousHitList().size());
    nAmbiguousViews += (nAmbiguousHitsW == 0) ? 0.f : 1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool AmbiguousRegionFeatureTool::GetViewAmbiguousHitVariables(const Algorithm *const pAlgorithm, const ProtoShowerMatch &protoShowerMatch,
    const HitType hitType, const CartesianVector &nuVertex3D, float &unaccountedHitEnergy)
{
    std::map<int, CaloHitList> ambiguousHitSpines;
    CaloHitList hitsToExcludeInEnergyCalcs; // to avoid double  counting
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertex3D, hitType));
    const ProtoShower &protoShower(hitType == TPC_VIEW_U
                                       ? protoShowerMatch.GetProtoShowerU()
                                       : (hitType == TPC_VIEW_V ? protoShowerMatch.GetProtoShowerV() : protoShowerMatch.GetProtoShowerW()));

    this->BuildAmbiguousSpines(pAlgorithm, hitType, protoShower, nuVertex2D, ambiguousHitSpines, hitsToExcludeInEnergyCalcs);

    if (ambiguousHitSpines.empty())
        return false;

    const CartesianVector &startDirection(protoShower.GetConnectionPathway().GetStartDirection());
    float startL(-std::numeric_limits<float>::max());
    float ambiguousHitEnergyMean(0.f);

    // Get average energy of the ambiguous hits
    for (const CaloHit *const pAmbiguousCaloHit : protoShower.GetAmbiguousHitList())
    {
        const float thisT(startDirection.GetCrossProduct(pAmbiguousCaloHit->GetPositionVector() - nuVertex2D).GetMagnitude());
        const float thisL(startDirection.GetDotProduct(pAmbiguousCaloHit->GetPositionVector() - nuVertex2D));

        if ((thisL > startL) && (thisT < m_maxTransverseDistance))
            startL = thisL;

        ambiguousHitEnergyMean += pAmbiguousCaloHit->GetElectromagneticEnergy() * 1000.f;
    }

    ambiguousHitEnergyMean /= protoShower.GetAmbiguousHitList().size();

    // Get mean energy of other pathways, avoiding the float counting hits
    float otherEnergyMeanSum(0.f);

    for (const auto &entry : ambiguousHitSpines)
    {
        int nOtherEnergyHits(0);
        float otherEnergyMean(0.f);

        for (const CaloHit *const pOtherCaloHit : entry.second)
        {
            if (std::find(hitsToExcludeInEnergyCalcs.begin(), hitsToExcludeInEnergyCalcs.end(), pOtherCaloHit) != hitsToExcludeInEnergyCalcs.end())
                continue;

            otherEnergyMean += pOtherCaloHit->GetElectromagneticEnergy() * 1000.f;
            ++nOtherEnergyHits;
        }

        if (nOtherEnergyHits == 0)
            continue;

        otherEnergyMean /= static_cast<float>(nOtherEnergyHits);
        otherEnergyMeanSum += otherEnergyMean;
    }

    // Get the spine mean energy, only consider a limited number of  non-ambiguous hits
    float spineEnergyMean(0.f);
    unsigned int nSpineEnergyHits(0);

    for (const CaloHit *const pSpineCaloHit : protoShower.GetSpineHitList())
    {
        if (std::find(protoShower.GetAmbiguousHitList().begin(), protoShower.GetAmbiguousHitList().end(), pSpineCaloHit) !=
            protoShower.GetAmbiguousHitList().end())
            continue;

        const float thisL(startDirection.GetDotProduct(pSpineCaloHit->GetPositionVector() - nuVertex2D));

        if ((thisL > startL) && (nSpineEnergyHits < m_maxSampleHits))
        {
            spineEnergyMean += pSpineCaloHit->GetElectromagneticEnergy() * 1000.f;
            ++nSpineEnergyHits;
        }
    }

    if (nSpineEnergyHits == 0)
        return false;

    spineEnergyMean /= static_cast<float>(nSpineEnergyHits);
    unaccountedHitEnergy = ambiguousHitEnergyMean - otherEnergyMeanSum - spineEnergyMean;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void AmbiguousRegionFeatureTool::BuildAmbiguousSpines(const Algorithm *const pAlgorithm, const HitType hitType, const ProtoShower &protoShower,
    const CartesianVector &nuVertex2D, std::map<int, CaloHitList> &ambiguousHitSpines, CaloHitList &hitsToExcludeInEnergyCalcs)
{
    const CaloHitList *pCaloHitList;

    if (this->GetHitListOfType(pAlgorithm, hitType, pCaloHitList) != STATUS_CODE_SUCCESS)
        return;

    std::map<int, CaloHitList> ambiguousHitSpinesTemp;

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (std::find(protoShower.GetAmbiguousHitList().begin(), protoShower.GetAmbiguousHitList().end(), pCaloHit) !=
            protoShower.GetAmbiguousHitList().end())
            continue;

        int count(0);

        // A hit can be in more than one spine
        for (unsigned int i = 0; i < protoShower.GetAmbiguousDirectionVector().size(); ++i)
        {
            const CartesianVector &significantDirection(protoShower.GetAmbiguousDirectionVector()[i]);
            const CartesianVector displacement(pCaloHit->GetPositionVector() - nuVertex2D);
            const float thisT(significantDirection.GetCrossProduct(displacement).GetMagnitude());
            const float thisL(significantDirection.GetDotProduct(displacement));

            if ((thisL > 0.f) && (thisT < m_maxTransverseDistance))
            {
                ++count;
                ambiguousHitSpinesTemp[i].push_back(pCaloHit);
            }

            if (count == 2)
                hitsToExcludeInEnergyCalcs.push_back(pCaloHit);
        }
    }

    if (ambiguousHitSpinesTemp.empty())
        return;

    // Make sure the pathways are continuous
    for (const auto &entry : ambiguousHitSpinesTemp)
    {
        CaloHitList continuousSpine(this->FindAmbiguousContinuousSpine(entry.second, protoShower.GetAmbiguousHitList(), nuVertex2D));

        if (continuousSpine.size() > 0)
            ambiguousHitSpines[entry.first] = continuousSpine;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AmbiguousRegionFeatureTool::GetHitListOfType(const Algorithm *const pAlgorithm, const HitType hitType, const CaloHitList *&pCaloHitList) const
{
    const std::string typeHitListName(hitType == TPC_VIEW_U ? m_caloHitListNameU : hitType == TPC_VIEW_V ? m_caloHitListNameV : m_caloHitListNameW);

    PANDORA_THROW_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*pAlgorithm, typeHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
        return STATUS_CODE_NOT_INITIALIZED;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList AmbiguousRegionFeatureTool::FindAmbiguousContinuousSpine(
    const CaloHitList &caloHitList, const CaloHitList &ambiguousHitList, const CartesianVector &nuVertex2D)
{
    CaloHitList continuousHitList;

    CaloHitVector caloHitVector(caloHitList.begin(), caloHitList.end());
    std::sort(caloHitVector.begin(), caloHitVector.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(nuVertex2D));

    for (unsigned int i = 0; i < caloHitVector.size(); ++i)
    {
        CaloHitList connectedHitList;
        connectedHitList.push_back(caloHitVector[i]);

        if (LArClusterHelper::GetClosestDistance(connectedHitList.front()->GetPositionVector(), ambiguousHitList) > m_maxHitSeparation)
            continue;

        bool found(true);

        while (found)
        {
            found = false;

            for (unsigned int j = (i + 1); j < caloHitVector.size(); ++j)
            {
                const CaloHit *const pCaloHit(caloHitVector[j]);

                if (std::find(connectedHitList.begin(), connectedHitList.end(), pCaloHit) != connectedHitList.end())
                    continue;

                if (LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), connectedHitList) < m_maxHitSeparation)
                {
                    // to avoid ends of tracks
                    if (static_cast<float>(connectedHitList.size()) < static_cast<float>(caloHitVector.size() * m_maxTrackFraction))
                    {
                        found = true;
                        connectedHitList.push_back(pCaloHit);
                    }

                    break;
                }
            }
        }

        if (connectedHitList.size() >= 2)
        {
            continuousHitList.insert(continuousHitList.begin(), connectedHitList.begin(), connectedHitList.end());
            break;
        }
    }

    return continuousHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode AmbiguousRegionFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DefaultFloat", m_defaultFloat));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameU", m_caloHitListNameU));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameV", m_caloHitListNameV));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameW", m_caloHitListNameW));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxTransverseDistance", m_maxTransverseDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSampleHits", m_maxSampleHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxHitSeparation", m_maxHitSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxTrackFraction", m_maxTrackFraction));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
