/**
 *  @file   larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.cc
 *
 *  @brief  Implementation of the track shower id feature fool class
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArTrackShowerId/CutClusterCharacterisationAlgorithm.h"
#include "larpandoracontent/LArTrackShowerId/TrackShowerIdFeatureTool.h"

using namespace pandora;

namespace lar_content
{

PositionInCryostat LocatePointInCryostat_ICARUS(const float point_X)
{       
    PositionInCryostat position;
    const float CathodeMinX = 210.14; ///< Cathode position for ICARUS cryostat #1
    const float CathodeMaxX = 210.29;

    float RealCathodeMinX = 0., RealCathodeMaxX = 0.;
    if (point_X < 0.) ///< ICARUS cryostat #0 
    {                                                                                                                                                                         
        RealCathodeMinX = CathodeMaxX*(-1.);
        RealCathodeMaxX = CathodeMinX*(-1.);
    }
    else ///< ICARUS cryostat #1
    {                                                                                                                                                                                       
        RealCathodeMinX = CathodeMinX;
        RealCathodeMaxX = CathodeMaxX;
    }

    if (point_X < RealCathodeMinX)
    {
        position = PositionInCryostat::BelowCathode;
    }
    else if (point_X > RealCathodeMaxX)
    {
        position = PositionInCryostat::AboveCathode;
    }
    else
    {
        position = PositionInCryostat::WithinCathode;
    }

    return position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoDShowerFitFeatureTool::TwoDShowerFitFeatureTool() : m_slidingShowerFitWindow(3), m_slidingLinearFitWindow(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDShowerFitFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float ratio(-1.f);
    try
    {
        const TwoDSlidingFitResult slidingFitResultLarge(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const float straightLineLength =
            (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
        if (straightLineLength > std::numeric_limits<float>::epsilon())
            ratio = (CutClusterCharacterisationAlgorithm::GetShowerFitWidth(pAlgorithm, pCluster, m_slidingShowerFitWindow)) / straightLineLength;
    }
    catch (const StatusCodeException &)
    {
        ratio = -1.f;
    }
    featureVector.push_back(ratio);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDShowerFitFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder, const std::string &featureToolName,
    const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pCluster);

    if (featureMap.find(featureToolName + "_WidthLenRatio") != featureMap.end())
    {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_WidthLenRatio");
    featureMap[featureToolName + "_WidthLenRatio"] = toolFeatureVec[0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDShowerFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingShowerFitWindow", m_slidingShowerFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TwoDLinearFitFeatureTool::TwoDLinearFitFeatureTool() :
    m_slidingLinearFitWindow(3),
    m_slidingLinearFitWindowLarge(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDLinearFitFeatureTool::Run(LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float dTdLWidth(-1.f), straightLineLengthLarge(-1.f), diffWithStraightLineMean(-1.f), diffWithStraightLineSigma(-1.f),
        maxFitGapLength(-1.f), rmsSlidingLinearFit(-1.f);
    this->CalculateVariablesSlidingLinearFit(pCluster, straightLineLengthLarge, diffWithStraightLineMean, diffWithStraightLineSigma,
        dTdLWidth, maxFitGapLength, rmsSlidingLinearFit);

    if (straightLineLengthLarge > std::numeric_limits<float>::epsilon())
    {
        diffWithStraightLineMean /= straightLineLengthLarge;
        diffWithStraightLineSigma /= straightLineLengthLarge;
        dTdLWidth /= straightLineLengthLarge;
        maxFitGapLength /= straightLineLengthLarge;
        rmsSlidingLinearFit /= straightLineLengthLarge;
    }

    featureVector.push_back(straightLineLengthLarge);
    featureVector.push_back(diffWithStraightLineMean);
    featureVector.push_back(diffWithStraightLineSigma);
    featureVector.push_back(dTdLWidth);
    featureVector.push_back(maxFitGapLength);
    featureVector.push_back(rmsSlidingLinearFit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDLinearFitFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder, const std::string &featureToolName,
    const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pCluster);

    if (featureMap.find(featureToolName + "_StLineLenLarge") != featureMap.end() ||
        featureMap.find(featureToolName + "_DiffStLineMean") != featureMap.end() ||
        featureMap.find(featureToolName + "_DiffStLineSigma") != featureMap.end() ||
        featureMap.find(featureToolName + "_dTdLWidth") != featureMap.end() ||
        featureMap.find(featureToolName + "_MaxFitGapLen") != featureMap.end() ||
        featureMap.find(featureToolName + "_rmsSlidingLinFit") != featureMap.end())
    {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_StLineLenLarge");
    featureOrder.push_back(featureToolName + "_DiffStLineMean");
    featureOrder.push_back(featureToolName + "_DiffStLineSigma");
    featureOrder.push_back(featureToolName + "_dTdLWidth");
    featureOrder.push_back(featureToolName + "_MaxFitGapLen");
    featureOrder.push_back(featureToolName + "_rmsSlidingLinFit");

    featureMap[featureToolName + "_StLineLenLarge"] = toolFeatureVec[0];
    featureMap[featureToolName + "_DiffStLineMean"] = toolFeatureVec[1];
    featureMap[featureToolName + "_DiffStLineSigma"] = toolFeatureVec[2];
    featureMap[featureToolName + "_dTdLWidth"] = toolFeatureVec[3];
    featureMap[featureToolName + "_MaxFitGapLen"] = toolFeatureVec[4];
    featureMap[featureToolName + "_rmsSlidingLinFit"] = toolFeatureVec[5];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDLinearFitFeatureTool::CalculateVariablesSlidingLinearFit(const pandora::Cluster *const pCluster, float &straightLineLengthLarge,
    float &diffWithStraightLineMean, float &diffWithStraightLineSigma, float &dTdLWidth, float &maxFitGapLength, float &rmsSlidingLinearFit) const
{
    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult slidingFitResultLarge(
            pCluster, m_slidingLinearFitWindowLarge, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        if (slidingFitResult.GetLayerFitResultMap().empty())
            throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

        straightLineLengthLarge =
            (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
        rmsSlidingLinearFit = 0.f;

        FloatVector diffWithStraightLineVector;
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        CartesianVector previousFitPosition(slidingFitResult.GetGlobalMinLayerPosition());
        float dTdLMin(+std::numeric_limits<float>::max()), dTdLMax(-std::numeric_limits<float>::max());

        for (const auto &mapEntry : slidingFitResult.GetLayerFitResultMap())
        {
            const LayerFitResult &layerFitResult(mapEntry.second);
            rmsSlidingLinearFit += layerFitResult.GetRms();

            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(layerFitResult.GetL(), layerFitResult.GetFitT(), thisFitPosition);

            LayerFitResultMap::const_iterator iterLarge =
                slidingFitResultLarge.GetLayerFitResultMap().find(slidingFitResultLarge.GetLayer(layerFitResult.GetL()));

            if (slidingFitResultLarge.GetLayerFitResultMap().end() == iterLarge)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            diffWithStraightLineVector.push_back(static_cast<float>(std::fabs(layerFitResult.GetFitT() - iterLarge->second.GetFitT())));

            const float thisGapLength((thisFitPosition - previousFitPosition).GetMagnitude());
            const float minZ(std::min(thisFitPosition.GetZ(), previousFitPosition.GetZ()));
            const float maxZ(std::max(thisFitPosition.GetZ(), previousFitPosition.GetZ()));

            if ((maxZ - minZ) > std::numeric_limits<float>::epsilon())
            {
                const float gapZ(LArGeometryHelper::CalculateGapDeltaZ(this->GetPandora(), minZ, maxZ, hitType));
                const float correctedGapLength(thisGapLength * (1.f - gapZ / (maxZ - minZ)));

                if (correctedGapLength > maxFitGapLength)
                    maxFitGapLength = correctedGapLength;
            }

            dTdLMin = std::min(dTdLMin, static_cast<float>(layerFitResult.GetGradient()));
            dTdLMax = std::max(dTdLMax, static_cast<float>(layerFitResult.GetGradient()));
            previousFitPosition = thisFitPosition;
        }

        if (diffWithStraightLineVector.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        diffWithStraightLineMean = 0.f;
        diffWithStraightLineSigma = 0.f;

        for (const float diffWithStraightLine : diffWithStraightLineVector)
            diffWithStraightLineMean += diffWithStraightLine;

        diffWithStraightLineMean /= static_cast<float>(diffWithStraightLineVector.size());

        for (const float diffWithStraightLine : diffWithStraightLineVector)
            diffWithStraightLineSigma += (diffWithStraightLine - diffWithStraightLineMean) * (diffWithStraightLine - diffWithStraightLineMean);

        if (diffWithStraightLineSigma < 0.f)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        diffWithStraightLineSigma = std::sqrt(diffWithStraightLineSigma / static_cast<float>(diffWithStraightLineVector.size()));
        dTdLWidth = dTdLMax - dTdLMin;
    }
    catch (const StatusCodeException &)
    {
        straightLineLengthLarge = -1.f;
        diffWithStraightLineMean = -1.f;
        diffWithStraightLineSigma = -1.f;
        dTdLWidth = -1.f;
        maxFitGapLength = -1.f;
        rmsSlidingLinearFit = -1.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDLinearFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SlidingLinearFitWindowLarge", m_slidingLinearFitWindowLarge));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TwoDVertexDistanceFeatureTool::TwoDVertexDistanceFeatureTool() :
    m_slidingLinearFitWindow(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDVertexDistanceFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float straightLineLength(-1.f), ratio(-1.f);
    try
    {
        const TwoDSlidingFitResult slidingFitResultLarge(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        straightLineLength = (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
        if (straightLineLength > std::numeric_limits<float>::epsilon())
            ratio = (CutClusterCharacterisationAlgorithm::GetVertexDistance(pAlgorithm, pCluster)) / straightLineLength;
    }
    catch (const StatusCodeException &)
    {
        ratio = -1.f;
    }
    featureVector.push_back(ratio);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDVertexDistanceFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder,
    const std::string &featureToolName, const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pCluster);

    if (featureMap.find(featureToolName + "_DistLenRatio") != featureMap.end())
    {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_DistLenRatio");
    featureMap[featureToolName + "_DistLenRatio"] = toolFeatureVec[0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDVertexDistanceFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

PfoHierarchyFeatureTool::PfoHierarchyFeatureTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    CaloHitList parent3DHitList;
    LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, parent3DHitList);
    const unsigned int nParentHits3D(parent3DHitList.size());

    PfoList allDaughtersPfoList;
    LArPfoHelper::GetAllDownstreamPfos(pInputPfo, allDaughtersPfoList);
    const unsigned int nDaughterPfos(allDaughtersPfoList.empty() ? 0 : allDaughtersPfoList.size() - 1);

    unsigned int nDaughterHits3DTotal(0);

    if (nDaughterPfos > 0)
    {
        // ATTN This relies on knowing that the first pfo in allDaughtersPfoList is the input pfo
        allDaughtersPfoList.pop_front();

        for (const ParticleFlowObject *const pDaughterPfo : allDaughtersPfoList)
        {
            CaloHitList daughter3DHitList;
            LArPfoHelper::GetCaloHits(pDaughterPfo, TPC_3D, daughter3DHitList);
            nDaughterHits3DTotal += daughter3DHitList.size();
        }
    }

    const LArMvaHelper::MvaFeature nDaughters(static_cast<double>(nDaughterPfos));
    const LArMvaHelper::MvaFeature nDaughterHits3D(static_cast<double>(nDaughterHits3DTotal));
    const LArMvaHelper::MvaFeature daughterParentNHitsRatio(
        (nParentHits3D > 0) ? static_cast<double>(nDaughterHits3DTotal) / static_cast<double>(nParentHits3D) : 0.);

    featureVector.push_back(nDaughters);
    featureVector.push_back(nDaughterHits3D);
    featureVector.push_back(daughterParentNHitsRatio);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoHierarchyFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder, const std::string &featureToolName,
    const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pInputPfo);

    if (featureMap.find(featureToolName + "_NDaughters") != featureMap.end() ||
        featureMap.find(featureToolName + "_NDaughterHits3D") != featureMap.end() ||
        featureMap.find(featureToolName + "_DaughterParentHitRatio") != featureMap.end())
    {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_NDaughters");
    featureOrder.push_back(featureToolName + "_NDaughterHits3D");
    featureOrder.push_back(featureToolName + "_DaughterParentHitRatio");

    featureMap[featureToolName + "_NDaughters"] = toolFeatureVec[0];
    featureMap[featureToolName + "_NDaughterHits3D"] = toolFeatureVec[1];
    featureMap[featureToolName + "_DaughterParentHitRatio"] = toolFeatureVec[2];
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PfoHierarchyFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ConeChargeFeatureTool::ConeChargeFeatureTool() :
    m_conMinHits(3),
    m_minCharge(0.1f),
    m_conFracRange(0.2f),
    m_MoliereRadius(10.1f),
    m_MoliereRadiusFrac(0.2f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeChargeFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{   
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ClusterList clusterListW;
    LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_W, clusterListW);

    LArMvaHelper::MvaFeature haloTotalRatio, concentration, conicalness;

    if (!clusterListW.empty())
    {
        CaloHitList clusterCaloHitList;
        clusterListW.front()->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);

        const CartesianVector &pfoStart(clusterCaloHitList.front()->GetPositionVector());
        CartesianVector centroid(0.f, 0.f, 0.f);
        LArPcaHelper::EigenVectors eigenVecs;
        LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
        LArPcaHelper::RunPca(clusterCaloHitList, centroid, eigenValues, eigenVecs);

        float chargeCore(0.f), chargeHalo(0.f), chargeCon(0.f);
        this->CalculateChargeDistribution(clusterCaloHitList, pfoStart, eigenVecs[0], chargeCore, chargeHalo, chargeCon);
        haloTotalRatio = (chargeCore + chargeHalo > std::numeric_limits<float>::epsilon()) ? chargeHalo / (chargeCore + chargeHalo) : -1.f;
        concentration = (chargeCore + chargeHalo > std::numeric_limits<float>::epsilon()) ? chargeCon / (chargeCore + chargeHalo) : -1.f;
        const float pfoLength(std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pInputPfo)));
        conicalness = (pfoLength > std::numeric_limits<float>::epsilon())
            ? this->CalculateConicalness(clusterCaloHitList, pfoStart, eigenVecs[0], pfoLength)
            : 1.f;
    }

    featureVector.push_back(haloTotalRatio);
    featureVector.push_back(concentration);
    featureVector.push_back(conicalness);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeChargeFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder, const std::string &featureToolName,
    const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pInputPfo);

    if (featureMap.find(featureToolName + "_HaloTotalRatio") != featureMap.end() ||
        featureMap.find(featureToolName + "_Concentration") != featureMap.end() ||
        featureMap.find(featureToolName + "_Conicalness") != featureMap.end())
    {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_HaloTotalRatio");
    featureOrder.push_back(featureToolName + "_Concentration");
    featureOrder.push_back(featureToolName + "_Conicalness");

    featureMap[featureToolName + "_HaloTotalRatio"] = toolFeatureVec[0];
    featureMap[featureToolName + "_Concentration"] = toolFeatureVec[1];
    featureMap[featureToolName + "_Conicalness"] = toolFeatureVec[2];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeChargeFeatureTool::CalculateChargeDistribution(const CaloHitList &caloHitList, const CartesianVector &pfoStart,
    const CartesianVector &pfoDir, float &chargeCore, float &chargeHalo, float &chargeCon)
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const float distFromTrackFit(((pCaloHit->GetPositionVector() - pfoStart).GetCrossProduct(pfoDir)).GetMagnitude());

        if (distFromTrackFit < m_MoliereRadiusFrac * m_MoliereRadius)
            chargeCore += pCaloHit->GetInputEnergy();
        else
            chargeHalo += pCaloHit->GetInputEnergy();

        chargeCon += pCaloHit->GetInputEnergy() / std::max(1.E-2f, distFromTrackFit); /* Set 1.E-2f to prevent division by 0 and to set max histogram bin as 100 */
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ConeChargeFeatureTool::CalculateConicalness(
    const CaloHitList &caloHitList, const CartesianVector &pfoStart, const CartesianVector &pfoDir, const float pfoLength)
{
    float totalChargeStart(0.f), totalChargeEnd(0.f);
    float chargeConStart(0.f), chargeConEnd(0.f);
    unsigned int nHitsConStart(0), nHitsConEnd(0);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const float distFromTrackFit(((pCaloHit->GetPositionVector() - pfoStart).GetCrossProduct(pfoDir)).GetMagnitude());
        const float hitLength(std::fabs((pCaloHit->GetPositionVector() - pfoStart).GetDotProduct(pfoDir)));

        if (hitLength / pfoLength < m_conFracRange)
        {
            chargeConStart += distFromTrackFit * distFromTrackFit * pCaloHit->GetInputEnergy();
            ++nHitsConStart;
            totalChargeStart += pCaloHit->GetInputEnergy();
        }
        else if (1.f - hitLength / pfoLength < m_conFracRange)
        {
            chargeConEnd += distFromTrackFit * distFromTrackFit * pCaloHit->GetInputEnergy();
            ++nHitsConEnd;
            totalChargeEnd += pCaloHit->GetInputEnergy();
        }
    }

    float conicalness(1.f);

    if (nHitsConStart >= m_conMinHits && nHitsConEnd >= m_conMinHits && totalChargeEnd > m_minCharge &&
        std::sqrt(chargeConStart) > m_minCharge && totalChargeStart > m_minCharge)
        conicalness = (std::sqrt(chargeConEnd / chargeConStart)) / (totalChargeEnd / totalChargeStart);

    return conicalness;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeChargeFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConMinHits", m_conMinHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCharge", m_minCharge));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConFracRange", m_conFracRange));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MoliereRadius", m_MoliereRadius));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MoliereRadiusFrac", m_MoliereRadiusFrac));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ConeChargeFeatureTool_ICARUS::ConeChargeFeatureTool_ICARUS() :
    m_useICARUSCollectionPlane(false),
    m_conMinHits(3),
    m_minCharge(0.1f),
    m_conFracRange(0.2f),
    m_MoliereRadius(10.1f),
    m_MoliereRadiusFrac(0.2f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeChargeFeatureTool_ICARUS::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    LArMvaHelper::MvaFeature haloTotalRatio, concentration, conicalness;

    CaloHitList orderedCaloHitList;

    if (!m_useICARUSCollectionPlane)
    {
        ClusterList clusterListW;
        LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_W, clusterListW);

        if (!clusterListW.empty()) 
        {
            clusterListW.front()->GetOrderedCaloHitList().FillCaloHitList(orderedCaloHitList);

            const CartesianVector &pfoStart(orderedCaloHitList.front()->GetPositionVector());
            CartesianVector centroid(0.f, 0.f, 0.f);
            LArPcaHelper::EigenVectors eigenVecs;
            LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
            LArPcaHelper::RunPca(orderedCaloHitList, centroid, eigenValues, eigenVecs);

            float chargeCore(0.f), chargeHalo(0.f), chargeCon(0.f);
            this->CalculateChargeDistribution(orderedCaloHitList, pfoStart, eigenVecs[0], chargeCore, chargeHalo, chargeCon);
            haloTotalRatio = (chargeCore + chargeHalo > std::numeric_limits<float>::epsilon()) ? chargeHalo / (chargeCore + chargeHalo) : -1.f;
            concentration = (chargeCore + chargeHalo > std::numeric_limits<float>::epsilon()) ? chargeCon / (chargeCore + chargeHalo) : -1.f;
            const float pfoLength(std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pInputPfo)));
            conicalness = (pfoLength > std::numeric_limits<float>::epsilon())
                            ? this->CalculateConicalness(orderedCaloHitList, pfoStart, eigenVecs[0], pfoLength)
                            : 1.f;
        }
    }
    else
    {
        // W 
        ClusterList clusterListW;
        LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_W, clusterListW);
        float minX(9999.), maxX(9999.);
        if (!clusterListW.empty())
            clusterListW.front()->GetClusterSpanX(minX, maxX);

        // U 
        ClusterList clusterListU;
        LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_U, clusterListU);
        float minX_U(9999.), maxX_U(9999.);
        if (!clusterListU.empty())
            clusterListU.front()->GetClusterSpanX(minX_U, maxX_U);

        // V 
        ClusterList clusterListV;
        LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_V, clusterListV);
        float minX_V(9999.), maxX_V(9999.);
        if (!clusterListV.empty())
            clusterListV.front()->GetClusterSpanX(minX_V, maxX_V);

        // Checks                                                             
        if (clusterListU.empty() && clusterListV.empty()) 
        {
            this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListW.front(), orderedCaloHitList);
        }
        else
        {
            // Get drift span from well-defined views
            if (!clusterListU.empty() && !clusterListV.empty()) 
            {
                minX = std::min(minX_U, minX_V); 
                maxX = std::max(maxX_U, maxX_V);  
            }
            else if (!clusterListU.empty() && clusterListV.empty())
            {
                minX = minX_U; 
                maxX = maxX_U; 
            }
            else if (clusterListU.empty() && !clusterListV.empty())
            {
                minX = minX_V; 
                maxX = maxX_V; 
            }

            // Based on whether the PFP crosses the cathode, views are combined to obtain Collection in ICARUS      
            // PFP crosses the cathode                                                          
            if (LocatePointInCryostat_ICARUS(minX) == PositionInCryostat::BelowCathode &&
                LocatePointInCryostat_ICARUS(maxX) == PositionInCryostat::AboveCathode) 
            {
                if (clusterListU.empty() || clusterListV.empty()) 
                {
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListW.front(), orderedCaloHitList);
                }
                else 
                {                                                                                                                                                                                                                                                                                                         
                    CaloHitList orderedCaloHitListU;
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListU.front(), orderedCaloHitListU);
                                                                                                                                                     
                    CaloHitList orderedCaloHitListV;
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListV.front(), orderedCaloHitListV);

                    if (!orderedCaloHitListU.empty() && !orderedCaloHitListV.empty())
                    {                                                   
                        const VertexList *pVertexList(nullptr);
                        (void) PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);
                        const float pVertexX = pVertexList->front()->GetPosition().GetX();
                        
                        if (LocatePointInCryostat_ICARUS(pVertexX) == PositionInCryostat::AboveCathode) 
                        {
                            this->CombineCaloHitListsToHaveCollection(pAlgorithm, 
                                orderedCaloHitListV, orderedCaloHitListU, orderedCaloHitList); ///< V, then U (TPC 2/3)
                        }
                        else if (LocatePointInCryostat_ICARUS(pVertexX) == PositionInCryostat::BelowCathode)
                        {
                            this->CombineCaloHitListsToHaveCollection(pAlgorithm, 
                                orderedCaloHitListU, orderedCaloHitListV, orderedCaloHitList); ///< U, then V (TPC 0/1)
                        }
                        else 
                        {
                            this->CombineCaloHitListsToHaveCollection(pAlgorithm, 
                                orderedCaloHitListU, orderedCaloHitListV, orderedCaloHitList); ///< Ordering cannot be established a priori (vertex within cathode)
                        }
                    }                                                                                                                                                            
                    else     
                    {                                               
                        orderedCaloHitList = orderedCaloHitListU; ///< No vertex exists, either U or V are fine
                    }
                }                                
            }
            // PFP is contained within TPC 2/3
            else if (LocatePointInCryostat_ICARUS(minX) == PositionInCryostat::AboveCathode)
            {
                if (!clusterListV.empty())
                {
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListV.front(), orderedCaloHitList);
                }
                else
                {
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListU.front(), orderedCaloHitList); ///< Fall back to Induction-2
                }
            }
            // PFP is contained within TPC 0/1
            else if (LocatePointInCryostat_ICARUS(maxX) == PositionInCryostat::BelowCathode)
            {
                if (!clusterListU.empty())
                {
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListU.front(), orderedCaloHitList);
                }
                else
                {
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListV.front(), orderedCaloHitList); ///< Fall back to Induction-2
                }
            }
        }

        if (!orderedCaloHitList.empty()) 
        {
            const CartesianVector &pfoStart(orderedCaloHitList.front()->GetPositionVector());
            CartesianVector centroid(0.f, 0.f, 0.f);
            LArPcaHelper::EigenVectors eigenVecs;
            LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
            LArPcaHelper::RunPca(orderedCaloHitList, centroid, eigenValues, eigenVecs);

            float chargeCore(0.f), chargeHalo(0.f), chargeCon(0.f);
            this->CalculateChargeDistribution(orderedCaloHitList, pfoStart, eigenVecs[0], chargeCore, chargeHalo, chargeCon);
            haloTotalRatio = (chargeCore + chargeHalo > std::numeric_limits<float>::epsilon()) ? chargeHalo / (chargeCore + chargeHalo) : -1.f;
            concentration = (chargeCore + chargeHalo > std::numeric_limits<float>::epsilon()) ? chargeCon / (chargeCore + chargeHalo) : -1.f;
            const float pfoLength(std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pInputPfo)));
            conicalness = (pfoLength > std::numeric_limits<float>::epsilon())
                            ? this->CalculateConicalness(orderedCaloHitList, pfoStart, eigenVecs[0], pfoLength)
                            : 1.f;
        }
    }
    
    featureVector.push_back(haloTotalRatio);
    featureVector.push_back(concentration);
    featureVector.push_back(conicalness);
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeChargeFeatureTool_ICARUS::OrderCaloHitsByDistanceToVertex(
    const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, CaloHitList &caloHitList)
{
    const VertexList *pVertexList(nullptr);
    (void)PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (!pVertexList || pVertexList->empty())
        return;

    unsigned int nInteractionVertices(0);
    const Vertex *pInteractionVertex(nullptr);

    for (const Vertex *pVertex : *pVertexList)
    {
        if ((pVertex->GetVertexLabel() == VERTEX_INTERACTION) && (pVertex->GetVertexType() == VERTEX_3D))
        {
            ++nInteractionVertices;
            pInteractionVertex = pVertex;
        }
    }
    
    if (pInteractionVertex && (1 == nInteractionVertices))
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), pInteractionVertex->GetPosition(), hitType));

        CaloHitList clusterCaloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);

        clusterCaloHitList.sort(ConeChargeFeatureTool_ICARUS::VertexComparator(vertexPosition2D));
        caloHitList.insert(caloHitList.end(), clusterCaloHitList.begin(), clusterCaloHitList.end());
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

void ConeChargeFeatureTool_ICARUS::CombineCaloHitListsToHaveCollection(const pandora::Algorithm *const pAlgorithm, const pandora::CaloHitList &orderedCaloHitList1,
                                                                  const pandora::CaloHitList &orderedCaloHitList2, pandora::CaloHitList &mergedCaloHitList)
{

    const VertexList *pVertexList(nullptr);
    (void)PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);
    const CartesianVector vertexPosition2DU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_U));
    const CartesianVector vertexPosition2DV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_V));

    float hitVertexDistance = 0., maxHitVertexDistance1 = 0., minHitVertexDistance2 = 9999.;

    // First hit list
    for (CaloHitList::const_iterator hIter = orderedCaloHitList1.begin(); hIter != orderedCaloHitList1.end(); hIter++) 
    {
        const CaloHit *const pCaloHit = *hIter;
        const CartesianVector &hit(pCaloHit->GetPositionVector());

        const CartesianVector &vertexPosition = (pCaloHit->GetHitType() == TPC_VIEW_U) ? vertexPosition2DU : vertexPosition2DV;
        hitVertexDistance = 
            pow(vertexPosition.GetX() - hit.GetX(), 2) +
            pow(vertexPosition.GetY() - hit.GetY(), 2) +
            pow(vertexPosition.GetZ() - hit.GetZ(), 2);

        bool selectHitViewU = (pCaloHit->GetHitType() == TPC_VIEW_U) && 
            (LocatePointInCryostat_ICARUS(hit.GetX()) == PositionInCryostat::BelowCathode); ///< U is Collection "below cathode"
        bool selectHitViewV = pCaloHit->GetHitType() == TPC_VIEW_V && 
            (LocatePointInCryostat_ICARUS(hit.GetX()) == PositionInCryostat::AboveCathode); ///< V is Collection "below cathode"

        if (selectHitViewU || selectHitViewV)
        {
            mergedCaloHitList.push_back(pCaloHit);
            maxHitVertexDistance1 = std::max(maxHitVertexDistance1, hitVertexDistance);
        }
    }

    // Second hit list                                                                                                                                                                 
    for (CaloHitList::const_iterator hIter = orderedCaloHitList2.begin(); hIter != orderedCaloHitList2.end(); hIter++) 
    {
        const CaloHit *const pCaloHit = *hIter;
        const CartesianVector &hit(pCaloHit->GetPositionVector());

        const CartesianVector &vertexPosition = (pCaloHit->GetHitType() == TPC_VIEW_U) ? vertexPosition2DU : vertexPosition2DV;
        hitVertexDistance = 
            pow(vertexPosition.GetX() - hit.GetX(), 2) +
            pow(vertexPosition.GetY() - hit.GetY(), 2) +
            pow(vertexPosition.GetZ() - hit.GetZ(), 2);

        bool selectHitViewU = (pCaloHit->GetHitType() == TPC_VIEW_U) && 
            (LocatePointInCryostat_ICARUS(hit.GetX()) == PositionInCryostat::BelowCathode); ///< U is Collection "below cathode"
        bool selectHitViewV = pCaloHit->GetHitType() == TPC_VIEW_V && 
            (LocatePointInCryostat_ICARUS(hit.GetX()) == PositionInCryostat::AboveCathode); ///< V is Collection "below cathode"

        if (selectHitViewU || selectHitViewV)
        {
            mergedCaloHitList.push_back(pCaloHit);
            minHitVertexDistance2 = std::min(minHitVertexDistance2, hitVertexDistance);
        }
    }

    // Reorder if needed                                                                                                                   
    if (minHitVertexDistance2 < maxHitVertexDistance1)
    {
        mergedCaloHitList.sort(ConeChargeFeatureTool_ICARUS::DistanceToVertexComparator(vertexPosition2DU, vertexPosition2DV, TPC_VIEW_U, TPC_VIEW_V));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ConeChargeFeatureTool_ICARUS::VertexComparator::VertexComparator(const CartesianVector vertexPosition2D) : m_neutrinoVertex(vertexPosition2D)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConeChargeFeatureTool_ICARUS::VertexComparator::operator()(const CaloHit *const left, const CaloHit *const right) const
{
    const float distanceL((left->GetPositionVector() - m_neutrinoVertex).GetMagnitudeSquared());
    const float distanceR((right->GetPositionVector() - m_neutrinoVertex).GetMagnitudeSquared());
    return distanceL < distanceR;
}

//------------------------------------------------------------------------------------------------------------------------------------------

ConeChargeFeatureTool_ICARUS::DistanceToVertexComparator::DistanceToVertexComparator(const CartesianVector vertexPosition2D_A, const CartesianVector vertexPosition2D_B,
                                                                                const HitType hitType_A, const HitType hitType_B) :
  m_neutrinoVertex_A(vertexPosition2D_A), m_neutrinoVertex_B(vertexPosition2D_B), m_hitType_A(hitType_A), m_hitType_B(hitType_B)
{                                                                                                                                            
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ConeChargeFeatureTool_ICARUS::DistanceToVertexComparator::operator()(const CaloHit *const left, const CaloHit *const right) const
{
    float distanceL, distanceR;

    if (left->GetHitType() == m_hitType_A)
    {
        distanceL = (left->GetPositionVector() - m_neutrinoVertex_A).GetMagnitudeSquared();
    }
    else 
    {
        distanceL = (left->GetPositionVector() - m_neutrinoVertex_B).GetMagnitudeSquared();
    }

    if (right->GetHitType() == m_hitType_A)
    {
        distanceR = (right->GetPositionVector() - m_neutrinoVertex_A).GetMagnitudeSquared();
    }
    else 
    {
        distanceR = (right->GetPositionVector() - m_neutrinoVertex_B).GetMagnitudeSquared();
    }

    return distanceL < distanceR;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeChargeFeatureTool_ICARUS::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder, const std::string &featureToolName,
    const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pInputPfo);

    if (featureMap.find(featureToolName + "_HaloTotalRatio") != featureMap.end() ||
        featureMap.find(featureToolName + "_Concentration") != featureMap.end() ||
        featureMap.find(featureToolName + "_Conicalness") != featureMap.end())
    {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_HaloTotalRatio");
    featureOrder.push_back(featureToolName + "_Concentration");
    featureOrder.push_back(featureToolName + "_Conicalness");

    featureMap[featureToolName + "_HaloTotalRatio"] = toolFeatureVec[0];
    featureMap[featureToolName + "_Concentration"] = toolFeatureVec[1];
    featureMap[featureToolName + "_Conicalness"] = toolFeatureVec[2];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ConeChargeFeatureTool_ICARUS::CalculateChargeDistribution(const CaloHitList &caloHitList, const CartesianVector &pfoStart,
    const CartesianVector &pfoDir, float &chargeCore, float &chargeHalo, float &chargeCon)
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const float distFromTrackFit(((pCaloHit->GetPositionVector() - pfoStart).GetCrossProduct(pfoDir)).GetMagnitude());

        if (distFromTrackFit < m_MoliereRadiusFrac * m_MoliereRadius)
            chargeCore += pCaloHit->GetInputEnergy();
        else
            chargeHalo += pCaloHit->GetInputEnergy();

        chargeCon += pCaloHit->GetInputEnergy() / std::max(1.E-2f, distFromTrackFit); /* Set 1.E-2f to prevent division by 0 and to set max histogram bin as 100 */
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ConeChargeFeatureTool_ICARUS::CalculateConicalness(
    const CaloHitList &caloHitList, const CartesianVector &pfoStart, const CartesianVector &pfoDir, const float pfoLength)
{
    float totalChargeStart(0.f), totalChargeEnd(0.f);
    float chargeConStart(0.f), chargeConEnd(0.f);
    unsigned int nHitsConStart(0), nHitsConEnd(0);

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const float distFromTrackFit(((pCaloHit->GetPositionVector() - pfoStart).GetCrossProduct(pfoDir)).GetMagnitude());
        const float hitLength(std::fabs((pCaloHit->GetPositionVector() - pfoStart).GetDotProduct(pfoDir)));

        if (hitLength / pfoLength < m_conFracRange)
        {
            chargeConStart += distFromTrackFit * distFromTrackFit * pCaloHit->GetInputEnergy();
            ++nHitsConStart;
            totalChargeStart += pCaloHit->GetInputEnergy();
        }
        else if (1.f - hitLength / pfoLength < m_conFracRange)
        {
            chargeConEnd += distFromTrackFit * distFromTrackFit * pCaloHit->GetInputEnergy();
            ++nHitsConEnd;
            totalChargeEnd += pCaloHit->GetInputEnergy();
        }
    }

    float conicalness(1.f);

    if (nHitsConStart >= m_conMinHits && nHitsConEnd >= m_conMinHits && totalChargeEnd > m_minCharge &&
        std::sqrt(chargeConStart) > m_minCharge && totalChargeStart > m_minCharge)
        conicalness = (std::sqrt(chargeConEnd / chargeConStart)) / (totalChargeEnd / totalChargeStart);

    return conicalness;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ConeChargeFeatureTool_ICARUS::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseICARUSCollectionPlane", m_useICARUSCollectionPlane));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConMinHits", m_conMinHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCharge", m_minCharge));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ConFracRange", m_conFracRange));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MoliereRadius", m_MoliereRadius));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MoliereRadiusFrac", m_MoliereRadiusFrac));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDLinearFitFeatureTool::ThreeDLinearFitFeatureTool() : m_slidingLinearFitWindow(3), m_slidingLinearFitWindowLarge(10000)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLinearFitFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    ClusterList clusterList;
    LArPfoHelper::GetTwoDClusterList(pInputPfo, clusterList);
    float diffWithStraightLineMean(0.f), maxFitGapLength(0.f), rmsSlidingLinearFit(0.f);
    LArMvaHelper::MvaFeature length, diff, gap, rms;
    unsigned int nClustersUsed(0);

    for (const Cluster *const pCluster : clusterList)
    {
        float straightLineLengthLargeCluster(-1.f), diffWithStraightLineMeanCluster(-1.f), maxFitGapLengthCluster(-1.f),
            rmsSlidingLinearFitCluster(-1.f);

        this->CalculateVariablesSlidingLinearFit(
            pCluster, straightLineLengthLargeCluster, diffWithStraightLineMeanCluster, maxFitGapLengthCluster, rmsSlidingLinearFitCluster);

        if (straightLineLengthLargeCluster > std::numeric_limits<float>::epsilon())
        {
            diffWithStraightLineMeanCluster /= straightLineLengthLargeCluster;
            maxFitGapLengthCluster /= straightLineLengthLargeCluster;
            rmsSlidingLinearFitCluster /= straightLineLengthLargeCluster;

            diffWithStraightLineMean += diffWithStraightLineMeanCluster;
            maxFitGapLength += maxFitGapLengthCluster;
            rmsSlidingLinearFit += rmsSlidingLinearFitCluster;

            ++nClustersUsed;
        }
    }

    if (nClustersUsed > 0)
    {
        const float nClusters(static_cast<float>(nClustersUsed));
        length = std::sqrt(LArPfoHelper::GetThreeDLengthSquared(pInputPfo));
        diff = diffWithStraightLineMean / nClusters;
        gap = maxFitGapLength / nClusters;
        rms = rmsSlidingLinearFit / nClusters;
    }

    featureVector.push_back(length);
    featureVector.push_back(diff);
    featureVector.push_back(gap);
    featureVector.push_back(rms);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLinearFitFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder,
    const std::string &featureToolName, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pInputPfo);

    if (featureMap.find(featureToolName + "_Length") != featureMap.end() ||
        featureMap.find(featureToolName + "_DiffStraightLineMean") != featureMap.end() ||
        featureMap.find(featureToolName + "_MaxFitGapLength") != featureMap.end() ||
        featureMap.find(featureToolName + "_SlidingLinearFitRMS") != featureMap.end())
    {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_Length");
    featureOrder.push_back(featureToolName + "_DiffStraightLineMean");
    featureOrder.push_back(featureToolName + "_MaxFitGapLength");
    featureOrder.push_back(featureToolName + "_SlidingLinearFitRMS");

    featureMap[featureToolName + "_Length"] = toolFeatureVec[0];
    featureMap[featureToolName + "_DiffStraightLineMean"] = toolFeatureVec[1];
    featureMap[featureToolName + "_MaxFitGapLength"] = toolFeatureVec[2];
    featureMap[featureToolName + "_SlidingLinearFitRMS"] = toolFeatureVec[3];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDLinearFitFeatureTool::CalculateVariablesSlidingLinearFit(const pandora::Cluster *const pCluster, float &straightLineLengthLarge,
    float &diffWithStraightLineMean, float &maxFitGapLength, float &rmsSlidingLinearFit) const
{
    try
    {
        const TwoDSlidingFitResult slidingFitResult(pCluster, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        const TwoDSlidingFitResult slidingFitResultLarge(
            pCluster, m_slidingLinearFitWindowLarge, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

        if (slidingFitResult.GetLayerFitResultMap().empty())
            throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

        straightLineLengthLarge =
            (slidingFitResultLarge.GetGlobalMaxLayerPosition() - slidingFitResultLarge.GetGlobalMinLayerPosition()).GetMagnitude();
        rmsSlidingLinearFit = 0.f;

        FloatVector diffWithStraightLineVector;
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        CartesianVector previousFitPosition(slidingFitResult.GetGlobalMinLayerPosition());
        float dTdLMin(+std::numeric_limits<float>::max()), dTdLMax(-std::numeric_limits<float>::max());

        for (const auto &mapEntry : slidingFitResult.GetLayerFitResultMap())
        {
            const LayerFitResult &layerFitResult(mapEntry.second);
            rmsSlidingLinearFit += layerFitResult.GetRms();

            CartesianVector thisFitPosition(0.f, 0.f, 0.f);
            slidingFitResult.GetGlobalPosition(layerFitResult.GetL(), layerFitResult.GetFitT(), thisFitPosition);

            LayerFitResultMap::const_iterator iterLarge =
                slidingFitResultLarge.GetLayerFitResultMap().find(slidingFitResultLarge.GetLayer(layerFitResult.GetL()));

            if (slidingFitResultLarge.GetLayerFitResultMap().end() == iterLarge)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            diffWithStraightLineVector.push_back(static_cast<float>(std::fabs(layerFitResult.GetFitT() - iterLarge->second.GetFitT())));

            const float thisGapLength((thisFitPosition - previousFitPosition).GetMagnitude());
            const float minZ(std::min(thisFitPosition.GetZ(), previousFitPosition.GetZ()));
            const float maxZ(std::max(thisFitPosition.GetZ(), previousFitPosition.GetZ()));

            if ((maxZ - minZ) > std::numeric_limits<float>::epsilon())
            {
                const float gapZ(LArGeometryHelper::CalculateGapDeltaZ(this->GetPandora(), minZ, maxZ, hitType));
                const float correctedGapLength(thisGapLength * (1.f - gapZ / (maxZ - minZ)));

                if (correctedGapLength > maxFitGapLength)
                    maxFitGapLength = correctedGapLength;
            }
            else
            {
                maxFitGapLength = 0.f;
            }

            dTdLMin = std::min(dTdLMin, static_cast<float>(layerFitResult.GetGradient()));
            dTdLMax = std::max(dTdLMax, static_cast<float>(layerFitResult.GetGradient()));
            previousFitPosition = thisFitPosition;
        }

        if (diffWithStraightLineVector.empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        diffWithStraightLineMean = 0.f;

        for (const float diffWithStraightLine : diffWithStraightLineVector)
            diffWithStraightLineMean += diffWithStraightLine;

        diffWithStraightLineMean /= static_cast<float>(diffWithStraightLineVector.size());
    }
    catch (const StatusCodeException &)
    {
        straightLineLengthLarge = -1.f;
        diffWithStraightLineMean = -1.f;
        maxFitGapLength = -1.f;
        rmsSlidingLinearFit = -1.f;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDLinearFitFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingLinearFitWindow", m_slidingLinearFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "SlidingLinearFitWindowLarge", m_slidingLinearFitWindowLarge));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDVertexDistanceFeatureTool::ThreeDVertexDistanceFeatureTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDVertexDistanceFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    LArMvaHelper::MvaFeature vertexDistance;

    const VertexList *pVertexList(nullptr);
    (void)PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (!pVertexList || pVertexList->empty())
    {
        featureVector.push_back(vertexDistance);
        return;
    }

    unsigned int nInteractionVertices(0);
    const Vertex *pInteractionVertex(nullptr);

    for (const Vertex *pVertex : *pVertexList)
    {
        if ((pVertex->GetVertexLabel() == VERTEX_INTERACTION) && (pVertex->GetVertexType() == VERTEX_3D))
        {
            ++nInteractionVertices;
            pInteractionVertex = pVertex;
        }
    }

    if (pInteractionVertex && (1 == nInteractionVertices))
    {
        try
        {
            vertexDistance = (pInteractionVertex->GetPosition() - LArPfoHelper::GetVertex(pInputPfo)->GetPosition()).GetMagnitude();
        }
        catch (const StatusCodeException &)
        {
            CaloHitList threeDCaloHitList;
            LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, threeDCaloHitList);

            if (!threeDCaloHitList.empty())
                vertexDistance = (pInteractionVertex->GetPosition() - (threeDCaloHitList.front())->GetPositionVector()).GetMagnitude();
        }
    }

    featureVector.push_back(vertexDistance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDVertexDistanceFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder,
    const std::string &featureToolName, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pInputPfo);

    if (featureMap.find(featureToolName + "_VertexDistance") != featureMap.end())
    {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_VertexDistance");

    featureMap[featureToolName + "_VertexDistance"] = toolFeatureVec[0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDVertexDistanceFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDOpeningAngleFeatureTool::ThreeDOpeningAngleFeatureTool() :
    m_hitFraction(0.5f),
    m_defaultValue(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDOpeningAngleFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    // Need the 3D hits to calculate PCA components
    CaloHitList threeDCaloHitList;
    LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, threeDCaloHitList);

    LArMvaHelper::MvaFeature diffAngle;
    if (!threeDCaloHitList.empty())
    {
        CartesianPointVector pointVectorStart, pointVectorEnd;
        this->Divide3DCaloHitList(pAlgorithm, threeDCaloHitList, pointVectorStart, pointVectorEnd);

        // Able to calculate angles only if > 1 point provided
        if ((pointVectorStart.size() > 1) && (pointVectorEnd.size() > 1))
        {
            try
            {
                // Run the PCA analysis twice
                CartesianVector centroidStart(0.f, 0.f, 0.f), centroidEnd(0.f, 0.f, 0.f);
                LArPcaHelper::EigenVectors eigenVecsStart, eigenVecsEnd;
                LArPcaHelper::EigenValues eigenValuesStart(0.f, 0.f, 0.f), eigenValuesEnd(0.f, 0.f, 0.f);

                LArPcaHelper::RunPca(pointVectorStart, centroidStart, eigenValuesStart, eigenVecsStart);
                LArPcaHelper::RunPca(pointVectorEnd, centroidEnd, eigenValuesEnd, eigenVecsEnd);

                const float openingAngle(this->OpeningAngle(eigenVecsStart.at(0), eigenVecsStart.at(1), eigenValuesStart));
                const float closingAngle(this->OpeningAngle(eigenVecsEnd.at(0), eigenVecsEnd.at(1), eigenValuesEnd));
                diffAngle = std::fabs(openingAngle - closingAngle);
            }
            catch (const StatusCodeException &)
            {
            }
        }
        else
        {
            diffAngle = m_defaultValue;
        }
    }

    featureVector.push_back(diffAngle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDOpeningAngleFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder,
    const std::string &featureToolName, const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pInputPfo);

    if (featureMap.find(featureToolName + "_AngleDiff") != featureMap.end())
    {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_AngleDiff");

    featureMap[featureToolName + "_AngleDiff"] = toolFeatureVec[0];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDOpeningAngleFeatureTool::Divide3DCaloHitList(const Algorithm *const pAlgorithm, const CaloHitList &threeDCaloHitList,
    CartesianPointVector &pointVectorStart, CartesianPointVector &pointVectorEnd)
{
    const VertexList *pVertexList(nullptr);
    (void)PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (threeDCaloHitList.empty() || !pVertexList || pVertexList->empty())
        return;

    unsigned int nInteractionVertices(0);
    const Vertex *pInteractionVertex(nullptr);

    for (const Vertex *pVertex : *pVertexList)
    {
        if ((pVertex->GetVertexLabel() == VERTEX_INTERACTION) && (pVertex->GetVertexType() == VERTEX_3D))
        {
            ++nInteractionVertices;
            pInteractionVertex = pVertex;
        }
    }

    if (pInteractionVertex && (1 == nInteractionVertices))
    {
        // Order by distance to vertex, so first ones are closer to nuvertex
        CaloHitVector threeDCaloHitVector(threeDCaloHitList.begin(), threeDCaloHitList.end());
        std::sort(threeDCaloHitVector.begin(), threeDCaloHitVector.end(),
            ThreeDChargeFeatureTool::VertexComparator(pInteractionVertex->GetPosition()));

        unsigned int iHit(1);
        const unsigned int nHits(threeDCaloHitVector.size());

        for (const CaloHit *const pCaloHit : threeDCaloHitVector)
        {
            if (static_cast<float>(iHit) / static_cast<float>(nHits) <= m_hitFraction)
                pointVectorStart.push_back(pCaloHit->GetPositionVector());

            if (static_cast<float>(iHit) / static_cast<float>(nHits) >= 1.f - m_hitFraction)
                pointVectorEnd.push_back(pCaloHit->GetPositionVector());

            ++iHit;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ThreeDOpeningAngleFeatureTool::OpeningAngle(const CartesianVector &principal, const CartesianVector &secondary, const CartesianVector &eigenValues) const
{
    const float principalMagnitude(principal.GetMagnitude());
    const float secondaryMagnitude(secondary.GetMagnitude());

    if (std::fabs(principalMagnitude) < std::numeric_limits<float>::epsilon())
    {
        std::cout << "ThreeDOpeningAngleFeatureTool::OpeningAngle - The principal eigenvector is 0." << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
    else if (std::fabs(secondaryMagnitude) < std::numeric_limits<float>::epsilon())
    {
        return 0.f;
    }

    const float cosTheta(principal.GetDotProduct(secondary) / (principalMagnitude * secondaryMagnitude));

    if (cosTheta > 1.f)
    {
        std::cout << "PcaShowerParticleBuildingAlgorithm::OpeningAngle - cos(theta) reportedly greater than 1." << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    const float sinTheta(std::sqrt(1.f - cosTheta * cosTheta));

    if (eigenValues.GetX() < std::numeric_limits<float>::epsilon())
    {
        std::cout << "PcaShowerParticleBuildingAlgorithm::OpeningAngle - principal eigenvalue less than or equal to 0." << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }
    else if (eigenValues.GetY() < std::numeric_limits<float>::epsilon())
    {
        return 0.f;
    }

    return std::atan(std::sqrt(eigenValues.GetY()) * sinTheta / std::sqrt(eigenValues.GetX()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDOpeningAngleFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "HitFraction", m_hitFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DefaultValue", m_defaultValue));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDPCAFeatureTool::ThreeDPCAFeatureTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDPCAFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    LArMvaHelper::MvaFeature pca1, pca2;

    // Need the 3D hits to calculate PCA components
    CaloHitList threeDCaloHitList;
    LArPfoHelper::GetCaloHits(pInputPfo, TPC_3D, threeDCaloHitList);

    if (!threeDCaloHitList.empty())
    {
        try
        {
            CartesianVector centroid(0.f, 0.f, 0.f);
            LArPcaHelper::EigenVectors eigenVecs;
            LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);

            LArPcaHelper::RunPca(threeDCaloHitList, centroid, eigenValues, eigenVecs);
            const float principalEigenvalue(eigenValues.GetX()), secondaryEigenvalue(eigenValues.GetY()), tertiaryEigenvalue(eigenValues.GetZ());

            if (principalEigenvalue > std::numeric_limits<float>::epsilon())
            {
                pca1 = secondaryEigenvalue / principalEigenvalue;
                pca2 = tertiaryEigenvalue / principalEigenvalue;
            }
            else
            {
                // ATTN if n3dHits == 1 then principal, secondary, and tertiary eigenvalues are zero hence default to zero
                pca1 = 0.;
                pca2 = 0.;
            }
        }
        catch (const StatusCodeException &)
        {
        }
    }

    featureVector.push_back(pca1);
    featureVector.push_back(pca2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDPCAFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder, const std::string &featureToolName,
    const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pInputPfo);

    if (featureMap.find(featureToolName + "_SecondaryPCARatio") != featureMap.end() ||
        featureMap.find(featureToolName + "_TertiaryPCARatio") != featureMap.end())
    {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_SecondaryPCARatio");
    featureOrder.push_back(featureToolName + "_TertiaryPCARatio");

    featureMap[featureToolName + "_SecondaryPCARatio"] = toolFeatureVec[0];
    featureMap[featureToolName + "_TertiaryPCARatio"] = toolFeatureVec[1];
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDPCAFeatureTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDChargeFeatureTool::ThreeDChargeFeatureTool() :
    m_endChargeFraction(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float totalCharge(-1.f), chargeSigma(-1.f), chargeMean(-1.f), endCharge(-1.f);
    LArMvaHelper::MvaFeature charge1, charge2;

    ClusterList clusterListW;
    LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_W, clusterListW);

    if (!clusterListW.empty())
        this->CalculateChargeVariables(pAlgorithm, clusterListW.front(), totalCharge, chargeSigma, chargeMean, endCharge);

    if (chargeMean > std::numeric_limits<float>::epsilon())
        charge1 = chargeSigma / chargeMean;

    if (totalCharge > std::numeric_limits<float>::epsilon())
        charge2 = endCharge / totalCharge;

    featureVector.push_back(charge1);
    featureVector.push_back(charge2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder, const std::string &featureToolName,
    const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pInputPfo);

    if (featureMap.find(featureToolName + "_FractionalSpread") != featureMap.end() ||
        featureMap.find(featureToolName + "_EndFraction") != featureMap.end())
    {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_FractionalSpread");
    featureOrder.push_back(featureToolName + "_EndFraction");

    featureMap[featureToolName + "_FractionalSpread"] = toolFeatureVec[0];
    featureMap[featureToolName + "_EndFraction"] = toolFeatureVec[1];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::CalculateChargeVariables(const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster,
    float &totalCharge, float &chargeSigma, float &chargeMean, float &endCharge)
{
    totalCharge = 0.f;
    chargeSigma = 0.f;
    chargeMean = 0.f;
    endCharge = 0.f;

    CaloHitList orderedCaloHitList;
    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, pCluster, orderedCaloHitList);

    FloatVector chargeVector;
    unsigned int hitCounter(0);
    const unsigned int nTotalHits(orderedCaloHitList.size());

    for (const CaloHit *const pCaloHit : orderedCaloHitList)
    {
        ++hitCounter;
        const float pCaloHitCharge(pCaloHit->GetInputEnergy());

        if (pCaloHitCharge >= 0.f)
        {
            totalCharge += pCaloHitCharge;
            chargeVector.push_back(pCaloHitCharge);

            if (hitCounter >= std::floor(static_cast<float>(nTotalHits) * (1.f - m_endChargeFraction)))
                endCharge += pCaloHitCharge;
        }
    }

    if (!chargeVector.empty())
    {
        chargeMean = totalCharge / static_cast<float>(chargeVector.size());

        for (const float charge : chargeVector)
            chargeSigma += (charge - chargeMean) * (charge - chargeMean);

        chargeSigma = std::sqrt(chargeSigma / static_cast<float>(chargeVector.size()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool::OrderCaloHitsByDistanceToVertex(
    const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, CaloHitList &caloHitList)
{
    const VertexList *pVertexList(nullptr);
    (void)PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (!pVertexList || pVertexList->empty())
        return;

    unsigned int nInteractionVertices(0);
    const Vertex *pInteractionVertex(nullptr);

    for (const Vertex *pVertex : *pVertexList)
    {
        if ((pVertex->GetVertexLabel() == VERTEX_INTERACTION) && (pVertex->GetVertexType() == VERTEX_3D))
        {
            ++nInteractionVertices;
            pInteractionVertex = pVertex;
        }
    }

    if (pInteractionVertex && (1 == nInteractionVertices))
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), pInteractionVertex->GetPosition(), hitType));

        CaloHitList clusterCaloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);

        clusterCaloHitList.sort(ThreeDChargeFeatureTool::VertexComparator(vertexPosition2D));
        caloHitList.insert(caloHitList.end(), clusterCaloHitList.begin(), clusterCaloHitList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDChargeFeatureTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EndChargeFraction", m_endChargeFraction));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDChargeFeatureTool::VertexComparator::VertexComparator(const CartesianVector vertexPosition2D) :
    m_neutrinoVertex(vertexPosition2D)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDChargeFeatureTool::VertexComparator::operator()(const CaloHit *const left, const CaloHit *const right) const
{
    const float distanceL((left->GetPositionVector() - m_neutrinoVertex).GetMagnitudeSquared());
    const float distanceR((right->GetPositionVector() - m_neutrinoVertex).GetMagnitudeSquared());
    return distanceL < distanceR;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDChargeFeatureTool_ICARUS::ThreeDChargeFeatureTool_ICARUS() : m_useICARUSCollectionPlane(false), m_endChargeFraction(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool_ICARUS::Run(
    LArMvaHelper::MvaFeatureVector &featureVector, const Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
        std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    float totalCharge(-1.f), chargeSigma(-1.f), chargeMean(-1.f), endCharge(-1.f);
    LArMvaHelper::MvaFeature charge1, charge2;
                                                                                       
    CaloHitList orderedCaloHitList;

    if (!m_useICARUSCollectionPlane)
    {
        ClusterList clusterListW;
        LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_W, clusterListW);

        if (!clusterListW.empty()) 
        {
            this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListW.front(), orderedCaloHitList);
            this->CalculateChargeVariables(orderedCaloHitList, totalCharge, chargeSigma, chargeMean, endCharge);
        }
    }
    else
    {
        // W 
        ClusterList clusterListW;
        LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_W, clusterListW);
        float minX(9999.), maxX(9999.);
        if (!clusterListW.empty())
            clusterListW.front()->GetClusterSpanX(minX, maxX);

        // U 
        ClusterList clusterListU;
        LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_U, clusterListU);
        float minX_U(9999.), maxX_U(9999.);
        if (!clusterListU.empty())
            clusterListU.front()->GetClusterSpanX(minX_U, maxX_U);

        // V 
        ClusterList clusterListV;
        LArPfoHelper::GetClusters(pInputPfo, TPC_VIEW_V, clusterListV);
        float minX_V(9999.), maxX_V(9999.);
        if (!clusterListV.empty())
            clusterListV.front()->GetClusterSpanX(minX_V, maxX_V);

        // Checks                                                             
        if (clusterListU.empty() && clusterListV.empty()) 
        {
            this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListW.front(), orderedCaloHitList);
        }
        else
        {
            // Get drift span from well-defined views
            if (!clusterListU.empty() && !clusterListV.empty()) 
            {
                minX = std::min(minX_U, minX_V); 
                maxX = std::max(maxX_U, maxX_V);  
            }
            else if (!clusterListU.empty() && clusterListV.empty())
            {
                minX = minX_U; 
                maxX = maxX_U; 
            }
            else if (clusterListU.empty() && !clusterListV.empty())
            {
                minX = minX_V; 
                maxX = maxX_V; 
            }

            // Based on whether the PFP crosses the cathode, views are combined to obtain Collection in ICARUS      
            // PFP crosses the cathode                                                          
            if (LocatePointInCryostat_ICARUS(minX) == PositionInCryostat::BelowCathode &&
                LocatePointInCryostat_ICARUS(maxX) == PositionInCryostat::AboveCathode) 
            {
                if (clusterListU.empty() || clusterListV.empty()) 
                {
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListW.front(), orderedCaloHitList);
                }
                else 
                {                                                                                                                                                                                                                                                                                                         
                    CaloHitList orderedCaloHitListU;
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListU.front(), orderedCaloHitListU);
                                                                                                                                                     
                    CaloHitList orderedCaloHitListV;
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListV.front(), orderedCaloHitListV);

                    if (!orderedCaloHitListU.empty() && !orderedCaloHitListV.empty())
                    {                                                   
                        const VertexList *pVertexList(nullptr);
                        (void) PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);
                        const float pVertexX = pVertexList->front()->GetPosition().GetX();
                        
                        if (LocatePointInCryostat_ICARUS(pVertexX) == PositionInCryostat::AboveCathode) 
                        {
                            this->CombineCaloHitListsToHaveCollection(pAlgorithm, 
                                orderedCaloHitListV, orderedCaloHitListU, orderedCaloHitList); ///< V, then U (TPC 2/3)
                        }
                        else if (LocatePointInCryostat_ICARUS(pVertexX) == PositionInCryostat::BelowCathode)
                        {
                            this->CombineCaloHitListsToHaveCollection(pAlgorithm, 
                                orderedCaloHitListU, orderedCaloHitListV, orderedCaloHitList); ///< U, then V (TPC 0/1)
                        }
                        else 
                        {
                            this->CombineCaloHitListsToHaveCollection(pAlgorithm, 
                                orderedCaloHitListU, orderedCaloHitListV, orderedCaloHitList); ///< Ordering cannot be established a priori (vertex within cathode)
                        }
                    }                                                                                                                                                            
                    else     
                    {                                               
                        orderedCaloHitList = orderedCaloHitListU; ///< No vertex exists, either U or V are fine
                    }
                }                                
            }
            // PFP is contained within TPC 2/3
            else if (LocatePointInCryostat_ICARUS(minX) == PositionInCryostat::AboveCathode)
            {
                if (!clusterListV.empty())
                {
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListV.front(), orderedCaloHitList);
                }
                else
                {
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListU.front(), orderedCaloHitList); ///< Fall back to Induction-2
                }
            }
            // PFP is contained within TPC 0/1
            else if (LocatePointInCryostat_ICARUS(maxX) == PositionInCryostat::BelowCathode)
            {
                if (!clusterListU.empty())
                {
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListU.front(), orderedCaloHitList);
                }
                else
                {
                    this->OrderCaloHitsByDistanceToVertex(pAlgorithm, clusterListV.front(), orderedCaloHitList); ///< Fall back to Induction-2
                }
            }
        }

        // Compute charge variables
        if (!orderedCaloHitList.empty())
            this->CalculateChargeVariables(orderedCaloHitList, totalCharge, chargeSigma, chargeMean, endCharge);
    }

    if (chargeMean > std::numeric_limits<float>::epsilon())
        charge1 = chargeSigma / chargeMean;

    if (totalCharge > std::numeric_limits<float>::epsilon())
        charge2 = endCharge / totalCharge;

    featureVector.push_back(charge1);
    featureVector.push_back(charge2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool_ICARUS::Run(LArMvaHelper::MvaFeatureMap &featureMap, StringVector &featureOrder, const std::string &featureToolName,
    const pandora::Algorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pInputPfo)
{
    LArMvaHelper::MvaFeatureVector toolFeatureVec;
    this->Run(toolFeatureVec, pAlgorithm, pInputPfo);

    if (featureMap.find(featureToolName + "_FractionalSpread") != featureMap.end() ||
        featureMap.find(featureToolName + "_EndFraction") != featureMap.end())
    {
        std::cout << "Already wrote this feature into map! Not writing again." << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    featureOrder.push_back(featureToolName + "_FractionalSpread");
    featureOrder.push_back(featureToolName + "_EndFraction");

    featureMap[featureToolName + "_FractionalSpread"] = toolFeatureVec[0];
    featureMap[featureToolName + "_EndFraction"] = toolFeatureVec[1];
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool_ICARUS::CalculateChargeVariables(const CaloHitList &orderedCaloHitList, float &totalCharge, 
    float &chargeSigma, float &chargeMean, float &endCharge)
{
    totalCharge = 0.f;
    chargeSigma = 0.f;
    chargeMean = 0.f;
    endCharge = 0.f;

    FloatVector chargeVector;
    unsigned int hitCounter(0);
    const unsigned int nTotalHits(orderedCaloHitList.size());

    for (const CaloHit *const pCaloHit : orderedCaloHitList)
    {
        ++hitCounter;
        const float pCaloHitCharge(pCaloHit->GetInputEnergy());

        if (pCaloHitCharge >= 0.f)
        {
            totalCharge += pCaloHitCharge;
            chargeVector.push_back(pCaloHitCharge);

            if (hitCounter >= std::floor(static_cast<float>(nTotalHits) * (1.f - m_endChargeFraction)))
                endCharge += pCaloHitCharge;
        }
    }

    if (!chargeVector.empty())
    {
        chargeMean = totalCharge / static_cast<float>(chargeVector.size());

        for (const float charge : chargeVector)
            chargeSigma += (charge - chargeMean) * (charge - chargeMean);

        chargeSigma = std::sqrt(chargeSigma / static_cast<float>(chargeVector.size()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool_ICARUS::OrderCaloHitsByDistanceToVertex(
    const Algorithm *const pAlgorithm, const pandora::Cluster *const pCluster, CaloHitList &caloHitList)
{
    const VertexList *pVertexList(nullptr);
    (void)PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);

    if (!pVertexList || pVertexList->empty())
        return;

    unsigned int nInteractionVertices(0);
    const Vertex *pInteractionVertex(nullptr);

    for (const Vertex *pVertex : *pVertexList)
    {
        if ((pVertex->GetVertexLabel() == VERTEX_INTERACTION) && (pVertex->GetVertexType() == VERTEX_3D))
        {
            ++nInteractionVertices;
            pInteractionVertex = pVertex;
        }
    }
    
    if (pInteractionVertex && (1 == nInteractionVertices))
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        const CartesianVector vertexPosition2D(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), pInteractionVertex->GetPosition(), hitType));

        CaloHitList clusterCaloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);

        clusterCaloHitList.sort(ThreeDChargeFeatureTool_ICARUS::VertexComparator(vertexPosition2D));
        caloHitList.insert(caloHitList.end(), clusterCaloHitList.begin(), clusterCaloHitList.end());
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDChargeFeatureTool_ICARUS::CombineCaloHitListsToHaveCollection(const pandora::Algorithm *const pAlgorithm, const pandora::CaloHitList &orderedCaloHitList1,
                                                                  const pandora::CaloHitList &orderedCaloHitList2, pandora::CaloHitList &mergedCaloHitList)
 {
                                                                                                                                                                  
    const VertexList *pVertexList(nullptr);
    (void)PandoraContentApi::GetCurrentList(*pAlgorithm, pVertexList);
    const CartesianVector vertexPosition2DU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_U));
    const CartesianVector vertexPosition2DV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), pVertexList->front()->GetPosition(), TPC_VIEW_V));

    float hitVertexDistance = 0., maxHitVertexDistance1 = 0., minHitVertexDistance2 = 9999.;

    // First hit list
    for (CaloHitList::const_iterator hIter = orderedCaloHitList1.begin(); hIter != orderedCaloHitList1.end(); hIter++) 
    {
        const CaloHit *const pCaloHit = *hIter;
        const CartesianVector &hit(pCaloHit->GetPositionVector());

        const CartesianVector &vertexPosition = (pCaloHit->GetHitType() == TPC_VIEW_U) ? vertexPosition2DU : vertexPosition2DV;
        hitVertexDistance = 
            pow(vertexPosition.GetX() - hit.GetX(), 2) +
            pow(vertexPosition.GetY() - hit.GetY(), 2) +
            pow(vertexPosition.GetZ() - hit.GetZ(), 2);

        bool selectHitViewU = (pCaloHit->GetHitType() == TPC_VIEW_U) && 
            (LocatePointInCryostat_ICARUS(hit.GetX()) == PositionInCryostat::BelowCathode); ///< U is Collection "below cathode"
        bool selectHitViewV = pCaloHit->GetHitType() == TPC_VIEW_V && 
            (LocatePointInCryostat_ICARUS(hit.GetX()) == PositionInCryostat::AboveCathode); ///< V is Collection "below cathode"

        if (selectHitViewU || selectHitViewV)
        {
            mergedCaloHitList.push_back(pCaloHit);
            maxHitVertexDistance1 = std::max(maxHitVertexDistance1, hitVertexDistance);
        }
    }

    // Second hit list                                                                                                                                                                 
    for (CaloHitList::const_iterator hIter = orderedCaloHitList2.begin(); hIter != orderedCaloHitList2.end(); hIter++) 
    {
        const CaloHit *const pCaloHit = *hIter;
        const CartesianVector &hit(pCaloHit->GetPositionVector());

        const CartesianVector &vertexPosition = (pCaloHit->GetHitType() == TPC_VIEW_U) ? vertexPosition2DU : vertexPosition2DV;
        hitVertexDistance = 
            pow(vertexPosition.GetX() - hit.GetX(), 2) +
            pow(vertexPosition.GetY() - hit.GetY(), 2) +
            pow(vertexPosition.GetZ() - hit.GetZ(), 2);

        bool selectHitViewU = (pCaloHit->GetHitType() == TPC_VIEW_U) && 
            (LocatePointInCryostat_ICARUS(hit.GetX()) == PositionInCryostat::BelowCathode); ///< U is Collection "below cathode"
        bool selectHitViewV = pCaloHit->GetHitType() == TPC_VIEW_V && 
            (LocatePointInCryostat_ICARUS(hit.GetX()) == PositionInCryostat::AboveCathode); ///< V is Collection "below cathode"

        if (selectHitViewU || selectHitViewV)
        {
            mergedCaloHitList.push_back(pCaloHit);
            minHitVertexDistance2 = std::min(minHitVertexDistance2, hitVertexDistance);
        }
    }

    // Reorder if needed                                                                                                                   
    if (minHitVertexDistance2 < maxHitVertexDistance1)
    {
        mergedCaloHitList.sort(ThreeDChargeFeatureTool_ICARUS::DistanceToVertexComparator(vertexPosition2DU, vertexPosition2DV, TPC_VIEW_U, TPC_VIEW_V));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDChargeFeatureTool_ICARUS::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "UseICARUSCollectionPlane", m_useICARUSCollectionPlane));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "EndChargeFraction", m_endChargeFraction));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDChargeFeatureTool_ICARUS::VertexComparator::VertexComparator(const CartesianVector vertexPosition2D) : m_neutrinoVertex(vertexPosition2D)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDChargeFeatureTool_ICARUS::VertexComparator::operator()(const CaloHit *const left, const CaloHit *const right) const
{
    const float distanceL((left->GetPositionVector() - m_neutrinoVertex).GetMagnitudeSquared());
    const float distanceR((right->GetPositionVector() - m_neutrinoVertex).GetMagnitudeSquared());
    return distanceL < distanceR;
}

//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDChargeFeatureTool_ICARUS::DistanceToVertexComparator::DistanceToVertexComparator(const CartesianVector vertexPosition2D_A, const CartesianVector vertexPosition2D_B,
                                                                                const HitType hitType_A, const HitType hitType_B) :
  m_neutrinoVertex_A(vertexPosition2D_A), m_neutrinoVertex_B(vertexPosition2D_B), m_hitType_A(hitType_A), m_hitType_B(hitType_B)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDChargeFeatureTool_ICARUS::DistanceToVertexComparator::operator()(const CaloHit *const left, const CaloHit *const right) const
{
    float distanceL, distanceR;

    if (left->GetHitType() == m_hitType_A)
    {
        distanceL = (left->GetPositionVector() - m_neutrinoVertex_A).GetMagnitudeSquared();
    }
    else 
    {
        distanceL = (left->GetPositionVector() - m_neutrinoVertex_B).GetMagnitudeSquared();
    }

    if (right->GetHitType() == m_hitType_A)
    {
        distanceR = (right->GetPositionVector() - m_neutrinoVertex_A).GetMagnitudeSquared();
    }
    else 
    {
        distanceR = (right->GetPositionVector() - m_neutrinoVertex_B).GetMagnitudeSquared();
    }

    return distanceL < distanceR;
}

} // namespace lar_content
