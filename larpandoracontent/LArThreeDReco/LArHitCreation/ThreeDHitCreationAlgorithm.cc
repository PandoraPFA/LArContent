/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.cc
 *
 *  @brief  Implementation of the three dimensional hit creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoracontent/LArPlugins/LArTransformationPlugin.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

#include <algorithm>

using namespace pandora;

namespace lar_content
{

void ThreeDHitCreationAlgorithm::FilterCaloHitsByType(const CaloHitVector &inputCaloHitVector, const HitType hitType, CaloHitVector &outputCaloHitVector) const
{
    for (const CaloHit *const pCaloHit : inputCaloHitVector)
    {
        if (hitType == pCaloHit->GetHitType())
            outputCaloHitVector.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDHitCreationAlgorithm::Run()
{
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ThreeDHitCreationAlgorithm: unable to find pfo list " << m_inputPfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    CaloHitList allNewThreeDHits;

    PfoVector pfoVector(pPfoList->begin(), pPfoList->end());
    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pPfo : pfoVector)
    {
        ProtoHitVector protoHitVector;

        for (HitCreationBaseTool *const pHitCreationTool : m_algorithmToolVector)
        {
            CaloHitVector remainingTwoDHits;
            this->SeparateTwoDHits(pPfo, protoHitVector, remainingTwoDHits);

            if (remainingTwoDHits.empty())
                break;

            pHitCreationTool->Run(this, pPfo, remainingTwoDHits, protoHitVector);
        }

        if (LArPfoHelper::IsTrack(pPfo))
            this->IterativeTreatment(protoHitVector);

        CaloHitList newThreeDHits;
        this->CreateThreeDHits(protoHitVector, newThreeDHits);
        this->AddThreeDHitsToPfo(pPfo, newThreeDHits);

        allNewThreeDHits.insert(allNewThreeDHits.end(), newThreeDHits.begin(), newThreeDHits.end());
    }

    if (!allNewThreeDHits.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, allNewThreeDHits, m_outputCaloHitListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::SeparateTwoDHits(const ParticleFlowObject *const pPfo, const ProtoHitVector &protoHitVector, CaloHitVector &remainingHitVector) const
{
    ClusterList twoDClusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);
    CaloHitList remainingHitList;

    for (const Cluster *const pCluster : twoDClusterList)
    {
        if (TPC_3D == LArClusterHelper::GetClusterHitType(pCluster))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        pCluster->GetOrderedCaloHitList().FillCaloHitList(remainingHitList);
    }

    CaloHitSet remainingHitSet(remainingHitList.begin(), remainingHitList.end());

    for (const ProtoHit &protoHit : protoHitVector)
    {
        CaloHitSet::iterator eraseIter = remainingHitSet.find(protoHit.GetParentCaloHit2D());

        if (remainingHitSet.end() == eraseIter)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        remainingHitSet.erase(eraseIter);
    }

    remainingHitVector.insert(remainingHitVector.end(), remainingHitSet.begin(), remainingHitSet.end());
    std::sort(remainingHitVector.begin(), remainingHitVector.end(), LArClusterHelper::SortHitsByPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::IterativeTreatment(ProtoHitVector &protoHitVector) const
{
    const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const float sigmaUVW(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetSigmaUVW());
    const float sigmaFit(sigmaUVW); // TODO
    const float sigmaHit(sigmaUVW); // TODO
    const float sigma3DFit(sigmaUVW / 5.f); // TODO
    const unsigned int layerWindow(10); // m_slidingFitHalfWindow TODO

    double originalChi2(0.);
    CartesianPointVector currentPoints3D;
    const ProtoHitVector originalProtoHitVector(protoHitVector);

    for (const ProtoHit &protoHit : protoHitVector)
    {
        originalChi2 += protoHit.GetChi2();
        currentPoints3D.push_back(protoHit.GetPosition3D());
    }

    try
    {
        double originalChi2Fit(0.);
        const ThreeDSlidingFitResult slidingFitResult(&currentPoints3D, layerWindow, layerPitch);

        for (const ProtoHit &protoHit : protoHitVector)
        {
            CartesianVector pointOnFit(0.f, 0.f, 0.f);
            const double rL(slidingFitResult.GetLongitudinalDisplacement(protoHit.GetPosition3D()));

            if (STATUS_CODE_SUCCESS == slidingFitResult.GetGlobalFitPosition(rL, pointOnFit))
            {
                const double uFit(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->YZtoU(pointOnFit.GetY(), pointOnFit.GetZ()));
                const double vFit(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->YZtoV(pointOnFit.GetY(), pointOnFit.GetZ()));
                const double wFit(pointOnFit.GetZ());

                const double outputU(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->YZtoU(protoHit.GetPosition3D().GetY(), protoHit.GetPosition3D().GetZ()));
                const double outputV(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->YZtoV(protoHit.GetPosition3D().GetY(), protoHit.GetPosition3D().GetZ()));
                const double outputW(protoHit.GetPosition3D().GetZ());

                const double deltaUFit(uFit - outputU), deltaVFit(vFit - outputV), deltaWFit(wFit - outputW);
                originalChi2Fit += ((deltaUFit * deltaUFit) / (sigma3DFit * sigma3DFit)) + ((deltaVFit * deltaVFit) / (sigma3DFit * sigma3DFit)) + ((deltaWFit * deltaWFit) / (sigma3DFit * sigma3DFit));
            }
        }

        originalChi2 += originalChi2Fit;
        unsigned int nIterations(0);
std::cout << " originalChi2 " << originalChi2 << ", originalChi2Fit " << originalChi2Fit << ", nOriginalPoints3D " << currentPoints3D.size() << std::endl;
for (const CartesianVector &point : currentPoints3D) PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "OriginalPoint", BLUE, 1);
PandoraMonitoringApi::ViewEvent(this->GetPandora());
        while (nIterations++ < 10) // TODO
        {
            double newChi2(0.);
            CartesianPointVector newPoints3D;
            ProtoHitVector newProtoHitVector(protoHitVector);

            const ThreeDSlidingFitResult newSlidingFitResult(&currentPoints3D, layerWindow, layerPitch);

            for (ProtoHit &protoHit : newProtoHitVector)
            {
                CartesianVector pointOnFit(0.f, 0.f, 0.f);
                const double rL(newSlidingFitResult.GetLongitudinalDisplacement(protoHit.GetPosition3D()));

                if (STATUS_CODE_SUCCESS == newSlidingFitResult.GetGlobalFitPosition(rL, pointOnFit))
                {
PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &pointOnFit, "FitPoint", ORANGE, 1);
                    const CaloHit *const pCaloHit2D(protoHit.GetParentCaloHit2D());
                    const HitType hitType(pCaloHit2D->GetHitType());

                    const double uFit(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->YZtoU(pointOnFit.GetY(), pointOnFit.GetZ()));
                    const double vFit(LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->YZtoV(pointOnFit.GetY(), pointOnFit.GetZ()));
                    const double wFit(pointOnFit.GetZ());

                    const double sigmaU((TPC_VIEW_U == hitType) ? sigmaHit : sigmaFit);
                    const double sigmaV((TPC_VIEW_V == hitType) ? sigmaHit : sigmaFit);
                    const double sigmaW((TPC_VIEW_W == hitType) ? sigmaHit : sigmaFit);

                    CartesianVector position3D(0.f, 0.f, 0.f);
                    double chi2(std::numeric_limits<double>::max());
                    double u(std::numeric_limits<double>::max()), v(std::numeric_limits<double>::max()), w(std::numeric_limits<double>::max());

                    if (protoHit.GetNTrajectorySamples() == 2)
                    {
                        u = (TPC_VIEW_U == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_U == protoHit.GetFirstTrajectorySample().GetHitType()) ? protoHit.GetFirstTrajectorySample().GetPosition().GetZ() : protoHit.GetLastTrajectorySample().GetPosition().GetZ();
                        v = (TPC_VIEW_V == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_V == protoHit.GetFirstTrajectorySample().GetHitType()) ? protoHit.GetFirstTrajectorySample().GetPosition().GetZ() : protoHit.GetLastTrajectorySample().GetPosition().GetZ();
                        w = (TPC_VIEW_W == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_W == protoHit.GetFirstTrajectorySample().GetHitType()) ? protoHit.GetFirstTrajectorySample().GetPosition().GetZ() : protoHit.GetLastTrajectorySample().GetPosition().GetZ();
                    }
                    else if (protoHit.GetNTrajectorySamples() == 1)
                    {
                        u = LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->YZtoU(protoHit.GetPosition3D().GetY(), protoHit.GetPosition3D().GetZ());
                        v = LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->YZtoV(protoHit.GetPosition3D().GetY(), protoHit.GetPosition3D().GetZ());
                        w = protoHit.GetPosition3D().GetZ();
                    }
                    else
                    {
                        std::cout << "ThreeDHitCreationAlgorithm::IterativeTreatment - Unexpected number of trajectory samples" << std::endl;
                        throw StatusCodeException(STATUS_CODE_FAILURE);
                    }

                    double bestY(std::numeric_limits<double>::max()), bestZ(std::numeric_limits<double>::max());
                    LArGeometryHelper::GetLArTransformationPlugin(this->GetPandora())->GetMinChiSquaredYZ(u, v, w, sigmaU, sigmaV, sigmaW, uFit, vFit, wFit, sigma3DFit, bestY, bestZ, chi2);
                    position3D.SetValues(pCaloHit2D->GetPositionVector().GetX(), static_cast<float>(bestY), static_cast<float>(bestZ));

                    protoHit.SetPosition3D(position3D, chi2);
                }

                newChi2 += protoHit.GetChi2();
                newPoints3D.push_back(protoHit.GetPosition3D());
            }

std::cout << " Iteration " << nIterations << " newChi2 " << newChi2 << ", nCurrentPoints3D " << newProtoHitVector.size() << ", prevChi2 " << originalChi2 << ", prevNHits " << protoHitVector.size() << std::endl;
for (const CartesianVector &point : newPoints3D) PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "NewPoint", RED, 1);
PandoraMonitoringApi::ViewEvent(this->GetPandora());
            if (newChi2 > 1. * originalChi2) // TODO
            {
                std::cout << "Iteration terminated " << std::endl;
                break;
            }

            originalChi2 = newChi2;
            protoHitVector = newProtoHitVector;
            currentPoints3D = newPoints3D;
        }
    }
    catch (const StatusCodeException &)
    {
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::CreateThreeDHits(const ProtoHitVector &protoHitVector, CaloHitList &newThreeDHits) const
{
    for (const ProtoHit &protoHit : protoHitVector)
    {
        const CaloHit *pCaloHit3D(nullptr);
        this->CreateThreeDHit(protoHit, pCaloHit3D);

        if (!pCaloHit3D)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        newThreeDHits.push_back(pCaloHit3D);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::CreateThreeDHit(const ProtoHit &protoHit, const CaloHit *&pCaloHit3D) const
{ 
    if (!this->CheckThreeDHit(protoHit))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    PandoraContentApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = protoHit.GetPosition3D();
    parameters.m_hitType = TPC_3D;

    const CaloHit *const pCaloHit2D(protoHit.GetParentCaloHit2D());
    parameters.m_pParentAddress = static_cast<const void*>(pCaloHit2D);

    // TODO Check these parameters, especially new cell dimensions
    parameters.m_cellThickness = pCaloHit2D->GetCellThickness();
    parameters.m_cellGeometry = RECTANGULAR;
    parameters.m_cellSize0 = pCaloHit2D->GetCellLengthScale();
    parameters.m_cellSize1 = pCaloHit2D->GetCellLengthScale();
    parameters.m_cellNormalVector = pCaloHit2D->GetCellNormalVector();
    parameters.m_expectedDirection = pCaloHit2D->GetExpectedDirection();
    parameters.m_nCellRadiationLengths = pCaloHit2D->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pCaloHit2D->GetNCellInteractionLengths();
    parameters.m_time = pCaloHit2D->GetTime();
    parameters.m_inputEnergy = pCaloHit2D->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pCaloHit2D->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pCaloHit2D->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pCaloHit2D->GetHadronicEnergy();
    parameters.m_isDigital = pCaloHit2D->IsDigital();
    parameters.m_hitRegion = pCaloHit2D->GetHitRegion();
    parameters.m_layer = pCaloHit2D->GetLayer();
    parameters.m_isInOuterSamplingLayer = pCaloHit2D->IsInOuterSamplingLayer();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pCaloHit3D));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDHitCreationAlgorithm::CheckThreeDHit(const ProtoHit &protoHit) const
{
    try
    {
        // Check that corresponding pseudo layer is within range
        (void) PandoraContentApi::GetPlugins(*this)->GetPseudoLayerPlugin()->GetPseudoLayer(protoHit.GetPosition3D());
    }
    catch (StatusCodeException &)
    {
        return false;
    }

    // TODO Check against detector geometry
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::AddThreeDHitsToPfo(const ParticleFlowObject *const pPfo, const CaloHitList &caloHitList) const
{
    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    ClusterList threeDClusterList;
    LArPfoHelper::GetThreeDClusterList(pPfo, threeDClusterList);

    if (!threeDClusterList.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const ClusterList *pClusterList(nullptr); std::string clusterListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListName));

    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList.insert(parameters.m_caloHitList.end(), caloHitList.begin(), caloHitList.end());

    const Cluster *pCluster3D(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster3D));

    if (!pCluster3D || !pClusterList || pClusterList->empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_outputClusterListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo, pCluster3D));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

const CartesianVector &ThreeDHitCreationAlgorithm::ProtoHit::GetPosition3D() const
{
    if (!m_isPositionSet)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return m_position3D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

double ThreeDHitCreationAlgorithm::ProtoHit::GetChi2() const
{
    if (!m_isPositionSet)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return m_chi2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ThreeDHitCreationAlgorithm::TrajectorySample &ThreeDHitCreationAlgorithm::ProtoHit::GetFirstTrajectorySample() const
{
    if (m_trajectorySampleVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return m_trajectorySampleVector.front();
}

//------------------------------------------------------------------------------------------------------------------------------------------

const ThreeDHitCreationAlgorithm::TrajectorySample &ThreeDHitCreationAlgorithm::ProtoHit::GetLastTrajectorySample() const
{
    if (m_trajectorySampleVector.size() < 2)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    return m_trajectorySampleVector.back();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDHitCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "HitCreationTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        HitCreationBaseTool *const pHitCreationTool(dynamic_cast<HitCreationBaseTool*>(*iter));

        if (!pHitCreationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pHitCreationTool);
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListName", m_outputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
