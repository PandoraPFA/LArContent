/**
 *  @file   LArContent/src/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional hit creation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArCalculators/LArTransformationCalculator.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArThreeDHelper.h"

#include "LArObjects/LArTwoDSlidingFitResult.h"

#include "LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ThreeDHitCreationAlgorithm::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    for (PfoList::const_iterator iter = pPfoList->begin(), iterEnd = pPfoList->end(); iter != iterEnd; ++iter)
    {
        try
        {
            ParticleFlowObject *pPfo = *iter;

            Cluster *pClusterU(NULL), *pClusterV(NULL), *pClusterW(NULL);
            this->GetClusters(pPfo, pClusterU, pClusterV, pClusterW);

            TwoDSlidingFitResult slidingFitResultU, slidingFitResultV, slidingFitResultW;
            LArClusterHelper::LArTwoDSlidingFit(pClusterU, m_slidingFitWindow, slidingFitResultU);
            LArClusterHelper::LArTwoDSlidingFit(pClusterV, m_slidingFitWindow, slidingFitResultV);
            LArClusterHelper::LArTwoDSlidingFit(pClusterW, m_slidingFitWindow, slidingFitResultW);

            CaloHitList caloHitListU, caloHitListV, caloHitListW;
            pClusterU->GetOrderedCaloHitList().GetCaloHitList(caloHitListU);
            pClusterV->GetOrderedCaloHitList().GetCaloHitList(caloHitListV);
            pClusterW->GetOrderedCaloHitList().GetCaloHitList(caloHitListW);

            CaloHitList newThreeDHits, ommittedTwoDHits;
            this->CreateThreeDHits(caloHitListU, slidingFitResultV, slidingFitResultW, newThreeDHits, ommittedTwoDHits);
            this->CreateThreeDHits(caloHitListV, slidingFitResultU, slidingFitResultW, newThreeDHits, ommittedTwoDHits);
            this->CreateThreeDHits(caloHitListW, slidingFitResultU, slidingFitResultV, newThreeDHits, ommittedTwoDHits);

            Cluster *pCluster3D(NULL);
            this->CreateThreeDCluster(newThreeDHits, pCluster3D);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo, pCluster3D));

            CaloHitList extrapolatedHits;
            this->CreateExtrapolatedHits(ommittedTwoDHits, pCluster3D, extrapolatedHits);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster3D, &extrapolatedHits));

            newThreeDHits.insert(extrapolatedHits.begin(), extrapolatedHits.end());
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, newThreeDHits, m_outputCaloHitListName));
        }
        catch (StatusCodeException &statusCodeException)
        {
            std::cout << "ThreeDHitCreationAlgorithm: Unable to create 3D hits for a given cluster " << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::GetClusters(const ParticleFlowObject *const pPfo, Cluster *&pClusterU, Cluster *&pClusterV, Cluster *&pClusterW) const
{
    pClusterU = NULL; pClusterV = NULL; pClusterW = NULL;
    const ClusterList &clusterList(pPfo->GetClusterList());

    if (3 != clusterList.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster(*iter);
        const HitType hitType(LArThreeDHelper::GetClusterHitType(pCluster));

        if ((TPC_VIEW_U != hitType) && (TPC_VIEW_V != hitType) && (TPC_VIEW_W != hitType))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        Cluster *&pTargetCluster((TPC_VIEW_U == hitType) ? pClusterU : (TPC_VIEW_V == hitType) ? pClusterV : pClusterW);
        pTargetCluster = pCluster;
    }

    if ((NULL == pClusterU) || (NULL == pClusterV) || (NULL == pClusterW))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::CreateThreeDHits(const CaloHitList &inputTwoDHits, const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
    CaloHitList &newThreeDHits, CaloHitList &ommittedTwoDHits) const
{
    for (CaloHitList::const_iterator iter = inputTwoDHits.begin(), iterEnd = inputTwoDHits.end(); iter != iterEnd; ++iter)
    {
        try
        {
            CaloHit *pCaloHit(*iter);

            CartesianVector position3D(0.f, 0.f, 0.f);
            float chiSquared(std::numeric_limits<float>::max());
            this->GetPosition3D(pCaloHit, fitResult1, fitResult2, position3D, chiSquared);

            if (chiSquared > m_chiSquaredCut)
                throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

            CaloHit *pCaloHit3D(NULL);
            PandoraContentApi::CaloHit::Parameters parameters;
            parameters.m_positionVector = position3D;
            parameters.m_hitType = TPC_3D;
            parameters.m_pParentAddress = static_cast<void*>(pCaloHit);

            // TODO Check these parameters, especially new cell dimensions
            parameters.m_cellThickness = pCaloHit->GetCellThickness();
            parameters.m_cellSizeU = pCaloHit->GetCellLengthScale();
            parameters.m_cellSizeV = pCaloHit->GetCellLengthScale();
            parameters.m_cellNormalVector = pCaloHit->GetCellNormalVector();
            parameters.m_expectedDirection = pCaloHit->GetExpectedDirection();
            parameters.m_nCellRadiationLengths = pCaloHit->GetNCellRadiationLengths();
            parameters.m_nCellInteractionLengths = pCaloHit->GetNCellInteractionLengths();
            parameters.m_time = pCaloHit->GetTime();
            parameters.m_inputEnergy = pCaloHit->GetInputEnergy();
            parameters.m_mipEquivalentEnergy = pCaloHit->GetMipEquivalentEnergy();
            parameters.m_electromagneticEnergy = pCaloHit->GetElectromagneticEnergy();
            parameters.m_hadronicEnergy = pCaloHit->GetHadronicEnergy();
            parameters.m_isDigital = pCaloHit->IsDigital();
            parameters.m_detectorRegion = pCaloHit->GetDetectorRegion();
            parameters.m_layer = pCaloHit->GetLayer();
            parameters.m_isInOuterSamplingLayer = pCaloHit->IsInOuterSamplingLayer();

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pCaloHit3D));
            newThreeDHits.insert(pCaloHit3D);
        }
        catch (StatusCodeException &)
        {
            ommittedTwoDHits.insert(*iter);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::GetPosition3D(const CaloHit *const pCaloHit, const TwoDSlidingFitResult &fitResult1, const TwoDSlidingFitResult &fitResult2,
    CartesianVector &position3D, float &chiSquared) const
{
    const CartesianVector &position(pCaloHit->GetPositionVector());
    CartesianVector fitPosition1(0.f, 0.f, 0.f), fitPosition2(0.f, 0.f, 0.f);
    fitResult1.GetGlobalFitPosition(position.GetX(), true, fitPosition1);
    fitResult2.GetGlobalFitPosition(position.GetX(), true, fitPosition2);

    const HitType hitType(pCaloHit->GetHitType());
    const HitType hitType1(LArThreeDHelper::GetClusterHitType(fitResult1.GetCluster()));
    const HitType hitType2(LArThreeDHelper::GetClusterHitType(fitResult2.GetCluster()));

    const float sigmaHit(LArGeometryHelper::GetLArTransformationCalculator()->GetSigmaUVW());
    const float sigmaFit(m_sigmaFitMultiplier * sigmaHit);

    if (m_useChiSquaredApproach)
    {
        const float u((TPC_VIEW_U == hitType) ? position.GetZ() : (TPC_VIEW_U == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());
        const float v((TPC_VIEW_V == hitType) ? position.GetZ() : (TPC_VIEW_V == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());
        const float w((TPC_VIEW_W == hitType) ? position.GetZ() : (TPC_VIEW_W == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());

        const float sigmaU((TPC_VIEW_U == hitType) ? sigmaHit : sigmaFit);
        const float sigmaV((TPC_VIEW_V == hitType) ? sigmaHit : sigmaFit);
        const float sigmaW((TPC_VIEW_W == hitType) ? sigmaHit : sigmaFit);

        float bestY(std::numeric_limits<float>::max()), bestZ(std::numeric_limits<float>::max());
        LArGeometryHelper::GetLArTransformationCalculator()->GetMinChiSquaredYZ(u, v, w, sigmaU, sigmaV, sigmaW, bestY, bestZ, chiSquared);
        position3D.SetValues(position.GetX(), bestY, bestZ);
    }
    else
    {
        const LArTransformationCalculator::PositionAndType hitPositionAndType(position.GetZ(), hitType);
        const LArTransformationCalculator::PositionAndType fitPositionAndType1(fitPosition1.GetZ(), hitType1);
        const LArTransformationCalculator::PositionAndType fitPositionAndType2(fitPosition2.GetZ(), hitType2);

        float bestY(std::numeric_limits<float>::max()), bestZ(std::numeric_limits<float>::max());
        LArGeometryHelper::GetLArTransformationCalculator()->GetProjectedYZ(hitPositionAndType, fitPositionAndType1, fitPositionAndType2, sigmaHit, sigmaFit, bestY, bestZ, chiSquared);
        position3D.SetValues(position.GetX(), bestY, bestZ);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::CreateThreeDCluster(const CaloHitList &caloHitList, Cluster *&pCluster) const
{
    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    const ClusterList *pClusterList = NULL; std::string clusterListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pClusterList, clusterListName));

    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList = caloHitList;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));

    if (!pClusterList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_outputClusterListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::CreateExtrapolatedHits(const CaloHitList &omittedHits, Cluster *pCluster, CaloHitList &extrapolatedHits) const
{
    // TODO
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDHitCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListName", m_outputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));

    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    m_useChiSquaredApproach = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseChiSquaredApproach", m_useChiSquaredApproach));

    m_sigmaFitMultiplier = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SigmaFitMultiplier", m_sigmaFitMultiplier));

    m_chiSquaredCut = 5.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ChiSquaredCut", m_chiSquaredCut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
