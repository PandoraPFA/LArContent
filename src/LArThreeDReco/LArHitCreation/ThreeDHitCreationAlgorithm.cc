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

#include "LArObjects/LArTwoDSlidingFitResult.h"

#include "LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

using namespace pandora;

namespace lar
{

void ThreeDHitCreationAlgorithm::GetRemainingTwoDHits(const ParticleFlowObject *const pPfo, CaloHitList &remainingHits) const
{
    CaloHitList usedHits;
    this->SeparateTwoDHits(pPfo, usedHits, remainingHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::GetUsedTwoDHits(const ParticleFlowObject *const pPfo, CaloHitList &usedHits) const
{
    CaloHitList remainingHits;
    this->SeparateTwoDHits(pPfo, usedHits, remainingHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::SeparateTwoDHits(const ParticleFlowObject *const pPfo, CaloHitList &usedHits, CaloHitList &remainingHits) const
{
    ClusterList threeDClusterList;
    const ClusterList &pfoClusterList(pPfo->GetClusterList());

    for (ClusterList::const_iterator iter = pfoClusterList.begin(), iterEnd = pfoClusterList.end(); iter != iterEnd; ++iter)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(*iter));

        if ((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType))
        {
            (*iter)->GetOrderedCaloHitList().GetCaloHitList(remainingHits);
        }
        else if (TPC_3D == hitType)
        {
            threeDClusterList.insert(*iter);
        }
        else
        {
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
        }
    }

    for (ClusterList::const_iterator iter = threeDClusterList.begin(), iterEnd = threeDClusterList.end(); iter != iterEnd; ++iter)
    {
        CaloHitList localCaloHitList;
        (*iter)->GetOrderedCaloHitList().GetCaloHitList(localCaloHitList);

        for (CaloHitList::const_iterator hIter = localCaloHitList.begin(), hIterEnd = localCaloHitList.end(); hIter != hIterEnd; ++hIter)
        {
            CaloHit *pTargetCaloHit = static_cast<CaloHit*>(const_cast<void*>((*hIter)->GetParentCaloHitAddress()));
            CaloHitList::iterator eraseIter = remainingHits.find(pTargetCaloHit);

            if (remainingHits.end() == eraseIter)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            usedHits.insert(pTargetCaloHit);
            remainingHits.erase(eraseIter);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::FilterCaloHitsByType(const CaloHitList &inputCaloHitList, const HitType hitType, CaloHitList &outputCaloHitList) const
{
    for (CaloHitList::const_iterator iter = inputCaloHitList.begin(), iterEnd = inputCaloHitList.end(); iter != iterEnd; ++iter)
    {
        if (hitType == (*iter)->GetHitType())
            outputCaloHitList.insert(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::AddThreeDHitsToPfo(ParticleFlowObject *const pPfo, CaloHitList &caloHitList, Cluster *&pCluster3D) const
{
    try
    {
        pCluster3D = this->GetThreeDCluster(pPfo);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pCluster3D, &caloHitList));
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_NOT_FOUND != statusCodeException.GetStatusCode())
            throw statusCodeException;

        this->CreateThreeDCluster(caloHitList, pCluster3D);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfo, pCluster3D));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDHitCreationAlgorithm::CreateThreeDHit(CaloHit *pCaloHit2D, const CartesianVector &position3D, CaloHit *&pCaloHit3D) const
{
    PandoraContentApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = position3D;
    parameters.m_hitType = TPC_3D;
    parameters.m_pParentAddress = static_cast<void*>(pCaloHit2D);

    // TODO Check these parameters, especially new cell dimensions
    parameters.m_cellThickness = pCaloHit2D->GetCellThickness();
    parameters.m_cellSizeU = pCaloHit2D->GetCellLengthScale();
    parameters.m_cellSizeV = pCaloHit2D->GetCellLengthScale();
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
    parameters.m_detectorRegion = pCaloHit2D->GetDetectorRegion();
    parameters.m_layer = pCaloHit2D->GetLayer();
    parameters.m_isInOuterSamplingLayer = pCaloHit2D->IsInOuterSamplingLayer();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pCaloHit3D));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDHitCreationAlgorithm::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    if (NULL == pPfoList)
        return STATUS_CODE_SUCCESS;

    for (PfoList::const_iterator pIter = pPfoList->begin(), pIterEnd = pPfoList->end(); pIter != pIterEnd; ++pIter)
    {
        ParticleFlowObject *pPfo = *pIter;

        CaloHitList allNewThreeDHits;

        for (HitCreationToolList::const_iterator tIter = m_algorithmToolList.begin(), tIterEnd = m_algorithmToolList.end(); tIter != tIterEnd; ++tIter)
        {
            CaloHitList remainingTwoDHits, newThreeDHits;
            this->GetRemainingTwoDHits(pPfo, remainingTwoDHits);

            if (remainingTwoDHits.empty())
                break;

            (*tIter)->Run(this, pPfo, remainingTwoDHits, newThreeDHits);

            if (!newThreeDHits.empty())
            {
                Cluster *pCluster3D(NULL);
                this->AddThreeDHitsToPfo(pPfo, newThreeDHits, pCluster3D);
                allNewThreeDHits.insert(newThreeDHits.begin(), newThreeDHits.end());
            }
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, allNewThreeDHits, m_outputCaloHitListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

Cluster *ThreeDHitCreationAlgorithm::GetThreeDCluster(ParticleFlowObject *const pPfo) const
{
    Cluster *pCluster3D(NULL);
    const ClusterList &clusterList(pPfo->GetClusterList());

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if (TPC_3D != hitType)
            continue;

        if (NULL != pCluster3D)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pCluster3D = pCluster;
    }

    if (NULL == pCluster3D)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return pCluster3D;
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

    if (pClusterList->empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, m_outputClusterListName));
}


//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void HitCreationTool::GetBestPosition3D(const CaloHit *const pCaloHit2D, const HitType hitType1, const HitType hitType2,
    const CartesianPointList &fitPositionList1, const CartesianPointList &fitPositionList2, CartesianVector &position3D, float &chiSquared) const
{
    if (fitPositionList1.empty() && fitPositionList2.empty())
    {
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }
    else if (fitPositionList1.empty())
    {
        if (fitPositionList2.size() != 1)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        CartesianPointList::const_iterator iter2 = fitPositionList2.begin();
        if (fitPositionList2.end() == iter2)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const CartesianVector &fitPosition2 = *iter2;
        this->GetPosition3D(pCaloHit2D, hitType2, fitPosition2, position3D, chiSquared);
    }
    else if (fitPositionList2.empty())
    {
        if (fitPositionList1.size() != 1)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

        CartesianPointList::const_iterator iter1 = fitPositionList1.begin();
        if (fitPositionList1.end() == iter1)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        const CartesianVector &fitPosition1 = *iter1;
        this->GetPosition3D(pCaloHit2D, hitType1, fitPosition1, position3D, chiSquared);
    }
    else
    {
        chiSquared = std::numeric_limits<float>::max();

        for (CartesianPointList::const_iterator iter1 = fitPositionList1.begin(), iterEnd1 = fitPositionList1.end(); iter1 != iterEnd1; ++iter1)
        {
            const CartesianVector &fitPosition1 = *iter1;
            for (CartesianPointList::const_iterator iter2 = fitPositionList2.begin(), iterEnd2 = fitPositionList2.end(); iter2 != iterEnd2; ++iter2)
            {
                const CartesianVector &fitPosition2 = *iter2;

                CartesianVector thisPosition3D(0.f, 0.f, 0.f);
                float thisChiSquared(std::numeric_limits<float>::max());
                this->GetPosition3D(pCaloHit2D, hitType1, hitType2, fitPosition1, fitPosition2, thisPosition3D, thisChiSquared);

                if (thisChiSquared < chiSquared)
                {
                    chiSquared = thisChiSquared;
                    position3D = thisPosition3D;
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitCreationTool::GetPosition3D(const CaloHit *const pCaloHit2D, const HitType hitType1, const HitType hitType2,
    const CartesianVector &fitPosition1, const CartesianVector &fitPosition2, CartesianVector &position3D, float &chiSquared) const
{
    const float sigmaHit(LArGeometryHelper::GetLArTransformationCalculator()->GetSigmaUVW());
    const float sigmaFit(sigmaHit); // TODO: Input uncertainties into this method
    const HitType hitType(pCaloHit2D->GetHitType());

    if (m_useChiSquaredApproach)
    {
        const float u((TPC_VIEW_U == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_U == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());
        const float v((TPC_VIEW_V == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_V == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());
        const float w((TPC_VIEW_W == hitType) ? pCaloHit2D->GetPositionVector().GetZ() : (TPC_VIEW_W == hitType1) ? fitPosition1.GetZ() : fitPosition2.GetZ());

        const float sigmaU((TPC_VIEW_U == hitType) ? sigmaHit : sigmaFit);
        const float sigmaV((TPC_VIEW_V == hitType) ? sigmaHit : sigmaFit);
        const float sigmaW((TPC_VIEW_W == hitType) ? sigmaHit : sigmaFit);

        float bestY(std::numeric_limits<float>::max()), bestZ(std::numeric_limits<float>::max());
        LArGeometryHelper::GetLArTransformationCalculator()->GetMinChiSquaredYZ(u, v, w, sigmaU, sigmaV, sigmaW, bestY, bestZ, chiSquared);
        position3D.SetValues(pCaloHit2D->GetPositionVector().GetX(), bestY, bestZ);
    }
    else
    {
        const LArTransformationCalculator::PositionAndType hitPositionAndType(pCaloHit2D->GetPositionVector().GetZ(), hitType);
        const LArTransformationCalculator::PositionAndType fitPositionAndType1(fitPosition1.GetZ(), hitType1);
        const LArTransformationCalculator::PositionAndType fitPositionAndType2(fitPosition2.GetZ(), hitType2);

        float bestY(std::numeric_limits<float>::max()), bestZ(std::numeric_limits<float>::max());
        LArGeometryHelper::GetLArTransformationCalculator()->GetProjectedYZ(hitPositionAndType, fitPositionAndType1, fitPositionAndType2, sigmaHit, sigmaFit, bestY, bestZ, chiSquared);
        position3D.SetValues(pCaloHit2D->GetPositionVector().GetX(), bestY, bestZ);
    }

    if (m_useDeltaXCorrection)
    {
        const float deltaX1(pCaloHit2D->GetPositionVector().GetX() - fitPosition1.GetX());
        const float deltaX2(pCaloHit2D->GetPositionVector().GetX() - fitPosition2.GetX());
        chiSquared += ((deltaX1 * deltaX1) / (m_sigmaX * m_sigmaX)) + ((deltaX2 * deltaX2) / (m_sigmaX * m_sigmaX));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitCreationTool::GetPosition3D(const CaloHit *const pCaloHit2D, const HitType hitType, const CartesianVector &fitPosition,
    CartesianVector &position3D, float &chiSquared) const
{
    if (pCaloHit2D->GetHitType() == hitType)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    LArGeometryHelper::MergeTwoPositions3D(pCaloHit2D->GetHitType(), hitType, pCaloHit2D->GetPositionVector(), fitPosition,
        position3D, chiSquared);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDHitCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "HitCreationTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        HitCreationTool *pHitCreationTool(dynamic_cast<HitCreationTool*>(*iter));

        if (NULL == pHitCreationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pHitCreationTool);
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListName", m_outputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitCreationTool::ReadSettings(const pandora::TiXmlHandle xmlHandle)
{
    m_useChiSquaredApproach = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseChiSquaredApproach", m_useChiSquaredApproach));

    m_useDeltaXCorrection = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseDeltaXCorrection", m_useDeltaXCorrection));

    m_sigmaX = 1.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SigmaX", m_sigmaX));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
