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

void ThreeDHitCreationAlgorithm::GetUnusedTwoDHits(const ParticleFlowObject *const pPfo, CaloHitList &caloHitList) const
{
    ClusterList threeDClusterList;
    const ClusterList &pfoClusterList(pPfo->GetClusterList());

    for (ClusterList::const_iterator iter = pfoClusterList.begin(), iterEnd = pfoClusterList.end(); iter != iterEnd; ++iter)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(*iter));

        if ((TPC_VIEW_U == hitType) || (TPC_VIEW_V == hitType) || (TPC_VIEW_W == hitType))
        {
            (*iter)->GetOrderedCaloHitList().GetCaloHitList(caloHitList);
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
            CaloHitList::iterator eraseIter = caloHitList.find(pTargetCaloHit);

            if (caloHitList.end() == eraseIter)
                throw StatusCodeException(STATUS_CODE_FAILURE);

            caloHitList.erase(eraseIter);
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

        CaloHitList allNewThreeDHits, inputTwoDHits;
        this->GetUnusedTwoDHits(pPfo, inputTwoDHits);

        for (HitCreationToolList::const_iterator tIter = m_algorithmToolList.begin(), tIterEnd = m_algorithmToolList.end(); tIter != tIterEnd; ++tIter)
        {
            CaloHitList newThreeDHits, omittedTwoDHits;
            (*tIter)->Run(this, pPfo, inputTwoDHits, newThreeDHits, omittedTwoDHits);

            if (!newThreeDHits.empty())
            {
                Cluster *pCluster3D(NULL);
                this->AddThreeDHitsToPfo(pPfo, newThreeDHits, pCluster3D);
                inputTwoDHits = omittedTwoDHits;
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

} // namespace lar
