/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.cc
 *
 *  @brief  Implementation of the three dimensional hit creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/HitCreationBaseTool.h"
#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
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
    ClusterList twoDClusterList, threeDClusterList;
    LArPfoHelper::GetTwoDClusterList(pPfo, twoDClusterList);
    LArPfoHelper::GetThreeDClusterList(pPfo, threeDClusterList);

    for (ClusterList::const_iterator iter = twoDClusterList.begin(), iterEnd = twoDClusterList.end(); iter != iterEnd; ++iter)
    {
        if (TPC_3D == LArClusterHelper::GetClusterHitType(*iter))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        (*iter)->GetOrderedCaloHitList().GetCaloHitList(remainingHits);
    }

    for (ClusterList::const_iterator iter = threeDClusterList.begin(), iterEnd = threeDClusterList.end(); iter != iterEnd; ++iter)
    {
        CaloHitList localCaloHitList;
        (*iter)->GetOrderedCaloHitList().GetCaloHitList(localCaloHitList);

        for (CaloHitList::const_iterator hIter = localCaloHitList.begin(), hIterEnd = localCaloHitList.end(); hIter != hIterEnd; ++hIter)
        {
            const CaloHit *const pTargetCaloHit = static_cast<const CaloHit*>((*hIter)->GetParentCaloHitAddress());
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

void ThreeDHitCreationAlgorithm::AddThreeDHitsToPfo(const ParticleFlowObject *const pPfo, CaloHitList &caloHitList, const Cluster *&pCluster3D) const
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

void ThreeDHitCreationAlgorithm::CreateThreeDHit(const CaloHit *const pCaloHit2D, const CartesianVector &position3D, const CaloHit *&pCaloHit3D) const
{ 
    if (!this->CheckThreeDHit(position3D))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    PandoraContentApi::CaloHit::Parameters parameters;
    parameters.m_positionVector = position3D;
    parameters.m_hitType = TPC_3D;
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

bool ThreeDHitCreationAlgorithm::CheckThreeDHit(const CartesianVector &position3D) const
{
    // Check that corresponding pseudo layer is within range
    try
    {
        (void) PandoraContentApi::GetPlugins(*this)->GetPseudoLayerPlugin()->GetPseudoLayer(position3D);
    }
    catch (StatusCodeException &)
    {
        return false;
    }

    // TODO Check against detector geometry

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDHitCreationAlgorithm::Run()
{
    const PfoList *pPfoList = NULL;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ThreeDHitCreationAlgorithm: unable to find pfo list " << m_inputPfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    for (PfoList::const_iterator pIter = pPfoList->begin(), pIterEnd = pPfoList->end(); pIter != pIterEnd; ++pIter)
    {
        const ParticleFlowObject *const pPfo = *pIter;

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
                const Cluster *pCluster3D(NULL);
                this->AddThreeDHitsToPfo(pPfo, newThreeDHits, pCluster3D);
                allNewThreeDHits.insert(newThreeDHits.begin(), newThreeDHits.end());
            }
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, allNewThreeDHits, m_outputCaloHitListName));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *ThreeDHitCreationAlgorithm::GetThreeDCluster(const ParticleFlowObject *const pThreeDPfo) const
{
    const Cluster *pCluster3D(NULL);
    const ClusterList &clusterList(pThreeDPfo->GetClusterList());

    for (ClusterList::const_iterator iter = clusterList.begin(), iterEnd = clusterList.end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (TPC_3D != LArClusterHelper::GetClusterHitType(pCluster))
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

void ThreeDHitCreationAlgorithm::CreateThreeDCluster(const CaloHitList &caloHitList, const Cluster *&pCluster) const
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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "HitCreationTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        HitCreationBaseTool *const pHitCreationTool(dynamic_cast<HitCreationBaseTool*>(*iter));

        if (NULL == pHitCreationTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pHitCreationTool);
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListName", m_outputCaloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputClusterListName", m_outputClusterListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
