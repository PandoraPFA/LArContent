/**
 *  @file   larpandoracontent/LArShowerRefinement/ElectronInitialRegionRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the electron initial region refinement algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/PeakDirectionFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ProtoShowerMatchingTool.h"
#include "larpandoracontent/LArShowerRefinement/ShowerSpineFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ElectronInitialRegionRefinementAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ElectronInitialRegionRefinementAlgorithm::ElectronInitialRegionRefinementAlgorithm() :
    m_showerSlidingFitWindow(1000),
    m_maxCoincideneTransverseSeparation(5.f),
    m_minSpinePurity(0.7f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronInitialRegionRefinementAlgorithm::Run()
{
    PfoVector showerPfoVector;
    this->FillShowerPfoVector(showerPfoVector);

    if (showerPfoVector.empty())
        return STATUS_CODE_SUCCESS;

    for (const ParticleFlowObject *const pShowerPfo : showerPfoVector)
    {
        // Only consider significant showers
        CaloHitList caloHits3D;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHits3D);

        if (caloHits3D.size() < 50)
            continue;

        this->RefineShower(pShowerPfo);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::FillShowerPfoVector(PfoVector &showerPfoVector)
{
    const PfoList *pPfoList(nullptr);

    if (PandoraContentApi::GetList(*this, m_showerPfoListName, pPfoList) != STATUS_CODE_SUCCESS)
        return;

    if (!pPfoList || pPfoList->empty())
    {
        std::cout << "ElectronInitialRegionRefinementAlgorithm: unable to find shower pfo list " << m_showerPfoListName << std::endl;
        return;
    }

    showerPfoVector.insert(showerPfoVector.begin(), pPfoList->begin(), pPfoList->end());

    std::sort(showerPfoVector.begin(), showerPfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::RefineShower(const ParticleFlowObject *const pShowerPfo)
{
    CartesianVector nuVertex3D(0.f, 0.f, 0.f);

    if (this->GetNeutrinoVertex(nuVertex3D) != STATUS_CODE_SUCCESS)
        return;

    ElectronProtoShowerVector protoShowerVectorU, protoShowerVectorV, protoShowerVectorW;

    std::cout << "Building U protoShowers" << std::endl;
    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_U, protoShowerVectorU);
    std::cout << "Building V protoShowers" << std::endl;
    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_V, protoShowerVectorV);
    std::cout << "Building W protoShowers" << std::endl;
    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_W, protoShowerVectorW);

    ProtoShowerMatchVector protoShowerMatchVector;
    m_pProtoShowerMatchingTool->Run(this, protoShowerVectorU, protoShowerVectorV, protoShowerVectorW, protoShowerMatchVector);

    for (const ProtoShowerMatch &protoShowerMatch : protoShowerMatchVector)
    {
    // Add in ambiguous hits



    }



    // Work out if it an electron
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronInitialRegionRefinementAlgorithm::GetNeutrinoVertex(CartesianVector &nuVertex3D)
{
    const VertexList *pNuVertexList(nullptr);
    const StatusCode statusCode(PandoraContentApi::GetList(*this, m_neutrinoVertexListName, pNuVertexList));

    if (statusCode != STATUS_CODE_SUCCESS)
        return statusCode;

    if (!pNuVertexList || (pNuVertexList->size() != 1))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ElectronInitialRegionRefinementAlgorithm: unable to find vertex list " << m_neutrinoVertexListName << " if it does exist, it may have more than one nu vertex" << std::endl;

        return STATUS_CODE_NOT_INITIALIZED;
    }

    nuVertex3D = pNuVertexList->front()->GetPosition();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::BuildViewProtoShowers(const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertex3D, 
   HitType hitType, ElectronProtoShowerVector &protoShowerVector)
{
    const CaloHitList *pViewHitList(nullptr);

    if (this->GetHitListOfType(hitType, pViewHitList) != STATUS_CODE_SUCCESS)
        return;

    /*
    // Only consider significant showers
    CaloHitList caloHits3D;
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHits3D);

    if (caloHits3D.size() < 50)
        return false;
    */

    CartesianVector showerVertexPosition(0.f, 0.f, 0.f);
    try
    {
        showerVertexPosition = this->GetShowerVertex(pShowerPfo, hitType, nuVertex3D);
    }
    catch(...)
    {
        return;
    }

    CartesianPointVector peakDirectionVector;
    if (m_pPeakDirectionFinderTool->Run(pShowerPfo, nuVertex3D, pViewHitList, hitType, peakDirectionVector) != STATUS_CODE_SUCCESS)
        return;

    CaloHitList unavailableHitList;
    for (CartesianVector &peakDirection : peakDirectionVector)
    {
        CaloHitList showerSpineHitList;
        if (m_pShowerSpineFinderTool->Run(nuVertex3D, pViewHitList, hitType, peakDirection, unavailableHitList, showerSpineHitList) != STATUS_CODE_SUCCESS)
            continue;

        this->RefineShowerVertex(pShowerPfo, hitType, nuVertex3D, peakDirection, showerVertexPosition);

        // If the spine passes the shower vertex, does it live inside the shower?
        if (!this->IsSpineCoincident(pShowerPfo, nuVertex3D, hitType, showerVertexPosition, showerSpineHitList))
            continue;

        CartesianVector showerStartPosition(0.f, 0.f, 0.f);
        CartesianVector showerStartDirection(0.f, 0.f, 0.f);

        if (m_pShowerStartFinderTool->Run(pShowerPfo, peakDirection, hitType, showerSpineHitList, showerStartPosition, showerStartDirection) != STATUS_CODE_SUCCESS)
            continue;
        
        std::cout << "ISOBEL - you need to get clean up the protoShower class (and add everything else to it...)" << std::endl;
        // hits to add are at the bottom of the code

        CaloHitList viewShowerHitList;
        LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

        const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

        ElectronProtoShower protoShower(ShowerCore(showerStartPosition, showerStartDirection), 
            ConnectionPathway(nuVertex2D, peakDirection), 
            showerSpineHitList, false, CaloHitList(), CartesianPointVector(), CaloHitList());

        protoShowerVector.push_back(protoShower);

        unavailableHitList.insert(unavailableHitList.begin(), showerSpineHitList.begin(), showerSpineHitList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronInitialRegionRefinementAlgorithm::GetHitListOfType(const HitType hitType, const CaloHitList *&pCaloHitList)
{
    const std::string typeHitListName(hitType == TPC_VIEW_U ? m_caloHitListNameU : hitType == TPC_VIEW_V ? m_caloHitListNameV : m_caloHitListNameW);

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, typeHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ShowerStartRefinementBaseTool: unable to find calo hit list " << typeHitListName << std::endl;

        return STATUS_CODE_NOT_INITIALIZED;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector ElectronInitialRegionRefinementAlgorithm::GetShowerVertex(const ParticleFlowObject *const pShowerPfo, const HitType hitType, 
    const CartesianVector &nuVertex3D) const
{
    ClusterList viewCusterList;
    LArPfoHelper::GetClusters(pShowerPfo, hitType, viewCusterList);

    if (viewCusterList.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const TwoDSlidingShowerFitResult twoDShowerSlidingFit(viewCusterList.front(), m_showerSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const CartesianVector &minLayerPosition(twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMinLayerPosition());
    const CartesianVector &maxLayerPosition(twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMaxLayerPosition());
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    if ((nuVertex2D.GetZ() > minLayerPosition.GetZ()) && (nuVertex2D.GetZ() < maxLayerPosition.GetZ()))
        return nuVertex2D;

    const float minSeparation((nuVertex2D - minLayerPosition).GetMagnitudeSquared());
    const float maxSeparation((nuVertex2D - maxLayerPosition).GetMagnitudeSquared());
    CartesianVector showerVertexPosition(minSeparation < maxSeparation ? minLayerPosition : maxLayerPosition);

    return (minSeparation < maxSeparation ? minLayerPosition : maxLayerPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronInitialRegionRefinementAlgorithm::RefineShowerVertex(const ParticleFlowObject *const pShowerPfo, const HitType hitType, 
    const CartesianVector &nuVertex3D, const CartesianVector &peakDirection, CartesianVector &showerVertexPosition) const
{
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    if (!this->IsShowerConnected(showerVertexPosition, nuVertex2D, peakDirection))
    {
        CaloHitList viewShowerHitList;
        LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

        float minL(std::numeric_limits<float>::max());

        for (const CaloHit *const pCaloHit : viewShowerHitList)
        {
            const CartesianVector displacement(pCaloHit->GetPositionVector() - nuVertex2D);
            const float longitudinalSeparation(peakDirection.GetDotProduct(displacement));

            if ((longitudinalSeparation < (-1.f)) || (longitudinalSeparation > minL))
                continue;

            const float transverseSeparation((peakDirection.GetCrossProduct(displacement)).GetMagnitude());

            if (transverseSeparation < m_maxCoincideneTransverseSeparation)
            {
                showerVertexPosition = pCaloHit->GetPositionVector();
                minL = longitudinalSeparation;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronInitialRegionRefinementAlgorithm::IsSpineCoincident(const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertex3D, 
    const HitType hitType, const CartesianVector &showerVertex, const CaloHitList &showerSpineHitList) const
{
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    CaloHitList viewShowerHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

    CaloHitList postShowerVertexSpineHits;
    const float showerDistanceFromNuSquared((showerVertex - nuVertex2D).GetMagnitudeSquared());

    for (const CaloHit * const pSpineHit : showerSpineHitList)
    {
        const CartesianVector &hitPosition(pSpineHit->GetPositionVector());
        const float separationSquared((hitPosition - nuVertex2D).GetMagnitudeSquared());

        if (separationSquared > showerDistanceFromNuSquared)
            postShowerVertexSpineHits.push_back(pSpineHit);
    }

    if (postShowerVertexSpineHits.size() == 0)
        return true;

    // Check whether shower spine is pure
    const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(postShowerVertexSpineHits, viewShowerHitList));
    const float spinePurity(static_cast<float>(sharedHitList.size()) / static_cast<float>(postShowerVertexSpineHits.size()));

    if (spinePurity < m_minSpinePurity)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronInitialRegionRefinementAlgorithm::IsShowerConnected(const CartesianVector &showerVertexPosition, const CartesianVector &nuVertex2D, 
    const CartesianVector &peakDirection) const
{
    CartesianVector displacement(showerVertexPosition - nuVertex2D);

    const float longitudinalSeparation(peakDirection.GetDotProduct(displacement));

    if (longitudinalSeparation < (-1.f))
        return false;

    const float transverseSeparation((peakDirection.GetCrossProduct(displacement)).GetMagnitude());

    if (transverseSeparation > m_maxCoincideneTransverseSeparation)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronInitialRegionRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoVertexListName", m_neutrinoVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameU", m_caloHitListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameV", m_caloHitListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameW", m_caloHitListNameW));

    AlgorithmTool *pAlgorithmTool1(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "PeakDirectionFinder", pAlgorithmTool1));
    m_pPeakDirectionFinderTool = dynamic_cast<PeakDirectionFinderTool *>(pAlgorithmTool1);

    if (!m_pPeakDirectionFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool2(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ShowerSpineFinder", pAlgorithmTool2));
    m_pShowerSpineFinderTool = dynamic_cast<ShowerSpineFinderTool *>(pAlgorithmTool2);

    if (!m_pShowerSpineFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool3(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ShowerStartFinder", pAlgorithmTool3));
    m_pShowerStartFinderTool = dynamic_cast<ShowerStartFinderTool *>(pAlgorithmTool3);

    if (!m_pShowerSpineFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool4(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ProtoShowerMatching", pAlgorithmTool4));
    m_pProtoShowerMatchingTool = dynamic_cast<ProtoShowerMatchingTool *>(pAlgorithmTool4);

    if (!m_pProtoShowerMatchingTool)
        return STATUS_CODE_INVALID_PARAMETER;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ShowerSlidingFitWindow", m_showerSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxCoincideneTransverseSeparation", m_maxCoincideneTransverseSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinSpinePurity", m_minSpinePurity));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content



        /////////////////////////
        // doesnt have to be here
        /*
        const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));
        CaloHitList viewShowerHitList;
        LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

        ElectronProtoShower protoShower(ShowerCore(showerStartPosition, showerStartDirection), 
            ConnectionPathway(nuVertex2D, peakDirection), 
            showerSpineHitList, false, CaloHitList(), CartesianPointVector(), CaloHitList());

        const float showerVertexL(std::max(peakDirection.GetDotProduct(showerStartPosition - nuVertex2D), 
            peakDirection.GetDotProduct(showerVertexPosition - nuVertex2D)));

        for (const CaloHit *const pCaloHit : showerSpineHitList)
        {
            if (std::find(viewShowerHitList.begin(), viewShowerHitList.end(), pCaloHit) != viewShowerHitList.end())
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const float showerL(peakDirection.GetDotProduct(hitPosition - nuVertex2D));

            if ((showerL > 0.f) && (showerL < showerVertexL))
                protoShower.m_hitsToAdd.push_back(pCaloHit);
        }
        */
        //////////////////////////////
