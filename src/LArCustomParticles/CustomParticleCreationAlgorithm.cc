/**
 *  @file   LArContent/src/LArCustomParticles/CustomParticleCreationAlgorithm.cc
 *
 *  @brief  Implementation of the 3D particle creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
            #include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArPfoHelper.h"

#include "LArCustomParticles/CustomParticleCreationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode CustomParticleCreationAlgorithm::Run()
{
    // Get input Pfo List
    const PfoList *pPfoList(NULL);

    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_pfoListName, pPfoList))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CustomParticleCreationAlgorithm: cannot find pfo list " << m_pfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    // Get input Vertex List
    const VertexList *pVertexList(NULL);

    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_vertexListName, pVertexList))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CustomParticleCreationAlgorithm: cannot find vertex list " << m_vertexListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    // Create temporary lists
    const PfoList *pTempPfoList = NULL; std::string tempPfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTempPfoList,
        tempPfoListName));

    const VertexList *pTempVertexList = NULL; std::string tempVertexListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pTempVertexList,
        tempVertexListName));

    // Loop over input Pfos
    PfoList pfoList(pPfoList->begin(), pPfoList->end());
    VertexList vertexList(pVertexList->begin(), pVertexList->end());

    for (PfoList::const_iterator iter = pfoList.begin(), iterEnd = pfoList.end(); iter != iterEnd; ++iter)
    {
        const ParticleFlowObject *const pInputPfo = *iter;

        if (pInputPfo->GetVertexList().empty())
            continue;

        const Vertex *const pInputVertex = LArPfoHelper::GetVertex(pInputPfo);

        if (0 == vertexList.count(pInputVertex))
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // Build a new pfo and vertex from the old pfo
        const ParticleFlowObject *pOutputPfo(NULL);

        this->CreatePfo(pInputPfo, pOutputPfo);

        if (NULL == pOutputPfo)
            continue;

        if (pOutputPfo->GetVertexList().empty())
            throw StatusCodeException(STATUS_CODE_FAILURE);

        // Transfer clusters and hierarchy information to new pfo, and delete old pfo and vertex
        ClusterList clusterList(pInputPfo->GetClusterList().begin(), pInputPfo->GetClusterList().end());
        PfoList parentList(pInputPfo->GetParentPfoList().begin(), pInputPfo->GetParentPfoList().end());
        PfoList daughterList(pInputPfo->GetDaughterPfoList().begin(), pInputPfo->GetDaughterPfoList().end());

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Pfo>(*this, pInputPfo, m_pfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete<Vertex>(*this, pInputVertex, m_vertexListName));

        for (ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo<Cluster>(*this, pOutputPfo, *cIter));
        }

        for (PfoList::const_iterator pIter = parentList.begin(), pIterEnd = parentList.end(); pIter != pIterEnd; ++pIter)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, *pIter, pOutputPfo));
        }

        for (PfoList::const_iterator dIter = daughterList.begin(), dIterEnd = daughterList.end(); dIter != dIterEnd; ++dIter)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SetPfoParentDaughterRelationship(*this, pOutputPfo, *dIter));
        }
    }

    if (!pTempPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_pfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_pfoListName));
    }

    if (!pTempVertexList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_vertexListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_vertexListName));
    }

const PfoList *pPfoList1(nullptr);
if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, m_pfoListName, pPfoList1))
{
    PfoVector pfoVector1(pPfoList1->begin(), pPfoList1->end());
    std::sort(pfoVector1.begin(), pfoVector1.end(), LArPfoHelper::SortByNHits);

    for (const Pfo *const pPfo1 : pfoVector1)
    {
        if ((pPfo1->GetParentPfoList().empty()) || ((1 == pPfo1->GetParentPfoList().size()) && LArPfoHelper::IsNeutrino(*(pPfo1->GetParentPfoList().begin()))))
            this->Print(pPfo1);
    }
}
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CustomParticleCreationAlgorithm::Print(const pandora::ParticleFlowObject *const pPfo1) const
{
    std::cout << "Alg " << this->GetType() << "---Pfo, nParents " << pPfo1->GetParentPfoList().size() << ", nDaughters " << pPfo1->GetDaughterPfoList().size() << std::endl;

    ClusterVector clusterVector1(pPfo1->GetClusterList().begin(), pPfo1->GetClusterList().end());
    std::sort(clusterVector1.begin(), clusterVector1.end(), LArClusterHelper::SortByNHits);
    for (const Cluster *const pCluster1 : clusterVector1)
        std::cout << "---PfoCluster " << this->GetType() << ", " << pCluster1->GetNCaloHits() << ", E " << pCluster1->GetHadronicEnergy()
         << " il " << pCluster1->GetInnerPseudoLayer() << " oc " << pCluster1->GetOrderedCaloHitList().size() << " span " << (pCluster1->GetOuterPseudoLayer() - pCluster1->GetInnerPseudoLayer()) << std::endl;

    for (const Vertex *const pVertex1 : pPfo1->GetVertexList())    
        std::cout << "---PfoVertex " << pVertex1->GetVertexLabel() << ", " << pVertex1->GetVertexType() << ", " << pVertex1->GetPosition() << std::endl;

    PfoVector pfoDaughters(pPfo1->GetDaughterPfoList().begin(), pPfo1->GetDaughterPfoList().end());
    std::sort(pfoDaughters.begin(), pfoDaughters.end(), LArPfoHelper::SortByNHits);

    for (const Pfo *const pDaughter : pfoDaughters)
        this->Print(pDaughter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CustomParticleCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
