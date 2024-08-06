/**
 *  @file   larpandoracontent/LArCheating/CheatingPfoCreationAlgorithm.cc
 *
 *  @brief  Implementation of the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingPfoCreationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"


using namespace pandora;

namespace lar_content
{

CheatingPfoCreationAlgorithm::CheatingPfoCreationAlgorithm() :
    m_collapseToPrimaryMCParticles(false),
    m_useOnlyAvailableClusters(true),
    m_addVertices(true),
    m_replaceCurrentVertexList(false),
    m_minGoodHitTypes(0),
    m_nHitsForGoodHitType(10)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingPfoCreationAlgorithm::Run()
{
    LArMCParticleHelper::MCRelationMap mcPrimaryMap;

    if (m_collapseToPrimaryMCParticles)
    {
        const MCParticleList *pMCParticleList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

        LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcPrimaryMap);
    }

    MCParticleToClusterListMap mcParticleToClusterListMap;

    for (const std::string &clusterListName : m_inputClusterListNames)
    {
        const ClusterList *pClusterList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, clusterListName, pClusterList))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "CheatingPfoCreationAlgorithm - Could not access cluster list with name " << clusterListName << std::endl;

            continue;
        }

        this->GetMCParticleToClusterListMap(pClusterList, mcPrimaryMap, mcParticleToClusterListMap);
    }

    this->CreatePfos(mcParticleToClusterListMap);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingPfoCreationAlgorithm::GetMCParticleToClusterListMap(const ClusterList *const pClusterList,
    const LArMCParticleHelper::MCRelationMap &mcPrimaryMap, MCParticleToClusterListMap &mcParticleToClusterListMap) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {
	try
        {
            if (m_useOnlyAvailableClusters && !PandoraContentApi::IsAvailable(*this, pCluster))
                continue;

            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCluster));

            if (m_collapseToPrimaryMCParticles)
            {
                LArMCParticleHelper::MCRelationMap::const_iterator primaryIter = mcPrimaryMap.find(pMCParticle);

                if (mcPrimaryMap.end() == primaryIter)
                    throw StatusCodeException(STATUS_CODE_NOT_FOUND);

                pMCParticle = primaryIter->second;
            }

            if (!m_particleIdList.empty() && !m_particleIdList.count(pMCParticle->GetParticleId()))
                continue;

            mcParticleToClusterListMap[pMCParticle].push_back(pCluster);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingPfoCreationAlgorithm::CreatePfos(const MCParticleToClusterListMap &mcParticleToClusterListMap) const
{
    if (mcParticleToClusterListMap.empty())
        return;

    const PfoList *pPfoList(nullptr);
    std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    const VertexList *pVertexList(nullptr);
    std::string vertexListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, vertexListName));

    MCParticleVector mcParticleVector;
    for (const auto &mapEntry : mcParticleToClusterListMap)
        mcParticleVector.push_back(mapEntry.first);
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcParticleVector)
    {
        const ClusterList &clusterList(mcParticleToClusterListMap.at(pMCParticle));

        if (clusterList.empty())
            continue;

        if (this->GetNHitTypesAboveThreshold(clusterList, m_nHitsForGoodHitType) < m_minGoodHitTypes)
            continue;

        try
        {
            PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
            pfoParameters.m_particleId = pMCParticle->GetParticleId();
            pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
            pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
            pfoParameters.m_energy = pMCParticle->GetEnergy();
            pfoParameters.m_momentum = pMCParticle->GetMomentum();
            pfoParameters.m_clusterList.insert(pfoParameters.m_clusterList.end(), clusterList.begin(), clusterList.end());

	    const ParticleFlowObject *pPfo(nullptr);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));
            
            
            if (m_addVertices && pMCParticle->GetEnergy() > 0)
            {
		try
		{
                    PandoraContentApi::Vertex::Parameters parameters;
	            parameters.m_position = pMCParticle->GetVertex();
                    parameters.m_vertexLabel = VERTEX_START;
                    parameters.m_vertexType = VERTEX_3D;

		    const Vertex *pVertex(nullptr);
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo<Vertex>(*this, pPfo, pVertex));
		}
		catch (StatusCodeException &)
		{
	            std::cout << "Background - Has no vertex" << std::endl;
		}

            }
        }
        catch (const StatusCodeException &)
        {
            std::cout << "CheatingPfoCreationAlgorithm: Could not create PFO for MCParticle with pdg code " << pMCParticle->GetParticleId()
                      << std::endl;
        }
    }

    if (!pPfoList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
    }

    if (!pVertexList->empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));

        if (m_replaceCurrentVertexList)
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int CheatingPfoCreationAlgorithm::GetNHitTypesAboveThreshold(const ClusterList &clusterList, const unsigned int nHitsThreshold) const
{
    HitTypeMap hitTypeMap;

    for (const Cluster *const pCluster : clusterList)
    {
        hitTypeMap[LArClusterHelper::GetClusterHitType(pCluster)] += pCluster->GetNCaloHits();
    }

    unsigned int nGoodViews(0);

    for (const HitTypeMap::value_type &mapEntry : hitTypeMap)
    {
        if (mapEntry.second > nHitsThreshold)
            ++nGoodViews;
    }

    return nGoodViews;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingPfoCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CollapseToPrimaryMCParticles", m_collapseToPrimaryMCParticles));

    if (m_collapseToPrimaryMCParticles)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "UseOnlyAvailableClusters", m_useOnlyAvailableClusters));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "AddVertices", m_addVertices));

    if (m_addVertices)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
            XmlHelper::ReadValue(xmlHandle, "ReplaceCurrentVertexList", m_replaceCurrentVertexList));
    }

    IntVector particleIdVector;
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ParticleIdList", particleIdVector));

    m_particleIdList.insert(particleIdVector.begin(), particleIdVector.end());

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinGoodHitTypes", m_minGoodHitTypes));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NHitsForGoodHitType", m_nHitsForGoodHitType));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
