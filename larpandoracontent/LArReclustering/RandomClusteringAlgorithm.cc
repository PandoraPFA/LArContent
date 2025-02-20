/**
 *  @file   larpandoracontent/LArReclustering/RandomClusteringAlgorithm.cc
 *
 *  @brief  Implementation file for the random clustering algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArReclustering/RandomClusteringAlgorithm.h"

#include <vector>
#include <algorithm>
#include <random>

using namespace pandora;

namespace lar_content
{

RandomClusteringAlgorithm::RandomClusteringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

RandomClusteringAlgorithm::~RandomClusteringAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RandomClusteringAlgorithm::Run()
{
    const CaloHitList *pCaloHits {nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::GetCurrentList(*this, pCaloHits));
    std::string hitListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::GetCurrentListName<CaloHit>(*this, hitListName));
    std::cout << "Random clustering alg - current calo hit list (" << hitListName << ") has " << pCaloHits->size() << " calo hits\n";

    if (pCaloHits->size() < m_numNewClusters)
        return STATUS_CODE_SUCCESS;

    std::vector<unsigned int> idxs;
    for (unsigned int i = 0; i < pCaloHits->size(); i++)
        idxs.push_back(i);
    std::shuffle(idxs.begin(), idxs.end(), std::default_random_engine(0));
    std::vector<unsigned int> splitIdxs;
    for (unsigned int i = 0; i < m_numNewClusters - 1; i++)
        splitIdxs.push_back(idxs.at(i));
    std::sort(splitIdxs.begin(), splitIdxs.end());

    std::vector<unsigned int>::const_iterator splitItr {splitIdxs.cbegin()};
    PandoraContentApi::Cluster::Parameters params;
    unsigned int cntr {0};
    std::cout << (splitItr == splitIdxs.cend()) << "\n";
    for (const CaloHit *const pCaloHit : *pCaloHits)
    {
        if (splitItr != splitIdxs.cend() && cntr++ > *splitItr)
        {
            const Cluster *pCluster {nullptr};
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                PandoraContentApi::Cluster::Create(*this, params, pCluster));
            std::cout << "Create with " << params.m_caloHitList.size() << " at split index " << *splitItr << "\n";
            params.m_caloHitList.clear();
            splitItr++;
        }
        params.m_caloHitList.push_back(pCaloHit);
    }
    const Cluster *pCluster {nullptr};
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::Cluster::Create(*this, params, pCluster));
    std::cout << "Create with " << params.m_caloHitList.size() << " at final split\n";
    params.m_caloHitList.clear();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode RandomClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        XmlHelper::ReadValue(xmlHandle, "NumNewClusters", m_numNewClusters));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
