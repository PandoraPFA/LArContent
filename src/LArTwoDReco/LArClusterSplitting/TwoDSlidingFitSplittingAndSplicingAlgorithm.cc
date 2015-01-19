/**
 *  @file   LArContent/src/LArTwoDReco/ClusterSplitting/TwoDSlidingFitSplittingAndSplicingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting and splicing algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArGeometryHelper.h"
#include "LArHelpers/LArPointingClusterHelper.h"

#include "LArTwoDReco/LArClusterSplitting/TwoDSlidingFitSplittingAndSplicingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TwoDSlidingFitSplittingAndSplicingAlgorithm::TwoDSlidingFitSplittingAndSplicingAlgorithm() :
    m_shortHalfWindowLayers(10),
    m_longHalfWindowLayers(20),
    m_minClusterLength(7.5f),
    m_vetoDisplacement(1.5f),
    m_runCosmicMode(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAndSplicingAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    TwoDSlidingFitResultMap branchSlidingFitResultMap, replacementSlidingFitResultMap;

    unsigned int nIterations(0);

    while (++nIterations < 100) // Protect against flip-flopping between two answers
    {
        // Get ordered list of candidate clusters
        ClusterVector clusterVector;
        this->GetListOfCleanClusters(pClusterList, clusterVector);

        // Calculate sliding fit results for branch clusters (use a soft sliding fit for these)
        this->BuildSlidingFitResultMap(clusterVector, m_shortHalfWindowLayers, branchSlidingFitResultMap);

        // Calculate sliding fit results for replacement clusters (use a hard linear fit for these)
        this->BuildSlidingFitResultMap(clusterVector, m_longHalfWindowLayers, replacementSlidingFitResultMap);

        // Compile a list of possible splits
        ClusterExtensionList splitList;

        if (m_runCosmicMode)
        {
            this->BuildClusterExtensionList(clusterVector, branchSlidingFitResultMap, replacementSlidingFitResultMap, splitList);
        }
        else
        {
            ClusterExtensionList intermediateList;
            this->BuildClusterExtensionList(clusterVector, branchSlidingFitResultMap, replacementSlidingFitResultMap, intermediateList);
            this->PruneClusterExtensionList(intermediateList, branchSlidingFitResultMap, replacementSlidingFitResultMap, splitList);
        }

        // Run splitting and extension
        if (STATUS_CODE_SUCCESS != this->RunSplitAndExtension(splitList, branchSlidingFitResultMap, replacementSlidingFitResultMap))
            break;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitSplittingAndSplicingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitSplittingAndSplicingAlgorithm::BuildSlidingFitResultMap(const ClusterVector &clusterVector, const unsigned int halfWindowLayers,
    TwoDSlidingFitResultMap &slidingFitResultMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        if (slidingFitResultMap.end() == slidingFitResultMap.find(*iter))
        {
            try
            {
                const TwoDSlidingFitResult slidingFitResult(*iter, halfWindowLayers, slidingFitPitch);

                if (!slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(*iter, slidingFitResult)).second)
                    throw StatusCodeException(STATUS_CODE_FAILURE);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitSplittingAndSplicingAlgorithm::BuildClusterExtensionList(const ClusterVector &clusterVector, 
    const TwoDSlidingFitResultMap &branchSlidingFitResultMap, const TwoDSlidingFitResultMap &replacementSlidingFitResultMap,
    ClusterExtensionList &clusterExtensionList) const
{
    // Loop over each possible pair of clusters
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterEndI = clusterVector.end(); iterI != iterEndI; ++iterI)
    {
        Cluster* pClusterI = *iterI;

        for (ClusterVector::const_iterator iterJ = iterI, iterEndJ = clusterVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            Cluster* pClusterJ = *iterJ;

            if (pClusterI == pClusterJ)
                continue;

            // Get the branch and replacement sliding fits for this pair of clusters
            TwoDSlidingFitResultMap::const_iterator iterBranchI = branchSlidingFitResultMap.find(*iterI);
            TwoDSlidingFitResultMap::const_iterator iterBranchJ = branchSlidingFitResultMap.find(*iterJ);

            TwoDSlidingFitResultMap::const_iterator iterReplacementI = replacementSlidingFitResultMap.find(*iterI);
            TwoDSlidingFitResultMap::const_iterator iterReplacementJ = replacementSlidingFitResultMap.find(*iterJ);

            if (branchSlidingFitResultMap.end() == iterBranchI || branchSlidingFitResultMap.end() == iterBranchJ ||
                replacementSlidingFitResultMap.end() == iterReplacementI || replacementSlidingFitResultMap.end() == iterReplacementJ)
            {
                // TODO May want to raise an exception under certain conditions
                continue;
            }

            const TwoDSlidingFitResult &branchSlidingFitI(iterBranchI->second);
            const TwoDSlidingFitResult &branchSlidingFitJ(iterBranchJ->second);

            const TwoDSlidingFitResult &replacementSlidingFitI(iterReplacementI->second);
            const TwoDSlidingFitResult &replacementSlidingFitJ(iterReplacementJ->second);

            // Search for a split in clusterI
            float branchChisqI(0.f);
            CartesianVector branchSplitPositionI(0.f, 0.f, 0.f);
            CartesianVector branchSplitDirectionI(0.f, 0.f, 0.f);
            CartesianVector replacementStartPositionJ(0.f, 0.f, 0.f);

            try
            {
                this->FindBestSplitPosition(branchSlidingFitI, replacementSlidingFitJ, replacementStartPositionJ, branchSplitPositionI, branchSplitDirectionI);
                branchChisqI = this->CalculateBranchChi2(pClusterI, branchSplitPositionI, branchSplitDirectionI);
            }
            catch (StatusCodeException &)
            {
            }

            // Search for a split in clusterJ
            float branchChisqJ(0.f);
            CartesianVector branchSplitPositionJ(0.f, 0.f, 0.f);
            CartesianVector branchSplitDirectionJ(0.f, 0.f, 0.f);
            CartesianVector replacementStartPositionI(0.f, 0.f, 0.f);

            try
            {
                this->FindBestSplitPosition(branchSlidingFitJ, replacementSlidingFitI, replacementStartPositionI, branchSplitPositionJ, branchSplitDirectionJ);
                branchChisqJ = this->CalculateBranchChi2(pClusterJ, branchSplitPositionJ, branchSplitDirectionJ);
            }
            catch (StatusCodeException &)
            {
            }

            // Re-calculate chi2 values if both clusters have a split
            if (branchChisqI > 0.f && branchChisqJ > 0.f)
            {
                const CartesianVector relativeDirection((branchSplitPositionJ - branchSplitPositionI).GetUnitVector());

                if (branchSplitDirectionI.GetDotProduct(relativeDirection) > 0.f &&
                    branchSplitDirectionJ.GetDotProduct(relativeDirection) < 0.f )
                {
                    try
                    {
                        const float newBranchChisqI(this->CalculateBranchChi2(pClusterI, branchSplitPositionI, relativeDirection));
                        const float newBranchChisqJ(this->CalculateBranchChi2(pClusterJ, branchSplitPositionJ, relativeDirection * -1.f));
                        branchChisqI = newBranchChisqI;
                        branchChisqJ = newBranchChisqJ;
                    }
                    catch (StatusCodeException &)
                    {
                    }
                }
            }

            // Select the overall best split position
            if (branchChisqI > branchChisqJ)
            {
                (void) clusterExtensionList.push_back(ClusterExtension(pClusterI, pClusterJ, replacementStartPositionJ, branchSplitPositionI, branchSplitDirectionI));
            }

            else if (branchChisqJ > branchChisqI)
            {
                (void) clusterExtensionList.push_back(ClusterExtension(pClusterJ, pClusterI, replacementStartPositionI, branchSplitPositionJ, branchSplitDirectionJ));
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitSplittingAndSplicingAlgorithm::PruneClusterExtensionList(const ClusterExtensionList &inputList,
    const TwoDSlidingFitResultMap &branchMap, const TwoDSlidingFitResultMap &replacementMap, ClusterExtensionList &outputList) const
{
    for (ClusterExtensionList::const_iterator eIter = inputList.begin(), eIterEnd = inputList.end(); eIter != eIterEnd; ++eIter)
    {
        const ClusterExtension &thisSplit = *eIter;

        Cluster* pBranchCluster = thisSplit.GetBranchCluster();
        Cluster* pReplacementCluster = thisSplit.GetReplacementCluster();
        const CartesianVector &branchVertex = thisSplit.GetBranchVertex();
        const CartesianVector &replacementVertex = thisSplit.GetReplacementVertex();

        const float distanceSquared((branchVertex - replacementVertex).GetMagnitudeSquared());
        const float vetoDistanceSquared(m_vetoDisplacement * m_vetoDisplacement);

        bool branchVeto(false), replacementVeto(false);

        // Veto the merge if another cluster is closer to the replacement vertex
        for (TwoDSlidingFitResultMap::const_iterator iter = branchMap.begin(), iterEnd = branchMap.end(); iter != iterEnd; ++iter)
        {
            const TwoDSlidingFitResult &slidingFit(iter->second);

            if (slidingFit.GetCluster() == pReplacementCluster || slidingFit.GetCluster() == pBranchCluster)
                continue;

            const float minDistanceSquared((replacementVertex - slidingFit.GetGlobalMinLayerPosition()).GetMagnitudeSquared());
            const float maxDistanceSquared((replacementVertex - slidingFit.GetGlobalMaxLayerPosition()).GetMagnitudeSquared());

            if (std::min(minDistanceSquared, maxDistanceSquared) < std::max(distanceSquared, vetoDistanceSquared))
            {
                branchVeto = true;
                break;
            }
        }

        // Veto the merge if another cluster is closer to the branch vertex
        for (TwoDSlidingFitResultMap::const_iterator iter = replacementMap.begin(), iterEnd = replacementMap.end(); iter != iterEnd; ++iter)
        {
            const TwoDSlidingFitResult &slidingFit(iter->second);

            if (slidingFit.GetCluster() == pReplacementCluster || slidingFit.GetCluster() == pBranchCluster)
                continue;

            const float minDistanceSquared((branchVertex - slidingFit.GetGlobalMinLayerPosition()).GetMagnitudeSquared());
            const float maxDistanceSquared((branchVertex - slidingFit.GetGlobalMaxLayerPosition()).GetMagnitudeSquared());

            if (std::min(minDistanceSquared, maxDistanceSquared) < std::max(distanceSquared, vetoDistanceSquared))
            {
                replacementVeto = true;
                break;
            }
        }

        if (branchVeto || replacementVeto)
            continue;

        outputList.push_back(thisSplit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TwoDSlidingFitSplittingAndSplicingAlgorithm::CalculateBranchChi2(const Cluster* const pCluster, const CartesianVector &splitPosition,
    const CartesianVector &splitDirection) const
{
    CaloHitList principalCaloHitList, branchCaloHitList;

    this->SplitBranchCluster(pCluster, splitPosition, splitDirection, principalCaloHitList, branchCaloHitList);

    float totalChi2(0.f);
    float totalHits(0.f);

    for (CaloHitList::const_iterator iter = branchCaloHitList.begin(), iterEnd = branchCaloHitList.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;

        const CartesianVector hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector projectedPosition(splitPosition + splitDirection * splitDirection.GetDotProduct(hitPosition - splitPosition));

        totalChi2 += (hitPosition - projectedPosition).GetMagnitudeSquared();
        totalHits += 1.f;
    }

    if (totalHits > 0.f)
        return std::sqrt(totalChi2/totalHits);

    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoDSlidingFitSplittingAndSplicingAlgorithm::SplitBranchCluster(const Cluster* const pCluster, const CartesianVector &splitPosition,
    const CartesianVector &splitDirection, CaloHitList &principalCaloHitList, CaloHitList &branchCaloHitList) const
{
    // Distribute hits in branch cluster between new principal and residual clusters
    CaloHitList caloHitsToDistribute;
    pCluster->GetOrderedCaloHitList().GetCaloHitList(caloHitsToDistribute);

    for (CaloHitList::const_iterator iter = caloHitsToDistribute.begin(), iterEnd = caloHitsToDistribute.end(); iter != iterEnd; ++iter)
    {
        CaloHit *pCaloHit = *iter;

        if (splitDirection.GetDotProduct((pCaloHit->GetPositionVector() - splitPosition)) > 0.f)
        {
            branchCaloHitList.insert(pCaloHit);
        }
        else
        {
            principalCaloHitList.insert(pCaloHit);
        }
    }

    if (branchCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAndSplicingAlgorithm::RunSplitAndExtension(const ClusterExtensionList &splitList,
    TwoDSlidingFitResultMap &branchResultMap, TwoDSlidingFitResultMap &replacementResultMap) const
{
    bool foundSplit(false);

    for (ClusterExtensionList::const_iterator iter = splitList.begin(), iterEnd = splitList.end(); iter != iterEnd; ++iter)
    {
        const ClusterExtension &thisSplit = *iter;

        Cluster* pBranchCluster = thisSplit.GetBranchCluster();
        Cluster* pReplacementCluster = thisSplit.GetReplacementCluster();
        const CartesianVector &branchSplitPosition = thisSplit.GetBranchVertex();
        const CartesianVector &branchSplitDirection = thisSplit.GetBranchDirection();

        TwoDSlidingFitResultMap::iterator iterBranch1 = branchResultMap.find(pBranchCluster);
        TwoDSlidingFitResultMap::iterator iterBranch2 = branchResultMap.find(pReplacementCluster);

        TwoDSlidingFitResultMap::iterator iterReplacement1 = replacementResultMap.find(pBranchCluster);
        TwoDSlidingFitResultMap::iterator iterReplacement2 = replacementResultMap.find(pReplacementCluster);

        if (branchResultMap.end() == iterBranch1 || branchResultMap.end() == iterBranch2 ||
            replacementResultMap.end() == iterReplacement1 || replacementResultMap.end() == iterReplacement2)
            continue;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->ReplaceBranch(pBranchCluster, pReplacementCluster,
                                                                              branchSplitPosition, branchSplitDirection));
        branchResultMap.erase(iterBranch1);
        branchResultMap.erase(iterBranch2);

        replacementResultMap.erase(iterReplacement1);
        replacementResultMap.erase(iterReplacement2);

        foundSplit = true;
    }

    if (foundSplit)
        return STATUS_CODE_SUCCESS;

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAndSplicingAlgorithm::ReplaceBranch(Cluster *const pBranchCluster, Cluster *const pReplacementCluster,
    const CartesianVector &branchSplitPosition, const CartesianVector &branchSplitDirection) const
{
    ClusterList clusterList;
    clusterList.insert(pBranchCluster);
    clusterList.insert(pReplacementCluster);

    std::string clusterListToSaveName, clusterListToDeleteName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName,
        clusterListToSaveName));

    // Entire replacement cluster goes into new principal cluster
    PandoraContentApi::Cluster::Parameters principalParameters;
    pReplacementCluster->GetOrderedCaloHitList().GetCaloHitList(principalParameters.m_caloHitList);

    // Distribute hits in branch cluster between new principal and residual clusters
    PandoraContentApi::Cluster::Parameters residualParameters;
    this->SplitBranchCluster(pBranchCluster, branchSplitPosition, branchSplitDirection, principalParameters.m_caloHitList, residualParameters.m_caloHitList);

    Cluster *pPrincipalCluster(NULL), *pResidualCluster(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, principalParameters, pPrincipalCluster));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, residualParameters, pResidualCluster));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoDSlidingFitSplittingAndSplicingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "ShortHalfWindow", m_shortHalfWindowLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LongHalfWindow", m_longHalfWindowLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VetoDisplacement", m_vetoDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CosmicMode", m_runCosmicMode));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
