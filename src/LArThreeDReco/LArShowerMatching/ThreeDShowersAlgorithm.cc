/**
 *  @file   LArContent/src/LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.cc
 * 
 *  @brief  Implementation of the three dimensional showers algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDReco/LArShowerMatching/ThreeDShowersAlgorithm.h"

using namespace pandora;

namespace lar
{

bool ThreeDShowersAlgorithm::SortByNMatchedSamplingPoints(const TensorType::Element &lhs, const TensorType::Element &rhs)
{
    if (lhs.GetOverlapResult().GetNMatchedSamplingPoints() != rhs.GetOverlapResult().GetNMatchedSamplingPoints())
        return (lhs.GetOverlapResult().GetNMatchedSamplingPoints() > rhs.GetOverlapResult().GetNMatchedSamplingPoints());

    if (lhs.GetOverlapResult().GetNSamplingPoints() != rhs.GetOverlapResult().GetNSamplingPoints())
        return (lhs.GetOverlapResult().GetNSamplingPoints() < rhs.GetOverlapResult().GetNSamplingPoints());

    return (lhs.GetOverlapResult().GetXOverlap().GetXOverlapSpan() < rhs.GetOverlapResult().GetXOverlap().GetXOverlapSpan());
}

//------------------------------------------------------------------------------------------------------------------------------------------

const TwoDSlidingShowerFitResult &ThreeDShowersAlgorithm::GetCachedSlidingFitResult(Cluster *const pCluster) const
{
    TwoDSlidingShowerFitResultMap::const_iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() == iter)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return iter->second;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::UpdateForNewCluster(Cluster *const pNewCluster)
{
    this->AddToSlidingFitCache(pNewCluster);
    ThreeDBaseAlgorithm<ShowerOverlapResult>::UpdateForNewCluster(pNewCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::UpdateUponDeletion(Cluster *const pDeletedCluster)
{
    this->RemoveFromSlidingFitCache(pDeletedCluster);
    ThreeDBaseAlgorithm<ShowerOverlapResult>::UpdateUponDeletion(pDeletedCluster);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::SelectInputClusters(const ClusterList *const pInputClusterList, ClusterList &selectedClusterList) const
{
    for (ClusterList::const_iterator iter = pInputClusterList->begin(), iterEnd = pInputClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (m_ignoreUnavailableClusters && !pCluster->IsAvailable())
            continue;

        if (pCluster->GetNCaloHits() < m_minClusterCaloHits)
            continue;

        if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLengthSquared)
            continue;

        selectedClusterList.insert(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::SetPfoParameters(const ProtoParticle &protoParticle, PandoraContentApi::ParticleFlowObject::Parameters &pfoParameters) const
{
    // TODO - correct these placeholder parameters
    pfoParameters.m_particleId = E_MINUS; // Shower
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList.insert(protoParticle.m_clusterListU.begin(), protoParticle.m_clusterListU.end());
    pfoParameters.m_clusterList.insert(protoParticle.m_clusterListV.begin(), protoParticle.m_clusterListV.end());
    pfoParameters.m_clusterList.insert(protoParticle.m_clusterListW.begin(), protoParticle.m_clusterListW.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::PreparationStep()
{
    ClusterList allClustersList;
    allClustersList.insert(this->m_clusterListU.begin(), this->m_clusterListU.end());
    allClustersList.insert(this->m_clusterListV.begin(), this->m_clusterListV.end());
    allClustersList.insert(this->m_clusterListW.begin(), this->m_clusterListW.end());

    for (ClusterList::const_iterator iter = allClustersList.begin(), iterEnd = allClustersList.end(); iter != iterEnd; ++iter)
    {
        this->AddToSlidingFitCache(*iter);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::TidyUp()
{
    m_slidingFitResultMap.clear();
    return ThreeDBaseAlgorithm<ShowerOverlapResult>::TidyUp();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::AddToSlidingFitCache(Cluster *const pCluster)
{
    const TwoDSlidingShowerFitResult slidingShowerFitResult(pCluster, m_slidingFitWindow);

    if (!m_slidingFitResultMap.insert(TwoDSlidingShowerFitResultMap::value_type(pCluster, slidingShowerFitResult)).second)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::RemoveFromSlidingFitCache(Cluster *const pCluster)
{
    TwoDSlidingShowerFitResultMap::iterator iter = m_slidingFitResultMap.find(pCluster);

    if (m_slidingFitResultMap.end() != iter)
        m_slidingFitResultMap.erase(iter);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    try
    {
        ShowerOverlapResult overlapResult;
        this->CalculateOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);

        if (overlapResult.IsInitialized())
            m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, overlapResult);
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (!(STATUS_CODE_NOT_FOUND == statusCodeException.GetStatusCode() || STATUS_CODE_NOT_INITIALIZED == statusCodeException.GetStatusCode()))
            throw statusCodeException;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW, ShowerOverlapResult &overlapResult)
{
    const TwoDSlidingShowerFitResult &fitResultU(this->GetCachedSlidingFitResult(pClusterU));
    const TwoDSlidingShowerFitResult &fitResultV(this->GetCachedSlidingFitResult(pClusterV));
    const TwoDSlidingShowerFitResult &fitResultW(this->GetCachedSlidingFitResult(pClusterW));

    const XSampling xSampling(fitResultU.GetShowerFitResult(), fitResultV.GetShowerFitResult(), fitResultW.GetShowerFitResult());

    ShowerPositionMapPair positionMapsU, positionMapsV, positionMapsW;
    this->GetShowerPositionMaps(fitResultU, fitResultV, fitResultW, xSampling, positionMapsU, positionMapsV, positionMapsW);

    unsigned int nSampledHitsU(0), nMatchedHitsU(0);
    this->GetBestHitOverlapFraction(pClusterU, xSampling, positionMapsU, nSampledHitsU, nMatchedHitsU);

    unsigned int nSampledHitsV(0), nMatchedHitsV(0);
    this->GetBestHitOverlapFraction(pClusterV, xSampling, positionMapsV, nSampledHitsV, nMatchedHitsV);

    unsigned int nSampledHitsW(0), nMatchedHitsW(0);
    this->GetBestHitOverlapFraction(pClusterW, xSampling, positionMapsW, nSampledHitsW, nMatchedHitsW);

    const unsigned int nMatchedHits(nMatchedHitsU + nMatchedHitsV + nMatchedHitsW);
    const unsigned int nSampledHits(nSampledHitsU + nSampledHitsV + nSampledHitsW);

    if (0 == nSampledHits)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const XOverlap xOverlapObject(xSampling.m_uMinX, xSampling.m_uMaxX, xSampling.m_vMinX, xSampling.m_vMaxX, xSampling.m_wMinX, xSampling.m_wMaxX, xSampling.m_xOverlapSpan);
    const ShowerOverlapResult showerOverlapResult(nMatchedHits, nSampledHits, xOverlapObject);

    if ((showerOverlapResult.GetMatchedFraction() < m_minShowerMatchedFraction) || (showerOverlapResult.GetNMatchedSamplingPoints() < m_minShowerMatchedPoints))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    overlapResult = showerOverlapResult;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::GetShowerPositionMaps(const TwoDSlidingShowerFitResult &fitResultU, const TwoDSlidingShowerFitResult &fitResultV,
    const TwoDSlidingShowerFitResult &fitResultW, const XSampling &xSampling, ShowerPositionMapPair &positionMapsU, ShowerPositionMapPair &positionMapsV,
    ShowerPositionMapPair &positionMapsW) const
{
    const unsigned int nPoints(static_cast<unsigned int>(xSampling.m_nPoints));
    
    for (unsigned n = 0; n <= nPoints; ++n)
    {
        const float x(xSampling.m_minX + (xSampling.m_maxX - xSampling.m_minX) * static_cast<float>(n) / static_cast<float>(nPoints));

        try
        {
            const int xBin(xSampling.GetBin(x));

            FloatVector uValues, vValues, wValues;
            fitResultU.GetShowerEdges(x, uValues);
            fitResultV.GetShowerEdges(x, vValues);
            fitResultW.GetShowerEdges(x, wValues);

            std::sort(uValues.begin(), uValues.end());
            std::sort(vValues.begin(), vValues.end());
            std::sort(wValues.begin(), wValues.end());

            if ((uValues.size() > 1) && (vValues.size() > 1))
            {
                const float uMin(uValues.front()), uMax(uValues.back());
                const float vMin(vValues.front()), vMax(vValues.back());
                const float uv2wMinMin(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_V, uMin, vMin));
                const float uv2wMaxMax(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_V, uMax, vMax));
                const float uv2wMinMax(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_V, uMin, vMax));
                const float uv2wMaxMin(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_V, uMax, vMin));
                positionMapsW.first.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, uv2wMinMin, uv2wMaxMax)));
                positionMapsW.second.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, uv2wMinMax, uv2wMaxMin)));
            }

            if ((uValues.size() > 1) && (wValues.size() > 1))
            {
                const float uMin(uValues.front()), uMax(uValues.back());
                const float wMin(wValues.front()), wMax(wValues.back());
                const float uw2vMinMin(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_W, uMin, wMin));
                const float uw2vMaxMax(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_W, uMax, wMax));
                const float uw2vMinMax(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_W, uMin, wMax));
                const float uw2vMaxMin(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_U, TPC_VIEW_W, uMax, wMin));
                positionMapsV.first.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, uw2vMinMin, uw2vMaxMax)));
                positionMapsV.second.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, uw2vMinMax, uw2vMaxMin)));
            }

            if ((vValues.size() > 1) && (wValues.size() > 1))
            {
                const float vMin(vValues.front()), vMax(vValues.back());
                const float wMin(wValues.front()), wMax(wValues.back());
                const float vw2uMinMin(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_V, TPC_VIEW_W, vMin, wMin));
                const float vw2uMaxMax(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_V, TPC_VIEW_W, vMax, wMax));
                const float vw2uMinMax(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_V, TPC_VIEW_W, vMin, wMax));
                const float vw2uMaxMin(LArGeometryHelper::MergeTwoPositions(TPC_VIEW_V, TPC_VIEW_W, vMax, wMin));
                positionMapsU.first.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, vw2uMinMin, vw2uMaxMax)));
                positionMapsU.second.insert(ShowerPositionMap::value_type(xBin, ShowerExtent(x, vw2uMinMax, vw2uMaxMin)));
            }
        }
        catch (StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::GetBestHitOverlapFraction(const Cluster *const pCluster, const XSampling &xSampling, const ShowerPositionMapPair &positionMaps,
    unsigned int &nSampledHits, unsigned int &nMatchedHits) const
{
    if ((xSampling.m_maxX - xSampling.m_minX) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    nSampledHits = 0; nMatchedHits = 0;
    unsigned int nMatchedHits1(0), nMatchedHits2(0);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hIter = iter->second->begin(), hIterEnd = iter->second->end(); hIter != hIterEnd; ++hIter)
        {
            CaloHit *pCaloHit = *hIter;
            const float x(pCaloHit->GetPositionVector().GetX());
            const float z(pCaloHit->GetPositionVector().GetZ());

            try
            {
                const int xBin(xSampling.GetBin(x));

                ++nSampledHits;

                ShowerPositionMap::const_iterator positionIter1 = positionMaps.first.find(xBin);
                ShowerPositionMap::const_iterator positionIter2 = positionMaps.second.find(xBin);

                if ((positionMaps.first.end() != positionIter1) && (z > positionIter1->second.GetLowEdgeZ()) && (z < positionIter1->second.GetHighEdgeZ()))
                    ++nMatchedHits1;

                if ((positionMaps.second.end() != positionIter2) && (z > positionIter2->second.GetLowEdgeZ()) && (z < positionIter2->second.GetHighEdgeZ()))
                    ++nMatchedHits2;
            }
            catch (StatusCodeException &)
            {
            }
        }
    }

    nMatchedHits = std::max(nMatchedHits1, nMatchedHits2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDShowersAlgorithm::ExamineTensor()
{
    unsigned int repeatCounter(0);

    for (TensorToolList::const_iterator iter = m_algorithmToolList.begin(), iterEnd = m_algorithmToolList.end(); iter != iterEnd; )
    {
        if ((*iter)->Run(this, m_overlapTensor))
        {
            iter = m_algorithmToolList.begin();

            if (++repeatCounter > m_nMaxTensorToolRepeats)
                break;
        }
        else
        {
            ++iter;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

ThreeDShowersAlgorithm::XSampling::XSampling(const TwoDSlidingFitResult &fitResultU, const TwoDSlidingFitResult &fitResultV,
    const TwoDSlidingFitResult &fitResultW)
{
    fitResultU.GetMinAndMaxX(m_uMinX, m_uMaxX);
    fitResultV.GetMinAndMaxX(m_vMinX, m_vMaxX);
    fitResultW.GetMinAndMaxX(m_wMinX, m_wMaxX);
    m_minX = std::max(m_uMinX, std::max(m_vMinX, m_wMinX));
    m_maxX = std::min(m_uMaxX, std::min(m_vMaxX, m_wMaxX));
    m_xOverlapSpan = (m_maxX - m_minX);
    m_nPoints = 1.f;

    if (m_xOverlapSpan < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const float nPointsU(std::fabs((m_xOverlapSpan / (m_uMaxX - m_uMinX)) * static_cast<float>(fitResultU.GetMaxLayer() - fitResultU.GetMinLayer())));
    const float nPointsV(std::fabs((m_xOverlapSpan / (m_vMaxX - m_vMinX)) * static_cast<float>(fitResultV.GetMaxLayer() - fitResultV.GetMinLayer())));
    const float nPointsW(std::fabs((m_xOverlapSpan / (m_wMaxX - m_wMinX)) * static_cast<float>(fitResultW.GetMaxLayer() - fitResultW.GetMinLayer())));

    m_nPoints = 1.f + ((nPointsU + nPointsV + nPointsW) / 3.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

int ThreeDShowersAlgorithm::XSampling::GetBin(const float x) const
{
    if (((x - m_minX) < -std::numeric_limits<float>::epsilon()) || ((x - m_maxX) > +std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    return static_cast<int>(0.5f + m_nPoints * (x - m_minX) / (m_maxX - m_minX));
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDShowersAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    AlgorithmToolList algorithmToolList;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle,
        "ShowerTools", algorithmToolList));

    for (AlgorithmToolList::const_iterator iter = algorithmToolList.begin(), iterEnd = algorithmToolList.end(); iter != iterEnd; ++iter)
    {
        ShowerTensorTool *pShowerTensorTool(dynamic_cast<ShowerTensorTool*>(*iter));

        if (NULL == pShowerTensorTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolList.push_back(pShowerTensorTool);
    }

    m_nMaxTensorToolRepeats = 5000;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NMaxTensorToolRepeats", m_nMaxTensorToolRepeats));

    m_slidingFitWindow = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    m_ignoreUnavailableClusters = true;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "IgnoreUnavailableClusters", m_ignoreUnavailableClusters));

    m_minClusterCaloHits = 5;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterCaloHits", m_minClusterCaloHits));

    float minClusterLength = 3.f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLength", minClusterLength));
    m_minClusterLengthSquared = minClusterLength * minClusterLength;

    m_minShowerMatchedFraction = 0.2f;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerMatchedFraction", m_minShowerMatchedFraction));

    m_minShowerMatchedPoints = 20;
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinShowerMatchedPoints", m_minShowerMatchedPoints));

    return ThreeDBaseAlgorithm<ShowerOverlapResult>::ReadSettings(xmlHandle);
}

} // namespace lar
