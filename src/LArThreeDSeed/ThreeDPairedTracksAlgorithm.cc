/**
 *  @file   LArContent/src/LArThreeDSeed/ThreeDPairedTracksAlgorithm.cc
 * 
 *  @brief  Implementation of the 3D seed finding algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArHelpers/LArClusterHelper.h"
#include "LArHelpers/LArGeometryHelper.h"

#include "LArThreeDSeed/ThreeDPairedTracksAlgorithm.h"

using namespace pandora;

namespace lar
{

StatusCode ThreeDPairedTracksAlgorithm::Run()
{
    
    m_overlapTensor.Clear();
    


    const ClusterList* pInputClusterListU = NULL;
    const ClusterList* pInputClusterListV = NULL;
    const ClusterList* pInputClusterListW = NULL;

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this,
	    m_inputClusterListNameU, pInputClusterListU));

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this,
            m_inputClusterListNameV, pInputClusterListV));
        
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetClusterList(*this,
            m_inputClusterListNameW, pInputClusterListW));

    if ((NULL == pInputClusterListU) || (NULL == pInputClusterListV) || (NULL == pInputClusterListW))
    {
        std::cout << "ThreeDPairedTracksAlgorithm: one or more input cluster lists unavailable" << std::endl;
        throw StatusCodeException(STATUS_CODE_SUCCESS);
    }


    ClusterVector availableClustersU, availableClustersV, availableClustersW;
    ClusterVector particleClustersU, particleClustersV, particleClustersW;

    for (ClusterList::const_iterator iter = pInputClusterListU->begin(), iterEnd = pInputClusterListU->end(); iter != iterEnd; ++iter)
    {
        if( (*iter)->IsAvailable() )
            availableClustersU.push_back(*iter);
        else
            particleClustersU.push_back(*iter);
    }

    for (ClusterList::const_iterator iter = pInputClusterListV->begin(), iterEnd = pInputClusterListV->end(); iter != iterEnd; ++iter)
    {
        if( (*iter)->IsAvailable() )
            availableClustersV.push_back(*iter);
        else
            particleClustersV.push_back(*iter);
    }

    for (ClusterList::const_iterator iter = pInputClusterListW->begin(), iterEnd = pInputClusterListW->end(); iter != iterEnd; ++iter)
    {
        if( (*iter)->IsAvailable() )
            availableClustersW.push_back(*iter);
        else
            particleClustersW.push_back(*iter);
    }
       
    std::sort(availableClustersU.begin(), availableClustersU.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(availableClustersV.begin(), availableClustersV.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(availableClustersW.begin(), availableClustersW.end(), LArClusterHelper::SortByNOccupiedLayers);
    
    std::sort(particleClustersU.begin(), particleClustersU.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(particleClustersV.begin(), particleClustersV.end(), LArClusterHelper::SortByNOccupiedLayers);
    std::sort(particleClustersW.begin(), particleClustersW.end(), LArClusterHelper::SortByNOccupiedLayers);


    
    this->CalculateOverlapResult(availableClustersU, availableClustersV, particleClustersW);

    
    this->CalculateOverlapResult(particleClustersU,  availableClustersV, availableClustersW);

    
    this->CalculateOverlapResult(availableClustersU, particleClustersV,  availableClustersW);
    



    while ( this->BuildNextParticle() )
    {

    }




    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDPairedTracksAlgorithm::CalculateOverlapResult(const ClusterVector& clusterVectorU, const ClusterVector& clusterVectorV, const ClusterVector& clusterVectorW)
{
    for (ClusterVector::const_iterator iterU = clusterVectorU.begin(), iterEndU = clusterVectorU.end(); iterU != iterEndU; ++iterU)
    {
        for (ClusterVector::const_iterator iterV = clusterVectorV.begin(), iterEndV = clusterVectorV.end(); iterV != iterEndV; ++iterV)
        {
            for (ClusterVector::const_iterator iterW = clusterVectorW.begin(), iterEndW = clusterVectorW.end(); iterW != iterEndW; ++iterW)
	    {
                this->CalculateOverlapResult(*iterU, *iterV, *iterW);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThreeDPairedTracksAlgorithm::CalculateOverlapResult(Cluster *pClusterU, Cluster *pClusterV, Cluster *pClusterW)
{
    if ( ( pClusterU->IsAvailable() == false && pClusterV->IsAvailable() == false ) ||
         ( pClusterV->IsAvailable() == false && pClusterW->IsAvailable() == false ) ||
         ( pClusterW->IsAvailable() == false && pClusterU->IsAvailable() == false ) )
        return;


    LArClusterHelper::TwoDSlidingFitResult slidingFitResultU, slidingFitResultV, slidingFitResultW;
    LArClusterHelper::LArTwoDSlidingFit(pClusterU, 10, slidingFitResultU);
    LArClusterHelper::LArTwoDSlidingFit(pClusterV, 10, slidingFitResultV);
    LArClusterHelper::LArTwoDSlidingFit(pClusterW, 10, slidingFitResultW);


    const float UminX( slidingFitResultU.GetGlobalMinLayerPosition().GetX() );
    const float UmaxX( slidingFitResultU.GetGlobalMaxLayerPosition().GetX() );
    const float VminX( slidingFitResultV.GetGlobalMinLayerPosition().GetX() );
    const float VmaxX( slidingFitResultV.GetGlobalMaxLayerPosition().GetX() );
    const float WminX( slidingFitResultW.GetGlobalMinLayerPosition().GetX() );
    const float WmaxX( slidingFitResultW.GetGlobalMaxLayerPosition().GetX() );

    const float minX( std::max( std::min(UminX,UmaxX), std::max( std::min(VminX,VmaxX), std::min(WminX,WmaxX) ) ) );
    const float maxX( std::min( std::max(UminX,UmaxX), std::min( std::max(VminX,VmaxX), std::max(WminX,WmaxX) ) ) );

    const float xSpanU(UmaxX - UminX);
    const float xSpanV(VmaxX - VminX);
    const float xSpanW(WmaxX - WminX);
    const float xOverlap(maxX - minX);

    const float nPointsU(std::fabs((xOverlap / xSpanU) * static_cast<float>(slidingFitResultU.GetMaxLayer() - slidingFitResultU.GetMinLayer())));
    const float nPointsV(std::fabs((xOverlap / xSpanV) * static_cast<float>(slidingFitResultV.GetMaxLayer() - slidingFitResultV.GetMinLayer())));
    const float nPointsW(std::fabs((xOverlap / xSpanW) * static_cast<float>(slidingFitResultW.GetMaxLayer() - slidingFitResultW.GetMinLayer())));

    const float xPitch(3.f * xOverlap / (nPointsU + nPointsV + nPointsW));

    // Chi2 calculations
    float pseudoChi2Sum(0.f);
    unsigned int nSamplingPoints(0), nMatchedSamplingPoints(0);

    for (float x = minX; x < maxX; x += xPitch)
    {
        try
	{
            CartesianVector fitUVector(0.f, 0.f, 0.f), fitVVector(0.f, 0.f, 0.f), fitWVector(0.f, 0.f, 0.f);
            slidingFitResultU.GetGlobalFitPosition(x, true, fitUVector);
            slidingFitResultV.GetGlobalFitPosition(x, true, fitVVector);
            slidingFitResultW.GetGlobalFitPosition(x, true, fitWVector);

            CartesianVector fitUDirection(0.f, 0.f, 0.f), fitVDirection(0.f, 0.f, 0.f), fitWDirection(0.f, 0.f, 0.f);
            slidingFitResultU.GetGlobalFitDirection(x, true, fitUDirection);
            slidingFitResultV.GetGlobalFitDirection(x, true, fitVDirection);
            slidingFitResultW.GetGlobalFitDirection(x, true, fitWDirection);

            const float u(fitUVector.GetZ()), v(fitVVector.GetZ()), w(fitWVector.GetZ());
            const float uv2w(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_V, u, v));
            const float uw2v(LArGeometryHelper::MergeTwoPositions(VIEW_U, VIEW_W, u, w));
            const float vw2u(LArGeometryHelper::MergeTwoPositions(VIEW_V, VIEW_W, v, w));

            ++nSamplingPoints;
            const float deltaU((vw2u - u) * fitUDirection.GetX());
            const float deltaV((uw2v - v) * fitVDirection.GetX());
            const float deltaW((uv2w - w) * fitWDirection.GetX());

            const float pseudoChi2(deltaW * deltaW + deltaV * deltaV + deltaU * deltaU);
            pseudoChi2Sum += pseudoChi2;

            if (pseudoChi2 < 3.0)
                ++nMatchedSamplingPoints;
	}
        catch (StatusCodeException &)
	{
	}
    }

    if ( nMatchedSamplingPoints<5 ) return;

    m_overlapTensor.SetOverlapResult(pClusterU, pClusterV, pClusterW, TrackOverlapResult(nMatchedSamplingPoints,nSamplingPoints,pseudoChi2Sum));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThreeDPairedTracksAlgorithm::BuildNextParticle()
{
    

    unsigned int bestSamplingPoints(5);
    float bestReducedChi2(0.f);

    Cluster *pBestClusterU(NULL), *pBestClusterV(NULL), *pBestClusterW(NULL);

    const ClusterList &clusterListU(m_overlapTensor.GetClusterListU());
    const ClusterList &clusterListV(m_overlapTensor.GetClusterListV());
    const ClusterList &clusterListW(m_overlapTensor.GetClusterListW());

    for (ClusterList::const_iterator iterU = clusterListU.begin(), iterUEnd = clusterListU.end(); iterU != iterUEnd; ++iterU)
    {
        for (ClusterList::const_iterator iterV = clusterListV.begin(), iterVEnd = clusterListV.end(); iterV != iterVEnd; ++iterV)
        {
            for (ClusterList::const_iterator iterW = clusterListW.begin(), iterWEnd = clusterListW.end(); iterW != iterWEnd; ++iterW)
            {
	        Cluster* pClusterU = *iterU;
                Cluster* pClusterV = *iterV;
                Cluster* pClusterW = *iterW;

                if ( ( pClusterU->IsAvailable() == false && pClusterV->IsAvailable() == false ) ||
                     ( pClusterV->IsAvailable() == false && pClusterW->IsAvailable() == false ) ||
                     ( pClusterW->IsAvailable() == false && pClusterU->IsAvailable() == false ) )
		    continue;

	        try{
                    const TrackOverlapResult &overlapResult(m_overlapTensor.GetOverlapResult(pClusterU, pClusterV, pClusterW));

                    if ( (overlapResult.GetNMatchedSamplingPoints() > bestSamplingPoints) ||
                         (overlapResult.GetNMatchedSamplingPoints() == bestSamplingPoints && overlapResult.GetReducedChi2() < bestReducedChi2) )
		    {
		        bestSamplingPoints = overlapResult.GetNMatchedSamplingPoints();
                        bestReducedChi2 = overlapResult.GetReducedChi2();

                        pBestClusterU = pClusterU;
                        pBestClusterV = pClusterV;
                        pBestClusterW = pClusterW;
		    }
		}
                catch (StatusCodeException &)
		{
		}
                
            }
        }
    }

    if ( pBestClusterU && pBestClusterV && pBestClusterW )
    {
// ClusterList tempListU; tempListU.insert(pBestClusterU);
// ClusterList tempListV; tempListV.insert(pBestClusterV);
// ClusterList tempListW; tempListW.insert(pBestClusterW); 
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempListU, "ClusterU", RED);
// PandoraMonitoringApi::VisualizeClusters(&tempListV, "ClusterV", GREEN);
// PandoraMonitoringApi::VisualizeClusters(&tempListW, "ClusterW", BLUE);
// PandoraMonitoringApi::ViewEvent();

        const PfoList *pPfoList = NULL; std::string pfoListName;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryPfoListAndSetCurrent(*this, pPfoList, pfoListName));

        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = 22;
        pfoParameters.m_charge = 0;
        pfoParameters.m_mass = 0.f;
        pfoParameters.m_energy = 0.f;
        pfoParameters.m_momentum = CartesianVector(0., 0., 0.);

        if ( pBestClusterU->IsAvailable() ) pfoParameters.m_clusterList.insert(pBestClusterU);
	if ( pBestClusterV->IsAvailable() ) pfoParameters.m_clusterList.insert(pBestClusterV);
	if ( pBestClusterW->IsAvailable() ) pfoParameters.m_clusterList.insert(pBestClusterW);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters));
    
        if (!pPfoList->empty())
	{
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SavePfoList(*this, m_outputPfoListName));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentPfoList(*this, m_outputPfoListName));
	}

        return true;
    }
                

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThreeDPairedTracksAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameU", m_inputClusterListNameU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameV", m_inputClusterListNameV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputClusterListNameW", m_inputClusterListNameW));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputPfoListName", m_outputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
