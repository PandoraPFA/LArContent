/**
 *  @file   larpandoradlcontent/LArSignalID/DlSignalAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning signal algorithm.
 *
 *  $Log: $
 */

#include <chrono>
#include <cmath>

#include <torch/script.h>
#include <torch/torch.h>

#include "Pandora/Pandora.h"

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArVertexHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArSignalId/DlSignalAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlSignalAlgorithm::DlSignalAlgorithm() :
    m_trainingMode{false},
    m_trainingOutputFile{""},
    m_event{-1},
    m_pass{1},
    m_height{256},
    m_width{256},
    m_driftStep{0.5f},
    m_visualise{false},
    m_writeTree{false},
    m_rng(static_cast<std::mt19937::result_type>(std::chrono::high_resolution_clock::now().time_since_epoch().count())),
    m_printOut{false},
    m_signalListNameU{""},
    m_signalListNameV{""},
    m_signalListNameW{""},
    m_signalListName2D{""},
    m_caloHitListName2D{""},
    m_backgroundListName{""},
    m_applyCheatedSeparation{false},
    m_simpleZoom{false},
    m_passOneTrustThreshold{0}
{
}

DlSignalAlgorithm::~DlSignalAlgorithm()
{
    if (m_writeTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "SignalAssessmentAlgorithm: Unable to write to ROOT tree" << std::endl;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSignalAlgorithm::Run()
{
    if (m_visualise)
    {
	PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    }
    
    if (!m_applyCheatedSeparation)
    {    
        if (m_trainingMode)	    
	    return this->PrepareTrainingSample();
        else
            return this->Infer();
    }
    else
        return this->CheatedSeparation();

    return STATUS_CODE_SUCCESS;
}

StatusCode DlSignalAlgorithm::PrepareTrainingSample()
{
    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetMCToHitsMap(mcToHitsMap));
    if (mcToHitsMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // Get boundaries for hits and make x dimension common
    std::map<int, float> wireMin, wireMax;
    std::map<int, bool> viewCalculated;
    float driftMin{std::numeric_limits<float>::max()}, driftMax{-std::numeric_limits<float>::max()};
    for (const std::string &listname : m_inputCaloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
	try
        {
            PandoraContentApi::GetList(*this, listname, pCaloHitList);
	}
	catch (const StatusCodeException &e)
	{
            continue;
	}
        if (!pCaloHitList || pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        float viewDriftMin{driftMin}, viewDriftMax{driftMax};
        try
        {
            this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
	    viewCalculated[view] = true;
        }
        catch (const StatusCodeException &e)
        {
            if (e.GetStatusCode() == STATUS_CODE_NOT_FOUND)
	        continue;
	    else
                throw;
        }

	driftMin = std::min(viewDriftMin, driftMin);
        driftMax = std::max(viewDriftMax, driftMax);
    }

    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        if (!pCaloHitList)
        {
            std::cout << "ERR: Could not find full CaloHitList - DlSignalAlgorithm unable to proceed" << std::endl;
            continue;
        }

        HitType view{pCaloHitList->front()->GetHitType()};
        float viewDriftMin{driftMin}, viewDriftMax{driftMax};
        if (viewCalculated[view] != true)
        {
            const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
            const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
            const float zSpan{pitch * (m_height - 1)};
            bool projected{false};

            if (viewCalculated[6] == true && viewCalculated[4] == true) //Check W and U - Project to V
            {
                //float x{(viewDriftMax - viewDriftMin) * 0.5};
                const double z1{(wireMax[6] - wireMin[6]) * 0.5};
                const double z2{(wireMax[4] - wireMin[4]) * 0.5};
                const double m_projectedCoordinate{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()->WUtoV(z1, z2)};
                std::cout << "Projected Coordinated: " << m_projectedCoordinate <<std::endl;
                wireMin[view] = m_projectedCoordinate - zSpan;
                wireMax[view] = m_projectedCoordinate + zSpan;
                projected = true;
            }
            if (viewCalculated[6] == true && viewCalculated[5] == true) //Check W and V - Project to U
            {
                const double z1{(wireMax[6] - wireMin[6]) * 0.5};
                const double z2{(wireMax[5] - wireMin[5]) * 0.5};
                this->GetPandora();
                const double m_projectedCoordinate{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()->VWtoU(z1, z2)};
                std::cout << "Projected Coordinated: " << m_projectedCoordinate <<std::endl;
                wireMin[view] = m_projectedCoordinate - zSpan;
                wireMax[view] = m_projectedCoordinate + zSpan;
                projected = true;
            }
            if (viewCalculated[4] == true && viewCalculated[5] == true) //Check U and V - Project to W
            {
                const double z1{(wireMax[4] - wireMin[4]) * 0.5};
                const double z2{(wireMax[5] - wireMin[5]) * 0.5};
                const double m_projectedCoordinate{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()->UVtoW(z1, z2)};
                std::cout << "Projected Coordinated: " << m_projectedCoordinate <<std::endl;
                wireMin[view] = m_projectedCoordinate - zSpan;
                wireMax[view] = m_projectedCoordinate + zSpan;
                projected = true;
            }
            if (projected == false)
            {
                try
                {
                    m_simpleZoom = true;
                    this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
                    m_simpleZoom = false;
                }
                catch (const StatusCodeException &e)
                {
                    std::cout << "ERR: Could not calculate zoom region - DlSignalAlgorithm unable to proceed" << std::endl;
                    if (e.GetStatusCode() == STATUS_CODE_NOT_FOUND)
                    {
                        continue;
                    }
                }
            }

        }
        driftMin = std::min(viewDriftMin, driftMin);
        driftMax = std::max(viewDriftMax, driftMax);
    }

    
    
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
        if (!(isU || isV || isW))
            return STATUS_CODE_NOT_ALLOWED;

        const std::string trainingFilename{m_trainingOutputFile + "_" + listname + ".csv"};
        unsigned long nHits{0};

        // Calo hits
        double xMin{driftMin}, xMax{driftMax}, zMin{wireMin[view]}, zMax{wireMax[view]};

        LArMvaHelper::MvaFeatureVector featureVector;
	featureVector.emplace_back(xMin);
        featureVector.emplace_back(xMax);
        featureVector.emplace_back(zMin);
        featureVector.emplace_back(zMax);

	for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            const float x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()}, adc{pCaloHit->GetMipEquivalentEnergy()};
            // If on a refinement pass, drop hits outside the region of interest
            if (m_pass > 1 && (x < xMin || x > xMax || z < zMin || z > zMax))
                continue;
	    featureVector.emplace_back(static_cast<double>(x));
            featureVector.emplace_back(static_cast<double>(z));
            featureVector.emplace_back(static_cast<double>(adc));
	    ++nHits;
  
	    const MCParticle *pMainMCParticle(nullptr);
	    try {pMainMCParticle = MCParticleHelper::GetMainMCParticle(pCaloHit);}
	    catch (const StatusCodeException &) {}

	    if (pMainMCParticle)
            {
                const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMainMCParticle));

                if (LArMCParticleHelper::IsNeutrino(pParentMCParticle))
                { 
                    if (pMainMCParticle->GetParticleId() == 11)
	            {
		        featureVector.emplace_back(11);
		    }
		    else if (pMainMCParticle->GetParticleId() == 22)
		    {
	                featureVector.emplace_back(22);
		    }
	        }
                else
                {
                    featureVector.emplace_back(0);
                }
            }
	    else
	    {
	        featureVector.emplace_back(0);
	    }
	}
	featureVector.insert(featureVector.begin() + 4, static_cast<double>(nHits));
        LArMvaHelper::ProduceTrainingExample(trainingFilename, true, featureVector);
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSignalAlgorithm::Infer()
{
    if (m_pass == 1)
        ++m_event;

    std::map<int, float> wireMin, wireMax;
    std::map<int, bool> viewCalculated;
    float driftMin{std::numeric_limits<float>::max()}, driftMax{-std::numeric_limits<float>::max()};
    for (const std::string &listname : m_inputCaloHitListNames)
    {
	const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));

        if (!pCaloHitList)
	    continue;

        if (pCaloHitList->size() < m_passOneTrustThreshold)
	    continue;

       	HitType view{pCaloHitList->front()->GetHitType()};
	float viewDriftMin{driftMin}, viewDriftMax{driftMax};
        try
	{
	    this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
	    viewCalculated[view] = true;
	}
	catch (const StatusCodeException &e)
	{
	    if (e.GetStatusCode() == STATUS_CODE_NOT_FOUND)
            {
                continue;
	    }
	}
	driftMin = std::min(viewDriftMin, driftMin);
        driftMax = std::max(viewDriftMax, driftMax);
    }

    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

	if (!pCaloHitList)
	{
	    std::cout << "ERR: Could not find full CaloHitList - DlSignalAlgorithm unable to proceed" << std::endl;
            continue;
	}
	
	HitType view{pCaloHitList->front()->GetHitType()};
	float viewDriftMin{driftMin}, viewDriftMax{driftMax};
	if (viewCalculated[view] != true)
	{
	    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
            const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
            const float zSpan{pitch * (m_height - 1)};
	    bool projected{false};

	    if (viewCalculated[6] == true && viewCalculated[4] == true) //Check W and U - Project to V
            {
                //float x{(viewDriftMax - viewDriftMin) * 0.5};
		const double z1{(wireMax[6] - wireMin[6]) * 0.5};
		const double z2{(wireMax[4] - wireMin[4]) * 0.5};
                const double m_projectedCoordinate{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()->WUtoV(z1, z2)};
		std::cout << "Projected Coordinated: " << m_projectedCoordinate <<std::endl;
                wireMin[view] = m_projectedCoordinate - zSpan;
                wireMax[view] = m_projectedCoordinate + zSpan;
		projected = true;
	    }
            if (viewCalculated[6] == true && viewCalculated[5] == true) //Check W and V - Project to U
	    {
		const double z1{(wireMax[6] - wireMin[6]) * 0.5};
                const double z2{(wireMax[5] - wireMin[5]) * 0.5};
		this->GetPandora();
                const double m_projectedCoordinate{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()->VWtoU(z1, z2)};
		std::cout << "Projected Coordinated: " << m_projectedCoordinate <<std::endl;
		wireMin[view] = m_projectedCoordinate - zSpan;
                wireMax[view] = m_projectedCoordinate + zSpan;
		projected = true;
	    }
	    if (viewCalculated[4] == true && viewCalculated[5] == true) //Check U and V - Project to W
	    {
		const double z1{(wireMax[4] - wireMin[4]) * 0.5};
                const double z2{(wireMax[5] - wireMin[5]) * 0.5};
                const double m_projectedCoordinate{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()->UVtoW(z1, z2)};
		std::cout << "Projected Coordinated: " << m_projectedCoordinate <<std::endl;
		wireMin[view] = m_projectedCoordinate - zSpan;
                wireMax[view] = m_projectedCoordinate + zSpan;
		projected = true;
            }
            
	    if (projected == false)
            {
                try
                {
	            m_simpleZoom = true;
                    this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
	            m_simpleZoom = false;
	        }
                catch (const StatusCodeException &e)
                {
                    std::cout << "ERR: Could not calculate zoom region - DlSignalAlgorithm unable to proceed" << std::endl;
	            if (e.GetStatusCode() == STATUS_CODE_NOT_FOUND)
                    {
                        continue;
                    }
                }
            }

	}
	driftMin = std::min(viewDriftMin, driftMin);
        driftMax = std::max(viewDriftMax, driftMax);
    }
 
	
    
    CaloHitList signalCandidatesU, signalCandidatesV, signalCandidatesW, signalCandidates2D, backgroundCaloHitList;
    for (const std::string &listName : m_caloHitListNames)
    {
	const CaloHitList *pCaloHitList{nullptr};
	PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));
	
	if (!pCaloHitList || pCaloHitList->empty())
   	    continue;
	
	HitType view;
	try
	{
	    view = pCaloHitList->front()->GetHitType();
	}
	catch (const StatusCodeException &)
	{
	    continue;
	}
       
        const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
        if (!isU && !isV && !isW)
            return STATUS_CODE_NOT_ALLOWED;

	LArDLHelper::TorchInput input;
        PixelMap pixelMap;
  	this->MakeNetworkInputFromHits(*pCaloHitList, view, driftMin, driftMax, wireMin[view], wireMax[view], input, pixelMap);

        // Run the input through the trained model
        LArDLHelper::TorchInputVector inputs;
        inputs.push_back(input);
        LArDLHelper::TorchOutput output;
        if (isU)
            LArDLHelper::Forward(m_modelU, inputs, output);
        else if (isV)
            LArDLHelper::Forward(m_modelV, inputs, output);
        else
            LArDLHelper::Forward(m_modelW, inputs, output);

        // we want the maximum value in the num_classes dimension (1) for every pixel
        auto classes{torch::argmax(output, 1)};
        // the argmax result is a 1 x height x width tensor where each element is a class id
        auto classesAccessor{classes.accessor<long, 3>()};
        std::map<int, bool> haveSeenMap;
	
	for (const auto &[pCaloHit, pixel] : pixelMap)
        {
	    //The ordering of the pixel is x coordinate, z coordinate 
            const auto cls{classesAccessor[0][pixel.second][pixel.first]};
	    if (m_printOut)
            {
		if (m_pass > 1 && cls == 3)    
                {
		    std::cout << "*Electron identified*" << std::endl;
	        }
	        if (m_pass > 1 && cls == 2)
	        {
	            std::cout << "*Photon identified*" << std::endl;    
		}
		if (m_pass < 2 && cls == 2)
	        {
                    std::cout << "*Signal Pixel identified*" << std::endl;
		}
	    }

            if (cls == 2 || cls == 3)
            {
		signalCandidates2D.emplace_back(pCaloHit);

		if (isU)
	        {
                    signalCandidatesU.emplace_back(pCaloHit);
		}
                else if (isV)
		{ 
	            signalCandidatesV.emplace_back(pCaloHit);
		}
                else
                {
		    signalCandidatesW.emplace_back(pCaloHit);
		}
            }
	    else
            {
	        backgroundCaloHitList.emplace_back(pCaloHit);
	    }


        }


	if (m_visualise)
        {
            if (isU)
	    {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                                         &signalCandidatesU, "candidate signal U", RED));
            }

            if (isV)
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                                         &signalCandidatesV, "candidate signal V", RED));
	    }

            if (isW)
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                                         &signalCandidatesW, "candidate signal W", RED));
	    }

            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }

#ifdef MONITORING
        if (m_visualise)
        {
            PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
            try
            {
		for (const CaloHit *pCaloHit : *pCaloHitList)
                {
                    const float x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()};
                    const MCParticle *pMainMCParticle(nullptr);
                    try {pMainMCParticle = MCParticleHelper::GetMainMCParticle(pCaloHit);}
                    catch (const StatusCodeException &) {}

		    if (pMainMCParticle)
                    {
                        const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMainMCParticle));

                        if (LArMCParticleHelper::IsNeutrino(pParentMCParticle))
                        {
                             const CartesianVector signalHit(x, 0.f, z);
                             PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &signalHit, "true signal", YELLOW, 2));
	  		}

                    }
		 }
            }
            catch (StatusCodeException &e)
            {
                std::cerr << "DlSignalAlgorithm: Warning. Couldn't find signal hits." << std::endl;
            }
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        }
#endif

    }
    if (m_printOut)
    {
        std::cout << "Printing U view length: " << signalCandidatesU.size() << std::endl;
        std::cout << "Printing V view length: " << signalCandidatesV.size() << std::endl;
        std::cout << "Printing W view length: " << signalCandidatesW.size() << std::endl;
        std::cout << "Printing 2D view length: " << signalCandidates2D.size() << std::endl;
        std::cout << "Printing background CaloHitList length: " << backgroundCaloHitList.size() << std::endl;
        std::cout << "Printing New CaloHitList Names: " << std::endl;
        std::cout << m_signalListNameU << " | " << m_signalListNameV << " | " << m_signalListNameW <<
	                         	" | " << m_signalListName2D << " | " << m_backgroundListName << std::endl;
    }
    
    if (signalCandidatesU.empty() || signalCandidatesV.empty() || signalCandidatesW.empty() || signalCandidates2D.empty())
    {
        std::cout << "Error: A CaloHitList is empty" << std::endl;
    }
    if (!signalCandidatesU.empty())    
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, signalCandidatesU, m_signalListNameU));
    
    if (!signalCandidatesV.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, signalCandidatesV, m_signalListNameV));
    
    if (!signalCandidatesW.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, signalCandidatesW, m_signalListNameW));
    
    if (!signalCandidates2D.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, signalCandidates2D, m_signalListName2D));
    
    if (!backgroundCaloHitList.empty())
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, backgroundCaloHitList, m_backgroundListName));
    
    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSignalAlgorithm::CheatedSeparation()
{
    CaloHitList signalCandidatesU, signalCandidatesV, signalCandidatesW, signalCandidates2D, backgroundCaloHitList;
    for (const std::string &listname : m_caloHitListNames)
    {

        const CaloHitList *pCaloHitList(nullptr);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;    
      
      	HitType view{pCaloHitList->front()->GetHitType()};
        const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
        if (!isU && !isV && !isW)
            return STATUS_CODE_NOT_ALLOWED;
       
	for (const CaloHit *pCaloHit : *pCaloHitList)
        {
            const MCParticle *pMainMCParticle(nullptr);
            try {pMainMCParticle = MCParticleHelper::GetMainMCParticle(pCaloHit);}
            catch (const StatusCodeException &) {}
 
	    if (pMainMCParticle)
            {
                const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMainMCParticle));

                if (LArMCParticleHelper::IsNeutrino(pParentMCParticle))
                {
                    if (std::abs(pMainMCParticle->GetParticleId()) == 11 || std::abs(pMainMCParticle->GetParticleId()) == 22)
                    {
                        signalCandidates2D.emplace_back(pCaloHit);

                        if (isU)
                        {
                            signalCandidatesU.emplace_back(pCaloHit);
                        }
                        else if (isV)
                        {
                            signalCandidatesV.emplace_back(pCaloHit);
                        }
                        else
                        {
                            signalCandidatesW.emplace_back(pCaloHit);
	                }
                    }
                }
	   	else
                {
                    backgroundCaloHitList.emplace_back(pCaloHit);
                }
	
	    }
	}

    

        if (m_visualise)
        {
            if (isU)
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                                         &signalCandidatesU, "candidate signal U", RED));
            }

            if (isV)
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                                         &signalCandidatesV, "candidate signal V", RED));
            }

            if (isW)
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                                         &signalCandidatesW, "candidate signal W", RED));
            }

                PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
                try
                {
                    for (const CaloHit *pCaloHit : *pCaloHitList)
                    {
                        const float x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()};
                        const MCParticle *pMainMCParticle(nullptr);
                        try {pMainMCParticle = MCParticleHelper::GetMainMCParticle(pCaloHit);}
                        catch (const StatusCodeException &) {}

                        if (pMainMCParticle)
                        {
                            const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMainMCParticle));

                            if (LArMCParticleHelper::IsNeutrino(pParentMCParticle))
                            {
                                const CartesianVector signalHit(x, 0.f, z);
                                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &signalHit, "true signal", YELLOW, 2));
                            }

                        }
                     }
                }
                catch (StatusCodeException &e)
                {
                    std::cerr << "DlSignalAlgorithm: Warning. Couldn't find signal hits." << std::endl;
                }
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

            }
    }

    if (signalCandidatesU.empty() || signalCandidatesV.empty() || signalCandidatesW.empty() || signalCandidates2D.empty())
    {
        std::cout << "Error: A CaloHitList is empty" << std::endl;
    }
    if (!signalCandidatesU.empty())                                                                                        
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, signalCandidatesU, m_signalListNameU));
    
    if (!signalCandidatesV.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, signalCandidatesV, m_signalListNameV));
    
    if (!signalCandidatesW.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, signalCandidatesW, m_signalListNameW));
    
    if (!signalCandidates2D.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, signalCandidates2D, m_signalListName2D));

    if (!backgroundCaloHitList.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, backgroundCaloHitList, m_backgroundListName));

    return STATUS_CODE_SUCCESS;

}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSignalAlgorithm::MakeNetworkInputFromHits(const CaloHitList &caloHits, const HitType view, const float xMin,
    const float xMax, const float zMin, const float zMax, LArDLHelper::TorchInput &networkInput, PixelMap &pixelMap) const
{
    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());

    // Determine the bin edges
    std::vector<double> xBinEdges(m_width + 1);
    std::vector<double> zBinEdges(m_height + 1);
    xBinEdges[0] = xMin - 0.5f * m_driftStep;
    const double dx = ((xMax + 0.5f * m_driftStep) - xBinEdges[0]) / m_width;
    for (int i = 1; i < m_width + 1; ++i)
        xBinEdges[i] = xBinEdges[i - 1] + dx;
    zBinEdges[0] = zMin - 0.5f * pitch;
    const double dz = ((zMax + 0.5f * pitch) - zBinEdges[0]) / m_height;
    for (int i = 1; i < m_height + 1; ++i)
        zBinEdges[i] = zBinEdges[i - 1] + dz;

    LArDLHelper::InitialiseInput({1, 1, m_height, m_width}, networkInput);
    auto accessor = networkInput.accessor<float, 4>();

    for (const CaloHit *pCaloHit : caloHits)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        if (m_pass > 1)
        {
            if (x < xMin || x > xMax || z < zMin || z > zMax)
                continue;
        }
        const float adc{pCaloHit->GetMipEquivalentEnergy()};
        const int pixelX{static_cast<int>(std::floor((x - xBinEdges[0]) / dx))};
        const int pixelZ{static_cast<int>(std::floor((z - zBinEdges[0]) / dz))};
        accessor[0][0][pixelZ][pixelX] += adc;
	pixelMap[pCaloHit] = std::make_pair(pixelX, pixelZ);
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSignalAlgorithm::GetMCToHitsMap(LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    const CaloHitList *pCaloHitList2D(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName2D, pCaloHitList2D));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    if (pMCParticleList->empty() || pCaloHitList2D->empty())
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_minPrimaryGoodHits = 3;
    parameters.m_minHitsForGoodView = 2;
    parameters.m_minPrimaryGoodViews = 2;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    LArMCParticleHelper::SelectReconstructableMCParticles(
        pMCParticleList, pCaloHitList2D, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcToHitsMap);

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSignalAlgorithm::CompleteMCHierarchy(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, MCParticleList &mcHierarchy) const
{
    try
    {
        for (const auto &[mc, hits] : mcToHitsMap)
        {
            (void)hits;
            mcHierarchy.push_back(mc);
            LArMCParticleHelper::GetAllAncestorMCParticles(mc, mcHierarchy);
        }
    }
    catch (const StatusCodeException &e)
    {
        return e.GetStatusCode();
    }

    // Move the neutrino to the front of the list
    auto pivot =
        std::find_if(mcHierarchy.begin(), mcHierarchy.end(), [](const MCParticle *mc) -> bool { return LArMCParticleHelper::IsNeutrino(mc); });
    (void)pivot;
    if (pivot != mcHierarchy.end())
        std::rotate(mcHierarchy.begin(), pivot, std::next(pivot));
    else
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void DlSignalAlgorithm::GetHitRegion(const CaloHitList &caloHitList, float &xMin, float &xMax, float &zMin, float &zMax) const
{
    if (m_pass == 1)
    {
        xMin = std::numeric_limits<float>::max();
        xMax = -std::numeric_limits<float>::max();
        zMin = std::numeric_limits<float>::max();
        zMax = -std::numeric_limits<float>::max();
    }
    // Find the range of x and z values in the view
    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        xMin = std::min(x, xMin);
        xMax = std::max(x, xMax);
        zMin = std::min(z, zMin);
        zMax = std::max(z, zMax);
    }


    //std::cout << "Returning Image Span | Pre Scaling | x | z | " << xMax - xMin << " | " << zMax - zMin << std::endl; 

    const HitType view{caloHitList.front()->GetHitType()};
    const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
    if (!(isU || isV || isW))
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());

    if (m_simpleZoom == false && m_pass > 1)
    {
	float xCount{0.f}, zCount{0.f};
	float nCount{0.f};

	for (const CaloHit *pCaloHit : caloHitList)
        {
            const float x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()};
            xCount += x;
            zCount += z;
	    ++nCount;
	}
	if (nCount == 0)
	    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        const CartesianVector &centre{ xCount / nCount, 0.f, zCount / nCount};

        // Get hit distribution left/right asymmetry
        int nHitsLeft{0}, nHitsRight{0};

        // Get hit distribution upstream/downstream asymmetry
        int nHitsUpstream{0}, nHitsDownstream{0};

	const float xCtr{centre.GetX()};
        const float zCtr{centre.GetZ()};

	for (const CaloHit *pCaloHit : caloHitList)
        {
            const float x{pCaloHit->GetPositionVector().GetX()}, z{pCaloHit->GetPositionVector().GetZ()};
            if (x <= xCtr)
                ++nHitsLeft;
            else
                ++nHitsRight;

	    if (z <= zCtr)
                ++nHitsUpstream;
            else
                ++nHitsDownstream;
	}

        const int nHitsTotal{nHitsLeft + nHitsRight};
        if (nHitsTotal == 0)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        const float xAsymmetry{nHitsLeft / static_cast<float>(nHitsTotal)};

        const int nHitsViewTotal{nHitsUpstream + nHitsDownstream};
        if (nHitsViewTotal == 0)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        const float zAsymmetry{nHitsUpstream / static_cast<float>(nHitsViewTotal)};

        const float xSpan{m_driftStep * (m_width - 1)};
        xMin = xCtr - xAsymmetry * xSpan;
        xMax = xMin + (m_driftStep * (m_width - 1));
        const float zSpan{pitch * (m_height - 1)};
        zMin = zCtr - zAsymmetry * zSpan;
        zMax = zMin + zSpan;
    }

    if (m_simpleZoom == true && m_pass > 1 )
    {
        float xPos{-std::numeric_limits<float>::max()}, zPos{-std::numeric_limits<float>::max()}, adcMax{0.f};
        int nCount{0}, tCount{0}, hCount{0}; 
        for (const CaloHit *pCaloHit : caloHitList)
        {
            const float xC{pCaloHit->GetPositionVector().GetX()}, zC{pCaloHit->GetPositionVector().GetZ()}, adc{pCaloHit->GetMipEquivalentEnergy()};; 
	    tCount += 1;
	    if (xC >= xMin && xC < xMax)
	    {
                nCount += 1;
                if (adc > adcMax)
	        {
                    hCount += 1;
		    adcMax = adc;
   		    xPos = xC;
		    zPos = zC;
                }
	    }
        }
	std::cout << "Total Hit within Drift limits " << nCount << "Total Hits Searched " << tCount << std::endl;
	std::cout << "Hits with greater Adc: " << hCount << std::endl;
        if (nCount == 0)
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);

	const float xSpan{m_driftStep * (m_width - 1)};
        xMin = xPos - xSpan;
        xMax = xPos + xSpan;
	const float zSpan{pitch * (m_height - 1)};
        zMin = zPos - zSpan;
        zMax = zPos + zSpan;
	std::cout << "It Worked! " << xPos << " " << zPos << "Max ADC: " << adcMax<< std::endl;
    }

    
    // Avoid unreasonable rescaling of very small hit regions, pixels are assumed to be 0.5cm in x and wire pitch in z
    // ATTN: Rescaling is to a size 1 pixel smaller than the intended image to ensure all hits fit within an imaged binned
    // to be one pixel wider than this
    const float xRange{xMax - xMin}, zRange{zMax - zMin};
    const float minXSpan{m_driftStep * (m_width - 1)};
    if (xRange < minXSpan)
    {
        const float padding{0.5f * (minXSpan - xRange)};
        xMin -= padding;
        xMax += padding;
    }
    const float minZSpan{pitch * (m_height - 1)};
    if (zRange < minZSpan)
    {
        const float padding{0.5f * (minZSpan - zRange)};
        zMin -= padding;
        zMax += padding;
    }
 
    //std::cout << "Returning Image Span | Post Scaling | x | z | " << xMax - xMin << " | " << zMax - zMin << std::endl;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlSignalAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Pass", m_pass));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageHeight", m_height));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ImageWidth", m_width));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DriftStep", m_driftStep));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SignalListNameU", m_signalListNameU));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SignalListNameV", m_signalListNameV));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SignalListNameW", m_signalListNameW));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SignalListName2D", m_signalListName2D));

    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
    }
    else
    {
        std::string modelName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameU", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_modelU);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameV", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_modelV);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameW", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_modelW);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));
        if (m_writeTree)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
        }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputCaloHitListNames", m_inputCaloHitListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName2D", m_caloHitListName2D));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PassOneTrustThreshold", m_passOneTrustThreshold));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PrintOut", m_printOut));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "BackgroundListName", m_backgroundListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ApplyCheatedSeparation", m_applyCheatedSeparation));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content
