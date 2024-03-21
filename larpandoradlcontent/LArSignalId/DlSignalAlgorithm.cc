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
    m_nClasses{0},
    m_height{256},
    m_width{256},
    m_driftStep{0.5f},
    m_visualise{false},
    m_writeTree{false},
    m_rng(static_cast<std::mt19937::result_type>(std::chrono::high_resolution_clock::now().time_since_epoch().count())),
    m_volumeType{"dune_fd_hd"},
    m_signalListNameU{""},
    m_signalListNameV{""},
    m_signalListNameW{""},
    m_signalListName2D{""}
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
    
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();

    return STATUS_CODE_SUCCESS;
}

StatusCode DlSignalAlgorithm::PrepareTrainingSample()
{
    LArMCParticleHelper::MCContributionMap mcToHitsMap;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->GetMCToHitsMap(mcToHitsMap));
    MCParticleList hierarchy;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->CompleteMCHierarchy(mcToHitsMap, hierarchy));

    // Get boundaries for hits and make x dimension common
    std::map<HitType, float> wireMin, wireMax;
    float driftMin{std::numeric_limits<float>::max()}, driftMax{-std::numeric_limits<float>::max()};
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        float viewDriftMin{driftMin}, viewDriftMax{driftMax};
        this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
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
                    featureVector.emplace_back(1);
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

    std::map<HitType, float> wireMin, wireMax;
    float driftMin{std::numeric_limits<float>::max()}, driftMax{-std::numeric_limits<float>::max()};
    for (const std::string &listname : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
        if (pCaloHitList->empty())
            continue;

        HitType view{pCaloHitList->front()->GetHitType()};
        float viewDriftMin{driftMin}, viewDriftMax{driftMax};
        this->GetHitRegion(*pCaloHitList, viewDriftMin, viewDriftMax, wireMin[view], wireMax[view]);
        driftMin = std::min(viewDriftMin, driftMin);
        driftMax = std::max(viewDriftMax, driftMax);
    }

    CartesianPointVector signalCandidatesPositionU, signalCandidatesPositionV, signalCandidatesPositionW;
    for (const std::string &listName : m_caloHitListNames)
    {
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listName, pCaloHitList));

        HitType view{pCaloHitList->front()->GetHitType()};
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
        CaloHitList signalCandidatesU, signalCandidatesV, signalCandidatesW, signalCandidates2D;
	
	for (const auto &[pCaloHit, pixel] : pixelMap)
        {
	    //The ordering of the pixel is x coordinate, z coordinate 
            const auto cls{classesAccessor[0][pixel.second][pixel.first]};
            if (cls == 2)
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
	       //add background list
	    }


        }

	PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, signalCandidatesU, m_signalListNameU));
	PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, signalCandidatesV, m_signalListNameV));
	PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, signalCandidatesW, m_signalListNameW));
	PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, signalCandidates2D, m_signalListName2D));

	if (m_visualise)
        {
            if (isU)
	    {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                                         &signalCandidatesU, "candidate signal U", GREEN));
            }

            if (isV)
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                                         &signalCandidatesV, "candidate signal V", GREEN));
	    }

            if (isW)
            {
                PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(),
                                         &signalCandidatesU, "candidate signal U", GREEN));
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
                             PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &signalHit, "true signal", BLUE, 3));
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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    LArMCParticleHelper::PrimaryParameters parameters;
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
    xMin = std::numeric_limits<float>::max();
    xMax = -std::numeric_limits<float>::max();
    zMin = std::numeric_limits<float>::max();
    zMax = -std::numeric_limits<float>::max();
    // Find the range of x and z values in the view
    for (const CaloHit *pCaloHit : caloHitList)
    {
        const float x{pCaloHit->GetPositionVector().GetX()};
        const float z{pCaloHit->GetPositionVector().GetZ()};
        xMin = std::min(x, xMin);
        xMax = std::max(x, xMax);
        zMin = std::min(z, zMin);
        zMax = std::max(z, zMax);
    }

    if (caloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const HitType view{caloHitList.front()->GetHitType()};
    const bool isU{view == TPC_VIEW_U}, isV{view == TPC_VIEW_V}, isW{view == TPC_VIEW_W};
    if (!(isU || isV || isW))
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // ATTN If wire w pitches vary between TPCs, exception will be raised in initialisation of lar pseudolayer plugin
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());

   // if (m_pass > 1)
   // {
   //     const VertexList *pVertexList(nullptr);
   //     PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputVertexListName, pVertexList));
   //     if (pVertexList->empty())
   //         throw StatusCodeException(STATUS_CODE_NOT_FOUND);
   //     const CartesianVector &vertex{pVertexList->front()->GetPosition()};

        // Get hit distribution left/right asymmetry
   //     int nHitsLeft{0}, nHitsRight{0};
   //     const double xVtx{vertex.GetX()};
   //     for (const std::string &listname : m_caloHitListNames)
   //     {
   //         const CaloHitList *pCaloHitList(nullptr);
   //         PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, listname, pCaloHitList));
   //         if (pCaloHitList->empty())
   //             continue;
   //         for (const CaloHit *const pCaloHit : *pCaloHitList)
   //         {
   //             const CartesianVector &pos{pCaloHit->GetPositionVector()};
   //             if (pos.GetX() <= xVtx)
   //                 ++nHitsLeft;
   //             else
   //                 ++nHitsRight;
   //         }
   //     }
   //     const int nHitsTotal{nHitsLeft + nHitsRight};
   //     if (nHitsTotal == 0)
   //         throw StatusCodeException(STATUS_CODE_NOT_FOUND);
   //     const float xAsymmetry{nHitsLeft / static_cast<float>(nHitsTotal)};

        // Vertices
   //     const LArTransformationPlugin *transform{this->GetPandora().GetPlugins()->GetLArTransformationPlugin()};
   //     double zVtx{0.};
   //     if (isW)
   //         zVtx += transform->YZtoW(vertex.GetY(), vertex.GetZ());
   //     else if (isV)
   //         zVtx += transform->YZtoV(vertex.GetY(), vertex.GetZ());
   //     else
   //         zVtx = transform->YZtoU(vertex.GetY(), vertex.GetZ());

        // Get hit distribution upstream/downstream asymmetry
   //     int nHitsUpstream{0}, nHitsDownstream{0};
   //     for (const CaloHit *const pCaloHit : caloHitList)
   //     {
   //         const CartesianVector &pos{pCaloHit->GetPositionVector()};
   //         if (pos.GetZ() <= zVtx)
   //             ++nHitsUpstream;
   //         else
   //             ++nHitsDownstream;
   //     }
   //     const int nHitsViewTotal{nHitsUpstream + nHitsDownstream};
   //     if (nHitsViewTotal == 0)
   //         throw StatusCodeException(STATUS_CODE_NOT_FOUND);
   //     const float zAsymmetry{nHitsUpstream / static_cast<float>(nHitsViewTotal)};

   //     const float xSpan{m_driftStep * (m_width - 1)};
   //     xMin = xVtx - xAsymmetry * xSpan;
   //     xMax = xMin + (m_driftStep * (m_width - 1));
   //     const float zSpan{pitch * (m_height - 1)};
   //     zMin = zVtx - zAsymmetry * zSpan;
   //     zMax = zMin + zSpan;
   // }

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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "DistanceThresholds", m_thresholds));
    m_nClasses = m_thresholds.size() - 1;
   // if (m_pass > 1)
   // {
   //     PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputSignalListNames", m_inputSignalListNames));
   // }

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
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VolumeType", m_volumeType));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content
