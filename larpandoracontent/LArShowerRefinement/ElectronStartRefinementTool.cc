/**
 *  @file   larpandoracontent/LArShowerRefinement/ElectronStartRefinementTool.cc
 *
 *  @brief  Implementation of the gamma start refinement tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "PandoraMonitoringApi.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArShowerRefinement/ElectronStartRefinementTool.h"

using namespace pandora;

int colourCount = 0;
std::vector<Color> colours = {BLUE, GREEN, VIOLET, DARKGREEN, DARKYELLOW, DARKVIOLET, DARKCYAN, DARKRED};


namespace lar_content
{

ElectronStartRefinementTool::ElectronStartRefinementTool() : 
    m_showerSlidingFitWindow(20),
    m_maxCoincideneTransverseSeparation(2.f),
    m_minSpinePurity(0.7f),
    m_maxAngularDeviation(5.f),
    m_maxXSeparation(5.f),
    m_maxSeparation(2.f),
    m_extendMode(true),
    m_moveVertexMode(true),
    m_truncateGammaMode(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

ElectronStartRefinementTool::~ElectronStartRefinementTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const CartesianVector &nuVertexPosition, const CaloHitList *const pCaloHitListU, const CaloHitList *const pCaloHitListV, const CaloHitList *const pCaloHitListW)
{

    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);


    // Only apply electron refinement algorithm to showers
    if (!LArPfoHelper::IsShower(pShowerPfo))
        return false;

    // Only consider significant showers
    CaloHitList caloHits3D;
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHits3D);

    if (caloHits3D.size() < 50)
        return false;

    CaloHitList usedHitListU, usedHitListV, usedHitListW;
    ElectronProtoShowerVector protoShowerVectorU, protoShowerVectorV, protoShowerVectorW;
    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_U, protoShowerVectorU, usedHitListU);
    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_V, protoShowerVectorV, usedHitListV);
    this->BuildProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, TPC_VIEW_W, protoShowerVectorW, usedHitListW);

    //bool isGamma(false);
    //float highestNHits(0.f);
    //float highestOpeningAngle(-std::numeric_limits<float>::max()); // tie-breaker
    //float showerVertexX(0.f),showerVertexY(0.f), showerVertexZ(0.f);

    int count(0);

    for (ElectronProtoShower &protoShowerU : protoShowerVectorU)
    {
        for (ElectronProtoShower &protoShowerV : protoShowerVectorV)
        {
            for (ElectronProtoShower &protoShowerW : protoShowerVectorW)
            {
                LArConnectionPathwayHelper::Consistency consistency(LArConnectionPathwayHelper::Consistency::POSITION);

                if (!this->ArePathwaysConsistent(pAlgorithm, nuVertexPosition, protoShowerU, protoShowerV, protoShowerW, consistency))
                {
                    //////////////////////////
                    //std::cout << "PATHWAYS ARE NOT CONSISTENT" << std::endl;
                    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
                    //////////////////////////
                    continue;
                }
                else
                {
                    ++count;
                    //////////////////////////////
                    /*
                    const CartesianVector projectedNuVertexPositionU(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, TPC_VIEW_U));
                    const CartesianVector projectedNuVertexPositionV(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, TPC_VIEW_V));
                    const CartesianVector projectedNuVertexPositionW(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, TPC_VIEW_W));

                    PfoList visualize;
                    visualize.push_back(pShowerPfo);
                    PandoraMonitoringApi::VisualizeParticleFlowObjects(pAlgorithm->GetPandora(), &visualize, "ShowerPfo", BLUE);
                    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerU.m_connectionPathway.m_pathwayHitList, "Shower Pathway U", BLACK);
                    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerU.m_spineHitList, "Shower Spine U", RED);
                    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &protoShowerU.m_showerCore.m_startPosition, "Shower start position U", BLACK, 2);
                    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerV.m_connectionPathway.m_pathwayHitList, "Shower Pathway V", BLACK);
                    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerV.m_spineHitList, "Shower Spine V", RED);
                    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &protoShowerV.m_showerCore.m_startPosition, "Shower start position V", BLACK, 2);
                    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerW.m_connectionPathway.m_pathwayHitList, "Shower Pathway W", BLACK);
                    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerW.m_spineHitList, "Shower Spine W", RED);
                    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &protoShowerW.m_showerCore.m_startPosition, "Shower start position W", BLACK, 2);
                    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPositionU, "Nu vertex position", GREEN, 2);
                    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPositionV, "Nu vertex position", GREEN, 2);
                    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPositionW, "Nu vertex position", GREEN, 2);
                    */
                    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
                    /////////////////////////////

                    LArConnectionPathwayHelper::ElectronTreeVariables electronTreeVariables;

                    LArConnectionPathwayHelper::FillElectronTreeVariables(pAlgorithm, pShowerPfo, protoShowerU, protoShowerV, protoShowerW, nuVertexPosition,
                        pCaloHitListU, pCaloHitListV, pCaloHitListW, consistency, electronTreeVariables);

                    /////////////////////////////
                    /*
                    std::cout << "electronTreeVariables.m_maxNPostShowerStartHits: " << electronTreeVariables.m_maxNPostShowerStartHits << std::endl;
                    std::cout << "electronTreeVariables.m_maxPostShowerStartOpeningAngle: " << electronTreeVariables.m_maxPostShowerStartOpeningAngle << std::endl;
                    std::cout << "electronTreeVariables.m_minPostShowerStartShowerStartMoliereRadius: " << electronTreeVariables.m_minPostShowerStartShowerStartMoliereRadius << std::endl;
                    std::cout << "electronTreeVariables.m_maxPostShowerStartShowerStartEnergyAsymmetry: " << electronTreeVariables.m_maxPostShowerStartShowerStartEnergyAsymmetry << std::endl;

                    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
                    */
                    /////////////////////////////

                    if (pAlgorithm->m_createTrainingTrees)
                    {
                        std::cout << "IN m_createTrainingTrees SECTION" << std::endl;
                        std::string treeString(pAlgorithm->IsElectron(pShowerPfo) ? "ElectronSignalTree" : "ElectronBackgroundTree");
                        pAlgorithm->FillTree(treeString, electronTreeVariables);
                        return false;
                    }

                    if (count == 1)
                        pAlgorithm->SetElectronTreeMetadata(pShowerPfo, electronTreeVariables);


                    return true;
                }
            }
        }
    }

                    /*
                    if (pAlgorithm->m_hybridMode)
                    {
                        std::cout << "IN HYBRID MODE SECTION" << std::endl;

                        if (!pAlgorithm->IsElectron(pShowerPfo))
                            continue;
                    }
                    else
                    {                    
                        std::cout << "IN REAL RECO MODE!!" << std::endl;
                        std::cout << "m_electronTMVACut: " << pAlgorithm->m_electronTMVACut << std::endl;

                        if (!pAlgorithm->TMVAIsElectron(electronTreeVariables, pShowerPfo, true))
                        {
                            if (pAlgorithm->TMVAIsGamma(electronTreeVariables, pShowerPfo) && this->IsSensibleShower(electronTreeVariables))
                            {
                                if ((std::fabs(electronTreeVariables.m_maxNPostShowerStartHits - highestNHits) == std::numeric_limits<float>::epsilon() && 
                                     (electronTreeVariables.m_maxPostShowerStartOpeningAngle > highestOpeningAngle)) || (electronTreeVariables.m_maxNPostShowerStartHits > highestNHits))
                                {
                                    isGamma = true;
                                    highestNHits = electronTreeVariables.m_maxNPostShowerStartHits;
                                    highestOpeningAngle = electronTreeVariables.m_maxPostShowerStartOpeningAngle;
                                    showerVertexX = electronTreeVariables.m_showerStartX;
                                    showerVertexY = electronTreeVariables.m_showerStartY;
                                    showerVertexZ = electronTreeVariables.m_showerStartZ;
                                }
                            }

                            continue;
                        }
                    }

                    // This is an electron so do electron things...

                    // if electron make sure this is the metadata saved
                    pAlgorithm->SetElectronTreeMetadata(pShowerPfo, electronTreeVariables);

                    if (m_extendMode)
                    {
                        if (this->IsShowerExtendable(pAlgorithm, protoShowerU, protoShowerV, protoShowerW))
                            this->ExtendShower(pAlgorithm, pShowerPfo, protoShowerU, protoShowerV, protoShowerW);
                    }

                    // Set vertex metadata parameters... Only referred to in the shift configuration
                    if (m_moveVertexMode)
                        pAlgorithm->SetElectronMetadata(nuVertexPosition, pShowerPfo);

                    return true;
                }
            }
        }
    }

    // Now remove if instead it looks to be a gamma
    if (m_truncateGammaMode && isGamma)
    {
        //ONLY MOVE VERTEX IF IT LOOKS LIKE THE VERTEX IS FURTHER AWAY FROM THE NU VERTEX THEN ORIGINAL SHOWER VERTEX?? AND IF IT IS CONNECTED I.E. DO WE NEED TO DO THIS?

        std::cout << "THIS IS A GAMMA" << std::endl;
        std::cout << "showerVertexX: " << showerVertexX << std::endl;
        std::cout << "showerVertexY: " << showerVertexY << std::endl;
        std::cout << "showerVertexZ: " << showerVertexZ << std::endl;

        ClusterList clusterList;
        LArPfoHelper::GetClusters(pShowerPfo, TPC_3D, clusterList);

        const float currentSeparation(LArClusterHelper::GetClosestDistance(nuVertexPosition, clusterList));

        const CartesianVector newVertex(showerVertexX, showerVertexY, showerVertexZ);
        const float newSeparation((newVertex - nuVertexPosition).GetMagnitude());

        if (newSeparation > currentSeparation)
            pAlgorithm->SetGammaVertex(newVertex, pShowerPfo);
    }
                    */

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::IsSensibleShower(const LArConnectionPathwayHelper::ElectronTreeVariables &electronTreeVariables)
{
    // Hits
    if (electronTreeVariables.m_maxNPostShowerStartHits < 200.f)
        return false;

    // Scattering Angle
    if (electronTreeVariables.m_maxPostShowerStartOpeningAngle < 2.f)
        return false;

    // Moliere Radius
    if (electronTreeVariables.m_minPostShowerStartShowerStartMoliereRadius < 2.f)
        return false;

    // Shower Energy Asymmetry
    if (electronTreeVariables.m_maxPostShowerStartShowerStartEnergyAsymmetry < 0.f)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::BuildProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertexPosition, const HitType hitType, ElectronProtoShowerVector &protoShowerVector, CaloHitList &usedHitList)
{
    /////////////////////////////////
    std::cout << "Investigating " << (hitType == TPC_VIEW_U ? "U" : hitType == TPC_VIEW_V ? "V" : "W") << " view..." << std::endl;
    /////////////////////////////////

    CaloHitList viewShowerHitList;
    const CaloHitList viewAllHitList(pAlgorithm->GetAllHitsOfType(hitType));
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

    const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertexPosition, hitType));

    CaloHitList viewROIHits;
    this->CollectHitsWithinROI(pAlgorithm, viewShowerHitList, viewAllHitList, projectedNuVertexPosition, viewROIHits);

    if (viewROIHits.empty())
        return;

    // Fill angular decomposition map
    AngularDecompositionMap angularDecompositionMap;
    this->FillAngularDecompositionMap(viewROIHits, projectedNuVertexPosition, angularDecompositionMap);

    if (angularDecompositionMap.empty())
        return;

    // Smooth angular decomposition map
    this->SmoothAngularDecompositionMap(angularDecompositionMap);

    // Find peak directions
    IntVector angularPeakVector;
    this->ObtainPeakVector(angularDecompositionMap, angularPeakVector);

    // Now, for future use, find the shower's vertex <- if layer positions span the vertex then put vertex as nu vertex?? <- might mean we are able to shift the vertex?
    CartesianVector showerVertexPosition(0.f, 0.f, 0.f);
    try
    {
        showerVertexPosition = this->GetShowerVertex(pShowerPfo, hitType, projectedNuVertexPosition);
    }
    catch(...)
    {
        return;
    }

    //////////////////////////////////////////
    // Draw neutrino vertex
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &projectedNuVertexPosition, "Nu Vertex", BLACK, 2);

    // Draw shower and shower vertex
    PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &viewShowerHitList, "Shower Hits", RED);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &showerVertexPosition, "Shower Vertex", RED, 2);

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    // Draw peak vector directions
    for (unsigned int i : angularPeakVector)
    {
        const float theta0XZ(static_cast<float>(i) * m_theta0XZBinSize);
        const CartesianVector peakDirection(std::cos(theta0XZ), 0.f, std::sin(theta0XZ));
        CartesianVector end(projectedNuVertexPosition + (peakDirection * 10.f));
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &projectedNuVertexPosition, &end, "PEAK DIRECTION", BLACK, 2, 1);
    }

    // Draw spine 'ambiguous' region
    //const CartesianVector searchLength(projectedNuVertexPosition + CartesianVector(m_growingFitInitialLength, 0.0, 0.0));
    //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &projectedNuVertexPosition, &searchLength, "Search Length", RED, 2, 1);

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    //////////////////////////////////////////



    // Investigate peaks
    IntVector investigatedPeaks;

    for (unsigned int i = 0; i < angularPeakVector.size(); ++i)
    {
        int bestTheta0XZBin(0);

        if (!this->FindBestAngularPeak(angularDecompositionMap, angularPeakVector, investigatedPeaks, bestTheta0XZBin))
            break;

        investigatedPeaks.push_back(bestTheta0XZBin);

        const float theta0XZ(bestTheta0XZBin * m_theta0XZBinSize);
        const CartesianVector peakDirection(std::cos(theta0XZ), 0.f, std::sin(theta0XZ));
        const bool isEndDownstream(peakDirection.GetZ() > 0.f);

        ///////////////////////////////////////////
        //const CartesianVector end(projectedNuVertexPosition + (peakDirection * 10.f));
        //PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, &end, "direction", BLACK, 2,2);
        ///////////////////////////////////////////

        // Find the shower spine
        CaloHitList showerSpineHitList;
        try
        {
            this->FindShowerSpine(pAlgorithm, viewAllHitList, projectedNuVertexPosition, peakDirection, usedHitList, showerSpineHitList);
        }
        catch (...)
        {
            ///////////////////////////////////////////  
            //std::cout << "SHOWER SPINE FIT FAILED" << std::endl;
            //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            ///////////////////////////////////////////  
            continue;
        }

        // Demand that spine is significant, be lenient here as some have small stubs and a gap
        if (showerSpineHitList.size() < 7) 
        {
            ///////////////////////////////////////////  
            //std::cout << "SHOWER NOT SIGNIFICANT" << std::endl;
            //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            ///////////////////////////////////////////  
            continue;
        }


        /////////////////////////////////////////////
        // Draw shower spine
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &showerSpineHitList, "shower spine", colours[colourCount]);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        /////////////////////////////////////////////


        // Does spine coincide with shower pfo
        if (!this->IsSpineCoincident(pAlgorithm, projectedNuVertexPosition, peakDirection, viewShowerHitList, showerVertexPosition, showerSpineHitList))
        {
            ///////////////////////////////////////////  
            //std::cout << "not coincident" << std::endl;
            //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            ///////////////////////////////////////////  
            continue;
        }

        // Now look for the shower start position, first obtain the longitudinal position of spine hits
        LongitudinalPositionMap longitudinalPositionMap;
        try
        {
            this->ObtainLongitudinalDecomposition(pAlgorithm, showerSpineHitList, longitudinalPositionMap);
        }
        catch(...)
        {
            /////////////////////////////////
            //std::cout << "the super fine fit failed probably a gap in the initial hits??" << std::endl;
            //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
            /////////////////////////////////
            continue;
        }

        // Obtain spine energy profile
        EnergySpectrumMap energySpectrumMap;
        this->GetEnergyDistribution(pAlgorithm, showerSpineHitList, longitudinalPositionMap, energySpectrumMap);

        // Postulate combining the two
        CaloHitList mergedHitList(viewShowerHitList);

        for (const CaloHit *const pCaloHit : showerSpineHitList)
        {
            if (std::find(viewShowerHitList.begin(), viewShowerHitList.end(), pCaloHit) == viewShowerHitList.end())
                mergedHitList.push_back(pCaloHit);
        }

        // I don't like the ProtoShowers being created here (do this so that i can push the vertex forwards if i am sure we have an electron!)
        ProtoShowerVector protoShowerVectorTemp;
        CartesianVector showerStartPosition(0.f, 0.f, 0.f);

        this->FindShowerStart(pAlgorithm, projectedNuVertexPosition, peakDirection, longitudinalPositionMap, energySpectrumMap, showerSpineHitList, 
             showerStartPosition, mergedHitList, isEndDownstream, protoShowerVectorTemp, false); 

        // add in the shower core bit (TO DO -  so then i can check??)
        CaloHitList hitsToAdd;
        const float showerVertexL(std::max(peakDirection.GetDotProduct(showerStartPosition - projectedNuVertexPosition), 
            peakDirection.GetDotProduct(showerVertexPosition - projectedNuVertexPosition)));

        for (const CaloHit *const pCaloHit : showerSpineHitList)
        {
            if (std::find(viewShowerHitList.begin(), viewShowerHitList.end(), pCaloHit) != viewShowerHitList.end())
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const float showerL(peakDirection.GetDotProduct(hitPosition - projectedNuVertexPosition));

            if ((showerL > 0.f) && (showerL < showerVertexL))
                hitsToAdd.push_back(pCaloHit);
        }

        usedHitList.insert(usedHitList.begin(), showerSpineHitList.begin(), showerSpineHitList.end());

        protoShowerVector.push_back(ElectronProtoShower(protoShowerVectorTemp.front().m_showerCore, protoShowerVectorTemp.front().m_connectionPathway, showerSpineHitList, 
            false, CaloHitList(), CartesianPointVector(), hitsToAdd));


        /////////////////////////////////////////////////////////
        // now identify the ambiguous hits in protoshower
        CartesianPointVector peakDirections;
        CaloHitList usedHitListForHelper(showerSpineHitList);
        this->BuildHelperProtoShowers(pAlgorithm, pShowerPfo, nuVertexPosition, hitType, peakDirections, usedHitListForHelper);

        this->RefineHitsToAdd(pAlgorithm, protoShowerVector.back(), nuVertexPosition, hitType, peakDirections);

        /////////////////////////////////
        /*
        PfoList visualize;
        visualize.push_back(pShowerPfo);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(pAlgorithm->GetPandora(), &visualize, "ShowerPfo", BLUE);
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &showerSpineHitList, "Shower Spine", BLACK);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "Shower start position", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, "Nu vertex position", GREEN, 2);

        for (const CaloHit *const pHitToAdd : hitsToAdd)
        {
            const CartesianVector &hitPosition(pHitToAdd->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "HIT TO ADD", RED, 2);
        }
        */
        /////////////////////////////////
    }

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::BuildHelperProtoShowers(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertexPosition, const HitType tpcView, CartesianPointVector &significantPeakDirections, CaloHitList &usedHitList)
{
    CaloHitList showerHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, tpcView, showerHitList);

    CaloHitList viewHitList;
    const CaloHitList allHitList(pAlgorithm->GetAllHitsOfType(tpcView));

    for (const CaloHit *const pCaloHit : allHitList)
    {
        if (std::find(showerHitList.begin(), showerHitList.end(), pCaloHit) != showerHitList.end())
            continue;

        viewHitList.push_back(pCaloHit);
    }

    const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, tpcView));

    AngularDecompositionMap angularDecompositionMap;
    this->FillAngularDecompositionMap(viewHitList, projectedNuVertexPosition, angularDecompositionMap);

    if (angularDecompositionMap.empty())
        return;

    this->SmoothAngularDecompositionMap(angularDecompositionMap);

    // Obtain distribution peaks
    IntVector angularPeakVector;
    this->ObtainPeakVector(angularDecompositionMap, angularPeakVector);

    // Investigate peaks, searching for pathways (peaks indexed by theta0XZ bin lower edge)
    IntVector investigatedPeaks;

    for (unsigned int i = 0; i < angularPeakVector.size(); ++i)
    {
        int bestTheta0XZBin(0);

        if (!this->FindBestAngularPeak(angularDecompositionMap, angularPeakVector, investigatedPeaks, bestTheta0XZBin))
            break;

        investigatedPeaks.push_back(bestTheta0XZBin);

        const float theta0XZ(bestTheta0XZBin * m_theta0XZBinSize);
        const CartesianVector peakDirection(std::cos(theta0XZ), 0.f, std::sin(theta0XZ));

        CaloHitList showerSpineHitList;
        try
        {
            this->FindShowerSpine(pAlgorithm, viewHitList, projectedNuVertexPosition, peakDirection, usedHitList, showerSpineHitList);
        }
        catch (...)
        {
            continue;
        }

        // Demand that spine is significant
        if (showerSpineHitList.size() < 4) // was 7
            continue;

        usedHitList.insert(usedHitList.begin(), showerSpineHitList.begin(), showerSpineHitList.end());

        significantPeakDirections.push_back(peakDirection);

        ///////////////////////////////////
        /*
        const CartesianVector end(projectedNuVertexPosition + (peakDirection * 10.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, &end, "peak direction", VIOLET, 2,2);
        
        for (const CaloHit *const pCaloHit : showerSpineHitList)
        {
            const CartesianVector hitPosition(pCaloHit->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "shower spine", BLACK, 2);
        }
        PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        */
        ///////////////////////////////////
    }

    ///////////////////////////////////
    //std::cout << "Peak directions..." << std::endl;
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    ///////////////////////////////////
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::AssignShowerHits(ShowerStartRefinementAlgorithm *const pAlgorithm, const HitType &hitType, ElectronProtoShowerVector &protoShowerVector)
{
    CaloHitList viewHitList;
    const CaloHitList allHitList(pAlgorithm->GetAllHitsOfType(hitType));

    CaloHitList connectionPathwayHits;
    for (const ProtoShower &protoShower : protoShowerVector)
        connectionPathwayHits.insert(connectionPathwayHits.begin(), protoShower.m_connectionPathway.m_pathwayHitList.begin(), 
            protoShower.m_connectionPathway.m_pathwayHitList.end());

    for (const CaloHit *const pCaloHit : allHitList)
    {
        if (std::find(connectionPathwayHits.begin(), connectionPathwayHits.end(), pCaloHit) != connectionPathwayHits.end())
            continue;

        int bestProtoShower(-1);
        float bestT(std::numeric_limits<float>::max());

        for (unsigned int index = 0; index < protoShowerVector.size(); ++index)
        {
            const ProtoShower &protoShower(protoShowerVector[index]);
            const CartesianVector &showerStartPosition(protoShower.m_showerCore.m_startPosition);
            const CartesianVector &showerStartDirection(protoShower.m_showerCore.m_startDirection);

            const float t(showerStartDirection.GetCrossProduct(pCaloHit->GetPositionVector() - showerStartPosition).GetMagnitude());
            const float l(showerStartDirection.GetDotProduct(pCaloHit->GetPositionVector() - showerStartPosition));

            if ((l > 0.f) && (t < m_molliereRadius) && (t < bestT))
            {
                bestT = t;
                bestProtoShower = index;
            }
        }

        if (bestProtoShower >= 0)
            protoShowerVector[bestProtoShower].m_showerCore.m_coreHitList.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::RefineHitsToAdd(ShowerStartRefinementAlgorithm *const pAlgorithm, ElectronProtoShower &protoShower,
    const CartesianVector &nuVertexPosition, const HitType hitType, const CartesianPointVector &significantPeakDirections)
{
    const CartesianVector projectedNuVertexPosition(LArGeometryHelper::ProjectPosition(pAlgorithm->GetPandora(), nuVertexPosition, hitType));
    const CartesianVector &peakDirection(protoShower.m_connectionPathway.m_startDirection);

    CaloHitList refinedHitList;

    for (const CaloHit *const pHitToAdd : protoShower.m_hitsToAdd)
    {
        bool found(false);

        const CartesianVector &hitPosition(pHitToAdd->GetPositionVector());
        const CartesianVector displacement(hitPosition - projectedNuVertexPosition);
        const float thisT((peakDirection.GetCrossProduct(displacement)).GetMagnitudeSquared());

        for (const CartesianVector &significantPeakDirection: significantPeakDirections)
        {
            if (peakDirection.GetOpeningAngle(significantPeakDirection) > (M_PI / 2))
                continue;

            const float otherT((significantPeakDirection.GetCrossProduct(displacement)).GetMagnitudeSquared());

            if ((otherT < thisT) || (otherT < 0.5f))
            {
                found = true;

                if (std::find(protoShower.m_ambiguousDirectionVector.begin(), protoShower.m_ambiguousDirectionVector.end(), significantPeakDirection) == 
                    protoShower.m_ambiguousDirectionVector.end())
                {
                    protoShower.m_ambiguousDirectionVector.push_back(significantPeakDirection);
                }
            }
        }

        found ? protoShower.m_ambiguousHitList.push_back(pHitToAdd) : refinedHitList.push_back(pHitToAdd);
    }

    //////////////////////////////
    /*
    std::cout << "BEFORE HIT REFINEMENT" << std::endl;
    for (const CaloHit *const pCaloHit :  protoShower.m_hitsToAdd)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "before hit refinement", RED, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    //////////////////////////////

    protoShower.m_hitsToAdd = this->FindContinuousPath(refinedHitList, projectedNuVertexPosition);

    //////////////////////////////
    /*
    std::cout << "AFTER HIT REFINEMENT" << std::endl;
    for (const CaloHit *const pCaloHit :  protoShower.m_hitsToAdd)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "after hit refinement", GREEN, 2);
    }
    */
    for (const CaloHit *const pCaloHit :  protoShower.m_ambiguousHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "after hit refinement", VIOLET, 2);
    }

    for (const CartesianVector &ambiguousDirection : protoShower.m_ambiguousDirectionVector)
    {
        CartesianVector jamEnd(projectedNuVertexPosition + (ambiguousDirection* 10.f));
        PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, &jamEnd, "ambiguous particle", VIOLET, 2, 1);
    }


    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
   
    //////////////////////////////
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList ElectronStartRefinementTool::FindContinuousPath(const CaloHitList &refinedHitList, const CartesianVector &projectedNuVertexPosition)
{
    CaloHitList continuousHitList;

    CaloHitVector refinedHitVector(refinedHitList.begin(), refinedHitList.end());
    std::sort(refinedHitVector.begin(), refinedHitVector.end(), LArConnectionPathwayHelper::SortByDistanceToPoint(projectedNuVertexPosition));

    unsigned int startIndex(refinedHitVector.size());

    for (unsigned int i = 0; i < refinedHitVector.size(); ++i)
    {
        CaloHitList connectedHitList;
        connectedHitList.push_back(refinedHitVector[i]);

        bool found(true);

        while(found)
        {
            found = false;

            for (unsigned int j = (i + 1); j < refinedHitVector.size(); ++j)
            {
                const CaloHit *const pCaloHit(refinedHitVector[j]);

                if (std::find(connectedHitList.begin(), connectedHitList.end(), pCaloHit) != connectedHitList.end())
                    continue;

                if (LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), connectedHitList) < 1.f)
                {
                    found = true;
                    connectedHitList.push_back(pCaloHit);
                    break;
                }
            }
        }

        if ((connectedHitList.size() >= 2) || (connectedHitList.size() == refinedHitVector.size()))
        {
            startIndex = i;
            break;
        }
    }

    for (unsigned int i = startIndex; i < refinedHitVector.size(); ++i)
        continuousHitList.push_back(refinedHitVector[i]);

    return continuousHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::CollectHitsWithinROI(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &showerHitList, 
    const CaloHitList &allHitList, const CartesianVector &projectedNuVertexPosition, CaloHitList &collectedHits)
{
    float lowestTheta(std::numeric_limits<float>::max()), highestTheta((-1.f) * std::numeric_limits<float>::max());

    this->GetAngularExtrema(pAlgorithm, showerHitList, projectedNuVertexPosition, lowestTheta, highestTheta);
    this->CollectHitsWithinExtrema(pAlgorithm, allHitList, projectedNuVertexPosition, lowestTheta, highestTheta, collectedHits);

    ///////////////////////////////
    /*
    for (const CaloHit *const pCaloHit : collectedHits)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ROI HIT", SPRING, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    ///////////////////////////////
}

//------------------------------------------------------------------------------------------------------------------------------------------
// i dont know if this is better off as a simple x thing? i.e. between nu vertex and shower edge?
void ElectronStartRefinementTool::GetAngularExtrema(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &viewHitList, 
    const CartesianVector &projectedNuVertexPosition, float &lowestTheta, float &highestTheta)
{
    lowestTheta = std::numeric_limits<float>::max();
    highestTheta = (-1.f) * std::numeric_limits<float>::max();

    const CartesianVector xAxis(1.f, 0.f, 0.f);

    for (const CaloHit *const pCaloHit : viewHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector displacementVector(hitPosition - projectedNuVertexPosition);

        float theta0XZ(xAxis.GetOpeningAngle(displacementVector));

        if (displacementVector.GetZ() < 0.f)
            theta0XZ += (3.14 - theta0XZ);

        if (theta0XZ < lowestTheta)
            lowestTheta = theta0XZ;

        if (theta0XZ > highestTheta)
            highestTheta = theta0XZ;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::CollectHitsWithinExtrema(ShowerStartRefinementAlgorithm *const pAlgorithm, const CaloHitList &viewHitList, 
    const CartesianVector &projectedNuVertexPosition, const float lowestTheta, const float highestTheta, CaloHitList &collectedHits)
{
    const CartesianVector xAxis(1.f, 0.f, 0.f);

    for (const CaloHit *const pCaloHit : viewHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        const CartesianVector displacementVector(hitPosition - projectedNuVertexPosition);

        float theta0XZ(xAxis.GetOpeningAngle(displacementVector));

        if (displacementVector.GetZ() < 0.f)
            theta0XZ += (3.14 - theta0XZ);

        if ((theta0XZ < lowestTheta) || (theta0XZ > highestTheta))
            continue;

        collectedHits.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector ElectronStartRefinementTool::GetShowerVertex(const ParticleFlowObject *const pShowerPfo, const HitType hitType, 
    const CartesianVector &projectedNuVertexPosition)
{
    const TwoDSlidingShowerFitResult twoDShowerSlidingFit(this->FitShower(pShowerPfo, hitType));

    const CartesianVector &minLayerPosition(twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMinLayerPosition());
    const CartesianVector &maxLayerPosition(twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMaxLayerPosition());

    if ((projectedNuVertexPosition.GetZ() > minLayerPosition.GetZ()) && (projectedNuVertexPosition.GetZ() < maxLayerPosition.GetZ()))
    {
        std::cout << "pfo spans nu vertex" << std::endl;
        return projectedNuVertexPosition;
    }

    const float minSeparation((projectedNuVertexPosition - minLayerPosition).GetMagnitudeSquared());
    const float maxSeparation((projectedNuVertexPosition - maxLayerPosition).GetMagnitudeSquared());
    const CartesianVector &showerVertexPosition(minSeparation < maxSeparation ? minLayerPosition : maxLayerPosition);

    return showerVertexPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoDSlidingShowerFitResult ElectronStartRefinementTool::FitShower(const ParticleFlowObject *const pShowerPfo, const HitType hitType)
{
    ClusterList viewCusterList;
    LArPfoHelper::GetClusters(pShowerPfo, hitType, viewCusterList);

    const TwoDSlidingShowerFitResult twoDShowerSlidingFit(viewCusterList.front(), m_showerSlidingFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    return twoDShowerSlidingFit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::IsSpineCoincident(ShowerStartRefinementAlgorithm *const pAlgorithm, const CartesianVector &projectedNuVertexPosition, 
    const CartesianVector &peakDirection, const CaloHitList &viewShowerHitList, const CartesianVector &showerVertexPosition, const CaloHitList &showerSpineHitList)
{
    // Check whether shower is an extension of found pathway
    CartesianVector showerVertex(showerVertexPosition);

    bool found(true);

    if (!this->IsShowerConnected(showerVertexPosition, projectedNuVertexPosition, peakDirection))
    {
        found = false;
        float minL(std::numeric_limits<float>::max());

        // sometimes the shower start can be in annoying place -.-
        for (const CaloHit *const pCaloHit : viewShowerHitList)
        {
            CartesianVector displacement(pCaloHit->GetPositionVector() - projectedNuVertexPosition);
            const float longitudinalSeparation(peakDirection.GetDotProduct(displacement));

            if ((longitudinalSeparation < (-1.f)) || (longitudinalSeparation > minL))
                continue;

            const float transverseSeparation((peakDirection.GetCrossProduct(displacement)).GetMagnitude());

            if (transverseSeparation < m_maxCoincideneTransverseSeparation)
            {
                found = true;
                showerVertex = pCaloHit->GetPositionVector();
                minL = longitudinalSeparation;
            }
        }
    }

    if (!found)
        return false;

    CaloHitList postShowerVertexSpineHits;
    const float showerDistanceFromNuSquared((showerVertex - projectedNuVertexPosition).GetMagnitudeSquared());


    for (const CaloHit * const pSpineHit : showerSpineHitList)
    {
        const CartesianVector &hitPosition(pSpineHit->GetPositionVector());
        const float separationSquared((hitPosition - projectedNuVertexPosition).GetMagnitudeSquared());

        if (separationSquared > showerDistanceFromNuSquared)
        {
            postShowerVertexSpineHits.push_back(pSpineHit);
            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "purity hit", colours[colourCount], 2);
        }
    }


    /////////////////////////////////
    
    //CartesianVector end(projectedNuVertexPosition + (peakDirection * 100.0));
    //PandoraMonitoringApi::AddLineToVisualization(pAlgorithm->GetPandora(), &projectedNuVertexPosition, &end, "DIRECTION", PINK, 2, 2);
    //PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerVertex, "SPINE VERTEX POSITION", RED, 2);
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    /////////////////////////////////


    // Check whether shower spine is pure
    const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(postShowerVertexSpineHits, viewShowerHitList));
    const float spinePurity(postShowerVertexSpineHits.size() == 0 ? 1.f : static_cast<float>(sharedHitList.size()) / static_cast<float>(postShowerVertexSpineHits.size()));

    //std::cout << "spinePurity: " << spinePurity << std::endl;

    ++colourCount;

    if (spinePurity < m_minSpinePurity)
    {
        /////////////////////////////////
        //std::cout <<  "spine purity issue" << std::endl;
        //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        /////////////////////////////////
        return false;
    }

    /////////////////////////////////
    //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    /////////////////////////////////

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::IsShowerConnected(const CartesianVector &showerVertexPosition, const CartesianVector &projectedNuVertexPosition, 
    const CartesianVector &peakDirection)
{
    CartesianVector displacement(showerVertexPosition - projectedNuVertexPosition);

    const float longitudinalSeparation(peakDirection.GetDotProduct(displacement));

    if (longitudinalSeparation < (-1.f))
        return false;

    const float transverseSeparation((peakDirection.GetCrossProduct(displacement)).GetMagnitude());

    if (transverseSeparation > m_maxCoincideneTransverseSeparation)
    {
        /////////////////////////////////
        //std::cout << "transverse separation issue" << std::endl;
        //PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
        /////////////////////////////////
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ElectronStartRefinementTool::IsShowerExtendable(ShowerStartRefinementAlgorithm *const pAlgorithm, const ElectronProtoShower &protoShowerU, 
    const ElectronProtoShower &protoShowerV, const ElectronProtoShower &protoShowerW)
{
    const CaloHitList &hitsToAddU(protoShowerU.m_hitsToAdd);
    const CaloHitList &hitsToAddV(protoShowerV.m_hitsToAdd);
    const CaloHitList &hitsToAddW(protoShowerW.m_hitsToAdd);

    // Return if no hits to add...
    if (hitsToAddU.empty() && hitsToAddV.empty() && hitsToAddW.empty())
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//DUDE what if i change this to be like, does the 3D direction match the 3D shower??? since the direction comparison is pretty tempremental? 
// could then even pick up two view issues? if that turns out to be a problem?

bool ElectronStartRefinementTool::ArePathwaysConsistent(ShowerStartRefinementAlgorithm *const pAlgorithm, const CartesianVector &nuVertexPosition,
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, LArConnectionPathwayHelper::Consistency &consistency)
{
    std::cout << "CHECKING IF SHOWERS ARE CONSISTENT" << std::endl;
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerU.m_connectionPathway.m_pathwayHitList, "Shower Pathway U", BLACK);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerU.m_spineHitList, "Shower Spine U", BLUE);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &protoShowerU.m_showerCore.m_startPosition, "Shower start position U", BLUE, 2);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerV.m_connectionPathway.m_pathwayHitList, "Shower Pathway V", BLACK);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerV.m_spineHitList, "Shower Spine V", BLUE);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &protoShowerV.m_showerCore.m_startPosition, "Shower start position V", BLUE, 2);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerW.m_connectionPathway.m_pathwayHitList, "Shower Pathway W", BLACK);
    PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &protoShowerW.m_spineHitList, "Shower Spine W", BLUE);
    PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &protoShowerW.m_showerCore.m_startPosition, "Shower start position W", BLUE, 2);
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());

    if (LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, m_maxXSeparation, m_maxSeparation))
    {
        consistency = LArConnectionPathwayHelper::Consistency::POSITION;
    }
    else if (LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, m_maxAngularDeviation))
    {
        consistency = LArConnectionPathwayHelper::Consistency::DIRECTION;
    }
    else
        return false;

    //std::cout << "isobel do not forget to take this out" << std::endl;
    //consistency = LArConnectionPathwayHelper::Consistency::X_PROJECTION;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::ExtendShower(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo,
    const ElectronProtoShower &protoShowerU, const ElectronProtoShower &protoShowerV, const ElectronProtoShower &protoShowerW)
{
    const CaloHitList &hitsToAddU(protoShowerU.m_hitsToAdd);
    const CaloHitList &hitsToAddV(protoShowerV.m_hitsToAdd);
    const CaloHitList &hitsToAddW(protoShowerW.m_hitsToAdd);

    ElectronProtoShowerVector showersToExtend;

    if (!hitsToAddU.empty())
        showersToExtend.push_back(protoShowerU);

    if (!hitsToAddV.empty())
        showersToExtend.push_back(protoShowerV);

    if (!hitsToAddW.empty())
        showersToExtend.push_back(protoShowerW);

    if (showersToExtend.size() == 1)
        this->ExtendShowerOneView(pAlgorithm, pShowerPfo, showersToExtend);
    else if (showersToExtend.size() == 2)
        this->ExtendShowerTwoView(pAlgorithm, pShowerPfo, showersToExtend);
    else if (showersToExtend.size() == 3)
        this->ExtendShowerThreeView(pAlgorithm, pShowerPfo, showersToExtend);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::ExtendShowerOneView(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ElectronProtoShowerVector &showersToExtend)
{
    const CaloHitList &hitsToAdd(showersToExtend.front().m_hitsToAdd);

    /////////////////////////////////
    /*
    std::cout << "ONE VIEW HIT ADDITION" << std::endl;
    for (const CaloHit *const pCaloHit : hitsToAdd)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ADDED HIT", BLUE, 2);
    }
    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::ExtendShowerTwoView(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ElectronProtoShowerVector &showersToExtend)
{
    const ElectronProtoShower &protoShower1(showersToExtend[0]), &protoShower2(showersToExtend[1]);

    const CaloHitList &hitsToAdd1(protoShower1.m_hitsToAdd);
    const CaloHitList &hitsToAdd2(protoShower2.m_hitsToAdd);

    /////////////////////////////////
    /*
    std::cout << "TWO VIEW HIT ADDITION" << std::endl;
    for (const CaloHit *const pCaloHit : hitsToAdd1)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ADDED HIT", BLUE, 2);
    }

    for (const CaloHit *const pCaloHit : hitsToAdd2)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ADDED HIT", BLUE, 2);
    }

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd1);
    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ElectronStartRefinementTool::ExtendShowerThreeView(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, 
    const ElectronProtoShowerVector &showersToExtend)
{
    const ElectronProtoShower &protoShower1(showersToExtend[0]), &protoShower2(showersToExtend[1]), &protoShower3(showersToExtend[2]);

    const CaloHitList &hitsToAdd1(protoShower1.m_hitsToAdd);
    const CaloHitList &hitsToAdd2(protoShower2.m_hitsToAdd);
    const CaloHitList &hitsToAdd3(protoShower3.m_hitsToAdd);

    /////////////////////////////////
    /*
    std::cout << "THREE VIEW HIT ADDITION" << std::endl;
    for (const CaloHit *const pCaloHit : hitsToAdd1)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ADDED HIT", BLUE, 2);
    }

    for (const CaloHit *const pCaloHit : hitsToAdd2)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ADDED HIT", BLUE, 2);
    }

    for (const CaloHit *const pCaloHit : hitsToAdd3)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &hitPosition, "ADDED HIT", BLUE, 2);
    }

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    */
    /////////////////////////////////

    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd1);
    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd2);
    pAlgorithm->AddElectronPathway(pShowerPfo, hitsToAdd3);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ElectronStartRefinementTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxCoincideneTransverseSeparation", m_maxCoincideneTransverseSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinSpinePurity", m_minSpinePurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxAngularDeviation", m_maxAngularDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxXSeparation", m_maxXSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxSeparation", m_maxSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ExtendMode", m_extendMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MoveVertexMode", m_moveVertexMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TruncateGammaMode", m_truncateGammaMode));

    ShowerStartRefinementBaseTool::ReadSettings(xmlHandle);

    return STATUS_CODE_SUCCESS;
}
}

//------------------------------------------------------------------------------------------------------------------------------------------


            /* // be a bit nicer
            CartesianVector showerExtremalPoint(0.f, 0.f, 0.f);
            showerExtremalPoint = isEndDownstream ? twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMaxLayerPosition() : 
                twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMinLayerPosition();

            CartesianVector displacement(showerStartPosition - showerExtremalPoint);

            const float l(peakDirection.GetDotProduct(displacement));

            if (l > -1.f)
            {
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerStartPosition, "Shower Start", BLACK, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(pAlgorithm->GetPandora(), &showerExtremalPoint, "Shower Extremal", GREEN, 2);
                std::cout << "l: " << l << std::endl;
                /////////////////////////////////
                std::cout << "NO SHOWER START FOUND" << std::endl;
                PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
                /////////////////////////////////
                continue;
            }
            */
