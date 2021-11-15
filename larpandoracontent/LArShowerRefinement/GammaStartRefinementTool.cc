/**
 *  @file   larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.cc
 *
 *  @brief  Implementation of the gamma start refinement tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArShowerRefinement/GammaStartRefinementTool.h"

using namespace pandora;

namespace lar_content
{

GammaStartRefinementTool::GammaStartRefinementTool() : m_counter(0)
{
}

GammaStartRefinementTool::~GammaStartRefinementTool()
{
    try
    {
        std::cout << "11111111111" << std::endl;
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ShowerDistribution", "ShowerDistribution.root", "UPDATE"));
        std::cout << "2222222" << std::endl;
    }
    catch (const StatusCodeException &)
    {
        std::cout << "BAD JAM" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool GammaStartRefinementTool::Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const ParticleFlowObject *const pShowerPfo, const CartesianVector &nuVertexPosition)
{
    std::cout << "AAAAAAAAA" << std::endl;

    // only apply gamma refinement algorithm to showers
    if (!LArPfoHelper::IsShower(pShowerPfo))
        return false;

    ////////////////////////////////
    // temporary 
    CaloHitList caloHits3D;
    LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHits3D);

    ClusterList clusters3D;
    LArPfoHelper::GetClusters(pShowerPfo, TPC_3D, clusters3D);

    if (caloHits3D.size() < 100)
        return false;
    ////////////////////////////////

    ++m_counter;
    float deviationAngleX = -999;
    float deviationAngleY = -999;
    float deviationAngle1 = -999; // transverse
    float deviationAngle2 = -999; // transverse
    float deviationAngle3 = -999; // main shower direction

    float xCoordinate = -999;
    float yCoordinate = -999;
    float zCoordinate = -999;
    float lCoordinate1 = -999;
    float tCoordinate1 = -999;
    float lCoordinate2 = -999;
    float tCoordinate2 = -999;
    float hitEnergy = -999;

    float lNu1 = -999;
    float tNu1 = -999;
    float lNu2 = -999;
    float tNu2 = -999;

    CartesianPointVector initialHitPositions3D;

    for (const CaloHit *const pCaloHit : caloHits3D)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        if ((hitPosition - nuVertexPosition).GetMagnitude() > 14.f)
            continue;

        initialHitPositions3D.push_back(pCaloHit->GetPositionVector());
    }

    const ThreeDSlidingFitResult threeDSlidingFit(&initialHitPositions3D, 10000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    //const CartesianVector &axisIntercept(threeDSlidingFit.GetAxisIntercept());
    const CartesianVector &axis3(threeDSlidingFit.GetAxisDirection());

    const TwoDSlidingFitResult plane1(threeDSlidingFit.GetFirstFitResult());
    //const CartesianVector &axisIntercept1(plane1.GetAxisIntercept());
    const CartesianVector &axis1(plane1.GetAxisDirection());

    const TwoDSlidingFitResult plane2(threeDSlidingFit.GetSecondFitResult());
    //const CartesianVector &axisIntercept2(plane2.GetAxisIntercept());
    const CartesianVector &axis2(plane2.GetAxisDirection());

    const CartesianVector xAxis(1.f, 0.f, 0.f);
    const CartesianVector yAxis(0.f, 1.f, 0.f);

    std::cout << "axis1: " << axis1 << std::endl;
    std::cout << "axis2: " << axis2 << std::endl;
    std::cout << "axis3: " << axis3 << std::endl;

    std::cout << "min layer direction1: " << plane1.GetGlobalMinLayerDirection() << std::endl;
    std::cout << "min layer direction2: " << plane2.GetGlobalMinLayerDirection() << std::endl;
    std::cout << "min layer direction3: " << threeDSlidingFit.GetGlobalMinLayerDirection() << std::endl;

    std::cout << "axis3 dot axis2: " << axis3.GetDotProduct(axis2) << std::endl;
    std::cout << "axis3 dot axis1: " << axis3.GetDotProduct(axis1) << std::endl;

    int highestTheta0XZ = 0, highestTheta0YZ = 0;

    for (const CaloHit *const pCaloHit : caloHits3D)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        if ((hitPosition - nuVertexPosition).GetMagnitude() > 14.f)
            continue;

        hitEnergy = pCaloHit->GetElectromagneticEnergy();

        xCoordinate = hitPosition.GetX();
        yCoordinate = hitPosition.GetY();
        zCoordinate = hitPosition.GetZ();

        plane1.GetLocalPosition(hitPosition, lCoordinate1, tCoordinate1);
        plane1.GetLocalPosition(nuVertexPosition, lNu1, tNu1);
        plane2.GetLocalPosition(hitPosition, lCoordinate2, tCoordinate2);
        plane2.GetLocalPosition(nuVertexPosition, lNu2, tNu2);

        std::cout << "lNu1: " << lNu1 << std::endl;
        std::cout << "tNu1: " << tNu1 << std::endl;
        std::cout << "lNu2: " << lNu2 << std::endl;
        std::cout << "tNu2: " << tNu2 << std::endl;

        deviationAngleX = hitPosition.GetOpeningAngle(xAxis); 
        deviationAngleY = hitPosition.GetOpeningAngle(yAxis); 
        //deviationAngle1 = (hitPosition - nuVertexPosition).GetOpeningAngle(axis1);
        //deviationAngle2 = (hitPosition - nuVertexPosition).GetOpeningAngle(axis2);

        deviationAngle1 = std::atan((tCoordinate1 - tNu1) / (lCoordinate1 - lNu1));
        deviationAngle2 = std::atan((tCoordinate2 - tNu2) / (lCoordinate2 - lNu2));
        deviationAngle3 = (hitPosition - nuVertexPosition).GetOpeningAngle(axis3);

        /*
        if (tCoordinate1 < 0.f)
            deviationAngle1 *= -1.f;
        */

        const int theta0XZFactor(std::floor(deviationAngleX / pAlgorithm->m_binSize));
        const int theta0YZFactor(std::floor(deviationAngleY / pAlgorithm->m_binSize));

        if (theta0XZFactor > highestTheta0XZ)
            highestTheta0XZ = theta0XZFactor;

        if (theta0YZFactor > highestTheta0YZ)
            highestTheta0YZ = theta0YZFactor;

        pAlgorithm->m_theta0XZMap[theta0XZFactor].push_back(pCaloHit);
        pAlgorithm->m_theta0YZMap[theta0YZFactor].push_back(pCaloHit);

        if (zCoordinate < 0.f)
            deviationAngleX = (3.14 - deviationAngleX);

        if (xCoordinate < 0.f)
            deviationAngleY = (3.14 - deviationAngleY);

        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "ShowerCounter", m_counter));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngleX", deviationAngleX));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngleY", deviationAngleY));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngle1", deviationAngle1));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngle2", deviationAngle2));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "DeviationAngle3", deviationAngle3));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "XCoordinate", xCoordinate));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "YCoordinate", yCoordinate));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "ZCoordinate", zCoordinate));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinate1", lCoordinate1));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "TCoordinate1", tCoordinate1));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "LCoordinate2", lCoordinate2));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "TCoordinate2", tCoordinate2));
        PANDORA_MONITORING_API(SetTreeVariable(pAlgorithm->GetPandora(), "ShowerDistribution", "HitEnergy", hitEnergy));
        PANDORA_MONITORING_API(FillTree(pAlgorithm->GetPandora(), "ShowerDistribution"));
    }


    // Get total shower charge
    float total = 0;
    for (const auto &entry : pAlgorithm->m_theta0XZMap)
    {
        for (const CaloHit *const pCaloHit : entry.second)
            total += pCaloHit->GetElectromagneticEnergy();
    }

    if (total == 0)
        return false;

    // Get number of significant bins
    int nSigBins0XZ(0), nSigBins0YZ(0);

    for (const auto &entry : pAlgorithm->m_theta0XZMap)
    {
        float weight = 0;

        for (const CaloHit *const pCaloHit : entry.second)
            weight += pCaloHit->GetElectromagneticEnergy();

        if ((weight / total) > 0.01)
            ++nSigBins0XZ;

        std::cout << entry.first * pAlgorithm->m_binSize << ": " << weight / total << std::endl;
    }

    std::cout << "nSigBins0XZ: " << nSigBins0XZ << std::endl;

    for (const auto &entry : pAlgorithm->m_theta0YZMap)
    {
        float weight = 0;

        for (const CaloHit *const pCaloHit : entry.second)
            weight += pCaloHit->GetElectromagneticEnergy();

        if ((weight / total) > 0.01)
            ++nSigBins0YZ;
    }

    std::cout << "nSigBins0YZ: " << nSigBins0YZ << std::endl;

    // merge significant bins
    /*
    bool previous(false);
    std::vector<CaloHitList> jam;

    for (int i = 0; i <= highestTheta0XZ; ++i)
    {
        if (pAlgorithm->m_theta0XZMap.find(i) == pAlgorithm->m_theta0XZMap.end())
        {
            previous = false;
            continue;
        }

        float weight = 0;

        for (const CaloHit *const pCaloHit : pAlgorithm->m_theta0XZMap.at(i))
            weight += pCaloHit->GetElectromagneticEnergy();

        if ((weight / total) < (1.f / static_cast<float>(nSigBins0XZ)))
        {
            previous = false;
            continue;
        }

        std::cout << "i: " << i << std::endl;
        std::cout << "weight: " << weight << std::endl;

        if (previous)
        {
            std::cout << "previous = true" << std::endl; 
            jam.back().insert(jam.back().begin(), pAlgorithm->m_theta0XZMap.at(i).begin(), pAlgorithm->m_theta0XZMap.at(i).end());
        }
        else
        {
            std::cout << "previous = false" << std::endl;
            jam.push_back(pAlgorithm->m_theta0XZMap.at(i));
        }

        previous = true;
    }
    */

    bool previous(false);
    std::vector<CaloHitList> jam;

    for (int i = 0; i <= highestTheta0YZ; ++i)
    {
        if (pAlgorithm->m_theta0YZMap.find(i) == pAlgorithm->m_theta0YZMap.end())
        {
            previous = false;
            continue;
        }

        float weight = 0;

        for (const CaloHit *const pCaloHit : pAlgorithm->m_theta0YZMap.at(i))
            weight += pCaloHit->GetElectromagneticEnergy();

        if ((weight / total) < (1.f / static_cast<float>(nSigBins0XZ)))
        {
            previous = false;
            continue;
        }

        std::cout << "i: " << i << std::endl;
        std::cout << "weight: " << weight << std::endl;

        if (previous)
        {
            std::cout << "previous = true" << std::endl; 
            jam.back().insert(jam.back().begin(), pAlgorithm->m_theta0YZMap.at(i).begin(), pAlgorithm->m_theta0YZMap.at(i).end());
        }
        else
        {
            std::cout << "previous = false" << std::endl;
            jam.push_back(pAlgorithm->m_theta0YZMap.at(i));
        }

        previous = true;
    }



    std::cout << "number of subsections: " << jam.size() << std::endl;

    /*
    int highestBin0XZ(0);
    float highestWeight0XZ(-std::numeric_limits<float>::max());

    for (const auto &entry : pAlgorithm->m_theta0XZMap)
    {
        std::cout << "theta0XZ weight: " << static_cast<float>(entry.second.size()) / total << std::endl;

        float weight = 0.0;

        for (const CaloHit *const pCaloHit : entry.second)
            weight += pCaloHit->GetElectromagneticEnergy();

        if (weight > highestWeight0XZ)
        {
            highestWeight0XZ = entry.second.size();
            highestBin0XZ = entry.first;
        }
    }

    //int highestBin0YZ(0);
    //float highestWeight0YZ(-std::numeric_limits<float>::max());

    for (const auto &entry : pAlgorithm->m_theta0YZMap)
    {
        std::cout << "theta0YZ weight: " << static_cast<float>(entry.second.size()) / total << std::endl;

    }

    CaloHitList visualise;
    int highestHits = 0;
    const CaloHitList &hits0XZ(pAlgorithm->m_theta0XZMap[highestBin0XZ]);

    for (auto &entry : pAlgorithm->m_theta0YZMap)
    {
        const CaloHitList &sharedHits(LArMCParticleHelper::GetSharedHits(hits0XZ, entry.second));
        const int nSharedHits(sharedHits.size());

        if (nSharedHits > highestHits)
        {
            std::cout << "nSharedHits: " << nSharedHits << std::endl;
            highestHits = nSharedHits;
            visualise = sharedHits;
        }
    }


    for (const CaloHit *const pCaloHit : caloHits3D)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        if ((hitPosition - nuVertexPosition).GetMagnitude() < 14.f)
            continue;

        deviationAngleX = hitPosition.GetOpeningAngle(xAxis); 
        deviationAngleY = hitPosition.GetOpeningAngle(yAxis); 

        if ((deviationAngleX < highestBin0XZ) || (deviationAngleX > (highestBin0XZ + pAlgorithm->m_binSize)))
            continue;

        if ((deviationAngleY < highestBin0YZ) || (deviationAngleY > (highestBin0YZ + pAlgorithm->m_binSize)))
            continue;

        visualise.push_back(pCaloHit);
    }
    */

    for (const CaloHitList &visualise : jam)
        PandoraMonitoringApi::VisualizeCaloHits(pAlgorithm->GetPandora(), &visualise, "subsection", RED);

    PandoraMonitoringApi::ViewEvent(pAlgorithm->GetPandora());
    


    std::cout << "BBBBBBBBBB" << std::endl;

    // ISOBEL - SHOULD I ALSO HAVE A HIT CUT?
    if (!this->HasPathToNuVertex(pShowerPfo, nuVertexPosition))
        return false;

    std::cout << "CCCCCCCCCCCC" << std::endl;

    ProtoShowerVector protoShowerVector;
    this->BuildProtoShowers(pShowerPfo, protoShowerVector);

    std::cout << "DDDDDDDDDDD" << std::endl;

    if (protoShowerVector.empty())
        return false;

    std::cout << "EEEEEEEEEEEE" << std::endl;

    // pfo splitting alg? (can be quite simple if we are able to build the protoShowers sufficiently...)
    // set ProtoShower parent pfo address to match!
    if (protoShowerVector.size() != 1)
        return false;

    for (const ProtoShower &protoShower : protoShowerVector)
    {
        if (this->IsElectronPathway(protoShower))
            continue;

        this->RemoveConnectionPathway(protoShower);

        // change metadata to say that we think that this gamma is a gamma (so that we don't extend it in the future tool)
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void GammaStartRefinementTool::RemoveConnectionPathway(const ProtoShower &protoShower)
{
    std::cout << protoShower.m_showerCore.m_startPosition.GetX() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode GammaStartRefinementTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
