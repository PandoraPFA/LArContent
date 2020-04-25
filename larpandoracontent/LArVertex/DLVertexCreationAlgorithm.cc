/**
 *  @file   larpandoracontent/LArVertex/DLVertexCreationAlgorithm.cc
 * 
 *  @brief  Implementation of the DLVertexCreation algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h" 

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"

#include "larpandoracontent/LArVertex/DLVertexCreationAlgorithm.h"

#include <torch/script.h>

using namespace pandora;

namespace lar_content
{

DLVertexCreationAlgorithm::DLVertexCreationAlgorithm() :
    m_filePathEnvironmentVariable("FW_SEARCH_PATH"),
    m_numClusterCaloHitsPar(5),
    m_npixels(128),
    m_imgLenVec({-1.0f, 50.0f, 40.0f}),
    m_trainingSetMode(false),
    m_trainingImgLenVecIndex(0),
    m_lenBuffer(10.0),
    m_numViews(3),
    m_trainingDataFileName("data"),
    m_trainingLabelsFileName("labels"),
    m_vertexXCorrection(0.f),
    m_hitWidthZ(0.5),
    m_mcParticleListName("Input"),
    m_caloHitListName("CaloHitList2D"),
    m_interactionTypeFileName("interactionType"),
    m_vecPModule({nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr})
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
DLVertexCreationAlgorithm::~DLVertexCreationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
StatusCode DLVertexCreationAlgorithm::Initialize()
{
    // Load the model file. 
    const StringVector view{"W", "V", "U"};

    const unsigned int numModels(m_trainingSetMode ? m_trainingImgLenVecIndex : m_imgLenVec.size());
    for (int iModel=0; iModel<numModels; iModel++)
    {
        for (int iView=0; iView<m_numViews; iView++)
        {
            const std::string fileNameString(m_modelFileNamePrefix+std::to_string(iModel)+view[iView]+".pt");
            const std::string fullFileNameString(LArFileHelper::FindFileInPath(fileNameString, m_filePathEnvironmentVariable));

            try
            {
                m_vecPModule[m_numViews*iModel+iView] = torch::jit::load(fullFileNameString);
            }
            catch (const c10::Error &e)
            {
                std::cout << "Error loading the PyTorch module" << std::endl;
                return STATUS_CODE_FAILURE;
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
StatusCode DLVertexCreationAlgorithm::Run()
{
    torch::manual_seed(0);

    const ClusterList *pClusterList0, *pClusterList1, *pClusterList2;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListNames[0], pClusterList0));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListNames[1], pClusterList1));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputClusterListNames[2], pClusterList2));

    const HitType hitType0(LArClusterHelper::GetClusterHitType(pClusterList0->front()));
    const HitType hitType1(LArClusterHelper::GetClusterHitType(pClusterList1->front()));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(pClusterList2->front()));

    std::vector<unsigned int> allHitsCountVec{0, 0, 0};
    if (!(this->EventViewCheck(pClusterList0, allHitsCountVec[0]))) return STATUS_CODE_SUCCESS;
    if (!(this->EventViewCheck(pClusterList1, allHitsCountVec[1]))) return STATUS_CODE_SUCCESS;
    if (!(this->EventViewCheck(pClusterList2, allHitsCountVec[2]))) return STATUS_CODE_SUCCESS;
    CartesianVector MCVertexPosition(0.f, 0.f, 0.f);
    int interactionType(-1);
    if (m_trainingSetMode)
    {
        MCVertexPosition = this->GetMCVertexPosition();
        if (this->DetectorCheck(MCVertexPosition)) return STATUS_CODE_SUCCESS;
        interactionType = this->GetInteractionType();
    }

    std::stringstream ssBuf[6];
    const CartesianVector positionInput(0.f, 0.f, 0.f); unsigned int vertReconCount(0);
    CartesianVector position0(this->GetDLVertexForView(pClusterList0,hitType0,positionInput,0,vertReconCount,ssBuf,MCVertexPosition,allHitsCountVec[0]));
    CartesianVector position1(this->GetDLVertexForView(pClusterList1,hitType1,positionInput,0,vertReconCount,ssBuf,MCVertexPosition,allHitsCountVec[1]));
    CartesianVector position2(this->GetDLVertexForView(pClusterList2,hitType2,positionInput,0,vertReconCount,ssBuf,MCVertexPosition,allHitsCountVec[2]));
    if (m_trainingSetMode && (m_trainingImgLenVecIndex==0)) 
    {
        if (vertReconCount == ((1+m_trainingImgLenVecIndex)*m_numViews))
            this->WriteTrainingFiles(ssBuf, interactionType);
        return STATUS_CODE_SUCCESS;
    }

    CartesianVector position3D(0.f, 0.f, 0.f); float chiSquared(0.f);
    for (int imgLenVecIndex=1; imgLenVecIndex<m_imgLenVec.size(); imgLenVecIndex++)
    {
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), hitType0, hitType1, hitType2,
            position0, position1, position2, position3D, chiSquared);

        const CartesianVector position3D0(LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, hitType0));
        const CartesianVector position3D1(LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, hitType1));
        const CartesianVector position3D2(LArGeometryHelper::ProjectPosition(this->GetPandora(), position3D, hitType2));

        position0=this->GetDLVertexForView(pClusterList0,hitType0,position3D0,imgLenVecIndex,vertReconCount,ssBuf,MCVertexPosition,allHitsCountVec[0]);
        position1=this->GetDLVertexForView(pClusterList1,hitType1,position3D1,imgLenVecIndex,vertReconCount,ssBuf,MCVertexPosition,allHitsCountVec[1]);
        position2=this->GetDLVertexForView(pClusterList2,hitType2,position3D2,imgLenVecIndex,vertReconCount,ssBuf,MCVertexPosition,allHitsCountVec[2]);
        if (m_trainingSetMode && (m_trainingImgLenVecIndex==imgLenVecIndex))
        {
            if (vertReconCount == ((1+m_trainingImgLenVecIndex)*m_numViews))
                this->WriteTrainingFiles(ssBuf, interactionType);
            return STATUS_CODE_SUCCESS;
        }
    }

    if (vertReconCount == m_numViews*m_imgLenVec.size())
    {
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), hitType0, hitType1, hitType2,
            position0, position1, position2, position3D, chiSquared);

        const VertexList *pVertexList(nullptr); std::string temporaryListName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pVertexList, temporaryListName));

        PandoraContentApi::Vertex::Parameters parameters;
        parameters.m_position = position3D;
        parameters.m_vertexLabel = VERTEX_INTERACTION;
        parameters.m_vertexType = VERTEX_3D;
        const Vertex *pVertex(nullptr);

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
    }

    return STATUS_CODE_SUCCESS;
}

//---------------------------------------------------------------------------------------------------------------------------------
CartesianVector DLVertexCreationAlgorithm::GetDLVertexForView(const ClusterList *const pClusterList, const HitType &view,
    const CartesianVector &positionInput, const int imgLenVecIndex, unsigned int &vertReconCount, 
    std::stringstream ssBuf[6], const CartesianVector &MCVertexPosition, const unsigned int allHitsCount)
{
    const double length(m_imgLenVec[imgLenVecIndex]);
    const bool useScaledImg(length <= 0);
    DoubleVector xVec, zVec, sigmaVec, heightVec;
    int iHits(0);
    double minx(std::numeric_limits<double>::max()), minz(std::numeric_limits<double>::max());
    double maxx(-std::numeric_limits<double>::max()), maxz(-std::numeric_limits<double>::max());

    xVec.resize(allHitsCount, 0.f) ; zVec.resize(allHitsCount, 0.f)  ;
    sigmaVec.resize(allHitsCount, 0.f); heightVec.resize(allHitsCount, 0.f);

    for (const Cluster *const pCluster : *pClusterList) 
    {
        const OrderedCaloHitList &orderedCaloHitList1(pCluster->GetOrderedCaloHitList());
        for (OrderedCaloHitList::const_iterator iter1 = orderedCaloHitList1.begin(); iter1 != orderedCaloHitList1.end(); ++iter1)
        {
            for (CaloHitList::const_iterator hitIter1 = iter1->second->begin(); hitIter1 != iter1->second->end(); ++hitIter1)
            {
                const CartesianVector &positionVector1((*hitIter1)->GetPositionVector());
                xVec[iHits] = positionVector1.GetX(); zVec[iHits] = positionVector1.GetZ();
                sigmaVec[iHits] = (*hitIter1)->GetCellSize1(); heightVec[iHits] = (*hitIter1)->GetInputEnergy();
                if( (pCluster->GetNCaloHits()>=m_numClusterCaloHitsPar) && useScaledImg )
                {
                    minx = std::min(minx,xVec[iHits]); minz = std::min(minz,zVec[iHits]);
                    maxx = std::max(maxx,xVec[iHits]); maxz = std::max(maxz,zVec[iHits]);
                }
                iHits++;
            }
        }
    }

    double nstepx(0), nstepz(0), lengthX(0), lengthZ(0), scalImgLen(0);
    TwoDImage out2dVec(m_npixels, DoubleVector(m_npixels, 0));

    if (useScaledImg)
    {
        minx -= m_lenBuffer; minz -= m_lenBuffer;
        maxx += m_lenBuffer; maxz += m_lenBuffer;
        lengthX = maxx - minx;
        lengthZ = maxz - minz;
        scalImgLen = std::max(lengthX, lengthZ);
        nstepx = scalImgLen / m_npixels;
        nstepz = scalImgLen / m_npixels;
    }
    else
    {
        minx = positionInput.GetX() - length/2.0;
        minz = positionInput.GetZ() - length/2.0;
        nstepx = length / m_npixels;
        nstepz = length / m_npixels;
    }

    for (int i=0; i<allHitsCount; i++)
    {
         const int ZPixPosit = (int)((zVec[i]-minz)/nstepz);
         const int XPixPosit = (int)((xVec[i]-minx)/nstepx);
         if (XPixPosit>m_npixels||XPixPosit<0||ZPixPosit>m_npixels||ZPixPosit<0) continue;

         for (int j=0; j<m_npixels; j++)
         {
             double tempval(0), inputvalue1(0), dist(0);
             tempval = heightVec[i]*(0.5*std::erfc(-((minx+(j+1)*nstepx)-xVec[i])/std::sqrt(2*sigmaVec[i]*sigmaVec[i]))
                 - 0.5*std::erfc(-((minx+j*nstepx)-xVec[i])/std::sqrt(2*sigmaVec[i]*sigmaVec[i])));

             for (int k=m_npixels-1; k>-1; k--)
             {
                 dist = (minz+(k+1)*nstepz)-(zVec[i]+m_hitWidthZ/2.0);

                 if (dist>nstepz)
                     inputvalue1=0;
                 else if (dist<nstepz&&dist>0)
                     inputvalue1=std::min((nstepz-dist),m_hitWidthZ)*tempval;
                 else if (dist<0 && std::fabs(dist)<nstepz)
                     {if (std::fabs(dist)>m_hitWidthZ) inputvalue1=0; else inputvalue1=std::min((m_hitWidthZ-std::fabs(dist)),nstepz)*tempval;}
                 else if (dist<0 && std::fabs(dist)>nstepz)
                     {if (std::fabs(dist)>m_hitWidthZ) inputvalue1=0; else inputvalue1=std::min((m_hitWidthZ-std::fabs(dist)),nstepz)*tempval;}

                 out2dVec[j][k] += inputvalue1;
             }
         }
    }

    if (m_trainingSetMode && (m_trainingImgLenVecIndex==imgLenVecIndex))
    {
        vertReconCount += (this->CreateTrainingFiles(out2dVec,view,minx,nstepx,minz,nstepz,ssBuf,MCVertexPosition));
        return(CartesianVector(0.f, 0.f, 0.f));
    }
    else
    {
        const CartesianVector pixelPosition(this->DeepLearning(out2dVec,view,imgLenVecIndex));

        const double recoX(minx+(pixelPosition.GetX())*nstepx);
        const double recoZ(minz+(pixelPosition.GetZ())*nstepz);
        const CartesianVector position(recoX, 0.f, recoZ);
        vertReconCount += 1;
        return(position);
    }
}

//---------------------------------------------------------------------------------------------------------------------------------
CartesianVector DLVertexCreationAlgorithm::DeepLearning(const TwoDImage &out2dVec, const HitType &view,
    const int imgLenVecIndex)
{
    //torch::NoGradGuard no_grad;

    /* Get the index for model */
    int index(0);
    if (view == TPC_VIEW_W) index = m_numViews*imgLenVecIndex+0;
    if (view == TPC_VIEW_V) index = m_numViews*imgLenVecIndex+1;
    if (view == TPC_VIEW_U) index = m_numViews*imgLenVecIndex+2;

    /* Convert image to Torch "Tensor" */
    torch::Tensor input = torch::zeros({1, 1, m_npixels, m_npixels}, torch::TensorOptions().dtype(torch::kFloat32));
    auto accessor = input.accessor<float, 4>();

    for (int i=m_npixels-1; i>-1; i--)
        for (int j=0; j<m_npixels; j++)
            accessor[0][0][i][j] = out2dVec[j][i]/100.0;

    std::vector<torch::jit::IValue> inputs;
    inputs.push_back(input);

    /* Use the Model on the "Tensor" */
    at::Tensor output = m_vecPModule[index]->forward(inputs).toTensor();

    /* Get predicted vertex. */
    const auto &outputAccessor(output.accessor<float, 2>());
    const CartesianVector position(outputAccessor[0][0], 0.f, outputAccessor[0][1]);

    return(position);
}

//------------------------------------------------------------------------------------------------------------------------------------------
int DLVertexCreationAlgorithm::CreateTrainingFiles(const TwoDImage &out2dVec, const HitType &view, const float minx, const float nstepx,
    const float minz, const float nstepz, std::stringstream ssBuf[6], const CartesianVector &MCVertexPosition) const
{
    int index(0);
    if (view == TPC_VIEW_W) index = 0;
    if (view == TPC_VIEW_V) index = 1;
    if (view == TPC_VIEW_U) index = 2;
    const CartesianVector VertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), MCVertexPosition, view));

    const double MCXPixPosit((VertexPosition.GetX()-minx)/nstepx);
    const double MCZPixPosit((VertexPosition.GetZ()-minz)/nstepz);
    if ((int)MCXPixPosit>m_npixels||(int)MCXPixPosit<0||(int)MCZPixPosit>m_npixels||(int)MCZPixPosit<0)
        return(0);

    for (int i=m_npixels-1; i>-1; i--)
    {
        for (int j=0; j<m_npixels; j++)
        {
            ssBuf[index] << out2dVec[j][i] << " ";
        }
            ssBuf[index] << "\n";
    }
    ssBuf[m_numViews+index] << MCXPixPosit << " " << MCZPixPosit << "\n";

    return(1);
}

//------------------------------------------------------------------------------------------------------------------------------------------
void DLVertexCreationAlgorithm::WriteTrainingFiles(std::stringstream ssBuf[6], const int interactionType) const
{
    const StringVector view{"W", "V", "U"};
    for (int iView=0; iView<m_numViews; iView++)
    {
        std::ofstream outputStream1(m_trainingDataFileName+view[iView]+".txt", std::ios::app);
        std::ofstream outputStream2(m_trainingLabelsFileName+view[iView]+".txt", std::ios::app);
        
        outputStream1 << ssBuf[iView].rdbuf();
        outputStream2 << ssBuf[m_numViews+iView].rdbuf();

        outputStream1.close();
        outputStream2.close();
    }

    std::ofstream outputStream3(m_interactionTypeFileName, std::ios::app);
    outputStream3 << interactionType << std::endl;
    outputStream3.close();
}

//------------------------------------------------------------------------------------------------------------------------------------------
bool DLVertexCreationAlgorithm::DetectorCheck(const pandora::CartesianVector &position3D) const
{
    const LArTPCMap &larTPCMap(this->GetPandora().GetGeometry()->GetLArTPCMap());
    const LArTPC *const pFirstLArTPC(larTPCMap.begin()->second);

    float parentMinX(pFirstLArTPC->GetCenterX() - 0.5f * pFirstLArTPC->GetWidthX());
    float parentMaxX(pFirstLArTPC->GetCenterX() + 0.5f * pFirstLArTPC->GetWidthX());
    float parentMinY(pFirstLArTPC->GetCenterY() - 0.5f * pFirstLArTPC->GetWidthY());
    float parentMaxY(pFirstLArTPC->GetCenterY() + 0.5f * pFirstLArTPC->GetWidthY());
    float parentMinZ(pFirstLArTPC->GetCenterZ() - 0.5f * pFirstLArTPC->GetWidthZ());
    float parentMaxZ(pFirstLArTPC->GetCenterZ() + 0.5f * pFirstLArTPC->GetWidthZ());

    for (const LArTPCMap::value_type &mapEntry : larTPCMap)
    {
        const LArTPC *const pLArTPC(mapEntry.second);
        parentMinX = std::min(parentMinX, pLArTPC->GetCenterX() - 0.5f * pLArTPC->GetWidthX());
        parentMaxX = std::max(parentMaxX, pLArTPC->GetCenterX() + 0.5f * pLArTPC->GetWidthX());
        parentMinY = std::min(parentMinY, pLArTPC->GetCenterY() - 0.5f * pLArTPC->GetWidthY());
        parentMaxY = std::max(parentMaxY, pLArTPC->GetCenterY() + 0.5f * pLArTPC->GetWidthY());
        parentMinZ = std::min(parentMinZ, pLArTPC->GetCenterZ() - 0.5f * pLArTPC->GetWidthZ());
        parentMaxZ = std::max(parentMaxZ, pLArTPC->GetCenterZ() + 0.5f * pLArTPC->GetWidthZ());
    }

    bool volumeCheck(true), gapCheck(true);
    volumeCheck = volumeCheck && ( (position3D.GetX()<parentMaxX) && (position3D.GetX()>parentMinX) );
    volumeCheck = volumeCheck && ( (position3D.GetY()<parentMaxY) && (position3D.GetY()>parentMinY) );
    volumeCheck = volumeCheck && ( (position3D.GetZ()<parentMaxZ) && (position3D.GetZ()>parentMinZ) );

    gapCheck = gapCheck && LArGeometryHelper::IsInGap3D(this->GetPandora(), position3D, TPC_VIEW_W, 0.f);
    gapCheck = gapCheck && LArGeometryHelper::IsInGap3D(this->GetPandora(), position3D, TPC_VIEW_V, 0.f);
    gapCheck = gapCheck && LArGeometryHelper::IsInGap3D(this->GetPandora(), position3D, TPC_VIEW_U, 0.f);

    return((!volumeCheck) || gapCheck);
}

//------------------------------------------------------------------------------------------------------------------------------------------
bool DLVertexCreationAlgorithm::EventViewCheck(const pandora::ClusterList *const pClusterList, unsigned int &allHitsCount) const
{
    unsigned int filtHitsCount(0);

    for (const Cluster *const pCluster : *pClusterList) 
    {
        allHitsCount += pCluster->GetNCaloHits();
        if (pCluster->GetNCaloHits() >= m_numClusterCaloHitsPar) filtHitsCount += pCluster->GetNCaloHits();
    }

    return(!(allHitsCount==0 || filtHitsCount==0));
}

//------------------------------------------------------------------------------------------------------------------------------------------
pandora::CartesianVector DLVertexCreationAlgorithm::GetMCVertexPosition() const
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    CartesianVector MCVertexPosition(0.f, 0.f, 0.f);
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (!LArMCParticleHelper::IsNeutrino(pMCParticle))
            continue;

        MCVertexPosition.SetValues(pMCParticle->GetEndpoint().GetX() + m_vertexXCorrection, pMCParticle->GetEndpoint().GetY(), 
            pMCParticle->GetEndpoint().GetZ());
    }

    return(MCVertexPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------
int DLVertexCreationAlgorithm::GetInteractionType() const
{
    // Extract input collections
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    // ATTN Assumes single neutrino is parent of all neutrino-induced mc particles
    LArMCParticleHelper::MCContributionMap nuMCParticlesToGoodHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, LArMCParticleHelper::PrimaryParameters(),
        LArMCParticleHelper::IsBeamNeutrinoFinalState, nuMCParticlesToGoodHitsMap);

    MCParticleList mcPrimaryList;
    for (const auto &mapEntry : nuMCParticlesToGoodHitsMap) mcPrimaryList.push_back(mapEntry.first);
    mcPrimaryList.sort(LArMCParticleHelper::SortByMomentum);

    const LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(mcPrimaryList));
    return((int)(interactionType));
}

//------------------------------------------------------------------------------------------------------------------------------------------
StatusCode DLVertexCreationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    // Read settings from xml file here
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "ModelFileNamePrefix", m_modelFileNamePrefix));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FilePathEnvironmentVariable", m_filePathEnvironmentVariable));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "NumClusterCaloHitsPar", m_numClusterCaloHitsPar));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Npixels", m_npixels));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "ImgLenVec", m_imgLenVec));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingSetMode", m_trainingSetMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingImgLenVecIndex", m_trainingImgLenVecIndex));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LenBuffer", m_lenBuffer));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingDataFileName", m_trainingDataFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrainingLabelsFileName", m_trainingLabelsFileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VertexXCorrection", m_vertexXCorrection));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HitWidthZ", m_hitWidthZ));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InteractionTypeFileName", m_interactionTypeFileName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
