#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArObjects/LArPfoObjects.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "ShowerMergingPfoMopUpAlgorithm.h"

using namespace pandora;

namespace lar_content
{


StatusCode ShowerMergingPfoMopUpAlgorithm::Run()
{

    PfoVector sortedPfos;
    this->GetSortedPfos(sortedPfos);
   
    std::set < const ParticleFlowObject * > deletedPfos;


    int i = 0;
    for (const ParticleFlowObject *const pPfo1 : sortedPfos)
    {  
       if (deletedPfos.find (pPfo1) != deletedPfos.end())
           continue;

       int j =0;
       i++;
 

        for (const ParticleFlowObject *const pPfo2 : sortedPfos)
        {   

           if (deletedPfos.find (pPfo2) != deletedPfos.end())
               continue;

            j++;
            if (!pPfo1 || !pPfo2 || (pPfo1 == pPfo2))
                continue;
            
            if (pPfo1->GetParentPfoList().empty() || pPfo2->GetParentPfoList().empty())
                continue;

            const bool isAligned{this->AreDirectionsAligned(pPfo1, pPfo2)};
            bool invert = false;
            const bool isSameVertex{this->HaveSameVertex(pPfo1, pPfo2, invert)};
           

            if (isAligned && isSameVertex)
            {   
                

               if (!invert)
               {
                   for ( const Pfo* const daughter : pPfo2->GetDaughterPfoList() )
                       if (daughter == pPfo1)
                        invert=true;
               }
               

                if (invert)
                { 
                    this->MergeAndDeletePfos(pPfo2, pPfo1);

                }

                else 
                {
                    this->MergeAndDeletePfos(pPfo1, pPfo2);
                }
                    
                    deletedPfos.insert(pPfo1);
                    deletedPfos.insert(pPfo2);
                    break;

            }
        }
    }

    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------

void ShowerMergingPfoMopUpAlgorithm::GetSortedPfos(PfoVector &sortedPfos) const
{
    
    for (const std::string &listName : m_inputPfoListNames)
    {   
        const PfoList *pPfoList(nullptr);
        
        if (STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, listName, pPfoList))
            sortedPfos.insert(sortedPfos.end(), pPfoList->begin(), pPfoList->end());
    }
    std::sort(sortedPfos.begin(), sortedPfos.end(), LArPfoHelper::SortByNHits);

}

//------------------------------------------------------------------------------------------------------------


bool ShowerMergingPfoMopUpAlgorithm::AreDirectionsAligned(const ParticleFlowObject *const pPfo1, const ParticleFlowObject *const pPfo2) const
{
    if (pPfo1->GetVertexList().empty() || pPfo2->GetVertexList().empty())
        return false;

    const Vertex *const pVertex1 = LArPfoHelper::GetVertex(pPfo1);
    const Vertex *const pVertex2 = LArPfoHelper::GetVertex(pPfo2);

    if (!pVertex1 || !pVertex2)
        return false;

    const LArShowerPCA pca1 = LArPfoHelper::GetPrincipalComponents(pPfo1,pVertex1);
    const LArShowerPCA pca2 = LArPfoHelper::GetPrincipalComponents(pPfo2,pVertex2);

    const CartesianVector dir1 = pca1.GetPrimaryAxis();
    const CartesianVector dir2 = pca2.GetPrimaryAxis();

    const float cosAngle = dir1.GetUnitVector().GetDotProduct(dir2.GetUnitVector());
    return (cosAngle > m_alignmentAngle);

}

//------------------------------------------------------------------------------------------------------------


bool ShowerMergingPfoMopUpAlgorithm::HaveSameVertex(const ParticleFlowObject *const pPfo1, const ParticleFlowObject *const pPfo2, bool&invert) const
{
    if (pPfo1->GetVertexList().empty() || pPfo2->GetVertexList().empty())
             return false;

    const Pfo *const pNu = LArPfoHelper::GetParentNeutrino(pPfo1);
    const Vertex *const pNuVertex = LArPfoHelper::GetVertex(pNu);

    const Vertex *const pVertex1 = LArPfoHelper::GetVertex(pPfo1);
    const Vertex *const pVertex2 = LArPfoHelper::GetVertex(pPfo2);

    const CartesianVector &nuVtx = pNuVertex->GetPosition();

    const CartesianVector vtx1 = pVertex1->GetPosition();
    const CartesianVector vtx2 = pVertex2->GetPosition();

    const float distNuVtxToVtx1 = (nuVtx - vtx1).GetMagnitude();
    const float distNuVtxToVtx2 = (nuVtx - vtx2).GetMagnitude();
    
    if  (distNuVtxToVtx1 > 3.f && distNuVtxToVtx2 > 3.f)
      return false;

    const Pfo * pStub = nullptr;
    const Pfo * pShw = nullptr;
   
    const CartesianVector StubVtx = distNuVtxToVtx1 < distNuVtxToVtx2 ? vtx1:vtx2;
    const CartesianVector ShwVtx = distNuVtxToVtx1 < distNuVtxToVtx2 ? vtx2:vtx1;

    if ( distNuVtxToVtx1 < distNuVtxToVtx2) 
        { pStub = pPfo1;
          pShw = pPfo2;
          invert = false;
        }
     else 
        { pStub = pPfo2;
          pShw = pPfo1;
          invert = true;
        }

    if (!LArPfoHelper::IsShower(pShw))
        return false;

    CaloHitList StubHits, ShwHits;
    LArPfoHelper::GetCaloHits(pStub, TPC_3D, StubHits);
    LArPfoHelper::GetCaloHits(pShw, TPC_3D, ShwHits);

    
    float largestDist = -1;
    const CaloHit  *pStubEndPoint = nullptr;

    for ( const CaloHit * const pCaloHit: StubHits)
       {
         const CartesianVector &pos = pCaloHit->GetPositionVector();
         const float endPointDist = ( pos - StubVtx).GetMagnitudeSquared(); 
         if ( endPointDist > largestDist)
            { 
              largestDist = endPointDist;
              pStubEndPoint = pCaloHit;
       
            }
        }
    
    const float distNuVtxToStubEnd = ( nuVtx - pStubEndPoint->GetPositionVector()).GetMagnitude();
    const float distNuVtxToShwStart = ( nuVtx - ShwVtx ).GetMagnitude();

    if ( distNuVtxToShwStart < distNuVtxToStubEnd)
      return false;


    const float StubShwSep = ( ShwVtx - pStubEndPoint->GetPositionVector() ).GetMagnitudeSquared(); 
    

    return StubShwSep<m_stubShowerSeperation;

}

//------------------------------------------------------------------------------------------------------------

StatusCode ShowerMergingPfoMopUpAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputPfoListNames", m_inputPfoListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "AlignmentAngle", m_alignmentAngle));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "StubShowerSeperation", m_stubShowerSeperation));   
    
    return PfoMopUpBaseAlgorithm::ReadSettings(xmlHandle);
}

}
