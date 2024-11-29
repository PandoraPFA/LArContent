/**
 *  @file   larpandoracontent/LArReclustering/SimplePCAThreeDClusteringTool.cc
 *
 *  @brief  Implementation file for the simple PCA-based clustering tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"
#include "larpandoracontent/LArReclustering/SimplePCAThreeDClusteringTool.h"

using namespace pandora;

namespace lar_content
{

SimplePCAThreeDClusteringTool::SimplePCAThreeDClusteringTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SimplePCAThreeDClusteringTool::Run(const Algorithm *const /*pAlgorithm*/, const CaloHitList &inputCaloHitList, std::vector<CaloHitList> &outputCaloHitListsVector)
{
    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(inputCaloHitList, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

    // Now define ortho directions
    const CartesianVector orthoDirection1(eigenVecs.at(1).GetZ() > 0.f ? eigenVecs.at(1) : eigenVecs.at(1) * -1.f);

    //Lists for hits that have a positive and negative projection on the secondary axis
    CaloHitList posCaloHitList, negCaloHitList;

    for (const CaloHit *const pCaloHit3D : inputCaloHitList)
    {
        const CartesianVector pCaloHit3DPosition = pCaloHit3D->GetPositionVector();

        if((pCaloHit3DPosition-centroid).GetDotProduct(orthoDirection1)<0)
        {
            negCaloHitList.push_back(pCaloHit3D);
        }
		else
		{ 
            posCaloHitList.push_back(pCaloHit3D);
		}
    }

    outputCaloHitListsVector.push_back(posCaloHitList);
    outputCaloHitListsVector.push_back(negCaloHitList);

    //Now, run PCA independently on pos and neg hit lists, to refine split

    //Get pos list PCA
    CartesianVector centroidPos(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecsPos;
    LArPcaHelper::EigenValues eigenValuesPos(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(posCaloHitList, centroidPos, eigenValuesPos, eigenVecsPos);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirectionPos(eigenVecsPos.at(0).GetZ() > 0.f ? eigenVecsPos.at(0) : eigenVecsPos.at(0) * -1.f);

    //Get neg list PCA
    CartesianVector centroidNeg(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecsNeg;
    LArPcaHelper::EigenValues eigenValuesNeg(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(negCaloHitList, centroidNeg, eigenValuesNeg, eigenVecsNeg);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirectionNeg(eigenVecsNeg.at(0).GetZ() > 0.f ? eigenVecsNeg.at(0) : eigenVecsNeg.at(0) * -1.f);

    //Get intersection point of the two new principal axes
    CartesianVector intersectionPoint(0.f, 0.f, 0.f);
    float displacementPos(0.f),displacementNeg(0.f);
    LArPointingClusterHelper::GetIntersection(centroidPos,axisDirectionPos,centroidNeg,axisDirectionNeg,intersectionPoint,displacementPos,displacementNeg);

    //Clear pos and neg lists
    posCaloHitList.clear();
    negCaloHitList.clear();

    //Loop over original hit list, check whether it is within smaller cone of pos or neg axis, attach to relevant list
    for (const CaloHit *const pCaloHit3D : inputCaloHitList)
    {
        const float cosConeAxisPos = axisDirectionPos.GetCosOpeningAngle(pCaloHit3D->GetPositionVector()-intersectionPoint);
        const float cosConeAxisNeg = axisDirectionNeg.GetCosOpeningAngle(pCaloHit3D->GetPositionVector()-intersectionPoint);
        if(cosConeAxisPos>cosConeAxisNeg)
            posCaloHitList.push_back(pCaloHit3D);
        else negCaloHitList.push_back(pCaloHit3D);
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimplePCAThreeDClusteringTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
