/**
 *  @file   larpandoracontent/LArReclustering/SimplePCAThreeDClusteringTool.cc
 *
 *  @brief  Implementation file for the simple PCA-based clustering tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArReclustering/SimplePCAThreeDClusteringTool.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

using namespace pandora;

namespace lar_content
{

SimplePCAThreeDClusteringTool::SimplePCAThreeDClusteringTool()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<std::reference_wrapper<pandora::CaloHitList>> SimplePCAThreeDClusteringTool::Run(const Algorithm *const /*pAlgorithm*/, std::reference_wrapper<pandora::CaloHitList> &inputCaloHitList)
{
    std::vector<std::reference_wrapper<pandora::CaloHitList>> newCaloHitListsVector;

    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(inputCaloHitList.get(), centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : inputCaloHitList.get())
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));

    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

    // Now define ortho directions
    const CartesianVector seedDirection((axisDirection.GetX() < std::min(axisDirection.GetY(), axisDirection.GetZ())) ? CartesianVector(1.f, 0.f, 0.f) :
        (axisDirection.GetY() < std::min(axisDirection.GetX(), axisDirection.GetZ())) ? CartesianVector(0.f, 1.f, 0.f) : CartesianVector(0.f, 0.f, 1.f));
    const CartesianVector orthoDirection1(seedDirection.GetCrossProduct(axisDirection).GetUnitVector());

    //Lists for hits that have a positive and negative projection on the secondary axis
    CaloHitList posCaloHitList, negCaloHitList;

    for (const CaloHit *const pCaloHit : inputCaloHitList.get())
    {
        const CartesianVector pCaloHitPosition = pCaloHit->GetPositionVector();

        if((pCaloHitPosition-centroid).GetDotProduct(centroid+orthoDirection1)<0)
        {
            negCaloHitList.push_back(pCaloHit);
        }
		else
		{ 
            posCaloHitList.push_back(pCaloHit);
		}
    }

    newCaloHitListsVector.push_back(posCaloHitList);
    newCaloHitListsVector.push_back(negCaloHitList);

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
    for (const CaloHit *const pCaloHit3D : inputCaloHitList.get())
    {
        const float cosConeAxisPos = axisDirectionPos.GetCosOpeningAngle(pCaloHit3D->GetPositionVector()-intersectionPoint);
        const float cosConeAxisNeg = axisDirectionNeg.GetCosOpeningAngle(pCaloHit3D->GetPositionVector()-intersectionPoint);
        if(cosConeAxisPos>cosConeAxisNeg)
            posCaloHitList.push_back(pCaloHit3D);
        else negCaloHitList.push_back(pCaloHit3D);
    }

    return newCaloHitListsVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SimplePCAThreeDClusteringTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
