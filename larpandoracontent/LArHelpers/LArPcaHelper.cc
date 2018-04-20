/**
 *  @file   larpandoracontent/LArHelpers/LArPcaHelper.cc
 *
 *  @brief  Implementation of the principal curve analysis helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include <Eigen/Dense>

using namespace pandora;

namespace lar_content
{

void LArPcaHelper::RunPca(const CaloHitList &caloHitList, CartesianVector &centroid, EigenValues &outputEigenValues, EigenVectors &outputEigenVectors)
{
    CartesianPointDoublePairVector cartesianPointDoublePairVector;

    for (const CaloHit *const pCaloHit : caloHitList)
        cartesianPointDoublePairVector.push_back(std::make_pair(pCaloHit->GetPositionVector(), 1.0));

    return LArPcaHelper::RunPca(cartesianPointDoublePairVector, centroid, outputEigenValues, outputEigenVectors);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPcaHelper::RunPca(const CartesianPointVector &pointVector, CartesianVector &centroid, EigenValues &outputEigenValues, EigenVectors &outputEigenVectors)
{
    CartesianPointDoublePairVector cartesianPointDoublePairVector;

    for (const CartesianVector &point : pointVector)
        cartesianPointDoublePairVector.push_back(std::make_pair(point, 1.0));

    return LArPcaHelper::RunPca(cartesianPointDoublePairVector, centroid, outputEigenValues, outputEigenVectors);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPcaHelper::RunPca(const CartesianPointDoublePairVector &pointVector, CartesianVector &centroid, EigenValues &outputEigenValues, EigenVectors &outputEigenVectors)
{
    // The steps are:
    // 1) do a mean normalization of the input vec points
    // 2) compute the covariance matrix
    // 3) run the SVD
    // 4) extract the eigen vectors and values

    // Run through the point vector and get the mean position of all points
    if (pointVector.empty())
    {
        std::cout << "LArPcaHelper::RunPca - no three dimensional hits provided" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    double meanPosition[3] = {0., 0., 0.};

    for (const auto &iter : pointVector)
    {
        const CartesianVector &point(iter.first);
        meanPosition[0] += point.GetX();
        meanPosition[1] += point.GetY();
        meanPosition[2] += point.GetZ();
    }

    const double nThreeDHitsAsDouble(static_cast<double>(pointVector.size()));
    meanPosition[0] /= nThreeDHitsAsDouble;
    meanPosition[1] /= nThreeDHitsAsDouble;
    meanPosition[2] /= nThreeDHitsAsDouble;
    centroid = CartesianVector(meanPosition[0], meanPosition[1], meanPosition[2]);

    // Define elements of our covariance matrix
    double xi2(0.);
    double xiyi(0.);
    double xizi(0.);
    double yi2(0.);
    double yizi(0.);
    double zi2(0.);
    double weightSum(0.);

    for (const auto &iter : pointVector)
    {
        const CartesianVector &point(iter.first);
        const double weight(iter.second);
        const double x((point.GetX() - meanPosition[0]) * weight);
        const double y((point.GetY() - meanPosition[1]) * weight);
        const double z((point.GetZ() - meanPosition[2]) * weight);

        xi2  += x * x;
        xiyi += x * y;
        xizi += x * z;
        yi2  += y * y;
        yizi += y * z;
        zi2  += z * z;
        weightSum += weight * weight;
    }

    // Using Eigen package
    Eigen::Matrix3f sig;

    sig <<  xi2, xiyi, xizi,
           xiyi,  yi2, yizi,
           xizi, yizi,  zi2;

    if (std::fabs(weightSum) < std::numeric_limits<double>::epsilon())
    {
        std::cout << "LArPcaHelper::RunPca - weight of three dimensional hits = " << weightSum << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    sig *= 1. / weightSum;

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);

    if (eigenMat.info() != Eigen::ComputationInfo::Success)
    {
        std::cout << "LArPcaHelper::RunPca - decomposition failure, nThreeDHits = " << pointVector.size() << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    typedef std::pair<float,size_t> EigenValColPair;
    typedef std::vector<EigenValColPair> EigenValColVector;

    EigenValColVector eigenValColVector;
    const auto &resultEigenMat(eigenMat.eigenvalues());
    eigenValColVector.emplace_back(resultEigenMat(0), 0);
    eigenValColVector.emplace_back(resultEigenMat(1), 1);
    eigenValColVector.emplace_back(resultEigenMat(2), 2);

    std::sort(eigenValColVector.begin(), eigenValColVector.end(), [](const EigenValColPair &left, const EigenValColPair &right){return left.first > right.first;});

    // Get the eigen values
    outputEigenValues = CartesianVector(eigenValColVector.at(0).first, eigenValColVector.at(1).first, eigenValColVector.at(2).first);

    // Get the principal axes
    const Eigen::Matrix3f &eigenVecs(eigenMat.eigenvectors());

    for (const EigenValColPair &pair : eigenValColVector)
        outputEigenVectors.emplace_back(eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second));
}

} // namespace lar_content
