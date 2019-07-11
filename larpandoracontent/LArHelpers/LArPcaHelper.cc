/**
 *  @file   larpandoracontent/LArHelpers/LArPcaHelper.cc
 *
 *  @brief  Implementation of the principal curve analysis helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArObjectHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include <Eigen_tf/Dense>

using namespace pandora;

namespace lar_content
{

template <typename T>
void LArPcaHelper::RunPca(const T &t, CartesianVector &centroid, EigenValues &outputEigenValues, EigenVectors &outputEigenVectors)
{
    WeightedPointVector weightedPointVector;

    for (const auto &point : t)
        weightedPointVector.push_back(std::make_pair(LArObjectHelper::TypeAdaptor::GetPosition(point), 1.));

    return LArPcaHelper::RunPca(weightedPointVector, centroid, outputEigenValues, outputEigenVectors);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArPcaHelper::RunPca(const WeightedPointVector &pointVector, CartesianVector &centroid, EigenValues &outputEigenValues, EigenVectors &outputEigenVectors)
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
    double sumWeight(0.);

    for (const WeightedPoint &weightedPoint : pointVector)
    {
        const CartesianVector &point(weightedPoint.first);
        const double weight(weightedPoint.second);

        if (weight < 0.)
        {
            std::cout << "LArPcaHelper::RunPca - negative weight found" << std::endl;
            throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
        }

        meanPosition[0] += static_cast<double>(point.GetX()) * weight;
        meanPosition[1] += static_cast<double>(point.GetY()) * weight;
        meanPosition[2] += static_cast<double>(point.GetZ()) * weight;
        sumWeight += weight;
    }

    if (std::fabs(sumWeight) < std::numeric_limits<double>::epsilon())
    {
        std::cout << "LArPcaHelper::RunPca - sum of weights is zero" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    meanPosition[0] /= sumWeight;
    meanPosition[1] /= sumWeight;
    meanPosition[2] /= sumWeight;
    centroid = CartesianVector(meanPosition[0], meanPosition[1], meanPosition[2]);

    // Define elements of our covariance matrix
    double xi2(0.);
    double xiyi(0.);
    double xizi(0.);
    double yi2(0.);
    double yizi(0.);
    double zi2(0.);

    for (const WeightedPoint &weightedPoint : pointVector)
    {
        const CartesianVector &point(weightedPoint.first);
        const double weight(weightedPoint.second);
        const double x(static_cast<double>((point.GetX()) - meanPosition[0]));
        const double y(static_cast<double>((point.GetY()) - meanPosition[1]));
        const double z(static_cast<double>((point.GetZ()) - meanPosition[2]));

        xi2  += x * x * weight;
        xiyi += x * y * weight;
        xizi += x * z * weight;
        yi2  += y * y * weight;
        yizi += y * z * weight;
        zi2  += z * z * weight;
    }

    // Using Eigen_tf package
    Eigen_tf::Matrix3f sig;

    sig <<  xi2, xiyi, xizi,
           xiyi,  yi2, yizi,
           xizi, yizi,  zi2;

    sig *= 1. / sumWeight;

    Eigen_tf::SelfAdjointEigenSolver<Eigen_tf::Matrix3f> eigenMat(sig);

    if (eigenMat.info() != Eigen_tf::ComputationInfo::Success)
    {
        std::cout << "LArPcaHelper::RunPca - decomposition failure, nThreeDHits = " << pointVector.size() << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    typedef std::pair<float, size_t> EigenValColPair;
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
    const Eigen_tf::Matrix3f &eigenVecs(eigenMat.eigenvectors());

    for (const EigenValColPair &pair : eigenValColVector)
        outputEigenVectors.emplace_back(eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template void LArPcaHelper::RunPca(const CartesianPointVector &, CartesianVector &, EigenValues &, EigenVectors &);
template void LArPcaHelper::RunPca(const CaloHitList &, CartesianVector &, EigenValues &, EigenVectors &);

} // namespace lar_content
