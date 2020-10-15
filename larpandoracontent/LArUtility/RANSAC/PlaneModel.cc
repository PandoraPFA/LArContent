/**
 *  @file   larpandoracontent/LArUtility/RANSAC/PlaneModel.cc
 *
 *  @brief  Implementation file for the PlaneModel, to be used in RANSAC.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArUtility/RANSAC/AbstractModel.h"
#include "larpandoracontent/LArUtility/RANSAC/PlaneModel.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/RANSACMethod.h"

#include <Eigen/Core>
#include <Eigen/SVD>

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

double PlaneModel::ComputeDistanceMeasure(const SharedParameter param) const
{
    const auto currentPoint = std::dynamic_pointer_cast<RANSACHit>(param);
    if(currentPoint == nullptr)
        throw std::runtime_error("PlaneModel::ComputeDistanceMeasure() - Passed parameter are not of type RANSACHit.");

    const RANSACHit hit = *currentPoint;
    const Eigen::Vector3f currentPos = hit.GetVector() - m_origin;

    const Eigen::Vector3f b = currentPos.dot(m_direction) * m_direction;
    const double distance = (currentPos - b).norm();

    return distance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector PlaneModel::GetDirection() const
{
    return pandora::CartesianVector(m_direction[0], m_direction[1], m_direction[2]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::CartesianVector PlaneModel::GetOrigin() const
{
    return pandora::CartesianVector(m_origin[0], m_origin[1], m_origin[2]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PlaneModel::Initialize(const ParameterVector &inputParams)
{
    if(inputParams.size() != 3)
        throw std::runtime_error("PlaneModel - Number of input parameters does not match minimum number required for this model.");

    Eigen::Matrix3f m;
    Eigen::Vector3f totals = {0.0, 0.0, 0.0};

    for (auto param : inputParams)
    {
        auto currentPoint = std::dynamic_pointer_cast<RANSACHit>(param);

        if(currentPoint == nullptr)
            throw std::runtime_error("PlaneModel - inputParams type mismatch. It is not a RANSACHit.");

        totals += (*currentPoint).GetVector();
    }

    m_origin = totals / inputParams.size();

    for (unsigned int i = 0; i < inputParams.size(); ++i)
    {
        auto currentPoint = *std::dynamic_pointer_cast<RANSACHit>(inputParams[i]);
        Eigen::Vector3f shiftedPoint = Eigen::Vector3f(currentPoint.GetVector() - m_origin);
        m.row(i) = shiftedPoint;
    }

    const Eigen::JacobiSVD<Eigen::MatrixXf> svd(m, Eigen::ComputeThinV);
    m_direction = svd.matrixV().col(0); // TODO: Check this, returns matrix not vector.
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::pair<double, ParameterVector> PlaneModel::Evaluate(const ParameterVector &paramsToEval, double threshold)
{
    ParameterVector inliers;
    const float totalParams = paramsToEval.size();

    for(auto& param : paramsToEval)
    {
        if(ComputeDistanceMeasure(param) < threshold)
            inliers.push_back(param);
    }

    // TODO: This could actually take into account the tool -> i.e. weight against the iffy tools.
    const double inlierFraction = inliers.size() / totalParams;
    return std::make_pair(inlierFraction, inliers);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PlaneModel::operator=(PlaneModel &other)
{
    const pandora::CartesianVector origin = other.GetOrigin();
    const Eigen::Vector3f newOrigin(origin.GetX(), origin.GetY(), origin.GetZ());
    m_origin = newOrigin;

    const pandora::CartesianVector direction = other.GetDirection();
    const Eigen::Vector3f newDirection(direction.GetX(), direction.GetY(), direction.GetZ());
    m_direction = newDirection;
}

} // namespace lar_content
