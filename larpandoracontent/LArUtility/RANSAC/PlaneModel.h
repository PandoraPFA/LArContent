/**
 *  @file   larpandoracontent/LArUtility/RANSAC/PlaneModel.h
 *
 *  @brief  Header file for the PlaneModel, to be used in RANSAC.
 *
 *  $Log: $
 */
#ifndef LAR_PLANE_MODEL_RANSAC_H
#define LAR_PLANE_MODEL_RANSAC_H 1

#include "larpandoracontent/LArUtility/RANSAC/AbstractModel.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include <Eigen/Core>
#include <Eigen/SVD>

namespace lar_content
{

/**
 *  @brief  Class that implements a PlaneModel, to be fit using RANSAC.
 *          This is a model that fits a plane through 3D points, requiring at least 3 3D points to be built.
 *          The final model will be the plane that minimises the distance between the plane and the points.
 */
class PlaneModel: public AbstractModel<3>
{
private:

    Eigen::Vector3f m_direction;
    Eigen::Vector3f m_origin;

public:

    /**
     *  @brief  Project point on to the current plane and work out the distance
     *          from the plane.
     *
     *  @param param  The parameter to compare to the current plane.
     *
     *  @return double The distance of the given 3D point from the plane.
     */
    virtual double ComputeDistanceMeasure(const SharedParameter param) const override;

    /**
     *  @brief  Default constructor.
     *
     *  @param inputParams A vector of all the input 3D points, as SharedParameters.
     */
    PlaneModel(ParameterVector inputParams) { Initialize(inputParams); };

    /**
     *  @brief  Destructor
     */
    virtual ~PlaneModel() {};

    /**
     *  @brief  Get the direction of the current plane.
     *
     *  @return pandora::CartesianVector A vector corresponding to the 3D direction on the plane.
     */
    pandora::CartesianVector GetDirection() const;

    /**
     *  @brief  Get the origin of the current plane.
     *
     *  @return pandora::CartesianVector A vector corresponding to the 3D origin on the plane.
     */
    pandora::CartesianVector GetOrigin() const;

    /**
     *  @brief  Initialise the current model with the given parameters. This will produce the plane from the given 3D points.
     */
    virtual void Initialize(const ParameterVector &inputParams) override;

    /**
     *  @brief  Evaluate the current plane against a vector of points.
     *
     *  @param paramsToEval All of the 3D points to compare to the current plane.
     *  @param threshold A distance threshold value, for how close a point must be to the plane to be considered fitted.
     *
     *  @return std::pair<double, ParameterVector> A pair, containing the fraction of fitted points, and all the points
     *                                             that fit.
     */
    virtual std::pair<double, ParameterVector> Evaluate(const ParameterVector &paramsToEval, double threshold) override;

    /**
     *  @brief  Assignment operator, to set a model to the given model.
     *
     *  @param other The other PlaneModel to assign the current model to.
     */
    void operator=(PlaneModel &other);
};

} // namespace lar_content

#endif // LAR_PLANE_MODEL_RANSAC_H
