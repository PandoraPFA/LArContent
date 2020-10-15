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
 */
class PlaneModel: public AbstractModel<3>
{
private:

    Eigen::Vector3f m_direction;
    Eigen::Vector3f m_origin;

public:

    /**
     *  @brief  Project point on to the current line and work out the distance
     *          from the line.
     *
     *  @param param  The parameter to compare to the current line.
     */
    virtual double ComputeDistanceMeasure(const SharedParameter param) const override;

    PlaneModel(ParameterVector inputParams) { Initialize(inputParams); };
    virtual ~PlaneModel() {};

    pandora::CartesianVector GetDirection() const;
    pandora::CartesianVector GetOrigin() const;

    virtual void Initialize(const ParameterVector &inputParams) override;
    virtual std::pair<double, ParameterVector> Evaluate(const ParameterVector &paramsToEval, double threshold) override;
    void operator=(PlaneModel &other);
};

} // namespace lar_content

#endif // LAR_PLANE_MODEL_RANSAC_H
