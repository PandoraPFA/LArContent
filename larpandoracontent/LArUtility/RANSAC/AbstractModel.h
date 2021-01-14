// Code is from https://github.com/drsrinathsridhar/GRANSAC
// Original license file can be found in larpandoracontent/LArUtility/RANSAC/LICENSE.

#ifndef LAR_RANSAC_TEMPLATES_H
#define LAR_RANSAC_TEMPLATES_H 1

#include <iostream>
#include <stdexcept>
#include <vector>
#include <array>
#include <memory>

namespace lar_content
{
    /**
     *  @brief  Each abstract model is made of abstract parameters. Could be anything from a point (that make a 2D line or 3D plane
     *          or image correspondences) to a line
     */
    class AbstractParameter
    {
    public:
        /**
         *  @brief  Dummy Destructor
         */
        virtual ~AbstractParameter(void) {};
    };

    typedef std::shared_ptr<AbstractParameter> SharedParameter;
    typedef std::vector<SharedParameter> ParameterVector;

    /**
     *  @brief  Abstract model type for generic RANSAC model fitting
     */
    template <int t_NumParams>
    class AbstractModel
    {
    public:

        /**
         *  @brief  Initialise the current model.
         *
         *  @param  inputParams The required parameters, to set the size and limits of the model.
         */
        virtual void Initialize(const ParameterVector &inputParams) = 0;

        /**
         *  @brief  Evaluate the current model, finding which given parameters fit it to the given threshold.
         *
         *  @param  evaluateParams The parameters to evaluate the model for.
         *  @param  threshold The threshold that the current parameter must be close to the model to be considered part of the model.
         *
         *  @return std::pair<double, ParameterVector> A pair, containing the fraction of fitted parameters, and all the parameters
         *                                             that fit.
         */
        virtual std::pair<double, ParameterVector> Evaluate(const ParameterVector &evaluateParams, double threshold) = 0;

        /**
         *  @brief  Compute the distance that the given parameter is away from the current model.
         *
         *  @param  param The parameter to evaluate against the current model.
         *
         *  @return double A measure of how far the given parameter is away from the current model.
         */
        virtual double ComputeDistanceMeasure(const SharedParameter param) const = 0;
    };
} // namespace lar_content

#endif // LAR_RANSAC_TEMPLATES_H
