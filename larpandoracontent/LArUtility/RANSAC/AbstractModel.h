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
    // Each abstract model is made of abstract parameters
    // Could be anything from a point (that make a 2D line or 3D plane or image correspondences) to a line
    class AbstractParameter
    {
    public:
        virtual ~AbstractParameter(void) {}; // To make this polymorphic we add dummy destructor
    };

    typedef std::shared_ptr<AbstractParameter> SharedParameter;
    typedef std::vector<SharedParameter> ParameterVector;

    // Abstract model type for generic RANSAC model fitting
    template <int t_NumParams> /* Minimum number of parameters required to define this model*/
    class AbstractModel
    {
    protected:
        std::array<SharedParameter, t_NumParams> m_MinModelParams;

    public:
        virtual void Initialize(const ParameterVector &inputParams) = 0;
        virtual std::pair<double, ParameterVector> Evaluate(const ParameterVector &evaluateParams, double threshold) = 0;
        virtual double ComputeDistanceMeasure(SharedParameter param) = 0;

        virtual std::array<SharedParameter, t_NumParams> GetModelParams(void) { return m_MinModelParams; };
    };
} // namespace lar_content

#endif // LAR_RANSAC_TEMPLATES_H
