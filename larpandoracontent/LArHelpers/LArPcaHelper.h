/**
 *  @file   larpandoracontent/LArHelpers/LArPcaHelper.h
 *
 *  @brief  Header file for the principal curve analysis helper class.
 *
 *  $Log: $
 */
#ifndef LAR_PCA_HELPER_H
#define LAR_PCA_HELPER_H 1

#include "Objects/CartesianVector.h"

namespace lar_content
{

/**
 *  @brief  LArPcaHelper class
 */
class LArPcaHelper
{
public:
    typedef pandora::CartesianVector EigenValues;
    typedef std::vector<pandora::CartesianVector> EigenVectors;
    typedef std::pair<const pandora::CartesianVector, double> WeightedPoint;
    typedef std::vector<WeightedPoint> WeightedPointVector;

    /**
     *  @brief  Run principal component analysis using input calo hits (TPC_VIEW_U,V,W or TPC_3D; all treated as 3D points)
     *
     *  @param  t the input information
     *  @param  centroid to receive the centroid position
     *  @param  outputEigenValues to receive the eigen values
     *  @param  outputEigenVectors to receive the eigen vectors
     */
    template <typename T>
    static void RunPca(const T &t, pandora::CartesianVector &centroid, EigenValues &outputEigenValues, EigenVectors &outputEigenVectors);

    /**
     *  @brief  Run principal component analysis using weighted input Cartesian vectors (TPC_VIEW_U,V,W or TPC_3D; all treated as 3D points)
     *
     *  @param  pointVector a vector of pairs of positions and weights
     *  @param  centroid to receive the centroid position
     *  @param  outputEigenValues to receive the eigen values
     *  @param  outputEigenVectors to receive the eigen vectors
     */
    static void RunPca(const WeightedPointVector &pointVector, pandora::CartesianVector &centroid, EigenValues &outputEigenValues,
        EigenVectors &outputEigenVectors);
};

} // namespace lar_content

#endif // #ifndef LAR_PCA_HELPER_H
