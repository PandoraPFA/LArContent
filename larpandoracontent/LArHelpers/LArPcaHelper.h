/**
 *  @file   larpandoracontent/LArHelpers/LArLArPcaHelper.h
 *
 *  @brief  Header file for the principal curve analysis helper class.
 *
 *  $Log: $
 */
#ifndef LAR_PCA_HELPER_H
#define LAR_PCA_HELPER_H 1


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

    /**
     *  @brief  Run principal component analysis using input calo hits (TPC_VIEW_U,V,W or TPC_3D; all treated as 3D points)
     *
     *  @param  caloHitList the calo hit list
     *  @param  centroid to receive the centroid position
     *  @param  outputEigenValues to receive the eigen values
     *  @param  outputEigenVectors to receive the eigen vectors
     */
    static void RunPca(const pandora::CaloHitList &caloHitList, pandora::CartesianVector &centroid, EigenValues &outputEigenValues,
        EigenVectors &outputEigenVectors);

    /**
     *  @brief  Run principal component analysis using input Cartesian vectors (TPC_VIEW_U,V,W or TPC_3D; all treated as 3D points)
     *
     *  @param  pointVector the vector of positions
     *  @param  centroid to receive the centroid position
     *  @param  outputEigenValues to receive the eigen values
     *  @param  outputEigenVectors to receive the eigen vectors
     */
    static void RunPca(const pandora::CartesianPointVector &pointVector, pandora::CartesianVector &centroid, EigenValues &outputEigenValues,
        EigenVectors &outputEigenVectors);
};

} // namespace lar_content

#endif // #ifndef LAR_PCA_HELPER_H