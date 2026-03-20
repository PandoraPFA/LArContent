/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/PlaneSolverAlgorithm.h
 *
 *  @brief  The Plane Solver algorithm takes the collection of 2D hits across all planes and constructs an optimal set of triplet and doublet
 *          relationships between them. These relationships are persisted via the EventContext object so that they are available to subsequent
 *          algorithms.
 *
 *  $Log: $
 */
#ifndef LAR_PLANE_SOLVER_ALGORITHM_H
#define LAR_PLANE_SOLVER_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  An algorithm to construct hit triplets and doublets representing candidate 3D hits
 */
class PlaneSolverAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PlaneSolverAlgorithm();

private:
    struct Volume
    {
        unsigned int m_tpc;
        unsigned int m_childVolume;

        bool operator==(const Volume& other) const
        {
            return m_tpc == other.m_tpc && m_childVolume == other.m_childVolume;
        }
    };

    struct VolumeHash
    {
        std::size_t operator()(const Volume& v) const noexcept
        {
            return std::hash<int>{}(v.m_tpc) ^ (std::hash<int>{}(v.m_childVolume) << 1);
        }
    };

    typedef std::unordered_map<pandora::HitType, pandora::CaloHitVector> PlaneToHitsMap;
    typedef std::unordered_map<Volume, PlaneToHitsMap, VolumeHash> VolumeToReadoutMap;
    typedef std::vector<std::vector<int>> IndexMatrix;
    typedef std::vector<std::vector<float>> CostMatrix;

    pandora::StatusCode Run();

    /**
     *  @brief  Fills the map of the 2D hits, keyed by volume
     *
     *  @param  caloHitList the list of 2D hits to be mapped
     */
    void FillHitMap(const pandora::CaloHitList &caloHitList);

    /**
     *  @brief  Solve for the optimal set of triplet and doublet relationships between the 2D hits
     */
    void Solve() const;

    /**
     *  @brief  Compute the cost matrix for the hits in a given slice, where the cost is based on the chi-squared value for triplets and doublets of hits
     *
     *  @param  planetoHitsMap the map of the 2D hits in the slice, keyed by view
     *  @param  unmatchedCost the cost to be assigned to unmatched hits
     *  @param  bestK the matrix to be filled with the index of the best matching hit in the W view for each UV pair
     *
     *  @return The cost matrix for the hits in the slice
     */
    CostMatrix ComputeCostMatrix(const PlaneToHitsMap &planeToHitsMap, const float unmatchedCost, IndexMatrix &bestK) const;

    /**
     *  @brief  Implementation of the Kunhne-Munkres (aka Hungarian) algorithm to solve the optimal matching between UV pairs and W hits.
     *          For each row, attempt to assign it and if there is a conflict, reroute via an augmenting path.
     *
     *  @param  cost the cost matrix for the hits in a volume
     *
     *  @return The optimal matching between UV pairs and W hits, where the vector elements align with the U hits and indicate matched index
     *          in the V hits (or -1 if unmatched)
     */
    pandora::IntVector KuhneMunkres(const CostMatrix& cost) const;

    /**
     *  @brief  Create a 3D hit from a double of U and V hits
     *
     *  @param  pCaloHitU the U view hit
     *  @param  pCaloHitV the V view hit
     *  @param[out] pCaloHit3D the 3D hit to be created
     */
    void CreateThreeDHit(const pandora::CaloHit *const pCaloHitU, const pandora::CaloHit *const pCaloHitV, const pandora::CaloHit *&pCaloHit3D) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName; ///< The name of the calo hit list containing all of the 2D hits
    VolumeToReadoutMap m_planeToReadoutMap; ///< A map from volume, to readout plane, to the 2D hits in that plane
};

} // namespace lar_content

#endif // #ifndef LAR_PLANE_SOLVER_ALGORITHM_H
