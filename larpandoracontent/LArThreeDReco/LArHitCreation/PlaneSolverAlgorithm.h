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

    struct UVPair
    {
        int m_uIndex;
        int m_vIndex;
    };

    struct Triplet
    {
        int m_uIndex;
        int m_vIndex;
        int m_wIndex;
    };

    typedef std::unordered_map<pandora::HitType, pandora::CaloHitVector> PlaneToHitsMap;
    typedef std::unordered_map<Volume, PlaneToHitsMap, VolumeHash> VolumeToReadoutMap;
    typedef std::vector<std::vector<int>> IndexMatrix;
    typedef std::vector<std::vector<float>> CostMatrix;
    typedef std::vector<UVPair> PairVector;
    typedef std::vector<Triplet> TripletVector;

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
     *  @brief  Compute the cost matrix for the hits, where the cost is based on the chi-squared value for hit triplets.
     *          This function computes the cost of matching a UV pair, using the W hits to provide a constraint, but allows for the W hits
     *          to be used more than once.
     *
     *  @param  planetoHitsMap the map of the 2D hits, keyed by view
     *  @param  unmatchedCost the cost to be assigned to unmatched hits
     *
     *  @return The cost matrix for the hits in the slice
     */
    CostMatrix ComputeCostMatrix(const PlaneToHitsMap &planeToHitsMap, const float unmatchedCost) const;

    /**
     *  @brief  Compute the cost matrix for the hits, where the cost is based on the chi-squared value for hit triplets.
     *          This function starts with the UV pairs from the first pass of the Kuhne-Munkres algorithm, and computes the cost of matching
     *          each UV pair to a unique W hit, constructing triplet matches where possible, but retaining the doublet matches otherwise.
     *
     *  @param  uVPairs the list of UV pairs to be considered for triplet matching
     *  @param  planetoHitsMap the map of the 2D hits, keyed by view
     *  @param  unmatchedCost the cost to be assigned to unmatched hits
     *
     *  @return The cost matrix for the hits in the slice
     */
    CostMatrix ComputeTripletCostMatrix(const PairVector &uvPairs, const PlaneToHitsMap &planeToHitsMap, const float unmatchedCost) const;

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
     *  @brief  Build the list of UV pairs based on the assignment of V hits to U hits from the Kuhne-Munkres algorithm
     *
     *  @param  assignment the output of the Kuhne-Munkres algorithm, indicating the index of the assigned V hit for each U hit (or -1 if unmatched)
     *  @param  nU the number of U hits
     *  @param  nV the number of V hits
     *
     *  @return The list of UV pairs
     */
    PairVector BuildUVPairs(const pandora::IntVector &assignment, int nU, int nV) const;

    /**
     *  @brief  Build the list of triplets based on the assignment of W hits to UV pairs from the second pass of the Kuhne-Munkres algorithm
     *
     *  @param  uvPairs the list of UV pairs to be considered for triplet matching
     *  @param  assignment the output of the Kuhne-Munkres algorithm, indicating the index of the assigned W hit for each UV pair (or -1 if unmatched)
     *  @param  nW the number of W hits
     *
     *  @return The list of triplets, where unmatched UV pairs are indicated with a W index of -1
     */
    TripletVector BuildTriplets(const PairVector& uvPairs, const pandora::IntVector& assignment, int nW) const;

    /**
     *  @brief  Create a 3D hit from a double of U and V hits
     *
     *  @param  pCaloHitU the U view hit
     *  @param  pCaloHitV the V view hit
     *  @param[out] pCaloHit3D the 3D hit to be created
     */
    void CreateThreeDHit(const pandora::CaloHit *const pCaloHitU, const pandora::CaloHit *const pCaloHitV, const pandora::CaloHit *&pCaloHit3D) const;

    /**
     *  @brief  Create a 3D hit from a triplet of U, V and W hits
     *
     *  @param  pCaloHitU the U view hit
     *  @param  pCaloHitV the V view hit
     *  @param  pCaloHitW the W view hit
     *  @param[out] pCaloHit3D the 3D hit to be created
     */
    void CreateThreeDHit(const pandora::CaloHit *const pCaloHitU, const pandora::CaloHit *const pCaloHitV, const pandora::CaloHit *const pCaloHitW,
        const pandora::CaloHit *&pCaloHit3D) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName; ///< The name of the calo hit list containing all of the 2D hits
    VolumeToReadoutMap m_planeToReadoutMap; ///< A map from volume, to readout plane, to the 2D hits in that plane
};

} // namespace lar_content

#endif // #ifndef LAR_PLANE_SOLVER_ALGORITHM_H
