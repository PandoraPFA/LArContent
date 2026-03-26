/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/PlaneSolverAlgorithm.cc
 *
 *  @brief  Implementation of the three dimensional hit creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArThreeDReco/LArHitCreation/PlaneSolverAlgorithm.h"

using namespace pandora;

namespace lar_content
{

PlaneSolverAlgorithm::PlaneSolverAlgorithm() :
    m_caloHitListName{"CaloHitList2D"}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PlaneSolverAlgorithm::Run()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    this->FillHitMap(*pCaloHitList);
    this->Solve();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PlaneSolverAlgorithm::FillHitMap(const CaloHitList &caloHitList)
{
    for (const auto *pCaloHit : caloHitList)
    {
        const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pCaloHit));
        if (!pLArCaloHit)
        {
            std::cout << "PlaneSolverAlgorithm - non LArCaloHit found in calo hit list: " << m_caloHitListName << std::endl;
            throw StatusCodeException(STATUS_CODE_FAILURE);
        }
        const unsigned int tpc{pLArCaloHit->GetLArTPCVolumeId()};
        const unsigned int childVolume{pLArCaloHit->GetDaughterVolumeId()};
        const Volume volume{tpc, childVolume};

        const HitType view(pLArCaloHit->GetHitType());
        m_planeToReadoutMap[volume][view].emplace_back(pLArCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PlaneSolverAlgorithm::Solve() const
{
    CaloHitList hits3D;
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    // Loop over each read out volume (e.g. an APA in a horizontal drift detector) and solve for the optimal set of triplet and doublet
    // relationships between the 2D hits in each unit.
    for (const auto &[_, readout] : m_planeToReadoutMap)
    {
        CostMatrix costMatrix{this->ComputeCostMatrix(readout, 100.f)};
        const IntVector assignment{this->KuhneMunkres(costMatrix)};
        int nHitsU{static_cast<int>(readout.at(TPC_VIEW_U).size())};
        int nHitsV{static_cast<int>(readout.at(TPC_VIEW_V).size())};
        const PairVector pairs{this->BuildUVPairs(assignment, nHitsU, nHitsV)};
        CostMatrix tripletCostMatrix{this->ComputeTripletCostMatrix(pairs, readout, 100.f)};
        const IntVector tripletAssignment{this->KuhneMunkres(tripletCostMatrix)};
        int nHitsW{static_cast<int>(readout.at(TPC_VIEW_W).size())};
        const TripletVector triplets{this->BuildTriplets(pairs, tripletAssignment, nHitsW)};
        //CaloHitList hitsU{readout.at(TPC_VIEW_U).begin(), readout.at(TPC_VIEW_U).end()};
        //CaloHitList hitsV{readout.at(TPC_VIEW_V).begin(), readout.at(TPC_VIEW_V).end()};
        //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hitsU, "U", RED));
        //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hitsV, "V", GREEN));
        int nMatchedU{0}, nMatchedV{0};
        std::cout << "UV pairs with non-unique W" << std::endl;
        std::cout << "-----------------------------" << std::endl;
        for (size_t i = 0; i < assignment.size(); ++i)
        {
            if (assignment[i] < 0)
                continue;
            const size_t j{static_cast<size_t>(assignment[i])};
            //if (i < hitsU.size())
            //{
            //    const CartesianVector u(readout.at(TPC_VIEW_U)[i]->GetPositionVector());
            //    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &u, "U", RED, 2));
            //}
            //if (j < hitsV.size())
            //{
            //    const CartesianVector v(readout.at(TPC_VIEW_V)[j]->GetPositionVector());
            //    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &v, "V", GREEN, 2));
            //}
            //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

            if (i < readout.at(TPC_VIEW_U).size() && j < readout.at(TPC_VIEW_V).size())
            {
                ++nMatchedU;
                ++nMatchedV;
                std::cout << "Assignment: i = " << i << "(" << readout.at(TPC_VIEW_U).size() << "), j = " << j << " (" <<
                    readout.at(TPC_VIEW_V).size() << ")" << std::endl;

                //const CaloHit *pHit3D{nullptr};
                //this->CreateThreeDHit(readout.at(TPC_VIEW_U)[i], readout.at(TPC_VIEW_V)[j], pHit3D);
                //hits3D.emplace_back(pHit3D);
            }
        }
        std::cout << "Matched " << nMatchedU << " out of " << nHitsU << " U hits, and " << nMatchedV << " out of " << nHitsV << " V hits. " <<
            "The W view has " << nHitsW << std::endl;
        std::cout << "-----------------------------" << std::endl;
        std::cout << "Triplets" << std::endl;
        std::cout << "-----------------------------" << std::endl;
        for (size_t i = 0; i < triplets.size(); ++i)
        {
            std::cout << "Triplet: i = " << triplets[i].m_uIndex << ", j = " << triplets[i].m_vIndex << ", k = " << triplets[i].m_wIndex << std::endl;
            const CaloHit *pHit3D{nullptr};
            if (triplets[i].m_wIndex >= 0)
                this->CreateThreeDHit(readout.at(TPC_VIEW_U)[triplets[i].m_uIndex], readout.at(TPC_VIEW_V)[triplets[i].m_vIndex],
                    readout.at(TPC_VIEW_W)[triplets[i].m_wIndex], pHit3D);
            else
                this->CreateThreeDHit(readout.at(TPC_VIEW_U)[triplets[i].m_uIndex], readout.at(TPC_VIEW_V)[triplets[i].m_vIndex], pHit3D);
            hits3D.emplace_back(pHit3D);
        }
        std::cout << "-----------------------------" << std::endl;
        std::cout << "Number of triplets: " << triplets.size() << std::endl;
    }
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, hits3D, "KMHitList"));
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &hits3D, "3D", BLACK));
    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

PlaneSolverAlgorithm::CostMatrix PlaneSolverAlgorithm::ComputeCostMatrix(const PlaneToHitsMap &planeToHitsMap, const float unmatchedCost) const
{
    const CaloHitVector &uHits(planeToHitsMap.at(TPC_VIEW_U));
    const CaloHitVector &vHits(planeToHitsMap.at(TPC_VIEW_V));
    const CaloHitVector &wHits(planeToHitsMap.at(TPC_VIEW_W));
    const int nU{static_cast<int>(uHits.size())};
    const int nV{static_cast<int>(vHits.size())};
    const int nW{static_cast<int>(wHits.size())};
    const int N{std::max(nU, nV)};

    // Initialize the cost matrix with the cost for unmatched hits.
    CostMatrix C(N, FloatVector(N, unmatchedCost));

    // Compute all chi-squared values
    for (int i = 0; i < nU; ++i)
    {
        const CartesianVector u(uHits[i]->GetPositionVector());
        const float dx_u(0.5f * uHits[i]->GetCellSize1());
        const float xMin_u(u.GetX() - dx_u);
        const float xMax_u(u.GetX() + dx_u);
        
        for (int j = 0; j < nV; ++j)
        {
            const CartesianVector v(vHits[j]->GetPositionVector());
            const float dx_v(0.5f * vHits[j]->GetCellSize1());
            const float xMin_v(v.GetX() - dx_v);
            const float xMax_v(v.GetX() + dx_v);
            // Check that the U and V views are compatible in x, otherwise skip the chi-squared calculation
            if ((xMax_u < xMin_v) || (xMin_u > xMax_v))
                continue;
            
            float bestChi2 = std::numeric_limits<float>::max();
            int bestK{-1};

            for (int k = 0; k < nW; ++k)
            {
                const CartesianVector w(wHits[k]->GetPositionVector());
                const float dx_w(0.5f * wHits[k]->GetCellSize1());
                const float xMin_w(w.GetX() - dx_w);
                const float xMax_w(w.GetX() + dx_w);
                // Check that W view is compatible with U and V in x, otherwise skip the chi-squared calculation
                if ((xMax_u < xMin_w) || (xMin_u > xMax_w) || (xMax_v < xMin_w) || (xMin_v > xMax_w))
                    continue;

                float chi2{LArGeometryHelper::CalculateChiSquared(this->GetPandora(), uHits[i], vHits[j], wHits[k])};
                if (chi2 < bestChi2)
                {
                    bestChi2 = chi2;
                    bestK = k;
                }
            }

            if (bestK >= 0)
                C[i][j] = bestChi2;
        }
    }

    return C;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PlaneSolverAlgorithm::CostMatrix PlaneSolverAlgorithm::ComputeTripletCostMatrix(const PairVector &uvPairs, const PlaneToHitsMap &planeToHitsMap,
    const float unmatchedCost) const
{
    const CaloHitVector &uHits(planeToHitsMap.at(TPC_VIEW_U));
    const CaloHitVector &vHits(planeToHitsMap.at(TPC_VIEW_V));
    const CaloHitVector &wHits(planeToHitsMap.at(TPC_VIEW_W));
    int nPairs{static_cast<int>(uvPairs.size())};
    int nW{static_cast<int>(wHits.size())};
    int N{std::max(nPairs, nW)};

    CostMatrix C(N, FloatVector(N, unmatchedCost));
    for (int p = 0; p < nPairs; ++p)
    {
        int i{uvPairs[p].m_uIndex};
        int j{uvPairs[p].m_vIndex};

        for (int k = 0; k < nW; ++k)
        {
            float chi2{LArGeometryHelper::CalculateChiSquared(this->GetPandora(), uHits[i], vHits[j], wHits[k])};
            C[p][k] = chi2;
        }
    }

    return C;
}

//------------------------------------------------------------------------------------------------------------------------------------------

IntVector PlaneSolverAlgorithm::KuhneMunkres(const CostMatrix& cost) const
{
    const int N{static_cast<int>(cost.size())};
    // Define the dual variables for rows (u) and columns (v), enforcing u[i] + v[j] <= cost[i][j]
    // We are essentially trying to maximise the sum of the duals over all rows/columns, while enforcing the constraint above
    FloatVector u(N + 1), v(N + 1);
    // Define the matching from a row to a column (p[col] = row) and the way matrix to reconstruct the augmenting path - indexing is
    // 1-based, zero is a dummy node
    IntVector p(N + 1), way(N + 1);

    // Iterate over the rows
    for (int i = 1; i <= N; ++i)
    {
        // Start from the dummy column zero and initialise lowest cost and element usage
        int j0 = 0;
        p[j0] = i;
        FloatVector minReducedCost(N + 1, std::numeric_limits<float>::max());
        IntVector used(N + 1, 0);

        // Build the augmenting path (Dijkstra-like)
        do
        {
            // Mark the current column (j0) and the current row (i0)
            used[j0] = 1;
            int i0{p[j0]}, j1{0};
            float delta{std::numeric_limits<float>::max()};

            // Evaluate all possible assignments - trying to find a path from the current row to the best column with the minimum reduced cost
            for (int j = 1; j <= N; ++j)
            {
                if (used[j])
                    continue;
                // Key step, computing the reduced cost to maintain feasibility
                // When an optimal edge between hits is found, the reduced cost will be zero
                float reducedCost{cost[i0 - 1][j - 1] - u[i0] - v[j]};
                if (reducedCost < minReducedCost[j])
                {
                    // Keep track of the best cost for reaching column j so far
                    minReducedCost[j] = reducedCost;
                    // ... and how we got there. We reach column j via column j0
                    way[j] = j0;
                }
                if (minReducedCost[j] < delta)
                {
                    // Keep track of the best column to expand to
                    delta = minReducedCost[j];
                    j1 = j;
                }
            }

            // Update the potentials (ensures reduced cost of zero for at least one edge)
            // Opposing deltas also ensure that u_i + v_j = (u_i + delta) + (v_j - delta), so the constraint still holds
            for (int j = 0; j <= N; ++j)
            {
                if (used[j])
                {
                    u[p[j]] += delta;
                    v[j] -= delta;
                }
                else
                {
                    // Reduce the slack for unmatched columns
                    minReducedCost[j] -= delta;
                }
            }

            // Step j0 forward - continue the search
            j0 = j1;
            // Stop when we get to an unmatched column, implying we've found an augmenting path. Otherwise, we have a clash and should look
            // to reassign the previous row (which happens because we've now set j0 to j1 for the next loop iteration)
        }
        while (p[j0] != 0);

        // Walk back through way to update the matching along the augmenting path
        do
        {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        }
        while (j0);
    }

    // p[j] = i means row i is matched to column j, so we need to invert this to get the assignment of columns to rows.
    // Unmatched columns will have p[j] = 0, so we assign -1 in the output.
    IntVector assignment(N, -1);
    for (int j = 1; j <= N; ++j)
        if (p[j] != 0)
            assignment[p[j] - 1] = j - 1;

    return assignment;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PlaneSolverAlgorithm::PairVector PlaneSolverAlgorithm::BuildUVPairs(const IntVector &assignment, int nU, int nV) const
{
    PairVector pairs;

    for (int i = 0; i < nU; ++i)
    {
        int j{assignment[i]};

        if (j < nV)
            pairs.push_back({i, j});
    }

    return pairs;
}

//------------------------------------------------------------------------------------------------------------------------------------------

PlaneSolverAlgorithm::TripletVector PlaneSolverAlgorithm::BuildTriplets(const PairVector& uvPairs, const IntVector& assignment, int nW) const
{
    TripletVector result;
    const int nPairs{static_cast<int>(uvPairs.size())};

    for (int p = 0; p < nPairs; ++p)
    {
        int k{assignment[p]};

        if (k < nW)
            result.push_back({uvPairs[p].m_uIndex, uvPairs[p].m_vIndex, k});
        else
            result.push_back({uvPairs[p].m_uIndex, uvPairs[p].m_vIndex, -1});
    }

    return result;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PlaneSolverAlgorithm::CreateThreeDHit(const CaloHit *const pCaloHitU, const CaloHit *const pCaloHitV, const CaloHit *&pCaloHit3D) const
{
    float chi2{0.f};
    CartesianVector pos3D(0, 0, 0);
    PandoraContentApi::CaloHit::Parameters parameters;
    LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, pCaloHitU->GetPositionVector(), pCaloHitV->GetPositionVector(),
        pos3D, chi2);
    parameters.m_positionVector = pos3D;
    parameters.m_hitType = TPC_3D;

    const CaloHit *const pCaloHit2D(pCaloHitU);
    parameters.m_pParentAddress = static_cast<const void *>(pCaloHit2D);
    parameters.m_cellThickness = pCaloHit2D->GetCellThickness();
    parameters.m_cellGeometry = RECTANGULAR;
    parameters.m_cellSize0 = pCaloHit2D->GetCellLengthScale();
    parameters.m_cellSize1 = pCaloHit2D->GetCellLengthScale();
    parameters.m_cellNormalVector = pCaloHit2D->GetCellNormalVector();
    parameters.m_expectedDirection = pCaloHit2D->GetExpectedDirection();
    parameters.m_nCellRadiationLengths = pCaloHit2D->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pCaloHit2D->GetNCellInteractionLengths();
    parameters.m_time = pCaloHit2D->GetTime();
    parameters.m_inputEnergy = pCaloHit2D->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pCaloHit2D->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pCaloHit2D->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pCaloHit2D->GetHadronicEnergy();
    parameters.m_isDigital = pCaloHit2D->IsDigital();
    parameters.m_hitRegion = pCaloHit2D->GetHitRegion();
    parameters.m_layer = pCaloHit2D->GetLayer();
    parameters.m_isInOuterSamplingLayer = pCaloHit2D->IsInOuterSamplingLayer();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pCaloHit3D));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PlaneSolverAlgorithm::CreateThreeDHit(const CaloHit *const pCaloHitU, const CaloHit *const pCaloHitV, const CaloHit *const pCaloHitW,
    const CaloHit *&pCaloHit3D) const
{
    float chi2{0.f};
    CartesianVector pos3D(0, 0, 0);
    PandoraContentApi::CaloHit::Parameters parameters;
    LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, pCaloHitU->GetPositionVector(),
        pCaloHitV->GetPositionVector(), pCaloHitW->GetPositionVector(), pos3D, chi2);
    parameters.m_positionVector = pos3D;
    parameters.m_hitType = TPC_3D;

    const CaloHit *const pCaloHit2D(pCaloHitU);
    parameters.m_pParentAddress = static_cast<const void *>(pCaloHit2D);
    parameters.m_cellThickness = pCaloHit2D->GetCellThickness();
    parameters.m_cellGeometry = RECTANGULAR;
    parameters.m_cellSize0 = pCaloHit2D->GetCellLengthScale();
    parameters.m_cellSize1 = pCaloHit2D->GetCellLengthScale();
    parameters.m_cellNormalVector = pCaloHit2D->GetCellNormalVector();
    parameters.m_expectedDirection = pCaloHit2D->GetExpectedDirection();
    parameters.m_nCellRadiationLengths = pCaloHit2D->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = pCaloHit2D->GetNCellInteractionLengths();
    parameters.m_time = pCaloHit2D->GetTime();
    parameters.m_inputEnergy = pCaloHit2D->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = pCaloHit2D->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = pCaloHit2D->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = pCaloHit2D->GetHadronicEnergy();
    parameters.m_isDigital = pCaloHit2D->IsDigital();
    parameters.m_hitRegion = pCaloHit2D->GetHitRegion();
    parameters.m_layer = pCaloHit2D->GetLayer();
    parameters.m_isInOuterSamplingLayer = pCaloHit2D->IsInOuterSamplingLayer();
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CaloHit::Create(*this, parameters, pCaloHit3D));
}


//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PlaneSolverAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
