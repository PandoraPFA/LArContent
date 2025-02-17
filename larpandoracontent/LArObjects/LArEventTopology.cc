/**
 *  @file   larpandoracontent/LArObjects/LArEventTopology.cc
 *
 *  @brief  Implementation of lar pfo objects.
 *
 *  $Log: $
 */

#include "Helpers/MCParticleHelper.h"
#include "Objects/CaloHit.h"
#include "Pandora/PdgTable.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArEventTopology.h"

#include <sstream>

using namespace pandora;

namespace lar_content
{

LArEventTopology::LArEventTopology(const CaloHitList &caloHitList2D) :
    m_pRoot{nullptr},
    m_pParticle{nullptr}
{
    for (const CaloHit *const pCaloHit : caloHitList2D)
    {
        try
        {
            const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            m_mcHitMap[pMCParticle].emplace_back(pCaloHit);
            if (!(m_pRoot && LArMCParticleHelper::IsNeutrino(m_pRoot)))
            {
                const MCParticle *pRoot{pMCParticle};
                while (!pRoot->GetParentList().empty())
                    pRoot = pRoot->GetParentList().front();
                m_pRoot = pRoot;
            }
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArEventTopology::~LArEventTopology()
{
    delete m_pParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArEventTopology::ConstructVisibleHierarchy()
{
    m_pParticle = new Particle(m_pRoot);
    this->ConstructVisibleHierarchy(m_pParticle, m_pRoot);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArEventTopology::GetVertices(CartesianPointVector &vertices) const
{
    m_pParticle->GetVertices(vertices);

    CartesianPointVector forRemoval;
    for (auto iter1 = vertices.begin(); iter1 != vertices.end(); ++iter1)
    {
        const CartesianVector &vertex1{*iter1};
        if (std::find(forRemoval.begin(), forRemoval.end(), vertex1) != forRemoval.end())
            continue;
        for (auto iter2 = std::next(iter1); iter2 != vertices.end(); ++iter2)
        {
            const CartesianVector &vertex2{*iter2};
            if (std::find(forRemoval.begin(), forRemoval.end(), vertex2) != forRemoval.end())
                continue;
            if (vertex1.GetDistanceSquared(vertex2) < 1)
                forRemoval.emplace_back(vertex2);
        }
    }

    for (auto iter = vertices.begin(); iter != vertices.end();)
    {
        if (std::find(forRemoval.begin(), forRemoval.end(), *iter) != forRemoval.end())
            iter = vertices.erase(iter);
        else
            ++iter;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArEventTopology::Print() const
{
    std::string str{"Neutrino hierarchy: " + std::to_string(m_pRoot->GetParticleId()) + " E: " + std::to_string(m_pRoot->GetEnergy()) +
        " GeV\n" + m_pParticle->Print(m_mcHitMap, "  ")};
    std::cout << str << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArEventTopology::ConstructVisibleHierarchy(Particle *pParticle, const MCParticle *const pRootMC)
{
    for (const MCParticle *const pMC : pRootMC->GetDaughterList())
    {
        if (m_mcHitMap.find(pMC) != m_mcHitMap.end())
        {
            Particle *pChild{new Particle(pMC, m_mcHitMap.at(pMC))};
            pParticle->AddChild(pChild);
            this->ConstructVisibleHierarchy(pChild, pMC);
        }
        else
        {
            this->ConstructVisibleHierarchy(pParticle, pMC);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArEventTopology::PruneHierarchy()
{
    m_pParticle->Parse(m_mcHitMap);
    m_pParticle->Prune();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArEventTopology::Particle::Particle(const MCParticle *const pRoot) :
    m_fold{false}
{
    m_particles.emplace_back(pRoot);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArEventTopology::Particle::Particle(const MCParticle *const pRoot, const CaloHitList &caloHitList) :
    m_caloHits{caloHitList},
    m_fold{false}
{
    m_particles.emplace_back(pRoot);
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArEventTopology::Particle::~Particle()
{
    for (Particle *pParticle : m_children)
        delete pParticle;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArEventTopology::Particle::AddChild(Particle *pChild)
{
    m_children.emplace_back(pChild);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArEventTopology::Particle::GetVertices(CartesianPointVector &vertices) const
{
    const MCParticle *const pMC{m_particles.front()};
    if (LArMCParticleHelper::IsNeutrino(pMC))
    {
        if (std::find(vertices.begin(), vertices.end(), pMC->GetVertex()) == vertices.end())
            vertices.emplace_back(pMC->GetVertex());
    }
    else
    {
        // Consider vertex
        const MCParticle *const pParentMC{pMC->GetParentList().front()};
        if (LArMCParticleHelper::IsNeutrino(pParentMC))
        {
            if (std::find(vertices.begin(), vertices.end(), pMC->GetVertex()) == vertices.end())
                vertices.emplace_back(pMC->GetVertex());
        }
        else
        {
            const CartesianVector &momParent{pParentMC->GetMomentum().GetUnitVector()};
            const CartesianVector &mom{pMC->GetMomentum().GetUnitVector()};
            if ((pMC->GetParticleId() != pParentMC->GetParticleId()) || (mom.GetDotProduct(momParent) < 0.996f))
            {
                // Visible photon vertices occur at the endpoint, not the vertex, so skip photons here
                if (pMC->GetParticleId() != PHOTON && std::find(vertices.begin(), vertices.end(), pMC->GetVertex()) == vertices.end())
                    vertices.emplace_back(pMC->GetVertex());
            }
        }
        // Consider end point
        int nTrackLikeChildren{0};
        for (const Particle *const pChild : m_children)
        {
            const MCParticle *const pChildMC{pChild->m_particles.front()};
            const int pdg{std::abs(pChildMC->GetParticleId())};
            if (!(pdg == PHOTON || pdg == E_MINUS))
                ++nTrackLikeChildren;
        }
        // If only shower children, or if more than 1 track, tag the endpoint as a vertex
        // If  we have exactly 1 track child, then downstream vertex checks will take care of this
        // Also note that the endpoint of a photon actually marks the start of the associated hit deposition
        if (nTrackLikeChildren != 1)
        {
            if (std::find(vertices.begin(), vertices.end(), pMC->GetEndpoint()) == vertices.end())
                vertices.emplace_back(pMC->GetEndpoint());
        }
    }

    for (const Particle *const pChild : m_children)
        pChild->GetVertices(vertices);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArEventTopology::Particle::Parse(const MCHitMap &mcHitMap)
{
    const MCParticle *const pMC{m_particles.front()};
    const CaloHitList &caloHitList{mcHitMap.find(pMC) != mcHitMap.end() ? mcHitMap.at(pMC) : CaloHitList()};
    int uHits{0}, vHits{0}, wHits{0};
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        switch (pCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                ++uHits;
                break;
            case TPC_VIEW_V:
                ++vHits;
                break;
            case TPC_VIEW_W:
                ++wHits;
                break;
            default:
                break;
        }
    }

    const int pdg{std::abs(pMC->GetParticleId())};
    const int target{pdg == PHOTON || pdg == E_MINUS ? 8 : 5};
    const int totalHits{uHits + vHits + wHits};
    int goodViews{uHits >= target ? 1 : 0};
    goodViews += vHits >= target ? 1 : 0;
    goodViews += wHits >= target ? 1 : 0;

    if (totalHits < (3 * target) || goodViews < 2)
        m_fold = true;

    for (Particle *const pChild : m_children)
        pChild->Parse(mcHitMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArEventTopology::Particle::Prune()
{
    for (Particle *const pChild : m_children)
        pChild->Prune();

    std::list<Particle> toPrune;
    for (auto iter = m_children.begin(); iter != m_children.end();)
    {
        if ((*iter)->m_fold && (*iter)->m_children.empty())
            iter = m_children.erase(iter);
        else
            ++iter;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string LArEventTopology::Particle::Print(const MCHitMap &mcHitMap, const std::string &indent) const
{
    const MCParticle *const pMC{m_particles.front()};
    const CaloHitList &caloHitList{mcHitMap.find(pMC) != mcHitMap.end() ? mcHitMap.at(pMC) : CaloHitList()};
    int uHits{0}, vHits{0}, wHits{0};
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        switch (pCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                ++uHits;
                break;
            case TPC_VIEW_V:
                ++vHits;
                break;
            case TPC_VIEW_W:
                ++wHits;
                break;
            default:
                break;
        }
    }

    std::ostringstream os;
    os << indent << pMC->GetParticleId() << ": E: " << pMC->GetEnergy() << " GeV";
    os << " (" << uHits << "," << vHits << "," << wHits << ") " << (m_fold ? "[fold]" : "") << "\n";
    for (const Particle *const pChild : m_children)
        os << pChild->Print(mcHitMap, indent + "  ");

    return os.str();
}

} // namespace lar_content
