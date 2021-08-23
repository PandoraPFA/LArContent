/**
 *  @file   larpandoracontent/LArHelpers/LArHierarchyHelper.cc
 *
 *  @brief  Implementation of the lar hierarchy helper class.
 *
 *  $Log: $
 */

#include "Pandora/PdgTable.h"
#include "Pandora/StatusCodes.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"

#include <numeric>

namespace lar_content
{

using namespace pandora;

LArHierarchyHelper::MCHierarchy::MCHierarchy(const ReconstructabilityCriteria &recoCriteria) :
    m_recoCriteria(recoCriteria),
    m_pNeutrino{nullptr}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::~MCHierarchy()
{
    for (const Node *pNode : m_rootNodes)
        delete pNode;
    m_rootNodes.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCHierarchy::FillHierarchy(
    const MCParticleList &mcParticleList, const CaloHitList &caloHitList, const bool foldToPrimaries, const bool foldToLeadingShowers)
{
    const auto predicate = [](const MCParticle *pMCParticle) { return std::abs(pMCParticle->GetParticleId()) == NEUTRON; };
    m_mcToHitsMap.clear();
    for (const CaloHit *pCaloHit : caloHitList)
    {
        try
        {
            const MCParticle *pMCParticle{MCParticleHelper::GetMainMCParticle(pCaloHit)};
            m_mcToHitsMap[pMCParticle].emplace_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }

    MCParticleSet primarySet;
    m_pNeutrino = LArHierarchyHelper::GetMCPrimaries(mcParticleList, primarySet);
    MCParticleList primaries(primarySet.begin(), primarySet.end());
    primaries.sort(LArMCParticleHelper::SortByMomentum);
    if (m_recoCriteria.m_removeNeutrons)
        primaries.erase(std::remove_if(primaries.begin(), primaries.end(), predicate), primaries.end());
    if (foldToPrimaries && !foldToLeadingShowers)
    {
        for (const MCParticle *pPrimary : primaries)
        {
            MCParticleList allParticles{pPrimary};
            if (!m_recoCriteria.m_removeNeutrons)
            {
                LArMCParticleHelper::GetAllDescendentMCParticles(pPrimary, allParticles);
            }
            else
            {
                // Collect track-like and shower-like particles together, but throw out neutrons and descendents
                MCParticleList dummy;
                LArMCParticleHelper::GetAllDescendentMCParticles(pPrimary, allParticles, allParticles, dummy);
            }
            CaloHitList allHits;
            for (const MCParticle *pMCParticle : allParticles)
            {
                // Not all MC particles will have hits
                if (m_mcToHitsMap.find(pMCParticle) != m_mcToHitsMap.end())
                {
                    const CaloHitList &caloHits(m_mcToHitsMap.at(pMCParticle));
                    allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
                }
            }
            m_rootNodes.emplace_back(new Node(*this, allParticles, allHits));
        }
    }
    else if (foldToPrimaries && foldToLeadingShowers)
    {
        for (const MCParticle *pPrimary : primaries)
        {
            MCParticleList allParticles{pPrimary}, showerParticles, neutrons;
            int pdg{std::abs(pPrimary->GetParticleId())};
            const bool isShower{pdg == E_MINUS || pdg == PHOTON};
            const bool isNeutron{pdg == NEUTRON};
            if (isShower || isNeutron)
            {
                if (!m_recoCriteria.m_removeNeutrons)
                {
                    LArMCParticleHelper::GetAllDescendentMCParticles(pPrimary, allParticles);
                }
                else
                {
                    // Throw away neutrons
                    MCParticleList dummy;
                    LArMCParticleHelper::GetAllDescendentMCParticles(pPrimary, allParticles, allParticles, dummy);
                }
            }
            else
            {
                LArMCParticleHelper::GetAllDescendentMCParticles(pPrimary, allParticles, showerParticles, neutrons);
            }
            CaloHitList allHits;
            for (const MCParticle *pMCParticle : allParticles)
            {
                // ATTN - Not all MC particles will have hits
                if (m_mcToHitsMap.find(pMCParticle) != m_mcToHitsMap.end())
                {
                    const CaloHitList &caloHits(m_mcToHitsMap.at(pMCParticle));
                    allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
                }
            }
            Node *pNode{new Node(*this, allParticles, allHits)};
            m_rootNodes.emplace_back(pNode);
            if (!showerParticles.empty())
            {
                // Collect up all descendent hits for each shower and add the nodes as a child of the root node
                for (const MCParticle *pChild : showerParticles)
                    pNode->FillFlat(pChild);
            }
            if (!m_recoCriteria.m_removeNeutrons && !neutrons.empty())
            {
                // Collect up all descendent hits for each neutron and add the nodes as a child of the root node
                for (const MCParticle *pChild : neutrons)
                    pNode->FillFlat(pChild);
            }
        }
    }
    else if (foldToLeadingShowers)
    {
        for (const MCParticle *pPrimary : primaries)
        {
            MCParticleList allParticles{pPrimary};
            int pdg{std::abs(pPrimary->GetParticleId())};
            const bool isShower{pdg == E_MINUS || pdg == PHOTON};
            const bool isNeutron{pdg == NEUTRON};
            if (isShower || (isNeutron && !m_recoCriteria.m_removeNeutrons))
                LArMCParticleHelper::GetAllDescendentMCParticles(pPrimary, allParticles);
            CaloHitList allHits;
            for (const MCParticle *pMCParticle : allParticles)
            {
                // ATTN - Not all MC particles will have hits
                if (m_mcToHitsMap.find(pMCParticle) != m_mcToHitsMap.end())
                {
                    const CaloHitList &caloHits(m_mcToHitsMap.at(pMCParticle));
                    allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
                }
            }
            Node *pNode{new Node(*this, allParticles, allHits)};
            m_rootNodes.emplace_back(pNode);
            if (!(isShower || isNeutron))
            {
                // Find the children of this particle and recursively add them to the hierarchy
                const MCParticleList &children{pPrimary->GetDaughterList()};
                for (const MCParticle *pChild : children)
                    pNode->FillHierarchy(pChild, foldToLeadingShowers);
            }
        }
    }
    else
    {
        for (const MCParticle *pPrimary : primaries)
        {
            MCParticleList allParticles{pPrimary};
            CaloHitList allHits;
            for (const MCParticle *pMCParticle : allParticles)
            {
                // ATTN - Not all MC particles will have hits
                if (m_mcToHitsMap.find(pMCParticle) != m_mcToHitsMap.end())
                {
                    const CaloHitList &caloHits(m_mcToHitsMap.at(pMCParticle));
                    allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
                }
            }
            Node *pNode{new Node(*this, allParticles, allHits)};
            m_rootNodes.emplace_back(pNode);
            // Find the children of this particle and recursively add them to the hierarchy
            const MCParticleList &children{pPrimary->GetDaughterList()};
            for (const MCParticle *pChild : children)
                pNode->FillHierarchy(pChild, foldToLeadingShowers);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCHierarchy::GetFlattenedNodes(NodeVector &nodeVector) const
{
    NodeList queue;
    for (const Node *pNode : m_rootNodes)
    {
        nodeVector.emplace_back(pNode);
        queue.emplace_back(pNode);
    }
    while (!queue.empty())
    {
        const NodeVector &children{queue.front()->GetChildren()};
        queue.pop_front();
        for (const Node *pChild : children)
        {
            nodeVector.emplace_back(pChild);
            queue.emplace_back(pChild);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string LArHierarchyHelper::MCHierarchy::ToString() const
{
    std::string str;
    for (const Node *pNode : m_rootNodes)
        str += pNode->ToString("") + "\n";

    return str;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::Node::Node(const MCHierarchy &hierarchy, const MCParticle *pMCParticle) :
    m_hierarchy(hierarchy),
    m_mainParticle(pMCParticle),
    m_pdg{0}
{
    if (pMCParticle)
    {
        m_pdg = pMCParticle->GetParticleId();
        m_mcParticles.emplace_back(pMCParticle);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::Node::Node(const MCHierarchy &hierarchy, const MCParticleList &mcParticleList, const CaloHitList &caloHitList) :
    m_hierarchy(hierarchy),
    m_mcParticles(mcParticleList),
    m_caloHits(caloHitList),
    m_mainParticle(nullptr),
    m_pdg{0}
{
    if (!mcParticleList.empty())
    {
        m_mainParticle = mcParticleList.front();
        m_pdg = m_mainParticle->GetParticleId();
    }
    m_mcParticles.sort(LArMCParticleHelper::SortByMomentum);
    m_caloHits.sort();
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::Node::~Node()
{
    m_mcParticles.clear();
    m_caloHits.clear();
    for (const Node *node : m_children)
        delete node;
    m_children.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCHierarchy::Node::FillHierarchy(const MCParticle *pRoot, const bool foldToLeadingShowers)
{
    MCParticleList allParticles{pRoot};
    const int pdg{std::abs(pRoot->GetParticleId())};
    const bool isShower{pdg == E_MINUS || pdg == PHOTON};
    const bool isNeutron{pdg == NEUTRON};
    if (foldToLeadingShowers && (isShower || (isNeutron && !m_hierarchy.m_recoCriteria.m_removeNeutrons)))
        LArMCParticleHelper::GetAllDescendentMCParticles(pRoot, allParticles);
    else if (m_hierarchy.m_recoCriteria.m_removeNeutrons && isNeutron)
        return;

    CaloHitList allHits;
    for (const MCParticle *pMCParticle : allParticles)
    {
        // ATTN - Not all MC particles will have hits
        if (m_hierarchy.m_mcToHitsMap.find(pMCParticle) != m_hierarchy.m_mcToHitsMap.end())
        {
            const CaloHitList &caloHits(m_hierarchy.m_mcToHitsMap.at(pMCParticle));
            allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
        }
    }

    if (!allParticles.empty())
    {
        Node *pNode{new Node(m_hierarchy, allParticles, allHits)};
        m_children.emplace_back(pNode);
        if (!foldToLeadingShowers || (foldToLeadingShowers && !(isShower || isNeutron)))
        {
            // Find the children of this particle and recursively add them to the hierarchy
            const MCParticleList &children{pRoot->GetDaughterList()};
            for (const MCParticle *pChild : children)
                pNode->FillHierarchy(pChild, foldToLeadingShowers);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCHierarchy::Node::FillFlat(const MCParticle *pRoot)
{
    MCParticleList allParticles{pRoot};
    if (!m_hierarchy.m_recoCriteria.m_removeNeutrons)
    {
        LArMCParticleHelper::GetAllDescendentMCParticles(pRoot, allParticles);
    }
    else
    {
        MCParticleList neutrons;
        LArMCParticleHelper::GetAllDescendentMCParticles(pRoot, allParticles, allParticles, neutrons);
    }
    CaloHitList allHits;
    for (const MCParticle *pMCParticle : allParticles)
    {
        // ATTN - Not all MC particles will have hits
        if (m_hierarchy.m_mcToHitsMap.find(pMCParticle) != m_hierarchy.m_mcToHitsMap.end())
        {
            const CaloHitList &caloHits(m_hierarchy.m_mcToHitsMap.at(pMCParticle));
            allHits.insert(allHits.begin(), caloHits.begin(), caloHits.end());
        }
    }
    if (!allParticles.empty())
    {
        Node *pNode{new Node(m_hierarchy, allParticles, allHits)};
        m_children.emplace_back(pNode);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHierarchyHelper::MCHierarchy::Node::IsReconstructable() const
{
    const bool enoughHits{m_caloHits.size() >= m_hierarchy.m_recoCriteria.m_minHits};
    if (!enoughHits)
        return false;
    bool enoughGoodViews{false};
    unsigned int nHitsU{0}, nHitsV{0}, nHitsW{0};
    for (const CaloHit *pCaloHit : m_caloHits)
    {
        switch (pCaloHit->GetHitType())
        {
            case TPC_VIEW_U:
                ++nHitsU;
                break;
            case TPC_VIEW_V:
                ++nHitsV;
                break;
            case TPC_VIEW_W:
                ++nHitsW;
                break;
            default:
                break;
        }
        unsigned int nGoodViews{0};
        if (nHitsU >= m_hierarchy.m_recoCriteria.m_minHitsForGoodView)
            ++nGoodViews;
        if (nHitsV >= m_hierarchy.m_recoCriteria.m_minHitsForGoodView)
            ++nGoodViews;
        if (nHitsW >= m_hierarchy.m_recoCriteria.m_minHitsForGoodView)
            ++nGoodViews;
        if (nGoodViews >= m_hierarchy.m_recoCriteria.m_minGoodViews)
        {
            enoughGoodViews = true;
            break;
        }
    }

    return enoughGoodViews;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHierarchyHelper::MCHierarchy::Node::IsTestBeamParticle() const
{
    if (m_mainParticle)
        return LArMCParticleHelper::IsBeamParticle(m_mainParticle);
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArHierarchyHelper::MCHierarchy::Node::IsCosmicRay() const
{
    if (m_mainParticle)
        return LArMCParticleHelper::IsCosmicRay(m_mainParticle);
    else
        return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string LArHierarchyHelper::MCHierarchy::Node::ToString(const std::string &prefix) const
{
    std::string str(prefix + "PDG: " + std::to_string(m_pdg) + " Energy: " + std::to_string(m_mainParticle ? m_mainParticle->GetEnergy() : 0) +
                    " Hits: " + std::to_string(m_caloHits.size()) + "\n");
    for (const Node *pChild : m_children)
        str += pChild->ToString(prefix + "   ");

    return str;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria::ReconstructabilityCriteria() :
    m_minHits{15},
    m_minHitsForGoodView{5},
    m_minGoodViews{2},
    m_removeNeutrons{true}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria::ReconstructabilityCriteria(const ReconstructabilityCriteria &obj) :
    m_minHits{obj.m_minHits},
    m_minHitsForGoodView{obj.m_minHitsForGoodView},
    m_minGoodViews{obj.m_minGoodViews},
    m_removeNeutrons{obj.m_removeNeutrons}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria::ReconstructabilityCriteria(
    const unsigned int minHits, const unsigned int minHitsForGoodView, const unsigned int minGoodViews, const bool removeNeutrons) :
    m_minHits{minHits},
    m_minHitsForGoodView{minHitsForGoodView},
    m_minGoodViews{minGoodViews},
    m_removeNeutrons{removeNeutrons}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::RecoHierarchy::RecoHierarchy() : m_pNeutrino{nullptr}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::RecoHierarchy::~RecoHierarchy()
{
    for (const Node *pNode : m_rootNodes)
        delete pNode;
    m_rootNodes.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::RecoHierarchy::FillHierarchy(const PfoList &pfoList, const bool foldToPrimaries, const bool foldToLeadingShowers)
{
    PfoSet primarySet;
    m_pNeutrino = LArHierarchyHelper::GetRecoPrimaries(pfoList, primarySet);
    PfoList primaries(primarySet.begin(), primarySet.end());
    primaries.sort(LArPfoHelper::SortByNHits);
    if (foldToPrimaries && !foldToLeadingShowers)
    {
        for (const ParticleFlowObject *pPrimary : primaries)
        {
            PfoList allParticles;
            // ATTN - pPrimary gets added to the list of downstream PFOs, not just the child PFOs
            LArPfoHelper::GetAllDownstreamPfos(pPrimary, allParticles);
            CaloHitList allHits;
            for (const ParticleFlowObject *pPfo : allParticles)
                LArPfoHelper::GetAllCaloHits(pPfo, allHits);
            m_rootNodes.emplace_back(new Node(*this, allParticles, allHits));
        }
    }
    else if (foldToPrimaries && foldToLeadingShowers)
    {
        for (const ParticleFlowObject *pPrimary : primaries)
        {
            PfoList allParticles, showerParticles;
            int pdg{std::abs(pPrimary->GetParticleId())};
            const bool isShower{pdg == E_MINUS};
            // ATTN - pPrimary gets added to the list of downstream PFOs, not just the child PFOs
            if (isShower)
                LArPfoHelper::GetAllDownstreamPfos(pPrimary, allParticles);
            else
                LArPfoHelper::GetAllDownstreamPfos(pPrimary, allParticles, showerParticles);
            CaloHitList allHits;
            for (const ParticleFlowObject *pPfo : allParticles)
                LArPfoHelper::GetAllCaloHits(pPfo, allHits);
            Node *pNode{new Node(*this, allParticles, allHits)};
            m_rootNodes.emplace_back(pNode);
            if (!showerParticles.empty())
            {
                // Collect up all descendent hits for each shower and add the nodes as a child of the root node
                for (const ParticleFlowObject *pChild : showerParticles)
                    pNode->FillFlat(pChild);
            }
        }
    }
    else if (foldToLeadingShowers)
    {
        for (const ParticleFlowObject *pPrimary : primaries)
        {
            PfoList allParticles;
            int pdg{std::abs(pPrimary->GetParticleId())};
            const bool isShower{pdg == E_MINUS};
            // ATTN - pPrimary gets added to the list of downstream PFOs, not just the child PFOs
            if (isShower)
                LArPfoHelper::GetAllDownstreamPfos(pPrimary, allParticles);
            else
                allParticles.emplace_back(pPrimary);

            CaloHitList allHits;
            for (const ParticleFlowObject *pPfo : allParticles)
                LArPfoHelper::GetAllCaloHits(pPfo, allHits);
            Node *pNode{new Node(*this, allParticles, allHits)};
            m_rootNodes.emplace_back(pNode);
            if (!isShower)
            {
                // Find the children of this particle and recursively add them to the hierarchy
                const PfoList &children{pPrimary->GetDaughterPfoList()};
                for (const ParticleFlowObject *pChild : children)
                    pNode->FillHierarchy(pChild, foldToLeadingShowers);
            }
        }
    }
    else
    {
        for (const ParticleFlowObject *pPrimary : primaries)
        {
            PfoList allParticles{pPrimary};
            CaloHitList allHits;
            for (const ParticleFlowObject *pPfo : allParticles)
                LArPfoHelper::GetAllCaloHits(pPfo, allHits);
            Node *pNode{new Node(*this, allParticles, allHits)};
            m_rootNodes.emplace_back(pNode);
            // Find the children of this particle and recursively add them to the hierarchy
            const PfoList &children{pPrimary->GetDaughterPfoList()};
            for (const ParticleFlowObject *pChild : children)
                pNode->FillHierarchy(pChild, foldToLeadingShowers);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::RecoHierarchy::GetFlattenedNodes(NodeVector &nodeVector) const
{
    NodeList queue;
    for (const Node *pNode : m_rootNodes)
    {
        nodeVector.emplace_back(pNode);
        queue.emplace_back(pNode);
    }
    while (!queue.empty())
    {
        const NodeVector &children{queue.front()->GetChildren()};
        queue.pop_front();
        for (const Node *pChild : children)
        {
            nodeVector.emplace_back(pChild);
            queue.emplace_back(pChild);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string LArHierarchyHelper::RecoHierarchy::ToString() const
{
    std::string str;
    for (const Node *pNode : m_rootNodes)
        str += pNode->ToString("") + "\n";

    return str;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::RecoHierarchy::Node::Node(const RecoHierarchy &hierarchy, const ParticleFlowObject *pPfo) :
    m_hierarchy(hierarchy),
    m_pdg{0}
{
    if (pPfo)
    {
        m_pdg = pPfo->GetParticleId();
        m_pfos.emplace_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::RecoHierarchy::Node::Node(const RecoHierarchy &hierarchy, const PfoList &pfoList, const CaloHitList &caloHitList) :
    m_hierarchy(hierarchy),
    m_pdg{0}
{
    if (!pfoList.empty())
        m_pdg = pfoList.front()->GetParticleId();
    m_pfos = pfoList;
    m_pfos.sort(LArPfoHelper::SortByNHits);
    m_caloHits = caloHitList;
    m_caloHits.sort();
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::RecoHierarchy::Node::~Node()
{
    m_pfos.clear();
    m_caloHits.clear();
    for (const Node *node : m_children)
        delete node;
    m_children.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::RecoHierarchy::Node::FillHierarchy(const ParticleFlowObject *pRoot, const bool foldToLeadingShowers)
{
    PfoList allParticles;
    int pdg{std::abs(pRoot->GetParticleId())};
    const bool isShower{pdg == E_MINUS};
    if (foldToLeadingShowers && isShower)
        LArPfoHelper::GetAllDownstreamPfos(pRoot, allParticles);
    else
        allParticles.emplace_back(pRoot);

    CaloHitList allHits;
    for (const ParticleFlowObject *pPfo : allParticles)
        LArPfoHelper::GetAllCaloHits(pPfo, allHits);
    Node *pNode{new Node(m_hierarchy, allParticles, allHits)};
    m_children.emplace_back(pNode);
    if (!foldToLeadingShowers || (foldToLeadingShowers && !isShower))
    {
        // Find the children of this particle and recursively add them to the hierarchy
        const PfoList &children{pRoot->GetDaughterPfoList()};
        for (const ParticleFlowObject *pChild : children)
            pNode->FillHierarchy(pChild, foldToLeadingShowers);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::RecoHierarchy::Node::FillFlat(const ParticleFlowObject *pRoot)
{
    PfoList allParticles;
    LArPfoHelper::GetAllDownstreamPfos(pRoot, allParticles);
    CaloHitList allHits;
    for (const ParticleFlowObject *pPfo : allParticles)
        LArPfoHelper::GetAllCaloHits(pPfo, allHits);
    Node *pNode{new Node(m_hierarchy, allParticles, allHits)};
    m_children.emplace_back(pNode);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const PfoList &LArHierarchyHelper::RecoHierarchy::Node::GetRecoParticles() const
{
    return m_pfos;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const CaloHitList &LArHierarchyHelper::RecoHierarchy::Node::GetCaloHits() const
{
    return m_caloHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArHierarchyHelper::RecoHierarchy::Node::GetParticleId() const
{
    return m_pdg;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const std::string LArHierarchyHelper::RecoHierarchy::Node::ToString(const std::string &prefix) const
{
    std::string str(prefix + "PDG: " + std::to_string(m_pdg) + " Hits: " + std::to_string(m_caloHits.size()) + "\n");
    for (const Node *pChild : m_children)
        str += pChild->ToString(prefix + "   ");

    return str;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MCMatches::MCMatches(const MCHierarchy::Node *pMCParticle) : m_pMCParticle{pMCParticle}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MCMatches::AddRecoMatch(const RecoHierarchy::Node *pReco, const int nSharedHits)
{
    m_recoNodes.emplace_back(pReco);
    m_sharedHits.emplace_back(nSharedHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArHierarchyHelper::MCMatches::GetSharedHits(const RecoHierarchy::Node *pReco) const
{
    auto iter{std::find(m_recoNodes.begin(), m_recoNodes.end(), pReco)};
    if (iter == m_recoNodes.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    int index = iter - m_recoNodes.begin();

    return static_cast<int>(m_sharedHits[index]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHierarchyHelper::MCMatches::GetPurity(const RecoHierarchy::Node *pReco, const bool adcWeighted) const
{
    auto iter{std::find(m_recoNodes.begin(), m_recoNodes.end(), pReco)};
    if (iter == m_recoNodes.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const CaloHitList &recoHits{pReco->GetCaloHits()};
    const CaloHitList &mcHits{m_pMCParticle->GetCaloHits()};
    CaloHitVector intersection;
    std::set_intersection(mcHits.begin(), mcHits.end(), recoHits.begin(), recoHits.end(), std::back_inserter(intersection));

    return this->GetPurity(intersection, recoHits, adcWeighted);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHierarchyHelper::MCMatches::GetPurity(const RecoHierarchy::Node *pReco, const HitType view, const bool adcWeighted) const
{
    (void)view;
    auto iter{std::find(m_recoNodes.begin(), m_recoNodes.end(), pReco)};
    if (iter == m_recoNodes.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CaloHitList recoHits;
    for (const CaloHit *pCaloHit : pReco->GetCaloHits())
        if (pCaloHit->GetHitType() == view)
            recoHits.emplace_back(pCaloHit);
    CaloHitList mcHits;
    for (const CaloHit *pCaloHit : m_pMCParticle->GetCaloHits())
        if (pCaloHit->GetHitType() == view)
            mcHits.emplace_back(pCaloHit);

    CaloHitVector intersection;
    std::set_intersection(mcHits.begin(), mcHits.end(), recoHits.begin(), recoHits.end(), std::back_inserter(intersection));

    return this->GetPurity(intersection, recoHits, adcWeighted);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHierarchyHelper::MCMatches::GetCompleteness(const RecoHierarchy::Node *pReco, const bool adcWeighted) const
{
    auto iter{std::find(m_recoNodes.begin(), m_recoNodes.end(), pReco)};
    if (iter == m_recoNodes.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const CaloHitList &recoHits{pReco->GetCaloHits()};
    const CaloHitList &mcHits{m_pMCParticle->GetCaloHits()};
    CaloHitVector intersection;
    std::set_intersection(mcHits.begin(), mcHits.end(), recoHits.begin(), recoHits.end(), std::back_inserter(intersection));

    return this->GetCompleteness(intersection, mcHits, adcWeighted);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHierarchyHelper::MCMatches::GetCompleteness(const RecoHierarchy::Node *pReco, const HitType view, const bool adcWeighted) const
{
    (void)view;
    auto iter{std::find(m_recoNodes.begin(), m_recoNodes.end(), pReco)};
    if (iter == m_recoNodes.end())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CaloHitList recoHits;
    for (const CaloHit *pCaloHit : pReco->GetCaloHits())
        if (pCaloHit->GetHitType() == view)
            recoHits.emplace_back(pCaloHit);
    CaloHitList mcHits;
    for (const CaloHit *pCaloHit : m_pMCParticle->GetCaloHits())
        if (pCaloHit->GetHitType() == view)
            mcHits.emplace_back(pCaloHit);

    CaloHitVector intersection;
    std::set_intersection(mcHits.begin(), mcHits.end(), recoHits.begin(), recoHits.end(), std::back_inserter(intersection));

    return this->GetCompleteness(intersection, mcHits, adcWeighted);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHierarchyHelper::MCMatches::GetPurity(const CaloHitVector &intersection, const CaloHitList &recoHits, const bool adcWeighted) const
{
    float purity{0.f};
    if (!intersection.empty())
    {
        if (adcWeighted)
        {
            float adcSum{0.f};
            for (const CaloHit *pCaloHit : recoHits)
                adcSum += pCaloHit->GetInputEnergy();
            if (adcSum > std::numeric_limits<float>::epsilon())
            {
                for (const CaloHit *pCaloHit : intersection)
                    purity += pCaloHit->GetInputEnergy();
                purity /= adcSum;
            }
        }
        else
        {
            purity = intersection.size() / static_cast<float>(recoHits.size());
        }
    }

    return purity;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArHierarchyHelper::MCMatches::GetCompleteness(const CaloHitVector &intersection, const CaloHitList &mcHits, const bool adcWeighted) const
{
    float completeness{0.f};
    if (!intersection.empty())
    {
        if (adcWeighted)
        {
            float adcSum{0.f};
            for (const CaloHit *pCaloHit : mcHits)
                adcSum += pCaloHit->GetInputEnergy();
            if (adcSum > std::numeric_limits<float>::epsilon())
            {
                for (const CaloHit *pCaloHit : intersection)
                    completeness += pCaloHit->GetInputEnergy();
                completeness /= adcSum;
            }
        }
        else
        {
            completeness = intersection.size() / static_cast<float>(mcHits.size());
        }
    }

    return completeness;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MatchInfo::QualityCuts::QualityCuts() : m_minPurity{0.5f}, m_minCompleteness{0.1f}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MatchInfo::QualityCuts::QualityCuts(const float minPurity, const float minCompleteness) :
    m_minPurity{minPurity},
    m_minCompleteness{minCompleteness}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MatchInfo::MatchInfo() : MatchInfo(QualityCuts())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArHierarchyHelper::MatchInfo::MatchInfo(const QualityCuts &qualityCuts) : m_qualityCuts{qualityCuts}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MatchInfo::Match(const MCHierarchy &mcHierarchy, const RecoHierarchy &recoHierarchy)
{
    MCHierarchy::NodeVector mcNodes;
    mcHierarchy.GetFlattenedNodes(mcNodes);
    RecoHierarchy::NodeVector recoNodes;
    recoHierarchy.GetFlattenedNodes(recoNodes);

    std::sort(mcNodes.begin(), mcNodes.end(),
        [](const MCHierarchy::Node *lhs, const MCHierarchy::Node *rhs) { return lhs->GetCaloHits().size() > rhs->GetCaloHits().size(); });
    std::sort(recoNodes.begin(), recoNodes.end(),
        [](const RecoHierarchy::Node *lhs, const RecoHierarchy::Node *rhs) { return lhs->GetCaloHits().size() > rhs->GetCaloHits().size(); });

    std::map<const MCHierarchy::Node *, MCMatches> mcToMatchMap;
    for (const RecoHierarchy::Node *pRecoNode : recoNodes)
    {
        const CaloHitList &recoHits{pRecoNode->GetCaloHits()};
        const MCHierarchy::Node *pBestNode{nullptr};
        size_t bestSharedHits{0};
        for (const MCHierarchy::Node *pMCNode : mcNodes)
        {
            if (!pMCNode->IsReconstructable())
                continue;
            const CaloHitList &mcHits{pMCNode->GetCaloHits()};
            CaloHitVector intersection;
            std::set_intersection(mcHits.begin(), mcHits.end(), recoHits.begin(), recoHits.end(), std::back_inserter(intersection));

            if (!intersection.empty())
            {
                const size_t sharedHits{intersection.size()};
                if (sharedHits > bestSharedHits)
                {
                    bestSharedHits = sharedHits;
                    pBestNode = pMCNode;
                }
            }
        }
        if (pBestNode)
        {
            auto iter{mcToMatchMap.find(pBestNode)};
            if (iter != mcToMatchMap.end())
            {
                MCMatches &match(iter->second);
                match.AddRecoMatch(pRecoNode, static_cast<int>(bestSharedHits));
            }
            else
            {
                MCMatches match(pBestNode);
                match.AddRecoMatch(pRecoNode, static_cast<int>(bestSharedHits));
                mcToMatchMap.insert(std::make_pair(pBestNode, match));
            }
        }
        else
        {
            m_unmatchedReco.emplace_back(pRecoNode);
        }
    }

    for (auto [pMCNode, matches] : mcToMatchMap)
    {
        const RecoHierarchy::NodeVector &nodeVector{matches.GetRecoMatches()};
        if (nodeVector.size() == 1)
        {
            const RecoHierarchy::Node *pRecoNode{nodeVector.front()};
            const float purity{matches.GetPurity(pRecoNode)};
            const float completeness{matches.GetCompleteness(pRecoNode)};
            if (purity >= m_qualityCuts.m_minPurity && completeness > m_qualityCuts.m_minCompleteness)
                m_goodMatches.emplace_back(matches);
            else
                m_subThresholdMatches.emplace_back(matches);
        }
        else
        {
            MCMatches aboveThresholdMatches(pMCNode), belowThresholdMatches(pMCNode);
            for (const RecoHierarchy::Node *pRecoNode : nodeVector)
            {
                const float purity{matches.GetPurity(pRecoNode)};
                const float completeness{matches.GetCompleteness(pRecoNode)};
                if (purity >= m_qualityCuts.m_minPurity && completeness > m_qualityCuts.m_minCompleteness)
                    aboveThresholdMatches.AddRecoMatch(pRecoNode, matches.GetSharedHits(pRecoNode));
                else
                    belowThresholdMatches.AddRecoMatch(pRecoNode, matches.GetSharedHits(pRecoNode));
            }
            const size_t nAboveThresholdMatches{aboveThresholdMatches.GetNRecoMatches()};
            if (nAboveThresholdMatches == 1)
                m_goodMatches.emplace_back(aboveThresholdMatches);
            else if (nAboveThresholdMatches > 1)
                m_aboveThresholdMatches.emplace_back(aboveThresholdMatches);
            if (belowThresholdMatches.GetNRecoMatches() > 0)
                m_subThresholdMatches.emplace_back(belowThresholdMatches);
        }
    }

    const auto predicate = [](const MCMatches &lhs, const MCMatches &rhs) {
        return lhs.GetMC()->GetCaloHits().size() > rhs.GetMC()->GetCaloHits().size();
    };
    std::sort(m_goodMatches.begin(), m_goodMatches.end(), predicate);
    std::sort(m_subThresholdMatches.begin(), m_subThresholdMatches.end(), predicate);

    for (const MCHierarchy::Node *pMCNode : mcNodes)
    {
        if (pMCNode->IsReconstructable() && mcToMatchMap.find(pMCNode) == mcToMatchMap.end())
            m_unmatchedMC.emplace_back(pMCNode);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArHierarchyHelper::MatchInfo::GetNMCNodes() const
{
    std::set<const MCHierarchy::Node *> mcNodeSet;
    for (const MCMatches &match : this->GetGoodMatches())
        mcNodeSet.insert(match.GetMC());
    for (const MCMatches &match : this->GetAboveThresholdMatches())
        mcNodeSet.insert(match.GetMC());
    for (const MCMatches &match : this->GetSubThresholdMatches())
        mcNodeSet.insert(match.GetMC());
    for (const MCMatches &match : this->GetUnmatchedMC())
        mcNodeSet.insert(match.GetMC());

    return static_cast<unsigned int>(mcNodeSet.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArHierarchyHelper::MatchInfo::GetNNeutrinoMCNodes() const
{
    std::set<const MCHierarchy::Node *> mcNodeSet;
    for (const MCMatches &match : this->GetGoodMatches())
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (!(pNode->IsCosmicRay() || pNode->IsTestBeamParticle()))
            mcNodeSet.insert(pNode);
    }
    for (const MCMatches &match : this->GetAboveThresholdMatches())
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (!(pNode->IsCosmicRay() || pNode->IsTestBeamParticle()))
            mcNodeSet.insert(match.GetMC());
    }
    for (const MCMatches &match : this->GetSubThresholdMatches())
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (!(pNode->IsCosmicRay() || pNode->IsTestBeamParticle()))
            mcNodeSet.insert(match.GetMC());
    }
    for (const MCMatches &match : this->GetUnmatchedMC())
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (!(pNode->IsCosmicRay() || pNode->IsTestBeamParticle()))
            mcNodeSet.insert(match.GetMC());
    }

    return static_cast<unsigned int>(mcNodeSet.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArHierarchyHelper::MatchInfo::GetNCosmicRayMCNodes() const
{
    std::set<const MCHierarchy::Node *> mcNodeSet;
    for (const MCMatches &match : this->GetGoodMatches())
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (pNode->IsCosmicRay())
            mcNodeSet.insert(pNode);
    }
    for (const MCMatches &match : this->GetAboveThresholdMatches())
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (pNode->IsCosmicRay())
            mcNodeSet.insert(match.GetMC());
    }
    for (const MCMatches &match : this->GetSubThresholdMatches())
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (pNode->IsCosmicRay())
            mcNodeSet.insert(match.GetMC());
    }
    for (const MCMatches &match : this->GetUnmatchedMC())
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (pNode->IsCosmicRay())
            mcNodeSet.insert(match.GetMC());
    }

    return static_cast<unsigned int>(mcNodeSet.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArHierarchyHelper::MatchInfo::GetNTestBeamMCNodes() const
{
    std::set<const MCHierarchy::Node *> mcNodeSet;
    for (const MCMatches &match : this->GetGoodMatches())
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (pNode->IsTestBeamParticle())
            mcNodeSet.insert(pNode);
    }
    for (const MCMatches &match : this->GetAboveThresholdMatches())
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (pNode->IsTestBeamParticle())
            mcNodeSet.insert(match.GetMC());
    }
    for (const MCMatches &match : this->GetSubThresholdMatches())
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (pNode->IsTestBeamParticle())
            mcNodeSet.insert(match.GetMC());
    }
    for (const MCMatches &match : this->GetUnmatchedMC())
    {
        const MCHierarchy::Node *pNode{match.GetMC()};
        if (pNode->IsCosmicRay())
            mcNodeSet.insert(match.GetMC());
    }

    return static_cast<unsigned int>(mcNodeSet.size());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MatchInfo::Print(const MCHierarchy &mcHierarchy) const
{
    unsigned int nNeutrinoMCParticles{this->GetNNeutrinoMCNodes()}, nNeutrinoRecoParticles{0}, nNeutrinoRecoBTParticles{0};
    unsigned int nCosmicMCParticles{this->GetNCosmicRayMCNodes()}, nCosmicRecoParticles{0}, nCosmicRecoBTParticles{0};
    unsigned int nTestBeamMCParticles{this->GetNTestBeamMCNodes()}, nTestBeamRecoParticles{0}, nTestBeamRecoBTParticles{0};
    std::cout << "=== Good matches ===" << std::endl;
    for (const MCMatches &match : this->GetGoodMatches())
    {
        const MCHierarchy::Node *pMCNode{match.GetMC()};
        const int pdg{pMCNode->GetParticleId()};
        const size_t mcHits{pMCNode->GetCaloHits().size()};
        const std::string tag{pMCNode->IsTestBeamParticle() ? "(Beam) " : pMCNode->IsCosmicRay() ? "(Cosmic) " : ""};
        std::cout << "MC " << tag << pdg << " hits " << mcHits << std::endl;
        const RecoHierarchy::NodeVector &nodeVector{match.GetRecoMatches()};

        for (const RecoHierarchy::Node *pRecoNode : nodeVector)
        {
            const unsigned int recoHits{static_cast<unsigned int>(pRecoNode->GetCaloHits().size())};
            const unsigned int sharedHits{match.GetSharedHits(pRecoNode)};
            const float purity{match.GetPurity(pRecoNode)};
            const float completeness{match.GetCompleteness(pRecoNode)};
            std::cout << "   Matched " << sharedHits << " out of " << recoHits << " with purity " << purity << " and completeness "
                      << completeness << std::endl;
        }
        if (pMCNode->IsTestBeamParticle())
            ++nTestBeamRecoParticles;
        else if (pMCNode->IsCosmicRay())
            ++nCosmicRecoParticles;
        else
            ++nNeutrinoRecoParticles;
    }

    std::cout << "=== Above threshold matches ===" << std::endl;
    for (const MCMatches &match : this->GetAboveThresholdMatches())
    {
        const MCHierarchy::Node *pMCNode{match.GetMC()};
        const int pdg{pMCNode->GetParticleId()};
        const size_t mcHits{pMCNode->GetCaloHits().size()};
        const std::string tag{pMCNode->IsTestBeamParticle() ? "(Beam) " : pMCNode->IsCosmicRay() ? "(Cosmic) " : ""};
        std::cout << "MC " << tag << pdg << " hits " << mcHits << std::endl;
        const RecoHierarchy::NodeVector &nodeVector{match.GetRecoMatches()};

        for (const RecoHierarchy::Node *pRecoNode : nodeVector)
        {
            const unsigned int recoHits{static_cast<unsigned int>(pRecoNode->GetCaloHits().size())};
            const unsigned int sharedHits{match.GetSharedHits(pRecoNode)};
            const float purity{match.GetPurity(pRecoNode)};
            const float completeness{match.GetCompleteness(pRecoNode)};
            std::cout << "   Matched " << sharedHits << " out of " << recoHits << " with purity " << purity << " and completeness "
                      << completeness << std::endl;
        }
        if (pMCNode->IsTestBeamParticle())
            ++nTestBeamRecoParticles;
        else if (pMCNode->IsCosmicRay())
            ++nCosmicRecoParticles;
        else
            ++nNeutrinoRecoParticles;
    }

    std::cout << "=== Below threshold matches ===" << std::endl;
    for (const MCMatches &match : this->GetSubThresholdMatches())
    {
        const MCHierarchy::Node *pMCNode{match.GetMC()};
        const int pdg{pMCNode->GetParticleId()};
        const size_t mcHits{pMCNode->GetCaloHits().size()};
        const std::string tag{pMCNode->IsTestBeamParticle() ? "(Beam) " : pMCNode->IsCosmicRay() ? "(Cosmic) " : ""};
        std::cout << "MC " << tag << pdg << " hits " << mcHits << std::endl;
        const RecoHierarchy::NodeVector &nodeVector{match.GetRecoMatches()};

        for (const RecoHierarchy::Node *pRecoNode : nodeVector)
        {
            const unsigned int recoHits{static_cast<unsigned int>(pRecoNode->GetCaloHits().size())};
            const unsigned int sharedHits{match.GetSharedHits(pRecoNode)};
            const float purity{match.GetPurity(pRecoNode)};
            const float completeness{match.GetCompleteness(pRecoNode)};
            std::cout << "   Matched (below threshold) " << sharedHits << " out of " << recoHits << " with purity " << purity
                      << " and completeness " << completeness << std::endl;
        }
        if (pMCNode->IsTestBeamParticle())
            ++nTestBeamRecoBTParticles;
        else if (pMCNode->IsCosmicRay())
            ++nCosmicRecoBTParticles;
        else
            ++nNeutrinoRecoBTParticles;
    }

    std::cout << "=== Unmatched ===" << std::endl;
    for (const MCHierarchy::Node *pMCNode : this->GetUnmatchedMC())
    {
        const int pdg{pMCNode->GetParticleId()};
        const size_t mcHits{pMCNode->GetCaloHits().size()};
        const std::string tag{pMCNode->IsTestBeamParticle() ? "(Beam) " : pMCNode->IsCosmicRay() ? "(Cosmic) " : ""};
        std::cout << "MC " << tag << pdg << " hits " << mcHits << std::endl;
        std::cout << "   Unmatched" << std::endl;
    }

    if (mcHierarchy.IsNeutrinoHierarchy())
    {
        std::cout << "Neutrino Interaction Summary:" << std::endl;
        std::cout << std::fixed << std::setprecision(1);
        if (nNeutrinoMCParticles)
        {
            std::cout << "Matched final state particles: " << nNeutrinoRecoParticles << " of " << nNeutrinoMCParticles << " : "
                      << (100 * nNeutrinoRecoParticles / static_cast<float>(nNeutrinoMCParticles)) << "%" << std::endl;
        }
        if (nCosmicMCParticles)
        {
            std::cout << "Matched cosmics: " << nCosmicRecoParticles << " of " << nCosmicMCParticles << " : "
                      << (100 * nCosmicRecoParticles / static_cast<float>(nCosmicMCParticles)) << "%" << std::endl;
        }
    }
    else if (mcHierarchy.IsTestBeamHierarchy())
    {
        std::cout << "Test Beam Interaction Summary:" << std::endl;
        std::cout << std::fixed << std::setprecision(1);
        if (nTestBeamMCParticles)
        {
            std::cout << "Matched test beam particles: " << nTestBeamRecoParticles << " of " << nTestBeamMCParticles << " : "
                      << (100 * nTestBeamRecoParticles / static_cast<float>(nTestBeamMCParticles)) << "%" << std::endl;
            std::cout << "Loosely matched test beam particles: " << (nTestBeamRecoParticles + nTestBeamRecoBTParticles) << " of " << nTestBeamMCParticles
                      << " : " << (100 * (nTestBeamRecoParticles + nTestBeamRecoBTParticles) / static_cast<float>(nTestBeamMCParticles))
                      << "%" << std::endl;
        }
        if (nCosmicMCParticles)
        {
            std::cout << "Matched cosmics: " << nCosmicRecoParticles << " of " << nCosmicMCParticles << " : "
                      << (100 * nCosmicRecoParticles / static_cast<float>(nCosmicMCParticles)) << "%" << std::endl;
            std::cout << "Loosely matched cosmics: " << (nCosmicRecoParticles + nCosmicRecoBTParticles) << " of " << nCosmicMCParticles << " : "
                      << (100 * (nCosmicRecoParticles + nCosmicRecoBTParticles) / static_cast<float>(nCosmicMCParticles)) << "%" << std::endl;
        }
    }
    if (!this->GetUnmatchedReco().empty())
        std::cout << "Unmatched reco: " << this->GetUnmatchedReco().size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::FillMCHierarchy(const MCParticleList &mcParticleList, const CaloHitList &caloHitList, const bool foldToPrimaries,
    const bool foldToLeadingShowers, MCHierarchy &hierarchy)
{
    hierarchy.FillHierarchy(mcParticleList, caloHitList, foldToPrimaries, foldToLeadingShowers);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::FillRecoHierarchy(const PfoList &pfoList, const bool foldToPrimaries, const bool foldToLeadingShowers, RecoHierarchy &hierarchy)
{
    hierarchy.FillHierarchy(pfoList, foldToPrimaries, foldToLeadingShowers);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArHierarchyHelper::MatchHierarchies(const MCHierarchy &mcHierarchy, const RecoHierarchy &recoHierarchy, MatchInfo &matchInfo)
{
    matchInfo.Match(mcHierarchy, recoHierarchy);
}

//------------------------------------------------------------------------------------------------------------------------------------------
// private

const MCParticle *LArHierarchyHelper::GetMCPrimaries(const MCParticleList &mcParticleList, MCParticleSet &primaries)
{
    const MCParticle *pRoot{nullptr};
    for (const MCParticle *pMCParticle : mcParticleList)
    {
        try
        {
            const MCParticle *const pPrimary{LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle)};
            if (!LArMCParticleHelper::IsTriggeredBeamParticle(pPrimary))
                primaries.insert(pPrimary);
            else
                primaries.insert(pPrimary);
        }
        catch (const StatusCodeException &)
        {
            if (LArMCParticleHelper::IsNeutrino(pMCParticle))
                pRoot = pMCParticle;
            else if (pMCParticle->GetParticleId() != 111 && pMCParticle->GetParticleId() < 1e9)
                std::cout << "LArHierarchyHelper::MCHierarchy::FillHierarchy: MC particle with PDG code " << pMCParticle->GetParticleId()
                          << " at address " << pMCParticle << " has no associated primary particle" << std::endl;
        }
    }

    return pRoot;
}

const ParticleFlowObject *LArHierarchyHelper::GetRecoPrimaries(const PfoList &pfoList, PfoSet &primaries)
{
    const ParticleFlowObject *pRoot{nullptr};
    PfoSet cosmicPfos;
    for (const ParticleFlowObject *pPfo : pfoList)
    {
        if (LArPfoHelper::IsNeutrino(pPfo))
        {
            pRoot = pPfo;
            break;
        }
        else
        {
            const ParticleFlowObject *const pParent{LArPfoHelper::GetParentPfo(pPfo)};
            if (pParent && LArPfoHelper::IsNeutrino(pParent))
            {
                pRoot = pParent;
                break;
            }
            else
            {
                // Should be in a test beam scenario
                const int tier{LArPfoHelper::GetHierarchyTier(pPfo)};
                if (tier == 0 && LArPfoHelper::IsTestBeam(pPfo))
                {
                    // Triggered beam particle
                    primaries.insert(pPfo);
                    continue;
                }
                if (tier > 1)
                    continue;
                if (!LArPfoHelper::IsTestBeam(pPfo))
                {
                    // Cosmic induced
                    cosmicPfos.insert(pPfo);
                }
            }
        }
    }
    if (pRoot && LArPfoHelper::IsNeutrino(pRoot))
    {
        const PfoList &children{pRoot->GetDaughterPfoList()};
        for (const ParticleFlowObject *pPrimary : children)
            primaries.insert(pPrimary);
    }
    else
    {
        for (const ParticleFlowObject *pPfo : cosmicPfos)
            primaries.insert(pPfo);
    }

    return pRoot;
}

} // namespace lar_content
