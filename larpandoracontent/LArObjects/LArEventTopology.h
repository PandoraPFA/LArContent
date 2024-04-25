/**
 *  @file   larpandoracontent/LArObjects/LArEventTopology.h
 *
 *  @brief  Header file for lar event topology.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_TOPOLOGY_H
#define LAR_EVENT_TOPOLOGY_H 1

#include "Objects/MCParticle.h"

#include <map>

namespace lar_content
{

/**
 *  @brief  LArEventTopology class
 */
class LArEventTopology
{
private:
    typedef std::map<const pandora::MCParticle *, pandora::CaloHitList> MCHitMap;

    class Particle
    {
    public:
        /**
         *  @brief  Default constructor
         */
        Particle(const pandora::MCParticle *const pRoot);

        /**
         *  @brief  Default constructor
         *
         *  @param  pRoot The MC particle
         *  @param  caloHitList The list of hits associated with this particle
         */
        Particle(const pandora::MCParticle *const pRoot, const pandora::CaloHitList &caloHitList);

        /**
         *  @brief  Destructor
         */
        virtual ~Particle();

        /**
         *  @brief  Add a child particle to this particle
         *
         *  @param  pChild The child particle to add
         */
        void AddChild(Particle *pChild);

        /**
         *  @brief  Extract the clear topological vertices from the event
         */
        void GetVertices(pandora::CartesianPointVector &vertices) const;

        /**
         *  @brief  Assess this child particle and its descendents to see if it is topologically significant either indepednently
         *          or within the context of downstream activity
         *
         *  @param  mcHitMap The map from MC particles to hits
         */
        void Parse(const MCHitMap &mcHitMap);

        /**
         *  @brief  Prune particles from the hierarchy if they are not deemed topologically significant
         */
        void Prune();

        /**
         *  @brief  Print the visible hierarchy
         *
         *  @param  mcHitMap The map from MC particles to hits
         *  @param  indent The level of indent for printing
         *
         *  @return The string representation of this particle and its children
         */
        const std::string Print(const MCHitMap &mcHitMap, const std::string &indent) const;

    private:
        pandora::MCParticleList m_particles;
        pandora::CaloHitList m_caloHits;
        std::list<Particle *> m_children;
        bool m_fold;
    };

public:
    /**
     *  @brief  Constructor
     *
     *  @param  caloHitList2D The collection of all hits across all views
     */
    LArEventTopology(const pandora::CaloHitList &caloHitList2D);

    /**
     *  @brief  Destructor
     */
    virtual ~LArEventTopology();

    /**
     *  @brief  Construct a particle hierarchy based on the key topological features in an event
     */
    void ConstructVisibleHierarchy();

    /**
     *  @brief  Extract the clear topological vertices from the event
     */
    void GetVertices(pandora::CartesianPointVector &vertices) const;

    /**
     *  @brief  Fold or remove particles that aren't substantive parts of the hierarchy
     */
    void PruneHierarchy();

    /**
     *  @brief  Print the visible hierarchy
     */
    void Print() const;

private:
    /**
     *  @brief  Construct a particle hierarhcy based on the key topological features in an event
     *
     *  @param  pParticle The parent particle relative to which the hierarchy should be constructed
     *  @param  pRootMC The underlying MC particle associated with the parent particle 
     */
    void ConstructVisibleHierarchy(Particle *pParticle, const pandora::MCParticle *const pRootMC);

    MCHitMap m_mcHitMap;
    const pandora::MCParticle *m_pRoot;
    Particle *m_pParticle;
};

} // namespace lar_content

#endif // #ifndef LAR_EVENT_TOPOLOGY_H
