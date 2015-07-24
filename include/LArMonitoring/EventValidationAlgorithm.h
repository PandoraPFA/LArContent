/**
 *  @file   LArContent/include/LArMonitoring/EventValidationAlgorithm.h
 *
 *  @brief  Header file for the event validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_VALIDATION_ALGORITHM_H
#define LAR_EVENT_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  EventValidationAlgorithm class
 */
class EventValidationAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Default constructor
     */
    EventValidationAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~EventValidationAlgorithm();

private:
    /**
     *  @brief SimpleMCPrimary class
     */
    class SimpleMCPrimary
    {
    public:
        /**
         *  @brief  Constructor
         */
        SimpleMCPrimary();

        int                                 m_pdgCode;                  ///< The pdg code
        unsigned int                        m_nMCHitsTotal;             ///< The total number of mc hits
        unsigned int                        m_nMCHitsU;                 ///< The number of u mc hits
        unsigned int                        m_nMCHitsV;                 ///< The number of v mc hits
        unsigned int                        m_nMCHitsW;                 ///< The number of w mc hits
        unsigned int                        m_nMatchedPfos;             ///< The number of matched pfos
        const pandora::MCParticle          *m_pPandoraAddress;          ///< The address of the Pandora mc primary
        // TODO other properties
    };

    typedef std::vector<SimpleMCPrimary> SimpleMCPrimaryList;

    /**
     *  @brief SimpleMatchedPfo class
     */
    class SimpleMatchedPfo
    {
    public:
        /**
         *  @brief  Constructor
         */
        SimpleMatchedPfo();

        int                                 m_pdgCode;                  ///< The pdg code
        unsigned int                        m_nPfoHitsTotal;            ///< The total number of pfo hits
        unsigned int                        m_nPfoHitsU;                ///< The number of u pfo hits
        unsigned int                        m_nPfoHitsV;                ///< The number of v pfo hits
        unsigned int                        m_nPfoHitsW;                ///< The number of w pfo hits
        unsigned int                        m_nMatchedHitsTotal;        ///< The total number of matched hits
        unsigned int                        m_nMatchedHitsU;            ///< The number of u matched hits
        unsigned int                        m_nMatchedHitsV;            ///< The number of v matched hits
        unsigned int                        m_nMatchedHitsW;            ///< The number of w matched hits
        const pandora::ParticleFlowObject  *m_pPandoraAddress;          ///< The address of the Pandora mc primary
        // TODO other properties
    };

    typedef std::vector<SimpleMatchedPfo> SimpleMatchedPfoList;

    /**
     *  @brief  Sort simple mc primaries by number of mc hits
     * 
     *  @param  lhs the left-hand side
     *  @param  rhs the right-hand side
     * 
     *  @return boolean
     */
    static bool SortSimpleMCPrimaries(const SimpleMCPrimary &lhs, const SimpleMCPrimary &rhs);

    /**
     *  @brief  Sort simple matched pfos by number of matched hits
     * 
     *  @param  lhs the left-hand side
     *  @param  rhs the right-hand side
     * 
     *  @return boolean
     */
    static bool SortSimpleMatchedPfos(const SimpleMatchedPfo &lhs, const SimpleMatchedPfo &rhs);

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string         m_caloHitListName;          ///< Name of input calo hit list
    std::string         m_mcParticleListName;       ///< Name of input MC particle list
    std::string         m_pfoListName;              ///< Name of input Pfo list
    std::string         m_fileName;                 ///< Name of output file
    std::string         m_treeName;                 ///< Name of output tree

    bool                m_extractNeutrinoDaughters; ///< Whether to treat each neutrino pfo daughter as a standalone top-level pfo
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *EventValidationAlgorithm::Factory::CreateAlgorithm() const
{
    return new EventValidationAlgorithm();
}

} // namespace lar_content

#endif // LAR_EVENT_VALIDATION_ALGORITHM_H
