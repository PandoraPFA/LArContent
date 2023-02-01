/**
 *  @file   larpandoracontent/LArHelpers/LArInteractionTypeHelper.h
 *
 *  @brief  Header file for the interaction type helper class.
 *
 *  $Log: $
 */
#ifndef LAR_INTERACTION_TYPE_HELPER_H
#define LAR_INTERACTION_TYPE_HELPER_H 1

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"

namespace lar_content
{

// clang-format off
#define INTERACTION_TYPE_TABLE(d)                         \
    d(CCQEL_MU                                          ) \
    d(CCQEL_MU_P                                        ) \
    d(CCQEL_MU_P_P                                      ) \
    d(CCQEL_MU_P_P_P                                    ) \
    d(CCQEL_MU_P_P_P_P                                  ) \
    d(CCQEL_MU_P_P_P_P_P                                ) \
    d(CCQEL_E                                           ) \
    d(CCQEL_E_P                                         ) \
    d(CCQEL_E_P_P                                       ) \
    d(CCQEL_E_P_P_P                                     ) \
    d(CCQEL_E_P_P_P_P                                   ) \
    d(CCQEL_E_P_P_P_P_P                                 ) \
    d(NCQEL_P                                           ) \
    d(NCQEL_P_P                                         ) \
    d(NCQEL_P_P_P                                       ) \
    d(NCQEL_P_P_P_P                                     ) \
    d(NCQEL_P_P_P_P_P                                   ) \
    d(CCRES_MU                                          ) \
    d(CCRES_MU_P                                        ) \
    d(CCRES_MU_P_P                                      ) \
    d(CCRES_MU_P_P_P                                    ) \
    d(CCRES_MU_P_P_P_P                                  ) \
    d(CCRES_MU_P_P_P_P_P                                ) \
    d(CCRES_MU_PIPLUS                                   ) \
    d(CCRES_MU_P_PIPLUS                                 ) \
    d(CCRES_MU_P_P_PIPLUS                               ) \
    d(CCRES_MU_P_P_P_PIPLUS                             ) \
    d(CCRES_MU_P_P_P_P_PIPLUS                           ) \
    d(CCRES_MU_P_P_P_P_P_PIPLUS                         ) \
    d(CCRES_MU_PHOTON                                   ) \
    d(CCRES_MU_P_PHOTON                                 ) \
    d(CCRES_MU_P_P_PHOTON                               ) \
    d(CCRES_MU_P_P_P_PHOTON                             ) \
    d(CCRES_MU_P_P_P_P_PHOTON                           ) \
    d(CCRES_MU_P_P_P_P_P_PHOTON                         ) \
    d(CCRES_MU_PIZERO                                   ) \
    d(CCRES_MU_P_PIZERO                                 ) \
    d(CCRES_MU_P_P_PIZERO                               ) \
    d(CCRES_MU_P_P_P_PIZERO                             ) \
    d(CCRES_MU_P_P_P_P_PIZERO                           ) \
    d(CCRES_MU_P_P_P_P_P_PIZERO                         ) \
    d(CCRES_E                                           ) \
    d(CCRES_E_P                                         ) \
    d(CCRES_E_P_P                                       ) \
    d(CCRES_E_P_P_P                                     ) \
    d(CCRES_E_P_P_P_P                                   ) \
    d(CCRES_E_P_P_P_P_P                                 ) \
    d(CCRES_E_PIPLUS                                    ) \
    d(CCRES_E_P_PIPLUS                                  ) \
    d(CCRES_E_P_P_PIPLUS                                ) \
    d(CCRES_E_P_P_P_PIPLUS                              ) \
    d(CCRES_E_P_P_P_P_PIPLUS                            ) \
    d(CCRES_E_P_P_P_P_P_PIPLUS                          ) \
    d(CCRES_E_PHOTON                                    ) \
    d(CCRES_E_P_PHOTON                                  ) \
    d(CCRES_E_P_P_PHOTON                                ) \
    d(CCRES_E_P_P_P_PHOTON                              ) \
    d(CCRES_E_P_P_P_P_PHOTON                            ) \
    d(CCRES_E_P_P_P_P_P_PHOTON                          ) \
    d(CCRES_E_PIZERO                                    ) \
    d(CCRES_E_P_PIZERO                                  ) \
    d(CCRES_E_P_P_PIZERO                                ) \
    d(CCRES_E_P_P_P_PIZERO                              ) \
    d(CCRES_E_P_P_P_P_PIZERO                            ) \
    d(CCRES_E_P_P_P_P_P_PIZERO                          ) \
    d(NCRES_P                                           ) \
    d(NCRES_P_P                                         ) \
    d(NCRES_P_P_P                                       ) \
    d(NCRES_P_P_P_P                                     ) \
    d(NCRES_P_P_P_P_P                                   ) \
    d(NCRES_PIPLUS                                      ) \
    d(NCRES_P_PIPLUS                                    ) \
    d(NCRES_P_P_PIPLUS                                  ) \
    d(NCRES_P_P_P_PIPLUS                                ) \
    d(NCRES_P_P_P_P_PIPLUS                              ) \
    d(NCRES_P_P_P_P_P_PIPLUS                            ) \
    d(NCRES_PIMINUS                                     ) \
    d(NCRES_P_PIMINUS                                   ) \
    d(NCRES_P_P_PIMINUS                                 ) \
    d(NCRES_P_P_P_PIMINUS                               ) \
    d(NCRES_P_P_P_P_PIMINUS                             ) \
    d(NCRES_P_P_P_P_P_PIMINUS                           ) \
    d(NCRES_PHOTON                                      ) \
    d(NCRES_P_PHOTON                                    ) \
    d(NCRES_P_P_PHOTON                                  ) \
    d(NCRES_P_P_P_PHOTON                                ) \
    d(NCRES_P_P_P_P_PHOTON                              ) \
    d(NCRES_P_P_P_P_P_PHOTON                            ) \
    d(NCRES_PIZERO                                      ) \
    d(NCRES_P_PIZERO                                    ) \
    d(NCRES_P_P_PIZERO                                  ) \
    d(NCRES_P_P_P_PIZERO                                ) \
    d(NCRES_P_P_P_P_PIZERO                              ) \
    d(NCRES_P_P_P_P_P_PIZERO                            ) \
    d(CCDIS_MU                                          ) \
    d(CCDIS_MU_P                                        ) \
    d(CCDIS_MU_P_P                                      ) \
    d(CCDIS_MU_P_P_P                                    ) \
    d(CCDIS_MU_P_P_P_P                                  ) \
    d(CCDIS_MU_P_P_P_P_P                                ) \
    d(CCDIS_MU_PIPLUS                                   ) \
    d(CCDIS_MU_P_PIPLUS                                 ) \
    d(CCDIS_MU_P_P_PIPLUS                               ) \
    d(CCDIS_MU_P_P_P_PIPLUS                             ) \
    d(CCDIS_MU_P_P_P_P_PIPLUS                           ) \
    d(CCDIS_MU_P_P_P_P_P_PIPLUS                         ) \
    d(CCDIS_MU_PHOTON                                   ) \
    d(CCDIS_MU_P_PHOTON                                 ) \
    d(CCDIS_MU_P_P_PHOTON                               ) \
    d(CCDIS_MU_P_P_P_PHOTON                             ) \
    d(CCDIS_MU_P_P_P_P_PHOTON                           ) \
    d(CCDIS_MU_P_P_P_P_P_PHOTON                         ) \
    d(CCDIS_MU_PIZERO                                   ) \
    d(CCDIS_MU_P_PIZERO                                 ) \
    d(CCDIS_MU_P_P_PIZERO                               ) \
    d(CCDIS_MU_P_P_P_PIZERO                             ) \
    d(CCDIS_MU_P_P_P_P_PIZERO                           ) \
    d(CCDIS_MU_P_P_P_P_P_PIZERO                         ) \
    d(NCDIS_P                                           ) \
    d(NCDIS_P_P                                         ) \
    d(NCDIS_P_P_P                                       ) \
    d(NCDIS_P_P_P_P                                     ) \
    d(NCDIS_P_P_P_P_P                                   ) \
    d(NCDIS_PIPLUS                                      ) \
    d(NCDIS_P_PIPLUS                                    ) \
    d(NCDIS_P_P_PIPLUS                                  ) \
    d(NCDIS_P_P_P_PIPLUS                                ) \
    d(NCDIS_P_P_P_P_PIPLUS                              ) \
    d(NCDIS_P_P_P_P_P_PIPLUS                            ) \
    d(NCDIS_PIMINUS                                     ) \
    d(NCDIS_P_PIMINUS                                   ) \
    d(NCDIS_P_P_PIMINUS                                 ) \
    d(NCDIS_P_P_P_PIMINUS                               ) \
    d(NCDIS_P_P_P_P_PIMINUS                             ) \
    d(NCDIS_P_P_P_P_P_PIMINUS                           ) \
    d(NCDIS_PHOTON                                      ) \
    d(NCDIS_P_PHOTON                                    ) \
    d(NCDIS_P_P_PHOTON                                  ) \
    d(NCDIS_P_P_P_PHOTON                                ) \
    d(NCDIS_P_P_P_P_PHOTON                              ) \
    d(NCDIS_P_P_P_P_P_PHOTON                            ) \
    d(NCDIS_PIZERO                                      ) \
    d(NCDIS_P_PIZERO                                    ) \
    d(NCDIS_P_P_PIZERO                                  ) \
    d(NCDIS_P_P_P_PIZERO                                ) \
    d(NCDIS_P_P_P_P_PIZERO                              ) \
    d(NCDIS_P_P_P_P_P_PIZERO                            ) \
    d(CCCOH                                             ) \
    d(NCCOH                                             ) \
    d(COSMIC_RAY_MU                                     ) \
    d(COSMIC_RAY_P                                      ) \
    d(COSMIC_RAY_E                                      ) \
    d(COSMIC_RAY_PHOTON                                 ) \
    d(COSMIC_RAY_OTHER                                  ) \
    d(BEAM_PARTICLE_MU                                  ) \
    d(BEAM_PARTICLE_P                                   ) \
    d(BEAM_PARTICLE_E                                   ) \
    d(BEAM_PARTICLE_PHOTON                              ) \
    d(BEAM_PARTICLE_PI_PLUS                             ) \
    d(BEAM_PARTICLE_PI_MINUS                            ) \
    d(BEAM_PARTICLE_KAON_PLUS                           ) \
    d(BEAM_PARTICLE_KAON_MINUS                          ) \
    d(BEAM_PARTICLE_OTHER                               ) \
    d(BEAM_PARTICLE_PI_PLUS_PI_PLUS                     ) \
    d(BEAM_PARTICLE_PI_PLUS_PI_PLUS_PHOTON              ) \
    d(BEAM_PARTICLE_PI_PLUS_PI_PLUS_PIZERO              ) \
    d(BEAM_PARTICLE_PI_PLUS_COMPLEX                     ) \
    d(BEAM_PARTICLE_PI_MINUS_PI_MINUS                   ) \
    d(BEAM_PARTICLE_PI_MINUS_PI_MINUS_PHOTON            ) \
    d(BEAM_PARTICLE_PI_MINUS_PI_MINUS_PIZERO            ) \
    d(BEAM_PARTICLE_PI_MINUS_COMPLEX                    ) \
    d(BEAM_PARTICLE_P_P                                 ) \
    d(BEAM_PARTICLE_P_P_PHOTON                          ) \
    d(BEAM_PARTICLE_P_P_PHOTON_PHOTON                   ) \
    d(BEAM_PARTICLE_P_P_PHOTON_PHOTON_PHOTON            ) \
    d(BEAM_PARTICLE_P_P_PHOTON_PHOTON_PHOTON_PHOTON     ) \
    d(BEAM_PARTICLE_P_P_P                               ) \
    d(BEAM_PARTICLE_P_P_P_PHOTON                        ) \
    d(BEAM_PARTICLE_P_P_P_PHOTON_PHOTON                 ) \
    d(BEAM_PARTICLE_P_P_P_PHOTON_PHOTON_PHOTON          ) \
    d(BEAM_PARTICLE_P_P_P_P                             ) \
    d(BEAM_PARTICLE_P_P_P_P_PHOTON                      ) \
    d(BEAM_PARTICLE_P_P_P_P_PHOTON_PHOTON               ) \
    d(BEAM_PARTICLE_P_P_P_P_P                           ) \
    d(BEAM_PARTICLE_P_P_P_P_P_PHOTON                    ) \
    d(BEAM_PARTICLE_P_P_P_P_P_P                         ) \
    d(BEAM_PARTICLE_P_COMPLEX                           ) \
    d(BEAM_PARTICLE_MU_E                                ) \
    d(BEAM_PARTICLE_MU_COMPLEX                          ) \
    d(BEAM_PARTICLE_KAON_PLUS_MU                        ) \
    d(BEAM_PARTICLE_KAON_PLUS_KAON_PLUS_KAON0L_COMPLEX  ) \
    d(BEAM_PARTICLE_KAON_PLUS_KAON_PLUS_COMPLEX         ) \
    d(BEAM_PARTICLE_KAON_PLUS_COMPLEX                   ) \
    d(BEAM_PARTICLE_KAON_MINUS_MU                       ) \
    d(BEAM_PARTICLE_KAON_MINUS_KAON_MINUS_KAON0L_COMPLEX) \
    d(BEAM_PARTICLE_KAON_MINUS_KAON_MINUS_COMPLEX       ) \
    d(BEAM_PARTICLE_KAON_MINUS_COMPLEX                  ) \
    d(BEAM_PARTICLE_E_COMPLEX                           ) \
    d(BEAM_PARTICLE_COMPLEX_HIERARCHY                   ) \
    d(BEAM_PARTICLE_UNKNOWN_HIERARCHY                   ) \
    d(OTHER_INTERACTION                                 ) \
    d(ALL_INTERACTIONS                                  )

/**
 *  @brief  The particle type enum macro
 */
#define GET_INTERACTION_TYPE_ENTRY(a)                     \
    a,

/**
 *  @brief  The name switch statement macro
 */
#define GET_INTERACTION_TYPE_NAME_SWITCH(a)               \
    case a : return std::string(#a);
// clang-format on

class InteractionDescriptor;

/**
 *  @brief  LArInteractionTypeHelper class
 */
class LArInteractionTypeHelper
{
public:
    /**
     *  @brief  InteractionType enum
     */
    enum InteractionType : unsigned int
    {
        INTERACTION_TYPE_TABLE(GET_INTERACTION_TYPE_ENTRY) UNKNOWN_INTERACTION_TYPE = 0
    };

    /**
     *  @brief  Interaction parameters
     */
    class InteractionParameters
    {
    public:
        /**
         *  @brief  Constructor
         */
        InteractionParameters();

        unsigned int m_nNonNeutrons;
        unsigned int m_nMuons;
        unsigned int m_nElectrons;
        unsigned int m_nPhotons;
        unsigned int m_nProtons;
        unsigned int m_nPiPlus;
        unsigned int m_nPiMinus;
        unsigned int m_nPiZero;
        unsigned int m_nKaonPlus;
        unsigned int m_nKaonMinus;
        unsigned int m_nKaon0L;
    };

    /**
     *  @brief  Get the interaction type of an event
     *
     *  @param  mcPrimaryList the list of primary mc particles
     *
     *  @return interaction type
     */
    static InteractionType GetInteractionType(const pandora::MCParticleList &mcPrimaryList);

    /**
     *  @brief  Get the interaction descriptor of an event
     *
     *  @param  mcPrimaryList the list of primary mc particles
     *
     *  @return interaction descriptor
     */
    static InteractionDescriptor GetInteractionDescriptor(const pandora::MCParticleList &mcPrimaryList);

    /**
     *  @brief  Get the test beam hierarchy interaction type of an event
     *
     *  @param  mcPrimaryList the list of primary mc particles
     *
     *  @return interaction type
     */
    static InteractionType GetTestBeamHierarchyInteractionType(const pandora::MCParticleList &mcPrimaryList);

    /**
     *  @brief  Set parameters describing the number and species of primary interaction products
     *
     *  @param  mcPrimaryList the list of primary mc particles
     *  @param  parameters the parameter block to set
     */
    static void SetInteractionParameters(const pandora::MCParticleList &mcPrimaryList, InteractionParameters &parameters);

    /**
     *  @brief  Get the interaction type of an event under a cosmic ray hypothesis
     *
     *  @param  mcPrimaryList the list of primary mc particles
     *  @param  parameters the parameter block to set
     *
     *  @return interaction type
     */
    static InteractionType CosmicRayHypothesis(const pandora::MCParticleList &mcPrimaryList, const InteractionParameters &parameters);

    /**
     *  @brief  Get a string representation of an interaction type
     *
     *  @param  interactionType the interaction type
     *
     *  @return string
     */
    static std::string ToString(const InteractionType interactionType);
};

class InteractionDescriptor
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  isCC Whether or not the interaction is a CC interaction
     *  @param  isQE Whether or not the interaction is a QE interaction
     *  @param  isRes Whether or not the interaction is a resonant interaction
     *  @param  isDIS Whether or not the interaction is a DIS interaction
     *  @param  isCoherent Whether or not the interaction is a coherent interaction
     *  @param  isNumu Whether or not the interaction is a numu interaction
     *  @param  isNue Whether or not the interaction is a nue interaction
     *  @param  nPiPlus Number of pi plus in the final state
     *  @param  nPiMinus Number of pi minus in the final state
     *  @param  nPhotons Number of photons in the final state
     *  @param  nProtons Number of protons in the final state
     */
    InteractionDescriptor(const bool isCC, const bool isQE, const bool isRes, const bool isDIS, const bool isCoherent, const bool isNumu,
        const bool isNue, const unsigned int nPiPlus, const unsigned int nPiMinus, const unsigned int nPhotons, const unsigned int nProtons);

    /**
     *  @brief  Whether or not the interaction is CC
     *
     *  @return Whether or not the interaction is CC
     */
    bool IsCC() const;

    /**
     *  @brief  Whether or not the interaction is QE
     *
     *  @return Whether or not the interaction is QE
     */
    bool IsQE() const;

    /**
     *  @brief  Whether or not the interaction is resonant
     *
     *  @return Whether or not the interaction is resonant
     */
    bool IsResonant() const;

    /**
     *  @brief  Whether or not the interaction is DIS
     *
     *  @return Whether or not the interaction is DIS
     */
    bool IsDIS() const;

    /**
     *  @brief  Whether or not the interaction is coherent
     *
     *  @return Whether or not the interaction is coherent
     */
    bool IsCoherent() const;

    /**
     *  @brief  Whether or not the interaction is muon neutrino
     *
     *  @return Whether or not the interaction is muon neutrino
     */
    bool IsMuonNeutrino() const;

    /**
     *  @brief  Whether or not the interaction is electron neutrino
     *
     *  @return Whether or not the interaction is electron neutrino
     */
    bool IsElectronNeutrino() const;

    /**
     *  @brief  Retrieve the number of pi zeros
     *
     *  @return Retrieve the number of pi zeros
     */
    unsigned int GetNumPiZero() const;

    /**
     *  @brief  Retrieve the number of pi plus
     *
     *  @return Retrieve the number of pi plus
     */
    unsigned int GetNumPiPlus() const;

    /**
     *  @brief  Retrieve the number of pi minus
     *
     *  @return Retrieve the number of pi minus
     */
    unsigned int GetNumPiMinus() const;

    /**
     *  @brief  Retrieve the number of photons
     *
     *  @return Retrieve the number of photons
     */
    unsigned int GetNumPhotons() const;

    /**
     *  @brief  Retrieve the number of protons
     *
     *  @return Retrieve the number of protons
     */
    unsigned int GetNumProtons() const;

    /**
     *  @brief  Retrieve a unique ID describing the event (that is, a given type of interaction with specific final state particles has a
     *          particular ID)
     *
     *  @return The unique identifier for the event
     */
    int GetUniqueId() const;

    /**
     *  @brief  Retrieve the string descriptor for the event
     *
     *  @return The string descriptor for the event
     */
    const std::string &ToString() const;

    static const int CC;
    static const int NC;
    static const int QE;
    static const int RES;
    static const int DIS;
    static const int COH;
    static const int OTH;
    static const int MU;
    static const int E;
    static const int PIZERO;
    static const int PIPLUS;
    static const int PIMINUS;
    static const int PHOTON;
    static const int NP;

private:
    const bool m_isCC;
    const bool m_isQE;
    const bool m_isResonant;
    const bool m_isDIS;
    const bool m_isCoherent;
    const bool m_isNumu;
    const bool m_isNue;
    const unsigned int m_nPiZero;
    const unsigned int m_nPiPlus;
    const unsigned int m_nPiMinus;
    const unsigned int m_nPhotons;
    const unsigned int m_nProtons;
    int m_id;
    std::string m_descriptor;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool InteractionDescriptor::IsCC() const
{
    return m_isCC;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool InteractionDescriptor::IsQE() const
{
    return m_isQE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool InteractionDescriptor::IsResonant() const
{
    return m_isResonant;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool InteractionDescriptor::IsDIS() const
{
    return m_isDIS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool InteractionDescriptor::IsCoherent() const
{
    return m_isCoherent;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool InteractionDescriptor::IsMuonNeutrino() const
{
    return m_isNumu;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool InteractionDescriptor::IsElectronNeutrino() const
{
    return m_isNue;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int InteractionDescriptor::GetNumPiZero() const
{
    return m_nPiZero;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int InteractionDescriptor::GetNumPiPlus() const
{
    return m_nPiPlus;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int InteractionDescriptor::GetNumPiMinus() const
{
    return m_nPiMinus;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int InteractionDescriptor::GetNumPhotons() const
{
    return m_nPhotons;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int InteractionDescriptor::GetNumProtons() const
{
    return m_nProtons;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int InteractionDescriptor::GetUniqueId() const
{
    return m_id;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &InteractionDescriptor::ToString() const
{
    return m_descriptor;
}

} // namespace lar_content

#endif // #ifndef LAR_INTERACTION_TYPE_HELPER_H
