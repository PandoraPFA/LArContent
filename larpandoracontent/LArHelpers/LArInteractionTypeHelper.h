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

// Specify (name)
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

/**
 *  @brief  LArInteractionTypeHelper class
 */
class LArInteractionTypeHelper
{
public:
    /**
     *  @brief   InteractionType enum
     */
    enum InteractionType : unsigned int
    {
        INTERACTION_TYPE_TABLE(GET_INTERACTION_TYPE_ENTRY)
        UNKNOWN_INTERACTION_TYPE = 0
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
     *  @brief  Get the test beam hierarchy interaction type of an event
     *
     *  @param  mcPrimaryList the list of primary mc particles
     *
     *  @return interaction type
     */
    static InteractionType GetTestBeamHierarchyInteractionType(const pandora::MCParticleList &mcPrimaryList);

    /**
     *  @brief  Get a string representation of an interaction type
     *
     *  @param  interactionType the interaction type
     *
     *  @return string
     */
    static std::string ToString(const InteractionType interactionType);
};

} // namespace lar_content

#endif // #ifndef LAR_INTERACTION_TYPE_HELPER_H

