/**
 *  @file   larpandoracontent/LArMonitoring/EventValidationBaseAlgorithm.h
 *
 *  @brief  Header file for the event validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_VALIDATION_BASE_ALGORITHM_H
#define LAR_EVENT_VALIDATION_BASE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

#include <map>

namespace lar_content
{

/**
 *  @brief  EventValidationBaseAlgorithm class
 */
class EventValidationBaseAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EventValidationBaseAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~EventValidationBaseAlgorithm();

protected:
   /**
     *  @brief  ValidationInfo class
     */
    class ValidationInfo
    {
    public:
        /**
         *  @brief  Get the all mc particle to hits map
         *
         *  @return the all mc particle to hits map
         */
        const LArMCParticleHelper::MCContributionMap &GetAllMCParticleToHitsMap() const;

        /**
         *  @brief  Get the target mc particle to hits map
         *
         *  @return the target mc particle to hits map
         */
        const LArMCParticleHelper::MCContributionMap &GetTargetMCParticleToHitsMap() const;

        /**
         *  @brief  Get the pfo to hits map
         *
         *  @return the pfo to hits map
         */
        const LArMCParticleHelper::PfoContributionMap &GetPfoToHitsMap() const;

        /**
         *  @brief  Get the mc to pfo hit sharing map
         *
         *  @return the mc to pfo hit sharing map
         */
        const LArMCParticleHelper::MCParticleToPfoHitSharingMap &GetMCToPfoHitSharingMap() const;

        /**
         *  @brief  Get the interpreted mc to pfo hit sharing map
         *
         *  @return the interpreted mc to pfo hit sharing map
         */
        const LArMCParticleHelper::MCParticleToPfoHitSharingMap &GetInterpretedMCToPfoHitSharingMap() const;

        /**
         *  @brief  Set the all mc particle to hits map
         *
         *  @param  allMCParticleToHitsMap the all mc particle to hits map
         */
        void SetAllMCParticleToHitsMap(const LArMCParticleHelper::MCContributionMap &allMCParticleToHitsMap);

        /**
         *  @brief  Set the target mc particle to hits map
         *
         *  @param  targetMCParticleToHitsMap the target mc particle to hits map
         */
        void SetTargetMCParticleToHitsMap(const LArMCParticleHelper::MCContributionMap &targetMCParticleToHitsMap);

        /**
         *  @brief  Set the pfo to hits map
         *
         *  @param  pfoToHitsMap the pfo to hits map
         */
        void SetPfoToHitsMap(const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap);

        /**
         *  @brief  Set the mc to pfo hit sharing map
         *
         *  @param  mcToPfoHitSharingMap the mc to pfo hit sharing map
         */
        void SetMCToPfoHitSharingMap(const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap);

        /**
         *  @brief  Set the interpreted mc to pfo hit sharing map
         *
         *  @param  interpretedMCToPfoHitSharingMap the interpreted mc to pfo hit sharing map
         */
        void SetInterpretedMCToPfoHitSharingMap(const LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap);

    private:
        LArMCParticleHelper::MCContributionMap              m_allMCParticleToHitsMap;               ///< The all mc particle to hits map
        LArMCParticleHelper::MCContributionMap              m_targetMCParticleToHitsMap;            ///< The target mc particle to hits map
        LArMCParticleHelper::PfoContributionMap             m_pfoToHitsMap;                         ///< The pfo to hits map
        LArMCParticleHelper::MCParticleToPfoHitSharingMap   m_mcToPfoHitSharingMap;                 ///< The mc to pfo hit sharing map
        LArMCParticleHelper::MCParticleToPfoHitSharingMap   m_interpretedMCToPfoHitSharingMap;      ///< The interpreted mc to pfo hit sharing map
    };

    /**
     *  @brief  Fill the validation info containers
     *
     *  @param  pMCParticleList the address of the mc particle list
     *  @param  pCaloHitList the address of the calo hit list
     *  @param  pPfoList the address of the pfo list
     *  @param  validationInfo to receive the validation info
     */
    virtual void FillValidationInfo(const pandora::MCParticleList *const pMCParticleList, const pandora::CaloHitList *const pCaloHitList,
        const pandora::PfoList *const pPfoList, ValidationInfo &validationInfo) const = 0;

    /**
     *  @brief  Print matching information in a provided validation info object, and write information to tree if configured to do so
     *
     *  @param  validationInfo the validation info
     *  @param  useInterpretedMatching whether to use the interpreted (rather than raw) matching information
     *  @param  printToScreen whether to print the information to screen
     *  @param  fillTree whether to write the information to tree
     */
    virtual void ProcessOutput(const ValidationInfo &validationInfo, const bool useInterpretedMatching, const bool printToScreen,
        const bool fillTree) const = 0;

    /**
     *  @brief  Apply an interpretative matching procedure to the comprehensive matches in the provided validation info object
     *
     *  @param  validationInfo the validation info
     *  @param  interpretedMCToPfoHitSharingMap the output, interpreted mc particle to pfo hit sharing map
     */
    void InterpretMatching(const ValidationInfo &validationInfo, LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const;

    /**
     *  @brief  Get the strongest pfo match (most matched hits) between an available mc primary and an available pfo
     *
     *  @param  validationInfo the validation info
     *  @param  mcPrimaryVector the mc primary vector
     *  @param  usedPfos the set of previously used pfos
     *  @param  interpretedMCToPfoHitSharingMap the output, interpreted mc particle to pfo hit sharing map
     *
     *  @return whether a strong match was identified
     */
    bool GetStrongestPfoMatch(const ValidationInfo &validationInfo, const pandora::MCParticleVector &mcPrimaryVector, pandora::PfoSet &usedPfos,
        LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const;

    /**
     *  @brief  Get the best matches for any pfos left-over after the strong matching procedure
     *
     *  @param  validationInfo the validation info
     *  @param  mcPrimaryVector the mc primary vector
     *  @param  usedPfos the set of previously used pfos
     *  @param  interpretedMCToPfoHitSharingMap the output, interpreted mc particle to pfo hit sharing map
     */
    void GetRemainingPfoMatches(const ValidationInfo &validationInfo, const pandora::MCParticleVector &mcPrimaryVector, const pandora::PfoSet &usedPfos,
        LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap) const;

    /**
     *  @brief  Whether a provided mc primary and pfo are deemed to be a good match
     *
     *  @param  trueHits the list of true hits
     *  @param  recoHits the list of reco hits
     *  @param  sharedHits the list of shared hits
     *
     *  @return boolean
     */
    bool IsGoodMatch(const pandora::CaloHitList &trueHits, const pandora::CaloHitList &recoHits, const pandora::CaloHitList &sharedHits) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool                    m_selectInputHits;              ///< Whether to use only hits passing mc-based quality (is "reconstructable") checks
    float                   m_minHitSharingFraction;        ///< Minimum fraction of energy deposited by selected primary in a single "good" hit
    float                   m_maxPhotonPropagation;         ///< Maximum distance travelled by photon, downstream of a track, in mc particle hierarchy

    int                     m_fileIdentifier;               ///< The input file identifier
    int                     m_eventNumber;                  ///< The event number

    std::string             m_treeName;                     ///< Name of output tree

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Print all/raw matching information to screen
     *
     *  @param  validationInfo the validation info
     */
    void PrintAllMatches(const ValidationInfo &validationInfo) const;

    /**
     *  @brief  Print interpreted matching information to screen
     *
     *  @param  validationInfo the validation info
     */
    void PrintInterpretedMatches(const ValidationInfo &validationInfo) const;

    /**
     *  @brief  Write interpreted matching information to tree
     *
     *  @param  validationInfo the validation info
     */
    void WriteInterpretedMatches(const ValidationInfo &validationInfo) const;

    std::string             m_caloHitListName;              ///< Name of input calo hit list
    std::string             m_mcParticleListName;           ///< Name of input MC particle list
    std::string             m_pfoListName;                  ///< Name of input Pfo list

    bool                    m_printAllToScreen;             ///< Whether to print all/raw matching details to screen
    bool                    m_printMatchingToScreen;        ///< Whether to print matching output to screen
    bool                    m_writeToTree;                  ///< Whether to write all/raw matching details to tree

    bool                    m_useSmallPrimaries;            ///< Whether to consider matches to mc primaries with fewer than m_matchingMinPrimaryHits
    unsigned int            m_matchingMinSharedHits;        ///< The minimum number of shared hits used in matching scheme
    float                   m_matchingMinCompleteness;      ///< The minimum particle completeness to declare a match
    float                   m_matchingMinPurity;            ///< The minimum particle purity to declare a match

    std::string             m_fileName;                     ///< Name of output file
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArMCParticleHelper::MCContributionMap &EventValidationBaseAlgorithm::ValidationInfo::GetAllMCParticleToHitsMap() const
{
    return m_allMCParticleToHitsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArMCParticleHelper::MCContributionMap &EventValidationBaseAlgorithm::ValidationInfo::GetTargetMCParticleToHitsMap() const
{
    return m_targetMCParticleToHitsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArMCParticleHelper::PfoContributionMap &EventValidationBaseAlgorithm::ValidationInfo::GetPfoToHitsMap() const
{
    return m_pfoToHitsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArMCParticleHelper::MCParticleToPfoHitSharingMap &EventValidationBaseAlgorithm::ValidationInfo::GetMCToPfoHitSharingMap() const
{
    return m_mcToPfoHitSharingMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const LArMCParticleHelper::MCParticleToPfoHitSharingMap &EventValidationBaseAlgorithm::ValidationInfo::GetInterpretedMCToPfoHitSharingMap() const
{
    return m_interpretedMCToPfoHitSharingMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationBaseAlgorithm::ValidationInfo::SetAllMCParticleToHitsMap(const LArMCParticleHelper::MCContributionMap &allMCParticleToHitsMap)
{
    m_allMCParticleToHitsMap = allMCParticleToHitsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationBaseAlgorithm::ValidationInfo::SetTargetMCParticleToHitsMap(const LArMCParticleHelper::MCContributionMap &targetMCParticleToHitsMap)
{
    m_targetMCParticleToHitsMap = targetMCParticleToHitsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationBaseAlgorithm::ValidationInfo::SetPfoToHitsMap(const LArMCParticleHelper::PfoContributionMap &pfoToHitsMap)
{
    m_pfoToHitsMap = pfoToHitsMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationBaseAlgorithm::ValidationInfo::SetMCToPfoHitSharingMap(const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap)
{
    m_mcToPfoHitSharingMap = mcToPfoHitSharingMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationBaseAlgorithm::ValidationInfo::SetInterpretedMCToPfoHitSharingMap(const LArMCParticleHelper::MCParticleToPfoHitSharingMap &interpretedMCToPfoHitSharingMap)
{
    m_interpretedMCToPfoHitSharingMap = interpretedMCToPfoHitSharingMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationBaseAlgorithm::PrintAllMatches(const ValidationInfo &validationInfo) const
{
    return this->ProcessOutput(validationInfo, false, true, false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationBaseAlgorithm::PrintInterpretedMatches(const ValidationInfo &validationInfo) const
{
    return this->ProcessOutput(validationInfo, true, true, false);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void EventValidationBaseAlgorithm::WriteInterpretedMatches(const ValidationInfo &validationInfo) const
{
    return this->ProcessOutput(validationInfo, true, false, true);
}

} // namespace lar_content

#endif // LAR_EVENT_VALIDATION_BASE_ALGORITHM_H
