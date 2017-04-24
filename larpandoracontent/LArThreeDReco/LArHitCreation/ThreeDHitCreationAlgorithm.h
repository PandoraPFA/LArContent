/**
 *  @file   larpandoracontent/LArThreeDReco/LArHitCreation/ThreeDHitCreationAlgorithm.h
 *
 *  @brief  Header file for the three dimensional hit creation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_HIT_CREATION_ALGORITHM_H
#define LAR_THREE_D_HIT_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include <vector>

namespace lar_content
{

class HitCreationBaseTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ThreeDHitCreationAlgorithm::Algorithm class
 */
class ThreeDHitCreationAlgorithm : public pandora::Algorithm
{
public:

    /**
     *  @brief  Trajectory samples record the results of sampling a particles in a particular view
     */
    class TrajectorySample
    {
    public:
        /**
         *  @brief  Constructor
         */
        TrajectorySample(const pandora::CartesianVector &position, const pandora::HitType hitType, const double sigma);

        /**
         *  @brief  Get the sampling position
         * 
         *  @return the sampling position
         */
        const pandora::CartesianVector &GetPosition() const;

        /**
         *  @brief  Get the sampling hit type
         * 
         *  @return the sampling hit type
         */
        pandora::HitType GetHitType() const;

        /**
         *  @brief  Get the sampling sigma
         * 
         *  @return the sampling sigma
         */
        double GetSigma() const;

    private:
        pandora::CartesianVector    m_position;             ///< The sampling position
        pandora::HitType            m_hitType;              ///< The sampling hit type
        double                      m_sigma;                ///< The sampling sigma
    };

    typedef std::vector<TrajectorySample> TrajectorySampleVector;

    /**
     *  @brief  Proto hits are temporary constructs to be used during iterative 3D hit procedure
     */
    class ProtoHit
    {
    public:
        /**
         *  @brief  Constructor
         * 
         *  @param  pParentCaloHit2D the address of the parent 2D calo hit
         */
        ProtoHit(const pandora::CaloHit *const pParentCaloHit2D);

        /**
         *  @brief  Get the address of the parent 2D calo hit
         * 
         *  @return the address of the parent 2D calo hit
         */
        const pandora::CaloHit *GetParentCaloHit2D() const;

        /**
         *  @brief  Whether the proto hit position is set
         * 
         *  @return boolean
         */
        bool IsPositionSet() const;

        /**
         *  @brief  Get the output 3D position
         * 
         *  @return the output 3D position, if set
         * 
         *  @throws StatusCodeException
         */
        const pandora::CartesianVector &GetPosition3D() const;

        /**
         *  @brief  Get the chi squared value
         * 
         *  @return the chi squared value, if set
         * 
         *  @throws StatusCodeException
         */
        double GetChi2() const;

        /**
         *  @brief  Get the number of trajectory samples
         * 
         *  @return the number of trajectory samples
         */
        unsigned int GetNTrajectorySamples() const;

        /**
         *  @brief  Get the first trajectory sample
         * 
         *  @return the first trajectory sample, if at least one sample is present
         * 
         *  @throws StatusCodeException
         */
        const TrajectorySample &GetFirstTrajectorySample() const;

        /**
         *  @brief  Get the last trajectory sample
         * 
         *  @return the last trajectory sample, if at least two samples are present
         * 
         *  @throws StatusCodeException
         */
        const TrajectorySample &GetLastTrajectorySample() const;

        /**
         *  @brief  Set position 3D
         * 
         *  @param  the output 3D position
         *  @param  the output chi squared value
         */
        void SetPosition3D(const pandora::CartesianVector &position3D, const double chi2);

        /**
         *  @brief  Add a trajectory sample
         * 
         *  @param  the trajectory sample
         */
        void AddTrajectorySample(const TrajectorySample &trajectorySample);

    private:
        const pandora::CaloHit     *m_pParentCaloHit2D;         ///< The address of the parent 2D calo hit
        bool                        m_isPositionSet;            ///< Whether the output 3D position has been set
        pandora::CartesianVector    m_position3D;               ///< The output 3D position
        double                      m_chi2;                     ///< The output chi squared value
        TrajectorySampleVector      m_trajectorySampleVector;   ///< The trajectory sample vector
    };

    typedef std::vector<ProtoHit> ProtoHitVector;

    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

    /**
     *  @brief  Get the subset of a provided calo hit vector corresponding to a specified hit type
     *
     *  @param  inputCaloHitVector the input calo hit vector
     *  @param  hitType the hit type to filter upon
     *  @param  outputCaloHitVector to receive the output calo hit vector
     */
    void FilterCaloHitsByType(const pandora::CaloHitVector &inputCaloHitVector, const pandora::HitType hitType, pandora::CaloHitVector &outputCaloHitVector) const;

private:
    pandora::StatusCode Run();

    /**
     *  @brief  Get the list of 2D calo hits in a pfo for which 3D hits have and have not been created
     *
     *  @param  pPfo the address of the pfo
     *  @param  protoHitVector the vector of proto hits, describing current state of 2D hit construction
     *  @param  remainingHitVector to receive the vector of two dimensional calo hits for which three dimensional hits have not been created
     */
    void SeparateTwoDHits(const pandora::ParticleFlowObject *const pPfo, const ProtoHitVector &protoHitVector, pandora::CaloHitVector &remainingHitVector) const;

    /**
     *  @brief  Create new three dimensional hits from two dimensional hits
     *
     *  @param  protoHitVector the input proto hit vector
     *  @param  newThreeDHits to receive the addresses of the new three dimensional calo hits
     *  @param  allNewThreeDHits to record the addresses of all new three dimensional calo hits
     */
    void CreateThreeDHits(const ProtoHitVector &protoHitVector, pandora::CaloHitList &newThreeDHits, pandora::CaloHitList &allNewThreeDHits) const;

    /**
     *  @brief  Create a new three dimensional hit from a two dimensional hit
     *
     *  @param  protoHit the proto hit containing all required information
     *  @param  pCaloHit3D to receive the address of the new three dimensional calo hit
     */
    void CreateThreeDHit(const ProtoHit &protoHit, const pandora::CaloHit *&pCaloHit3D) const;

    /**
     *  @brief  Check that a new three dimensional position is not unphysical
     *
     *  @param  protoHit the proto hit
     *
     *  @param  boolean  
     */
    bool CheckThreeDHit(const ProtoHit &protoHit) const;

    /**
     *  @brief  Add a specified list of three dimensional hits to a cluster in a pfo, creating the new cluster if required
     *
     *  @param  pPfo the address of the pfo
     *  @param  caloHitList the list of three dimensional hits
     */
    void AddThreeDHitsToPfo(const pandora::ParticleFlowObject *const pPfo, const pandora::CaloHitList &caloHitList) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<HitCreationBaseTool*> HitCreationToolVector;
    HitCreationToolVector   m_algorithmToolVector;      ///< The algorithm tool vector

    std::string             m_inputPfoListName;         ///< The name of the input pfo list
    std::string             m_outputCaloHitListName;    ///< The name of the output calo hit list
    std::string             m_outputClusterListName;    ///< The name of the output cluster list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ThreeDHitCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new ThreeDHitCreationAlgorithm();
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ThreeDHitCreationAlgorithm::TrajectorySample::TrajectorySample(const pandora::CartesianVector &position, const pandora::HitType hitType, const double sigma) :
    m_position(position),
    m_hitType(hitType),
    m_sigma(sigma)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ThreeDHitCreationAlgorithm::TrajectorySample::GetPosition() const
{
    return m_position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::HitType ThreeDHitCreationAlgorithm::TrajectorySample::GetHitType() const
{
    return m_hitType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double ThreeDHitCreationAlgorithm::TrajectorySample::GetSigma() const
{
    return m_sigma;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline ThreeDHitCreationAlgorithm::ProtoHit::ProtoHit(const pandora::CaloHit *const pParentCaloHit2D) :
    m_pParentCaloHit2D(pParentCaloHit2D),
    m_isPositionSet(false),
    m_position3D(0.f, 0.f, 0.f),
    m_chi2(std::numeric_limits<double>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHit *ThreeDHitCreationAlgorithm::ProtoHit::GetParentCaloHit2D() const
{
    return m_pParentCaloHit2D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool ThreeDHitCreationAlgorithm::ProtoHit::IsPositionSet() const
{
    return m_isPositionSet;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int ThreeDHitCreationAlgorithm::ProtoHit::GetNTrajectorySamples() const
{
    return m_trajectorySampleVector.size();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ThreeDHitCreationAlgorithm::ProtoHit::SetPosition3D(const pandora::CartesianVector &position3D, const double chi2)
{
    m_position3D = position3D;
    m_chi2 = chi2;
    m_isPositionSet = true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ThreeDHitCreationAlgorithm::ProtoHit::AddTrajectorySample(const TrajectorySample &trajectorySample)
{
    m_trajectorySampleVector.push_back(trajectorySample);
}

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_HIT_CREATION_ALGORITHM_H
