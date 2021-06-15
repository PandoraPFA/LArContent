/**
 *  @file   larpandoracontent/LArObjects/LArMCParticle.h
 *
 *  @brief  Header file for the lar mc particle class.
 *
 *  $Log: $
 */
#ifndef LAR_MC_PARTICLE_H
#define LAR_MC_PARTICLE_H 1

#include "Objects/MCParticle.h"

#include "Pandora/ObjectCreation.h"
#include "Pandora/PandoraObjectFactories.h"

#include "Persistency/BinaryFileReader.h"
#include "Persistency/BinaryFileWriter.h"
#include "Persistency/XmlFileReader.h"
#include "Persistency/XmlFileWriter.h"

namespace lar_content
{

// Enumeration maps onto G4 process IDs, plus an ID for the incident neutrino
enum MCProcess
{
    MC_PROC_INCIDENT_NU = -1,
    MC_PROC_UNKNOWN,
    MC_PROC_PRIMARY,
    MC_PROC_COMPT,
    MC_PROC_PHOT,
    MC_PROC_ANNIHIL,
    MC_PROC_E_IONI,
    MC_PROC_E_BREM,
    MC_PROC_CONV,
    MC_PROC_MU_IONI,
    MC_PROC_MU_MINUS_CAPTURE_AT_REST,
    MC_PROC_NEUTRON_INELASTIC,
    MC_PROC_N_CAPTURE,
    MC_PROC_HAD_ELASTIC,
    MC_PROC_DECAY,
    MC_PROC_COULOMB_SCAT,
    MC_PROC_UNUSED,
    MC_PROC_MU_BREM,
    MC_PROC_MU_PAIR_PROD,
    MC_PROC_PHOTON_INELASTIC,
    MC_PROC_HAD_IONI,
    MC_PROC_PROTON_INELASTIC,
    MC_PROC_PI_PLUS_INELASTIC,
    MC_PROC_CHIPS_NUCLEAR_CAPTURE_AT_REST,
    MC_PROC_PI_MINUS_INELASTIC
};

// Optional Features for the LArMCParticleFactory
enum LArMCParticleFeature {
    MC_3D_STEP_POSITIONS = 0,
    MC_3D_STEP_MOMENTAS
};

typedef std::set<LArMCParticleFeature> LArMCFeatureSet;

/**
 *  @brief  LAr mc particle parameters
 */
class LArMCParticleParameters : public object_creation::MCParticle::Parameters
{
public:
    pandora::InputInt                          m_nuanceCode;      ///< The nuance code
    pandora::InputInt                          m_process;         ///< The process creating the particle

    std::vector<pandora::InputCartesianVector> m_mcStepPositions; ///< The positions of the geant4 steps
    std::vector<pandora::InputCartesianVector> m_mcStepMomentas;  ///< The momenta of the geant4 steps
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LAr mc particle class
 */
class LArMCParticle : public object_creation::MCParticle::Object
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  parameters the lar mc particle parameters
     */
    LArMCParticle(const LArMCParticleParameters &parameters);

    /**
     *  @brief  Get the nuance code
     *
     *  @return the nuance code
     */
    int GetNuanceCode() const;

    /**
     *  @brief  Get the position of the geant4 steps
     *
     *  @return vector of the positions
     */
    std::vector<pandora::CartesianVector> GetMCStepPositions() const;

    /**
     *  @brief  Get the momenta of the geant4 steps
     *
     *  @return vector of the momenta
     */
    std::vector<pandora::CartesianVector> GetMCStepMomentas() const;

    /**
     *  @brief  Get the process
     *
     *  @return the process
     */
    MCProcess GetProcess() const;

private:
    int                                   m_nuanceCode;      ///< The nuance code
    int                                   m_process;         ///< The process that created the particle

    std::vector<pandora::CartesianVector> m_mcStepPositions; ///< The positions of the geant4 steps
    std::vector<pandora::CartesianVector> m_mcStepMomentas;  ///< The momenta of the geant4 steps
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArMCParticleFactory responsible for object creation
 */
class LArMCParticleFactory : public pandora::ObjectFactory<object_creation::MCParticle::Parameters, object_creation::MCParticle::Object>
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  version the LArMCParticle version
     */
    LArMCParticleFactory(const unsigned int version = 2);

    /**
     *  @brief  Constructor
     *
     *  @param  features optional features to load/write.
     *  @param  version the LArMCParticle version
     */
    LArMCParticleFactory(const LArMCFeatureSet &features, const unsigned int version = 2);

    /**
     *  @brief  Create new parameters instance on the heap (memory-management to be controlled by user)
     *
     *  @return the address of the new parameters instance
     */
    Parameters *NewParameters() const;

    /**
     *  @brief  Read any additional (derived class only) object parameters from file using the specified file reader
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  fileReader the file reader, used to extract any additional parameters from file
     */
    pandora::StatusCode Read(Parameters &parameters, pandora::FileReader &fileReader) const;

    /**
     *  @brief  Persist any additional (derived class only) object parameters using the specified file writer
     *
     *  @param  pObject the address of the object to persist
     *  @param  fileWriter the file writer
     */
    pandora::StatusCode Write(const Object *const pObject, pandora::FileWriter &fileWriter) const;

    /**
     *  @brief  Create an object with the given parameters
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  pObject to receive the address of the object created
     */
    pandora::StatusCode Create(const Parameters &parameters, const Object *&pObject) const;

    /**
     *  @brief  Parse string vector to optional features set
     *
     *  @param  features the optional features, as strings
     */
    static pandora::StatusCode ParseFeatures(const pandora::StringVector &parameters, LArMCFeatureSet &optionalFeatures);

private:

    unsigned int m_version; ///< The LArMCParticle version
    LArMCFeatureSet m_optionalFeatures; ///< Optional features to add to the LArMCParticle
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArMCParticle::LArMCParticle(const LArMCParticleParameters &parameters) :
    object_creation::MCParticle::Object(parameters),
    m_nuanceCode(parameters.m_nuanceCode.Get()),
    m_process(parameters.m_process.Get())
{
    for (auto const &mcStepPosition : parameters.m_mcStepPositions)
        m_mcStepPositions.emplace_back(mcStepPosition.Get());

    for (auto const &mcStepMomenta : parameters.m_mcStepMomentas)
        m_mcStepMomentas.emplace_back(mcStepMomenta.Get());
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline int LArMCParticle::GetNuanceCode() const
{
    return m_nuanceCode;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<pandora::CartesianVector> LArMCParticle::GetMCStepPositions() const
{
    return m_mcStepPositions;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::vector<pandora::CartesianVector> LArMCParticle::GetMCStepMomentas() const
{
    return m_mcStepMomentas;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline MCProcess LArMCParticle::GetProcess() const
{
    return MCProcess(m_process);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArMCParticleFactory::LArMCParticleFactory(const unsigned int version) : m_version(version)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArMCParticleFactory::LArMCParticleFactory(const LArMCFeatureSet &features, const unsigned int version) :
    m_version(version),
    m_optionalFeatures(features)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArMCParticleFactory::Parameters *LArMCParticleFactory::NewParameters() const
{
    return (new LArMCParticleParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArMCParticleFactory::ParseFeatures(const pandora::StringVector &parameters, LArMCFeatureSet &optionalFeatures)
{
    for (auto feature : parameters)
    {
        if (feature == "3D_Step_Positions")
        {
            optionalFeatures.insert(LArMCParticleFeature::MC_3D_STEP_POSITIONS);
        }
        else if (feature == "3D_Step_Momentas")
        {
            optionalFeatures.insert(LArMCParticleFeature::MC_3D_STEP_MOMENTAS);
        }
        else
        {
            std::cout << "Error: Invalid LArMCParticleFactory feature given: " << feature << std::endl;
            return pandora::STATUS_CODE_INVALID_PARAMETER;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArMCParticleFactory::Create(const Parameters &parameters, const Object *&pObject) const
{
    const LArMCParticleParameters &larMCParticleParameters(dynamic_cast<const LArMCParticleParameters &>(parameters));
    pObject = new LArMCParticle(larMCParticleParameters);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArMCParticleFactory::Read(Parameters &parameters, pandora::FileReader &fileReader) const
{
    // ATTN: To receive this call-back must have already set file reader mc particle factory to this factory
    int nuanceCode(0);
    int process(0);

    unsigned int nMCStepPositions(0);
    std::vector<pandora::InputCartesianVector> mcStepPositions;

    unsigned int nMCStepMomentas(0);
    std::vector<pandora::InputCartesianVector> mcStepMomentas;

    if (pandora::BINARY == fileReader.GetFileType())
    {
        pandora::BinaryFileReader &binaryFileReader(dynamic_cast<pandora::BinaryFileReader &>(fileReader));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(nuanceCode));

        if (m_version > 1)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(process));

        if (m_optionalFeatures.count(LArMCParticleFeature::MC_3D_STEP_POSITIONS) > 0)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(nMCStepPositions));
            for (unsigned int step = 0; step < nMCStepPositions; ++step)
            {
                pandora::CartesianVector mcStepPosition(0.f, 0.f, 0.f);
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(mcStepPosition));
                mcStepPositions.emplace_back(mcStepPosition);
            }
        }

        if (m_optionalFeatures.count(LArMCParticleFeature::MC_3D_STEP_MOMENTAS) > 0)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(nMCStepMomentas));
            for (unsigned int step = 0; step < nMCStepMomentas; ++step)
            {
                pandora::CartesianVector mcStepMomenta(0.f, 0.f, 0.f);
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(mcStepMomenta));
                mcStepMomentas.emplace_back(mcStepMomenta);
            }
        }
    }
    else if (pandora::XML == fileReader.GetFileType())
    {
        pandora::XmlFileReader &xmlFileReader(dynamic_cast<pandora::XmlFileReader &>(fileReader));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("NuanceCode", nuanceCode));

        if (m_version > 1)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("Process", process));

        if (m_optionalFeatures.count(LArMCParticleFeature::MC_3D_STEP_POSITIONS) > 0)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("NumberOfMCStepPositions", nMCStepPositions));
            for (unsigned int step = 0; step < nMCStepPositions; ++step)
            {
                pandora::CartesianVector mcStepPosition(0.f, 0.f, 0.f);
                std::stringstream mcStepPositionName;
                mcStepPositionName << "MCStepPosition" << step;
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable(mcStepPositionName.str(), mcStepPosition));
                mcStepPositions.emplace_back(mcStepPosition);
            }
        }

        if (m_optionalFeatures.count(LArMCParticleFeature::MC_3D_STEP_MOMENTAS) > 0)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("NumberOfMCStepMomentas", nMCStepMomentas));
            for (unsigned int step = 0; step < nMCStepMomentas; ++step)
            {
                pandora::CartesianVector mcStepMomenta(0.f, 0.f, 0.f);
                std::stringstream mcStepMomentaName;
                mcStepMomentaName << "MCStepMomenta" << step;
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable(mcStepMomentaName.str(), mcStepMomenta));
                mcStepMomentas.emplace_back(mcStepMomenta);
            }
        }
    }
    else
    {
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    LArMCParticleParameters &larMCParticleParameters(dynamic_cast<LArMCParticleParameters &>(parameters));
    larMCParticleParameters.m_nuanceCode = nuanceCode;
    larMCParticleParameters.m_process = process;

    larMCParticleParameters.m_mcStepPositions = mcStepPositions;
    larMCParticleParameters.m_mcStepMomentas = mcStepMomentas;

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArMCParticleFactory::Write(const Object *const pObject, pandora::FileWriter &fileWriter) const
{
    // ATTN: To receive this call-back must have already set file writer mc particle factory to this factory
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle *>(pObject));

    if (!pLArMCParticle)
        return pandora::STATUS_CODE_INVALID_PARAMETER;

    const std::vector<pandora::CartesianVector> &mcStepPositions(pLArMCParticle->GetMCStepPositions());
    const int nMCStepPositions(mcStepPositions.size());

    const std::vector<pandora::CartesianVector> &mcStepMomentas(pLArMCParticle->GetMCStepMomentas());
    const int nMCStepMomentas(mcStepMomentas.size());

    if (pandora::BINARY == fileWriter.GetFileType())
    {
        pandora::BinaryFileWriter &binaryFileWriter(dynamic_cast<pandora::BinaryFileWriter &>(fileWriter));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArMCParticle->GetNuanceCode()));

        if (m_version > 1)
            PANDORA_RETURN_RESULT_IF(
                pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(static_cast<int>(pLArMCParticle->GetProcess())));

        if (m_optionalFeatures.count(LArMCParticleFeature::MC_3D_STEP_POSITIONS) > 0)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(nMCStepPositions));
            for (auto const &mcStepPosition : mcStepPositions)
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(mcStepPosition));
        }

        if (m_optionalFeatures.count(LArMCParticleFeature::MC_3D_STEP_MOMENTAS) > 0)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(nMCStepMomentas));
            for (auto const &mcStepMomenta : mcStepMomentas)
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(mcStepMomenta));
        }

    }
    else if (pandora::XML == fileWriter.GetFileType())
    {
        pandora::XmlFileWriter &xmlFileWriter(dynamic_cast<pandora::XmlFileWriter &>(fileWriter));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("NuanceCode", pLArMCParticle->GetNuanceCode()));

        if (m_version > 1)
            PANDORA_RETURN_RESULT_IF(
                pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("Process", static_cast<int>(pLArMCParticle->GetProcess())));

        if (m_optionalFeatures.count(LArMCParticleFeature::MC_3D_STEP_POSITIONS) > 0)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("NumberOfMCStepPositions", nMCStepPositions));
            for (int step = 0; step < nMCStepPositions; ++step)
            {
                const pandora::CartesianVector &position(mcStepPositions[step]);
                std::stringstream mcStepPositionName;
                mcStepPositionName << "MCStepPosition" << step;
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable(mcStepPositionName.str(), position));
            }
        }

        if (m_optionalFeatures.count(LArMCParticleFeature::MC_3D_STEP_MOMENTAS) > 0)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("NumberOfMCStepMomentas", nMCStepMomentas));
            for (int step = 0; step < nMCStepMomentas; ++step)
            {
                const pandora::CartesianVector &momenta(mcStepMomentas[step]);
                std::stringstream mcStepMomentaName;
                mcStepMomentaName << "MCStepMomenta" << step;
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable(mcStepMomentaName.str(), momenta));
            }
        }
    }
    else
    {
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lar_content

#endif // #ifndef LAR_MC_PARTICLE_H
