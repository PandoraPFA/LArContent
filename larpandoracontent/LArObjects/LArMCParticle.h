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

/**
 *  @brief  LAr mc particle parameters
 */
class LArMCParticleParameters : public object_creation::MCParticle::Parameters
{
public:
    pandora::InputInt                          m_nuanceCode;        ///< The nuance code
    std::vector<pandora::InputCartesianVector> m_mcStepPositions;    ///< The positions of the geant4 steps
    std::vector<pandora::InputCartesianVector> m_mcStepMomentas;    ///< The momenta of the geant4 steps
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

private:
    int                                   m_nuanceCode;       ///< The nuance code
    std::vector<pandora::CartesianVector> m_mcStepPositions;  ///< The positions of the geant4 steps
    std::vector<pandora::CartesianVector> m_mcStepMomentas;   ///< The momenta of the geant4 steps
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArMCParticleFactory responsible for object creation
 */
class LArMCParticleFactory : public pandora::ObjectFactory<object_creation::MCParticle::Parameters, object_creation::MCParticle::Object>
{
public:
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
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArMCParticle::LArMCParticle(const LArMCParticleParameters &parameters) :
    object_creation::MCParticle::Object(parameters),
    m_nuanceCode(parameters.m_nuanceCode.Get())
{
    for (auto const &mcStepPosition : parameters.m_mcStepPositions)
    {
        m_mcStepPositions.push_back(mcStepPosition.Get());
    }

    for (auto const &mcStepMomenta : parameters.m_mcStepMomentas)
    {
        m_mcStepMomentas.push_back(mcStepMomenta.Get());
    }
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
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArMCParticleFactory::Parameters *LArMCParticleFactory::NewParameters() const
{
    return (new LArMCParticleParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArMCParticleFactory::Create(const Parameters &parameters, const Object *&pObject) const
{
    const LArMCParticleParameters &larMCParticleParameters(dynamic_cast<const LArMCParticleParameters&>(parameters));
    pObject = new LArMCParticle(larMCParticleParameters);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArMCParticleFactory::Read(Parameters &parameters, pandora::FileReader &fileReader) const
{
    // ATTN: To receive this call-back must have already set file reader mc particle factory to this factory
    int nuanceCode(0);

    int nMCStepPositions(0);
    std::vector<pandora::InputCartesianVector> mcStepPositions;

    int nMCStepMomentas(0);
    std::vector<pandora::InputCartesianVector> mcStepMomentas;

    if (pandora::BINARY == fileReader.GetFileType())
    {
        pandora::BinaryFileReader &binaryFileReader(dynamic_cast<pandora::BinaryFileReader&>(fileReader));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(nuanceCode));

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(nMCStepPositions));
        for (int i = 0; i < nMCStepPositions; ++i)
        {   
            pandora::CartesianVector mcStepPosition(0.0,0.0,0.0);
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(mcStepPosition));
            mcStepPositions.push_back(pandora::InputCartesianVector(mcStepPosition));
        }

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(nMCStepMomentas));
        for (int i = 0; i < nMCStepMomentas; ++i)
        {
            pandora::CartesianVector mcStepMomenta(0.0,0.0,0.0);
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(mcStepMomenta));
            mcStepMomentas.push_back(pandora::InputCartesianVector(mcStepMomenta));
        }
    }
    else if (pandora::XML == fileReader.GetFileType())
    {
        pandora::XmlFileReader &xmlFileReader(dynamic_cast<pandora::XmlFileReader&>(fileReader));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("NuanceCode", nuanceCode));

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("NumberOfMCStepPositions", nMCStepPositions));
        for (int i = 0; i < nMCStepPositions; ++i)
        {
            pandora::CartesianVector mcStepPosition(0.0,0.0,0.0);
            std::stringstream mcStepPositionName;
            mcStepPositionName << "MCStepPosition" << i;
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable(mcStepPositionName.str(), mcStepPosition));
            mcStepPositions.push_back(pandora::InputCartesianVector(mcStepPosition));
        }

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("NumberOfMCStepMomentas", nMCStepMomentas));
        for (int i = 0; i < nMCStepMomentas; ++i)
        {
            pandora::CartesianVector mcStepMomenta(0.0,0.0,0.0); 
            std::stringstream mcStepMomentaName;
            mcStepMomentaName << "MCStepMomenta" << i;
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable(mcStepMomentaName.str(), mcStepMomenta));
            mcStepMomentas.push_back(pandora::InputCartesianVector(mcStepMomenta));
        }
    }
    else
    {
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    LArMCParticleParameters &larMCParticleParameters(dynamic_cast<LArMCParticleParameters&>(parameters));
    larMCParticleParameters.m_nuanceCode = nuanceCode;
    larMCParticleParameters.m_mcStepPositions = mcStepPositions;
    larMCParticleParameters.m_mcStepMomentas = mcStepMomentas;

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArMCParticleFactory::Write(const Object *const pObject, pandora::FileWriter &fileWriter) const
{
    // ATTN: To receive this call-back must have already set file writer mc particle factory to this factory
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle*>(pObject));

    if (!pLArMCParticle)
        return pandora::STATUS_CODE_INVALID_PARAMETER;

    const std::vector<pandora::CartesianVector> &mcStepPositions = pLArMCParticle->GetMCStepPositions();
    int nMCStepPositions = mcStepPositions.size();
    std::cout << "LArMCParticleFactory::Write - nMCStepPositions = " << nMCStepPositions << std::endl;

    const std::vector<pandora::CartesianVector> &mcStepMomentas = pLArMCParticle->GetMCStepMomentas();
    int nMCStepMomentas = mcStepMomentas.size();
    std::cout << "LArMCParticleFactory::Write - nMCStepMomentas = " << nMCStepMomentas << std::endl;

    if (pandora::BINARY == fileWriter.GetFileType())
    {
        pandora::BinaryFileWriter &binaryFileWriter(dynamic_cast<pandora::BinaryFileWriter&>(fileWriter));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArMCParticle->GetNuanceCode()));

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(nMCStepPositions)); 
        for (auto const &mcStepPosition : mcStepPositions)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(mcStepPosition));
        }

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(nMCStepMomentas));
        for (auto const &mcStepMomenta : mcStepMomentas)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(mcStepMomenta));
        }
    }
    else if (pandora::XML == fileWriter.GetFileType())
    {
        pandora::XmlFileWriter &xmlFileWriter(dynamic_cast<pandora::XmlFileWriter&>(fileWriter));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("NuanceCode", pLArMCParticle->GetNuanceCode()));

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("NumberOfMCStepPositions", nMCStepPositions));
        int mcStepPositionCounter(0);
        for (auto const &mcStepPosition : mcStepPositions)
        {
            std::stringstream mcStepPositionName;
            mcStepPositionName << "MCStepPosition" << mcStepPositionCounter;
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable(mcStepPositionName.str(), mcStepPosition));
            ++mcStepPositionCounter;
        }

        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("NumberOfMCStepMomentas", nMCStepMomentas));
        int mcStepMomentasCounter(0);
        for (auto const &mcStepMomenta : mcStepMomentas)
        {
            std::stringstream mcStepMomentaName;
            mcStepMomentaName << "MCStepMomenta" << mcStepMomentasCounter;
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable(mcStepMomentaName.str(), mcStepMomenta));
            ++mcStepMomentasCounter;
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
