/**
 *  @file   larpandoracontent/LArObjects/LArCaloHit.h
 * 
 *  @brief  Header file for the lar calo hit class.
 * 
 *  $Log: $
 */
#ifndef LAR_CALO_HIT_H
#define LAR_CALO_HIT_H 1

#include "Objects/CaloHit.h"

#include "Pandora/ObjectCreation.h"
#include "Pandora/PandoraObjectFactories.h"

#include "Persistency/BinaryFileReader.h"
#include "Persistency/BinaryFileWriter.h"
#include "Persistency/XmlFileReader.h"
#include "Persistency/XmlFileWriter.h"

namespace lar_content
{

/**
 *  @brief  LAr calo hit parameters
 */
class LArCaloHitParameters : public object_creation::CaloHit::Parameters
{
public:
    pandora::InputFloat     m_showerProbability;    ///< The probability for electromagnetic shower hypothesis
    pandora::InputFloat     m_trackProbability;     ///< The probability for track hypothesis
    pandora::InputFloat     m_michelProbability;    ///< The probability for michel hypothesis
    pandora::InputFloat     m_otherProbability;     ///< The probability for other hypotheses
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LAr calo hit class
 */
class LArCaloHit : public object_creation::CaloHit::Object
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  parameters the lar calo hit parameters
     */
    LArCaloHit(const LArCaloHitParameters &parameters);

    /**
     *  @brief  Get the probability for electromagnetic shower hypothesis
     *
     *  @return the probability for electromagnetic shower hypothesis
     */
    float GetShowerProbability() const;

    /**
     *  @brief  Get the probability for track hypothesis
     *
     *  @return the probability for track hypothesis
     */
    float GetTrackProbability() const;

    /**
     *  @brief  Get the probability for michel hypothesis
     *
     *  @return the probability for michel hypothesis
     */
    float GetMichelProbability() const;

    /**
     *  @brief  Get the probability for other hypotheses
     * 
     *  @return the probability for other hypotheses
     */
    float GetOtherProbability() const;

private:
    float                   m_showerProbability;    ///< The probability for electromagnetic shower hypothesis
    float                   m_trackProbability;     ///< The probability for track hypothesis
    float                   m_michelProbability;    ///< The probability for michel hypothesis
    float                   m_otherProbability;     ///< The probability for other hypotheses
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArCaloHitFactory responsible for object creation
 */
class LArCaloHitFactory : public pandora::ObjectFactory<object_creation::CaloHit::Parameters, object_creation::CaloHit::Object>
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

inline LArCaloHit::LArCaloHit(const LArCaloHitParameters &parameters) :
    object_creation::CaloHit::Object(parameters),
    m_showerProbability(parameters.m_showerProbability.Get()),
    m_trackProbability(parameters.m_trackProbability.Get()),
    m_michelProbability(parameters.m_michelProbability.Get()),
    m_otherProbability(parameters.m_otherProbability.Get())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetShowerProbability() const
{
    return m_showerProbability;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetTrackProbability() const
{
    return m_trackProbability;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetMichelProbability() const
{
    return m_michelProbability;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetOtherProbability() const
{
    return m_otherProbability;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArCaloHitFactory::Parameters *LArCaloHitFactory::NewParameters() const
{
    return (new LArCaloHitParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArCaloHitFactory::Create(const Parameters &parameters, const Object *&pObject) const
{
    const LArCaloHitParameters &larCaloHitParameters(dynamic_cast<const LArCaloHitParameters&>(parameters));
    pObject = new LArCaloHit(larCaloHitParameters);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArCaloHitFactory::Read(Parameters &parameters, pandora::FileReader &fileReader) const
{
    // ATTN: To receive this call-back must have already set file reader mc particle factory to this factory
    float showerProbability(0.f);
    float trackProbability(0.f);
    float michelProbability(0.f);
    float otherProbability(0.f);

    if (pandora::BINARY == fileReader.GetFileType())
    {
        pandora::BinaryFileReader &binaryFileReader(dynamic_cast<pandora::BinaryFileReader&>(fileReader));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(showerProbability));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(trackProbability));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(michelProbability));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(otherProbability));
    }
    else if (pandora::XML == fileReader.GetFileType())
    {
        pandora::XmlFileReader &xmlFileReader(dynamic_cast<pandora::XmlFileReader&>(fileReader));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("ShowerProbability", showerProbability));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("TrackProbability", trackProbability));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("MichelProbability", michelProbability));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("OtherProbability", otherProbability));
    }
    else
    {
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    LArCaloHitParameters &larCaloHitParameters(dynamic_cast<LArCaloHitParameters&>(parameters));
    larCaloHitParameters.m_showerProbability = showerProbability;
    larCaloHitParameters.m_trackProbability = trackProbability;
    larCaloHitParameters.m_michelProbability = michelProbability;
    larCaloHitParameters.m_otherProbability = otherProbability;

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArCaloHitFactory::Write(const Object *const pObject, pandora::FileWriter &fileWriter) const
{
    // ATTN: To receive this call-back must have already set file writer mc particle factory to this factory
    const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit*>(pObject));

    if (!pLArCaloHit)
        return pandora::STATUS_CODE_INVALID_PARAMETER;

    if (pandora::BINARY == fileWriter.GetFileType())
    {
        pandora::BinaryFileWriter &binaryFileWriter(dynamic_cast<pandora::BinaryFileWriter&>(fileWriter));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetShowerProbability()));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetTrackProbability()));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetMichelProbability()));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetOtherProbability()));
    }
    else if (pandora::XML == fileWriter.GetFileType())
    {
        pandora::XmlFileWriter &xmlFileWriter(dynamic_cast<pandora::XmlFileWriter&>(fileWriter));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("ShowerProbability", pLArCaloHit->GetShowerProbability()));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("TrackProbability", pLArCaloHit->GetTrackProbability()));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("MichelProbability", pLArCaloHit->GetMichelProbability()));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("OtherProbability", pLArCaloHit->GetOtherProbability()));
    }
    else
    {
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lar_content

#endif // #ifndef LAR_CALO_HIT_H
