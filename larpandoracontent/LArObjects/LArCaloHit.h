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
    pandora::InputUInt m_larTPCVolumeId;   ///< The lar tpc volume id
    pandora::InputUInt m_daughterVolumeId; ///< The daughter volume id

    pandora::InputFloat m_pTrack;          ///< Probability the hit is track-like
    pandora::InputFloat m_pShower;         ///< Probability the hit is shower-like
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
     *  @brief  Get the lar tpc volume id
     *
     *  @return the lar tpc volume id
     */
    unsigned int GetLArTPCVolumeId() const;

    /**
     *  @brief  Get the daughter volume id
     *
     *  @return the daughter volume id
     */
    unsigned int GetDaughterVolumeId() const;

    /**
     *  @brief  Fill the parameters associated with this calo hit
     *
     *  @param  parameters the output parameters
     */
    void FillParameters(LArCaloHitParameters &parameters) const;

    /**
     *  @brief  Get the probability that the hit is track-like
     *
     *  @return the probability the hit is track-like
     */
    float GetTrackProbability() const;

    /**
     *  @brief  Get the probability that the hit is shower-like
     *
     *  @return the probability the hit is shower-like
     */
    float GetShowerProbability() const;

    /**
     *  @brief  Set the probability that the hit is track-like
     *
     *  @param  probability the probability the hit is track-like
     */
    void SetTrackProbability(const float probability);

    /**
     *  @brief  Set the probability that the hit is shower-like
     *
     *  @param  probability the probability the hit is shower-like
     */
    void SetShowerProbability(const float probability);

private:
    unsigned int m_larTPCVolumeId;   ///< The lar tpc volume id
    unsigned int m_daughterVolumeId; ///< The daughter volume id
    pandora::InputFloat m_pTrack;    ///< The probability that the hit is track-like
    pandora::InputFloat m_pShower;   ///< The probability that the hit is shower-like
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArCaloHitFactory responsible for object creation
 */
class LArCaloHitFactory : public pandora::ObjectFactory<object_creation::CaloHit::Parameters, object_creation::CaloHit::Object>
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  version the LArCaloHit version
     */
    LArCaloHitFactory(const unsigned int version = 1);

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

private:
    unsigned int m_version; ///< The LArCaloHit version
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArCaloHit::LArCaloHit(const LArCaloHitParameters &parameters) :
    object_creation::CaloHit::Object(parameters),
    m_larTPCVolumeId(parameters.m_larTPCVolumeId.Get()),
    m_daughterVolumeId(parameters.m_daughterVolumeId.IsInitialized() ? parameters.m_daughterVolumeId.Get() : 0),
    m_pTrack(parameters.m_pTrack.IsInitialized() ? parameters.m_pTrack.Get() : -1.),
    m_pShower(parameters.m_pShower.IsInitialized() ? parameters.m_pShower.Get() : -1.)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArCaloHit::GetLArTPCVolumeId() const
{
    return m_larTPCVolumeId;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArCaloHit::GetDaughterVolumeId() const
{
    return m_daughterVolumeId;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArCaloHit::FillParameters(LArCaloHitParameters &parameters) const
{
    parameters.m_positionVector = this->GetPositionVector();
    parameters.m_expectedDirection = this->GetExpectedDirection();
    parameters.m_cellNormalVector = this->GetCellNormalVector();
    parameters.m_cellGeometry = this->GetCellGeometry();
    parameters.m_cellSize0 = this->GetCellSize0();
    parameters.m_cellSize1 = this->GetCellSize1();
    parameters.m_cellThickness = this->GetCellThickness();
    parameters.m_nCellRadiationLengths = this->GetNCellRadiationLengths();
    parameters.m_nCellInteractionLengths = this->GetNCellInteractionLengths();
    parameters.m_time = this->GetTime();
    parameters.m_inputEnergy = this->GetInputEnergy();
    parameters.m_mipEquivalentEnergy = this->GetMipEquivalentEnergy();
    parameters.m_electromagneticEnergy = this->GetElectromagneticEnergy();
    parameters.m_hadronicEnergy = this->GetHadronicEnergy();
    parameters.m_isDigital = this->IsDigital();
    parameters.m_hitType = this->GetHitType();
    parameters.m_hitRegion = this->GetHitRegion();
    parameters.m_layer = this->GetLayer();
    parameters.m_isInOuterSamplingLayer = this->IsInOuterSamplingLayer();
    // ATTN Set the parent address to the original owner of the calo hit
    parameters.m_pParentAddress = static_cast<const void *>(this);
    parameters.m_larTPCVolumeId = this->GetLArTPCVolumeId();
    parameters.m_daughterVolumeId = this->GetDaughterVolumeId();
    parameters.m_pTrack = this->GetTrackProbability();
    parameters.m_pShower = this->GetShowerProbability();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetTrackProbability() const
{
    return m_pTrack.Get();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float LArCaloHit::GetShowerProbability() const
{
    return m_pShower.Get();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArCaloHit::SetTrackProbability(const float probability)
{
    if (probability >= 0.f)
        m_pTrack = probability;
    else
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArCaloHit::SetShowerProbability(const float probability)
{
    if (probability >= 0.f)
        m_pShower = probability;
    else
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LArCaloHitFactory::LArCaloHitFactory(const unsigned int version) :
    m_version(version)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArCaloHitFactory::Parameters *LArCaloHitFactory::NewParameters() const
{
    return (new LArCaloHitParameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArCaloHitFactory::Create(const Parameters &parameters, const Object *&pObject) const
{
    const LArCaloHitParameters &larCaloHitParameters(dynamic_cast<const LArCaloHitParameters &>(parameters));
    pObject = new LArCaloHit(larCaloHitParameters);

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArCaloHitFactory::Read(Parameters &parameters, pandora::FileReader &fileReader) const
{
    // ATTN: To receive this call-back must have already set file reader mc particle factory to this factory
    unsigned int larTPCVolumeId(std::numeric_limits<unsigned int>::max());
    unsigned int daughterVolumeId(0);

    if (pandora::BINARY == fileReader.GetFileType())
    {
        pandora::BinaryFileReader &binaryFileReader(dynamic_cast<pandora::BinaryFileReader &>(fileReader));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(larTPCVolumeId));
        if (m_version > 1)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(daughterVolumeId));
    }
    else if (pandora::XML == fileReader.GetFileType())
    {
        pandora::XmlFileReader &xmlFileReader(dynamic_cast<pandora::XmlFileReader &>(fileReader));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("LArTPCVolumeId", larTPCVolumeId));
        if (m_version > 1)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("DaughterVolumeId", daughterVolumeId));
    }
    else
    {
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    LArCaloHitParameters &larCaloHitParameters(dynamic_cast<LArCaloHitParameters &>(parameters));
    larCaloHitParameters.m_larTPCVolumeId = larTPCVolumeId;
    larCaloHitParameters.m_daughterVolumeId = daughterVolumeId;

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArCaloHitFactory::Write(const Object *const pObject, pandora::FileWriter &fileWriter) const
{
    // ATTN: To receive this call-back must have already set file writer mc particle factory to this factory
    const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pObject));

    if (!pLArCaloHit)
        return pandora::STATUS_CODE_INVALID_PARAMETER;

    if (pandora::BINARY == fileWriter.GetFileType())
    {
        pandora::BinaryFileWriter &binaryFileWriter(dynamic_cast<pandora::BinaryFileWriter &>(fileWriter));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetLArTPCVolumeId()));
        if (m_version > 1)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetDaughterVolumeId()));
    }
    else if (pandora::XML == fileWriter.GetFileType())
    {
        pandora::XmlFileWriter &xmlFileWriter(dynamic_cast<pandora::XmlFileWriter &>(fileWriter));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("LArTPCVolumeId", pLArCaloHit->GetLArTPCVolumeId()));
        if (m_version > 1)
            PANDORA_RETURN_RESULT_IF(
                pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("DaughterVolumeId", pLArCaloHit->GetDaughterVolumeId()));
    }
    else
    {
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lar_content

#endif // #ifndef LAR_CALO_HIT_H
