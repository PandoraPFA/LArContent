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
    pandora::InputUInt m_larTPCVolumeId;    ///< The lar tpc volume id
    pandora::InputUInt m_daughterVolumeId;  ///< The daughter volume id
    pandora::InputUInt m_plane;             ///< The plane id    
    pandora::InputUInt m_wireId;            ///< The wire id (only unique to daughter volume)
    pandora::InputUInt m_plane1;            ///< The first projection plane id
    pandora::InputUInt m_minIntersectWire1; ///< The id of the first (lowest wire id) intersecting wire (in proj plane 1)
    pandora::InputUInt m_maxIntersectWire1; ///< The id of the last (highest wire id) intersecting wire (in proj plane 1)
    pandora::InputUInt m_plane2;            ///< The second projection plane id    
    pandora::InputUInt m_minIntersectWire2; ///< The id of the first (lowest wire id) intersecting wire (in proj plane 2)
    pandora::InputUInt m_maxIntersectWire2; ///< The id of the last (highest wire id) intersecting wire (in proj plane 2)
    pandora::FloatVector m_hitScores;       ///< Hit scores
    pandora::StringVector m_hitScoreLabels; ///< Labels for the hit scores
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
     *  @brief  Get the plane id
     *
     *  @return the plane id
     */
    unsigned int GetPlane() const;
    
    /**
     *  @brief  Get the wire id
     *
     *  @return the wire id
     */
    unsigned int GetWireId() const;

    /**
     *  @brief  Get the first projection plane id
     *
     *  @return the first projection plane id
     */
    unsigned int GetPlane1() const;
    
    /**
     *  @brief  Get the id of the minimum (lowest wire id) intersecting wire in proj plane 1
     *
     *  @return the minimum intersecting wire id for proj plane 1
     */
    unsigned int GetMinIntersectWire1() const;

    /**
     *  @brief  Get the id of the maximum (highest wire id) intersecting wire in proj plane 1
     *
     *  @return the maximum intersecting wire id for proj plane 1
     */
    unsigned int GetMaxIntersectWire1() const;

    /**
     *  @brief  Get the second projection plane id
     *
     *  @return the second projection plane id
     */
    unsigned int GetPlane2() const;
    
    /**
     *  @brief  Get the id of the minimum (lowest wire id) intersecting wire in proj plane 2
     *
     *  @return the minimum intersecting wire id for proj plane 2
     */
    unsigned int GetMinIntersectWire2() const;

    /**
     *  @brief  Get the id of the maximum (highest wire id) intersecting wire in proj plane 2
     *
     *  @return the maximum intersecting wire id for proj plane 2
     */
    unsigned int GetMaxIntersectWire2() const;    
    
    /**
     *  @brief  Fill the parameters associated with this calo hit
     *
     *  @param  parameters the output parameters
     */
    void FillParameters(LArCaloHitParameters &parameters) const;

    /**
     *  @brief  Get scores for each hit
     *
     *  @return a vector of scores for each category
     */
    const pandora::FloatVector &GetHitScores() const;

    /**
     *  @brief  Get scores for each hit
     *
     *  @return a vector of score labels for each category
     */
    const pandora::StringVector &GetHitScoreLabels() const;

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
    unsigned int m_larTPCVolumeId;          ///< The lar tpc volume id
    unsigned int m_daughterVolumeId;        ///< The daughter volume id
    unsigned int m_plane;                   ///< The plane id
    unsigned int m_wireId;                  ///< The wire id (only unique to daughter volume)
    unsigned int m_plane1;                  ///< The first projection plane id        
    unsigned int m_minIntersectWire1;       ///< The id of the first (lowest wire id) intersecting wire (in proj plane 1)
    unsigned int m_maxIntersectWire1;       ///< The id of the last (highest wire id) intersecting wire (in proj plane 1)
    unsigned int m_plane2;                  ///< The second projection plane id    
    unsigned int m_minIntersectWire2;       ///< The id of the first (lowest wire id) intersecting wire (in proj plane 2)
    unsigned int m_maxIntersectWire2;       ///< The id of the last (highest wire id) intersecting wire (in proj plane 2)
    pandora::FloatVector m_hitScores;       ///< Hit scores
    pandora::StringVector m_hitScoreLabels; ///< Labels for the hit scores
    pandora::InputFloat m_pTrack;           ///< The probability that the hit is track-like
    pandora::InputFloat m_pShower;          ///< The probability that the hit is shower-like
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
    m_plane(parameters.m_plane.Get()),
    m_wireId(parameters.m_wireId.Get()),
    m_plane1(parameters.m_plane1.Get()),
    m_minIntersectWire1(parameters.m_minIntersectWire1.Get()),
    m_maxIntersectWire1(parameters.m_maxIntersectWire1.Get()),
    m_plane2(parameters.m_plane2.Get()),
    m_minIntersectWire2(parameters.m_minIntersectWire2.Get()),
    m_maxIntersectWire2(parameters.m_maxIntersectWire2.Get()),
    m_hitScores(parameters.m_hitScores),
    m_hitScoreLabels(parameters.m_hitScoreLabels)
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

inline unsigned int LArCaloHit::GetPlane() const
{
    return m_plane;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArCaloHit::GetWireId() const
{
    return m_wireId;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArCaloHit::GetPlane1() const
{
    return m_plane1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArCaloHit::GetMinIntersectWire1() const
{
    return m_minIntersectWire1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArCaloHit::GetMaxIntersectWire1() const
{
    return m_maxIntersectWire1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArCaloHit::GetPlane2() const
{
    return m_plane2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArCaloHit::GetMinIntersectWire2() const
{
    return m_minIntersectWire2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int LArCaloHit::GetMaxIntersectWire2() const
{
    return m_maxIntersectWire2;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::FloatVector &LArCaloHit::GetHitScores() const
{
    return m_hitScores;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::StringVector &LArCaloHit::GetHitScoreLabels() const
{
    return m_hitScoreLabels;
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
    parameters.m_plane = this->GetPlane();
    parameters.m_wireId = this->GetWireId();
    parameters.m_plane1 = this->GetPlane1();
    parameters.m_minIntersectWire1 = this->GetMinIntersectWire1();
    parameters.m_maxIntersectWire1 = this->GetMaxIntersectWire1();
    parameters.m_plane2 = this->GetPlane2();
    parameters.m_minIntersectWire2 = this->GetMinIntersectWire2();
    parameters.m_maxIntersectWire2 = this->GetMaxIntersectWire2();
    parameters.m_hitScores = this->GetHitScores();
    parameters.m_hitScoreLabels = this->GetHitScoreLabels();
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
    unsigned int nLabels(0);
    pandora::FloatVector hitScores;
    pandora::StringVector hitScoreLabels;
    unsigned int plane(0), plane1(0), plane2(0);
    unsigned int wireId(0);
    unsigned int minIntersectWire1(0), maxIntersectWire1(0);
    unsigned int minIntersectWire2(0), maxIntersectWire2(0);     

    if (pandora::BINARY == fileReader.GetFileType())
    {
        pandora::BinaryFileReader &binaryFileReader(dynamic_cast<pandora::BinaryFileReader &>(fileReader));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(larTPCVolumeId));
        if (m_version > 1)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(daughterVolumeId));
        if (m_version > 2)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(nLabels));
            hitScores.reserve(nLabels);
            hitScoreLabels.reserve(nLabels);
            for (unsigned int i = 0; i < nLabels; ++i)
            {
                float score;
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(score));
                hitScores.emplace_back(std::move(score));

                std::string label;
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(label));
                hitScoreLabels.emplace_back(std::move(label));
            }
        }
        if (m_version > 3)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(plane));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(wireId));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(plane1));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(minIntersectWire1));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(maxIntersectWire1));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(plane2));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(minIntersectWire2));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileReader.ReadVariable(maxIntersectWire2));
        }
        
    }
    else if (pandora::XML == fileReader.GetFileType())
    {
        pandora::XmlFileReader &xmlFileReader(dynamic_cast<pandora::XmlFileReader &>(fileReader));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("LArTPCVolumeId", larTPCVolumeId));
        if (m_version > 1)
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("DaughterVolumeId", daughterVolumeId));
        if (m_version > 2)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("NHitLabels", nLabels));
            hitScores.reserve(nLabels);
            hitScoreLabels.reserve(nLabels);
            for (unsigned int i = 0; i < nLabels; ++i)
            {
                float score;
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("HitScore", score));
                hitScores.emplace_back(std::move(score));

                std::string label;
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("HitScoreLabel", label));
                hitScoreLabels.emplace_back(std::move(label));
            }
        }
        if (m_version > 3)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("Plane", plane));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("WireID", wireId));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("Plane1", plane1));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("MinIntersectWire1", minIntersectWire1));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("MaxIntersectWire1", maxIntersectWire1));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("Plane2", plane2));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("MinIntersectWire2", minIntersectWire2));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileReader.ReadVariable("MaxIntersectWire2", maxIntersectWire2));
        }
    }
    else
    {
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    LArCaloHitParameters &larCaloHitParameters(dynamic_cast<LArCaloHitParameters &>(parameters));
    larCaloHitParameters.m_larTPCVolumeId = larTPCVolumeId;
    larCaloHitParameters.m_daughterVolumeId = daughterVolumeId;
    larCaloHitParameters.m_hitScores = std::move(hitScores);
    larCaloHitParameters.m_hitScoreLabels = std::move(hitScoreLabels);
    larCaloHitParameters.m_plane = plane;
    larCaloHitParameters.m_wireId = wireId;
    larCaloHitParameters.m_plane1 = plane1;
    larCaloHitParameters.m_minIntersectWire1 = minIntersectWire1;
    larCaloHitParameters.m_maxIntersectWire1 = maxIntersectWire1;
    larCaloHitParameters.m_plane2 = plane2;
    larCaloHitParameters.m_minIntersectWire2 = minIntersectWire2;
    larCaloHitParameters.m_maxIntersectWire2 = maxIntersectWire2;    

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
        if (m_version > 2)
        {
            const pandora::FloatVector &hitScores(pLArCaloHit->GetHitScores());
            const pandora::StringVector &hitScoreLabels(pLArCaloHit->GetHitScoreLabels());
            const unsigned int nLabels(hitScores.size());
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(nLabels));
            for (unsigned int i = 0; i < nLabels; ++i)
            {
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(hitScores.at(i)));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(hitScoreLabels.at(i)));
            }
        }
        if (m_version > 3)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetPlane()));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetWireId()));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetPlane1()));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetMinIntersectWire1()));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetMaxIntersectWire1()));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetPlane2()));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetMinIntersectWire2()));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, binaryFileWriter.WriteVariable(pLArCaloHit->GetMaxIntersectWire2()));
        }
    }
    else if (pandora::XML == fileWriter.GetFileType())
    {
        pandora::XmlFileWriter &xmlFileWriter(dynamic_cast<pandora::XmlFileWriter &>(fileWriter));
        PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("LArTPCVolumeId", pLArCaloHit->GetLArTPCVolumeId()));
        if (m_version > 1)
            PANDORA_RETURN_RESULT_IF(
                pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("DaughterVolumeId", pLArCaloHit->GetDaughterVolumeId()));
        if (m_version > 2)
        {
            const pandora::FloatVector &hitScores(pLArCaloHit->GetHitScores());
            const pandora::StringVector &hitScoreLabels(pLArCaloHit->GetHitScoreLabels());
            const unsigned int nLabels(hitScores.size());
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("NHitLabels", nLabels));
            for (unsigned int i = 0; i < nLabels; ++i)
            {
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("HitScore", hitScores.at(i)));
                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("HitScoreLabel", hitScoreLabels.at(i)));
            }
        }
        if (m_version > 3)
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("Plane", pLArCaloHit->GetPlane()));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("WireId", pLArCaloHit->GetWireId()));
             PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("Plane1", pLArCaloHit->GetPlane1()));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("MinIntersectWire1", pLArCaloHit->GetMinIntersectWire1()));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("MaxIntersectWire1", pLArCaloHit->GetMaxIntersectWire1()));
             PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("Plane2", pLArCaloHit->GetPlane2()));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("MinIntersectWire2", pLArCaloHit->GetMinIntersectWire2()));
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, xmlFileWriter.WriteVariable("MaxIntersectWire2", pLArCaloHit->GetMaxIntersectWire2()));            
        }
    }
    else
    {
        return pandora::STATUS_CODE_INVALID_PARAMETER;
    }

    return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lar_content

#endif // #ifndef LAR_CALO_HIT_H
