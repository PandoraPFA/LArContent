/**
 *  @file   LArContent/include/LArObjects/LAr3DTrackPfo.h
 * 
 *  @brief  Header file for the lar 3d track pfo class.
 * 
 *  $Log: $
 */
#ifndef LAR_3D_TRACK_PFO_H
#define LAR_3D_TRACK_PFO_H 1

#include "Objects/ParticleFlowObject.h"

#include "Pandora/ObjectFactory.h"

#include <string>

namespace lar_content
{

/**
 *  @brief  LAr3DTrackPfo - simply add an additional property to the pandora particle flow object
 */
class LAr3DTrackPfo : public pandora::ParticleFlowObject
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  parameters the parameters
     */
    LAr3DTrackPfo(const PandoraContentApi::ParticleFlowObject::Parameters &parameters);

    /**
     *  @brief  Set the additional property string
     * 
     *  @param  additionalProperty the additional property string
     */
    void SetAdditionalProperty(const std::string &additionalProperty);

    /**
     *  @brief  Get the additional property string
     * 
     *  @return the additional property string
     */
    const std::string &GetAdditionalProperty() const;

private:
    std::string     m_additionalProperty;       ///< The additional property string
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LAr3DTrackPfoFactory
 */
class LAr3DTrackPfoFactory : public pandora::ObjectFactory<PandoraContentApi::ParticleFlowObject::Parameters, pandora::ParticleFlowObject>
{
public:
    /**
     *  @brief  Create an object with the given parameters
     *
     *  @param  parameters the parameters to pass in constructor
     *  @param  pObject to receive the address of the object created
     */
    pandora::StatusCode Create(const PandoraContentApi::ParticleFlowObject::Parameters &parameters, const pandora::ParticleFlowObject *&pObject) const;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline LAr3DTrackPfo::LAr3DTrackPfo(const PandoraContentApi::ParticleFlowObject::Parameters &parameters) :
    pandora::ParticleFlowObject(parameters),
    m_additionalProperty(std::string())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LAr3DTrackPfo::SetAdditionalProperty(const std::string &additionalProperty)
{
    m_additionalProperty = additionalProperty;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const std::string &LAr3DTrackPfo::GetAdditionalProperty() const
{
    return m_additionalProperty;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LAr3DTrackPfoFactory::Create(const Parameters &parameters, const pandora::ParticleFlowObject *&pObject) const
{
    pObject = new LAr3DTrackPfo(parameters);

    return pandora::STATUS_CODE_SUCCESS;
}

} // namespace lar_content

#endif // #ifndef LAR_3D_TRACK_PFO_H
