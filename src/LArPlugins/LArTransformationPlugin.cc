/**
 *  @file   LArContent/src/LArPlugins/LArTranformationPlugin.cc
 * 
 *  @brief  Implementation of the lar transformation plugin class.
 * 
 *  $Log: $
 */

#include "LArPlugins/LArTransformationPlugin.h"

namespace lar_content
{

LArTransformationPlugin::LArTransformationPlugin(const float wireZPitch) :
    m_wireZPitch(wireZPitch)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArTransformationPlugin::GetWireZPitch() const
{
    return m_wireZPitch;
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

float LArTransformationPlugin::PUPVtoPW(const float pu, const float pv) const
{
    return this->UVtoW(pu, pv); 
}

//------------------------------------------------------------------------------------------------------------------------------------------ 
     
float LArTransformationPlugin::PVPWtoPU(const float pv, const float pw) const
{
    return this->VWtoU(pv, pw);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArTransformationPlugin::PWPUtoPV(const float pw, const float pu) const
{
    return this->WUtoV(pw, pu);
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

float LArTransformationPlugin::PUPVtoPY(const float pu, const float pv)  const
{
    return this->UVtoY(pu, pv);
}

//------------------------------------------------------------------------------------------------------------------------------------------
  
float LArTransformationPlugin::PUPVtoPZ(const float pu, const float pv) const
{
    return this->UVtoZ(pu, pv); 
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
float LArTransformationPlugin::PYPZtoPU(const float py, const float pz) const
{
    return this->YZtoU(py, pz);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArTransformationPlugin::PYPZtoPV(const float py, const float pz) const
{
    return this->YZtoV(py, pz);
}

} // namespace lar_content
