/**
 *  @file   larpandoracontent/LArPlugins/LArTranformationPlugin.cc
 * 
 *  @brief  Implementation of the lar transformation plugin class.
 * 
 *  $Log: $
 */

#include "larpandoracontent/LArPlugins/LArTransformationPlugin.h"

namespace lar_content
{

LArTransformationPlugin::~LArTransformationPlugin()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

double LArTransformationPlugin::PUPVtoPW(const double pu, const double pv) const
{
    return this->UVtoW(pu, pv); 
}

//------------------------------------------------------------------------------------------------------------------------------------------ 
     
double LArTransformationPlugin::PVPWtoPU(const double pv, const double pw) const
{
    return this->VWtoU(pv, pw);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArTransformationPlugin::PWPUtoPV(const double pw, const double pu) const
{
    return this->WUtoV(pw, pu);
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

double LArTransformationPlugin::PUPVtoPY(const double pu, const double pv)  const
{
    return this->UVtoY(pu, pv);
}

//------------------------------------------------------------------------------------------------------------------------------------------
  
double LArTransformationPlugin::PUPVtoPZ(const double pu, const double pv) const
{
    return this->UVtoZ(pu, pv); 
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
double LArTransformationPlugin::PYPZtoPU(const double py, const double pz) const
{
    return this->YZtoU(py, pz);
}

//------------------------------------------------------------------------------------------------------------------------------------------

double LArTransformationPlugin::PYPZtoPV(const double py, const double pz) const
{
    return this->YZtoV(py, pz);
}

} // namespace lar_content
