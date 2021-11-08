/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartRefinementBaseTool.h
 *
 *  @brief  Header file for the shower characterisation tool class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_START_REFINEMENT_BASE_TOOL_H
#define LAR_SHOWER_START_REFINEMENT_BASE_TOOL_H 1

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.h"

namespace lar_content
{

class ShowerCore
{
public:
    ShowerCore();

    void SetStartPosition(const pandora::CartesianVector &startPosition);

    void SetStartDirection(const pandora::CartesianVector &startDirection);

    // when you're not feeling lazy, change this to private 
    pandora::CartesianVector m_startPosition;
    pandora::CartesianVector m_startDirection;
    pandora::CaloHitList m_coreHitList;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline ShowerCore::ShowerCore() : m_startPosition(0.f, 0.f, 0.f), m_startDirection(0.f, 0.f, 0.f), m_coreHitList()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ShowerCore::SetStartPosition(const pandora::CartesianVector &startPosition)
{
    m_startPosition = startPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ShowerCore::SetStartDirection(const pandora::CartesianVector &startDirection)
{
    m_startDirection = startDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

/*
inline ShowerCore::SetStartDirection(const CartesianVector &startDirection)
{
    m_startDirection = startDirection;
}
*/

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class ConnectionPathway
{
public:
    ConnectionPathway();

    // when you're not feeling lazy, change this to private
    pandora::CartesianVector m_startPosition;
    pandora::CartesianVector m_startDirection;
    pandora::CaloHitList m_pathwayHitList;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline ConnectionPathway::ConnectionPathway() : m_startPosition(0.f, 0.f, 0.f), m_startDirection(0.f, 0.f, 0.f), m_pathwayHitList()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class ProtoShower
{
public: 
    ProtoShower();

    // when you're not feeling lazy, change this to private  
    ShowerCore m_showerCore;
    ConnectionPathway m_connectionPathway;
};

typedef std::vector<ProtoShower> ProtoShowerVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShower::ProtoShower() : m_showerCore(), m_connectionPathway()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class ShowerStartRefinementBaseTool : public pandora::AlgorithmTool
{
public:
    ShowerStartRefinementBaseTool();

    virtual bool Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertexPosition) = 0;

protected:
    bool HasPathToNuVertex(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &neutrinoVertex) const;

    void BuildProtoShowers(const pandora::ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const;

    void FindShowerCores(const pandora::ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const;

    void FindShowerStartPositions(const pandora::ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const;

    void FindConnectionPathways(const pandora::ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const;

    bool IsElectronPathway(const ProtoShower &protoShower);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_maxDistanceForConnection;

    //private:
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------



} // namespace lar_content

#endif // #ifndef LAR_SHOWER_START_REFINEMENT_BASE_TOOL_H
