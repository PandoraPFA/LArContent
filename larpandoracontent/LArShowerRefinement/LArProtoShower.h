/**
 *  @file   larpandoracontent/LArShowerRefinement/LArProtoShower.h
 *
 *  @brief  Header file for the ProtoShower class.
 *
 *  $Log: $
 */
#ifndef LAR_PROTO_SHOWER_H
#define LAR_PROTO_SHOWER_H 1

#include "Pandora/PandoraInternal.h"

#include "Objects/CartesianVector.h"

namespace lar_content
{

class ShowerCore
{
public:
    ShowerCore();

    ShowerCore(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection);

    pandora::CartesianVector m_startPosition;
    pandora::CartesianVector m_startDirection;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline ShowerCore::ShowerCore(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection) : 
    m_startPosition(startPosition), m_startDirection(startDirection)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ShowerCore::ShowerCore() : m_startPosition(0.f, 0.f, 0.f), m_startDirection(0.f, 0.f, 0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class ConnectionPathway
{
public:
    ConnectionPathway();

    ConnectionPathway(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection);

    pandora::CartesianVector m_startPosition;
    pandora::CartesianVector m_startDirection;
};

typedef std::vector<ConnectionPathway> ConnectionPathwayVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ConnectionPathway::ConnectionPathway(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection) : 
    m_startPosition(startPosition), m_startDirection(startDirection)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ConnectionPathway::ConnectionPathway() : m_startPosition(0.f, 0.f, 0.f), m_startDirection(0.f, 0.f, 0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class ProtoShower
{
public: 
    ProtoShower(const ShowerCore &showerCore, const ConnectionPathway &connectionPathway, const pandora::CaloHitList &spineHitList, 
        const pandora::CaloHitList &ambiguousHitList, const pandora::CartesianPointVector &ambiguousDirectionVector, const pandora::CaloHitList &hitsToAdd);

    ProtoShower(const ProtoShower &electronProtoShower);

    ShowerCore m_showerCore;
    ConnectionPathway m_connectionPathway;
    pandora::CaloHitList m_spineHitList;
    pandora::CaloHitList m_ambiguousHitList;
    pandora::CartesianPointVector m_ambiguousDirectionVector;
    pandora::CaloHitList m_hitsToAdd;
};

typedef std::vector<ProtoShower> ProtoShowerVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShower::ProtoShower(const ShowerCore &showerCore, const ConnectionPathway &connectionPathway, const pandora::CaloHitList &spineHitList, 
    const pandora::CaloHitList &ambiguousHitList, const pandora::CartesianPointVector &ambiguousDirectionVector, const pandora::CaloHitList &hitsToAdd) : 
    m_showerCore(showerCore), m_connectionPathway(connectionPathway), m_spineHitList(spineHitList), m_ambiguousHitList(ambiguousHitList), 
    m_ambiguousDirectionVector(ambiguousDirectionVector), m_hitsToAdd(hitsToAdd)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShower::ProtoShower(const ProtoShower &electronProtoShower)
{
    m_showerCore = ShowerCore(electronProtoShower.m_showerCore.m_startPosition, electronProtoShower.m_showerCore.m_startDirection);
    m_connectionPathway = ConnectionPathway(electronProtoShower.m_connectionPathway.m_startPosition, electronProtoShower.m_connectionPathway.m_startDirection);
    m_spineHitList = electronProtoShower.m_spineHitList;
    m_ambiguousHitList = electronProtoShower.m_ambiguousHitList;
    m_ambiguousDirectionVector = electronProtoShower.m_ambiguousDirectionVector;
    m_hitsToAdd = electronProtoShower.m_hitsToAdd;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

enum Consistency
{
    POSITION,
    DIRECTION,
    X_PROJECTION
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class ProtoShowerMatch
{
public: 
    ProtoShowerMatch();

    ProtoShowerMatch(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, 
        const Consistency consistencyType);

    // when you're not feeling lazy, change this to private  
    ProtoShower m_protoShowerU;
    ProtoShower m_protoShowerV;
    ProtoShower m_protoShowerW;

    Consistency m_consistencyType;
};

typedef std::vector<ProtoShowerMatch> ProtoShowerMatchVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShowerMatch::ProtoShowerMatch(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW,
    const Consistency consistencyType) : 
        m_protoShowerU(protoShowerU), m_protoShowerV(protoShowerV), m_protoShowerW(protoShowerW), m_consistencyType(consistencyType)
{
}

} // namespace lar_content

#endif // #ifndef LAR_PROTO_SHOWER_H
