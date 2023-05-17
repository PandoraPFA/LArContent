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

    void SetStartPosition(const pandora::CartesianVector &startPosition);
    void SetStartDirection(const pandora::CartesianVector &startDirection);

    // when you're not feeling lazy, change this to private 
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
//------------------------------------------------------------------------------------------------------------------------------------------

class ConnectionPathway
{
public:
    ConnectionPathway();

    ConnectionPathway(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection);

    // when you're not feeling lazy, change this to private
    pandora::CartesianVector m_startPosition;
    pandora::CartesianVector m_startDirection;
};

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
    ProtoShower();

    ProtoShower(const ShowerCore &showerCore, const ConnectionPathway &connectionPathway, const pandora::CaloHitList &spineHitList, const bool isHelper,
        const pandora::CaloHitList &ambiguousHitList, const pandora::CartesianPointVector &ambiguousDirectionVector);

    // when you're not feeling lazy, change this to private  
    ShowerCore m_showerCore;
    ConnectionPathway m_connectionPathway;
    pandora::CaloHitList m_spineHitList;
    bool m_isHelper;
    pandora::CaloHitList m_ambiguousHitList;
    pandora::CartesianPointVector m_ambiguousDirectionVector;
};

typedef std::vector<ProtoShower> ProtoShowerVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShower::ProtoShower(const ShowerCore &showerCore, const ConnectionPathway &connectionPathway, const pandora::CaloHitList &spineHitList, 
    const bool isHelper, const pandora::CaloHitList &ambiguousHitList, const pandora::CartesianPointVector &ambiguousDirectionVector) : 
        m_showerCore(showerCore), m_connectionPathway(connectionPathway), m_spineHitList(spineHitList), m_isHelper(isHelper), 
        m_ambiguousHitList(ambiguousHitList), m_ambiguousDirectionVector(ambiguousDirectionVector)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShower::ProtoShower() : m_showerCore(), m_connectionPathway(), m_spineHitList(), m_isHelper(false), m_ambiguousHitList(), m_ambiguousDirectionVector()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class ElectronProtoShower : public ProtoShower
{
public: 
    ElectronProtoShower(const ShowerCore &showerCore, const ConnectionPathway &connectionPathway, const pandora::CaloHitList &spineHitList, 
        const bool isHelper, const pandora::CaloHitList &ambiguousHitList, const pandora::CartesianPointVector &ambiguousDirectionVector, const pandora::CaloHitList &hitsToAdd);

    ElectronProtoShower(const ElectronProtoShower &electronProtoShower);

    pandora::CaloHitList m_hitsToAdd;
};

typedef std::vector<ElectronProtoShower> ElectronProtoShowerVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ElectronProtoShower::ElectronProtoShower(const ShowerCore &showerCore, const ConnectionPathway &connectionPathway, const pandora::CaloHitList &spineHitList, 
    const bool isHelper, const pandora::CaloHitList &ambiguousHitList, const pandora::CartesianPointVector &ambiguousDirectionVector, const pandora::CaloHitList &hitsToAdd) : 
        ProtoShower(showerCore, connectionPathway, spineHitList, isHelper, ambiguousHitList, ambiguousDirectionVector), m_hitsToAdd(hitsToAdd)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ElectronProtoShower::ElectronProtoShower(const ElectronProtoShower &electronProtoShower)
{
    m_showerCore = ShowerCore(electronProtoShower.m_showerCore.m_startPosition, electronProtoShower.m_showerCore.m_startDirection);
    m_connectionPathway = ConnectionPathway(electronProtoShower.m_connectionPathway.m_startPosition, electronProtoShower.m_connectionPathway.m_startDirection);
    m_spineHitList = electronProtoShower.m_spineHitList;
    m_isHelper = electronProtoShower.m_isHelper;
    m_ambiguousHitList = electronProtoShower.m_ambiguousHitList;
    m_ambiguousDirectionVector = electronProtoShower.m_ambiguousDirectionVector;

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

    ProtoShowerMatch(const ElectronProtoShower &protoShowerU, const ElectronProtoShower &protoShowerV, const ElectronProtoShower &protoShowerW, 
        const Consistency consistencyType);

    // when you're not feeling lazy, change this to private  
    ElectronProtoShower m_protoShowerU;
    ElectronProtoShower m_protoShowerV;
    ElectronProtoShower m_protoShowerW;

    Consistency m_consistencyType;
};

typedef std::vector<ProtoShowerMatch> ProtoShowerMatchVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShowerMatch::ProtoShowerMatch(const ElectronProtoShower &protoShowerU, const ElectronProtoShower &protoShowerV, const ElectronProtoShower &protoShowerW,
    const Consistency consistencyType) : 
        m_protoShowerU(protoShowerU), m_protoShowerV(protoShowerV), m_protoShowerW(protoShowerW), m_consistencyType(consistencyType)
{
}

} // namespace lar_content

#endif // #ifndef LAR_PROTO_SHOWER_H
