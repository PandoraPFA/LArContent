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

/**
 *  @brief  ShowerCore class
 */
class ShowerCore
{
public:
    /**
     *  @brief  Default constructor
     */
    ShowerCore();

    /**
     *  @brief  Constructor
     *
     *  @param  startPosition the 2D position at which the cascade looks to begin
     *  @param  startDirection the initial 2D direction of the shower cascade
     */
    ShowerCore(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection);

    pandora::CartesianVector m_startPosition;   ///< the 2D position at which the cascade looks to begin
    pandora::CartesianVector m_startDirection;  ///< the initial 2D direction of the shower cascade
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

/**
 *  @brief  ConnectionPathway class
 */
class ConnectionPathway
{
public:
    /**
     *  @brief  Default constructor
     */
    ConnectionPathway();

    /**
     *  @brief  Constructor
     *
     *  @param  startPosition the start position of the connection pathway
     *  @param  startDirection the initial direction of the connection pathway
     */
    ConnectionPathway(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection);

    pandora::CartesianVector m_startPosition;   ///< the start position of the connection pathway
    pandora::CartesianVector m_startDirection;  ///< the initial direction of the connection pathway
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

/**
 *  @brief  ProtoShower class
 */
class ProtoShower
{
public: 
    /**
     *  @brief  Constructor
     *
     *  @param  showerCore the ShowerCore object
     *  @param  connectionPathway the ConnectionPathway object
     *  @param  spineHitList the shower spine hit list
     *  @param  ambiguousHitList the list of ambiguous hits (those with shared energy deposits)
     *  @param  ambiguousDirectionVector the initial directions of the ambiguous hit owners
     *  @param  hitsToAdd the list of hits to add to an electron shower pfo
     */
    ProtoShower(const ShowerCore &showerCore, const ConnectionPathway &connectionPathway, const pandora::CaloHitList &spineHitList, 
        const pandora::CaloHitList &ambiguousHitList, const pandora::CartesianPointVector &ambiguousDirectionVector, const pandora::CaloHitList &hitsToAdd);

    /**
     *  @brief  Copy constructor
     *
     *  @param  protoShower the input ProtoShower object to copy
     */
    ProtoShower(const ProtoShower &protoShower);

    ShowerCore m_showerCore;                                   ///< the ShowerCore object
    ConnectionPathway m_connectionPathway;                     ///< the ConnectionPathway object
    pandora::CaloHitList m_spineHitList;                       ///< the shower spine hit list 
    pandora::CaloHitList m_ambiguousHitList;                   ///< the list of ambiguous hits (those with shared energy deposits)
    pandora::CartesianPointVector m_ambiguousDirectionVector;  ///< the initial directions of the ambiguous hit owners 
    pandora::CaloHitList m_hitsToAdd;                          ///< the list of hits to add to an electron shower pfo
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

inline ProtoShower::ProtoShower(const ProtoShower &protoShower)
{
    m_showerCore = ShowerCore(protoShower.m_showerCore.m_startPosition, protoShower.m_showerCore.m_startDirection);
    m_connectionPathway = ConnectionPathway(protoShower.m_connectionPathway.m_startPosition, protoShower.m_connectionPathway.m_startDirection);
    m_spineHitList = protoShower.m_spineHitList;
    m_ambiguousHitList = protoShower.m_ambiguousHitList;
    m_ambiguousDirectionVector = protoShower.m_ambiguousDirectionVector;
    m_hitsToAdd = protoShower.m_hitsToAdd;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Consistency enumeration
 */
enum Consistency
{
    POSITION,
    DIRECTION,
    X_PROJECTION
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ProtoShowerMatch class
 */
class ProtoShowerMatch
{
public: 
    /**
     *  @brief  Default constructor
     */
    ProtoShowerMatch();

    /**
     *  @brief  Constructor
     *
     *  @param  protoShowerU the U view ProtoShower
     *  @param  protoShowerV the V view ProtoShower
     *  @param  protoShowerW the W view ProtoShower
     *  @param  consistencyType the nature of the 2D->3D match
     */
    ProtoShowerMatch(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, 
        const Consistency consistencyType);

    ProtoShower m_protoShowerU;     ///< the U view ProtoShower   
    ProtoShower m_protoShowerV;     ///< the V view ProtoShower
    ProtoShower m_protoShowerW;     ///< the W view ProtoShower 
    Consistency m_consistencyType;  ///< the nature of the 2D->3D match 
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
