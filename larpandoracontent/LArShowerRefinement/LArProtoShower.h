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

#include "Objects/CaloHit.h"
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

    /**
     *  @brief  Get the start position of the shower core
     * 
     *  @return the start position
     */
    const pandora::CartesianVector &GetStartPosition() const;

    /**
     *  @brief  Get the start direction of the shower core
     * 
     *  @return the start direction
     */
    const pandora::CartesianVector &GetStartDirection() const;

private:
    pandora::CartesianVector m_startPosition;  ///< the 2D position at which the cascade looks to begin
    pandora::CartesianVector m_startDirection; ///< the initial 2D direction of the shower cascade
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline ShowerCore::ShowerCore(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection) :
    m_startPosition(startPosition),
    m_startDirection(startDirection)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ShowerCore::ShowerCore() : m_startPosition(0.f, 0.f, 0.f), m_startDirection(0.f, 0.f, 0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ShowerCore::GetStartPosition() const
{
    return m_startPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ShowerCore::GetStartDirection() const
{
    return m_startDirection;
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

    /**
     *  @brief  Get the start position of the connection pathway
     * 
     *  @return the start position
     */
    const pandora::CartesianVector &GetStartPosition() const;

    /**
     *  @brief  Get the start direction of the connection pathway
     * 
     *  @return the start direction
     */
    const pandora::CartesianVector &GetStartDirection() const;

private:
    pandora::CartesianVector m_startPosition;  ///< the start position of the connection pathway
    pandora::CartesianVector m_startDirection; ///< the initial direction of the connection pathway
};

typedef std::vector<ConnectionPathway> ConnectionPathwayVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ConnectionPathway::ConnectionPathway(const pandora::CartesianVector &startPosition, const pandora::CartesianVector &startDirection) :
    m_startPosition(startPosition),
    m_startDirection(startDirection)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ConnectionPathway::ConnectionPathway() : m_startPosition(0.f, 0.f, 0.f), m_startDirection(0.f, 0.f, 0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ConnectionPathway::GetStartPosition() const
{
    return m_startPosition;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianVector &ConnectionPathway::GetStartDirection() const
{
    return m_startDirection;
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
        const pandora::CaloHitList &ambiguousHitList, const pandora::CartesianPointVector &ambiguousDirectionVector,
        const pandora::CaloHitList &hitsToAdd);

    /**
     *  @brief  Copy constructor
     *
     *  @param  protoShower the input ProtoShower object to copy
     */
    ProtoShower(const ProtoShower &protoShower);

    /**
     *  @brief  Get the shower core
     * 
     *  @return the shower core
     */
    const ShowerCore &GetShowerCore() const;

    /**
     *  @brief  Get the connection pathway
     * 
     *  @return the connection pathway
     */
    const ConnectionPathway &GetConnectionPathway() const;

    /**
     *  @brief  Get the spine hit list
     * 
     *  @return the spine hit list
     */
    const pandora::CaloHitList &GetSpineHitList() const;

    /**
     *  @brief  Get the ambiguous hit list
     * 
     *  @return the ambiguous hit list
     */
    const pandora::CaloHitList &GetAmbiguousHitList() const;

    /**
     *  @brief  Add an ambiguous hit to the ambiguous hit list
     *
     *  @param  ambiguousHit the ambiguous hit to add
     */
    void AddAmbiguousHit(const pandora::CaloHit *const ambiguousHit);

    /**
     *  @brief  Get the ambiguous direction vector
     * 
     *  @return the ambiguous direction vector
     */
    const pandora::CartesianPointVector &GetAmbiguousDirectionVector() const;

    /**
     *  @brief  Add an ambiguous direction to the ambiguous direction vector
     *
     *  @param  ambiguousDirection the ambiguous direction to add
     */
    void AddAmbiguousDirection(const pandora::CartesianVector &ambiguousDirection);

    /**
     *  @brief  Get the hits to add list
     * 
     *  @return the hits to add list
     */
    const pandora::CaloHitList &GetHitsToAddList() const;

    /**
     *  @brief  Set the hits to add list
     * 
     *  @param  the hits to add list
     */
    void SetHitsToAddList(const pandora::CaloHitList &hitsToAddList);

    /**
     *  @brief  Add a hit to the hits to add list
     *
     *  @param  hitToAdd the hit to add
     */
    void AddHitToAdd(const pandora::CaloHit *const hitToAdd);

private:
    const ShowerCore m_showerCore;                            ///< the ShowerCore object
    const ConnectionPathway m_connectionPathway;              ///< the ConnectionPathway object
    pandora::CaloHitList m_spineHitList;                      ///< the shower spine hit list
    pandora::CaloHitList m_ambiguousHitList;                  ///< the list of ambiguous hits (those with shared energy deposits)
    pandora::CartesianPointVector m_ambiguousDirectionVector; ///< the initial directions of the ambiguous hit owners
    pandora::CaloHitList m_hitsToAdd;                         ///< the list of hits to add to an electron shower pfo
};

typedef std::vector<ProtoShower> ProtoShowerVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShower::ProtoShower(const ShowerCore &showerCore, const ConnectionPathway &connectionPathway,
    const pandora::CaloHitList &spineHitList, const pandora::CaloHitList &ambiguousHitList,
    const pandora::CartesianPointVector &ambiguousDirectionVector, const pandora::CaloHitList &hitsToAdd) :
    m_showerCore(showerCore),
    m_connectionPathway(connectionPathway),
    m_spineHitList(spineHitList),
    m_ambiguousHitList(ambiguousHitList),
    m_ambiguousDirectionVector(ambiguousDirectionVector),
    m_hitsToAdd(hitsToAdd)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShower::ProtoShower(const ProtoShower &protoShower) :
    m_showerCore(ShowerCore(protoShower.GetShowerCore().GetStartPosition(), protoShower.GetShowerCore().GetStartDirection())),
    m_connectionPathway(
        ConnectionPathway(protoShower.GetConnectionPathway().GetStartPosition(), protoShower.GetConnectionPathway().GetStartDirection())),
    m_spineHitList(protoShower.GetSpineHitList()),
    m_ambiguousHitList(protoShower.GetAmbiguousHitList()),
    m_ambiguousDirectionVector(protoShower.GetAmbiguousDirectionVector()),
    m_hitsToAdd(protoShower.GetHitsToAddList())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ShowerCore &ProtoShower::GetShowerCore() const
{
    return m_showerCore;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ConnectionPathway &ProtoShower::GetConnectionPathway() const
{
    return m_connectionPathway;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &ProtoShower::GetSpineHitList() const
{
    return m_spineHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &ProtoShower::GetAmbiguousHitList() const
{
    return m_ambiguousHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ProtoShower::AddAmbiguousHit(const pandora::CaloHit *const ambiguousHit)
{
    m_ambiguousHitList.push_back(ambiguousHit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CartesianPointVector &ProtoShower::GetAmbiguousDirectionVector() const
{
    return m_ambiguousDirectionVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ProtoShower::AddAmbiguousDirection(const pandora::CartesianVector &ambiguousDirection)
{
    m_ambiguousDirectionVector.push_back(ambiguousDirection);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::CaloHitList &ProtoShower::GetHitsToAddList() const
{
    return m_hitsToAdd;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ProtoShower::SetHitsToAddList(const pandora::CaloHitList &hitsToAddList)
{
    m_hitsToAdd = hitsToAddList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void ProtoShower::AddHitToAdd(const pandora::CaloHit *const hitToAdd)
{
    m_hitsToAdd.push_back(hitToAdd);
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
    ProtoShowerMatch(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const Consistency consistencyType);

    /**
     *  @brief  Get the U view ProtoShower
     * 
     *  @return the U view ProtoShower
     */
    const ProtoShower &GetProtoShowerU() const;

    /**
     *  @brief  Get the V view ProtoShower
     * 
     *  @return the V view ProtoShower
     */
    const ProtoShower &GetProtoShowerV() const;

    /**
     *  @brief  Get the W view ProtoShower
     * 
     *  @return the W view ProtoShower
     */
    const ProtoShower &GetProtoShowerW() const;

    /**
     *  @brief  Get a modifiable ProtoShower in a given view
     * 
     *  @param  hitType the 2D view
     *
     *  @return the modifiable ProtoShower in the specified view
     */
    ProtoShower &GetProtoShowerToModify(const pandora::HitType hitType);

    /**
     *  @brief  Get the consistency type
     * 
     *  @return the consistency type
     */
    const Consistency &GetConsistencyType() const;

private:
    ProtoShower m_protoShowerU;    ///< the U view ProtoShower
    ProtoShower m_protoShowerV;    ///< the V view ProtoShower
    ProtoShower m_protoShowerW;    ///< the W view ProtoShower
    Consistency m_consistencyType; ///< the nature of the 2D->3D match
};

typedef std::vector<ProtoShowerMatch> ProtoShowerMatchVector;

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShowerMatch::ProtoShowerMatch(
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, const Consistency consistencyType) :
    m_protoShowerU(protoShowerU),
    m_protoShowerV(protoShowerV),
    m_protoShowerW(protoShowerW),
    m_consistencyType(consistencyType)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ProtoShower &ProtoShowerMatch::GetProtoShowerU() const
{
    return m_protoShowerU;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ProtoShower &ProtoShowerMatch::GetProtoShowerV() const
{
    return m_protoShowerV;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const ProtoShower &ProtoShowerMatch::GetProtoShowerW() const
{
    return m_protoShowerW;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline ProtoShower &ProtoShowerMatch::GetProtoShowerToModify(const pandora::HitType hitType)
{
    return hitType == pandora::TPC_VIEW_U ? m_protoShowerU : (hitType == pandora::TPC_VIEW_V ? m_protoShowerV : m_protoShowerW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline const Consistency &ProtoShowerMatch::GetConsistencyType() const
{
    return m_consistencyType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_content

#endif // #ifndef LAR_PROTO_SHOWER_H
