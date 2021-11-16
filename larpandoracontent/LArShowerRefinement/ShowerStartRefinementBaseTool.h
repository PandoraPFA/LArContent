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

class AngularPeak
{
public:
    AngularPeak(int uPeakBin, int vPeakBin, int wPeakBin, float chiSquared, float binWeightSum);

    int m_uPeakBin;
    int m_vPeakBin;
    int m_wPeakBin;
    float m_chiSquared;
    float m_binWeightSum;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline AngularPeak::AngularPeak(int uPeakBin, int vPeakBin, int wPeakBin, float chiSquared, float binWeightSum)
{
    m_uPeakBin = uPeakBin;
    m_vPeakBin = vPeakBin;
    m_wPeakBin = wPeakBin;
    m_chiSquared = chiSquared;
    m_binWeightSum = binWeightSum;
}

typedef std::vector<AngularPeak> AngularPeakVector;

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class ShowerStartRefinementBaseTool : public pandora::AlgorithmTool
{
public:
    ShowerStartRefinementBaseTool();

    virtual bool Run(ShowerStartRefinementAlgorithm *const pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertexPosition) = 0;

    typedef std::map<int, float> DeviationAngleMap;
    typedef std::map<const pandora::CaloHit*, float> LongitudinalPositionMap;

protected:
    bool HasPathToNuVertex(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &neutrinoVertex) const;

    void FindShowerSpine(const ShowerStartRefinementAlgorithm *pAlgorithm, const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &neutrinoVertex,
                         const pandora::CartesianVector &initialDirection, const pandora::HitType hitType, pandora::CaloHitList &showerSpineHitList, LongitudinalPositionMap &longitudinalPositionMap);

    void GetHitsInBoundingBox(const pandora::CartesianVector &firstCorner, const pandora::CartesianVector &secondCorner, const pandora::CaloHitList &inputHitList,
        const float distanceToLine, pandora::CaloHitList &hitsInBoundingBox) const; 

    bool IsInBoundingBox(const float minX, const float maxX, const float minZ, const float maxZ, const pandora::CartesianVector &hitPosition) const;

    bool IsCloseToLine(const pandora::CartesianVector &hitPosition, const pandora::CartesianVector &lineStart, const pandora::CartesianVector &lineDirection, 
        const float distanceToLine) const;

    bool IsInLineSegment(const pandora::CartesianVector &lowerBoundary, const pandora::CartesianVector &upperBoundary, const pandora::CartesianVector &point) const;

    void BuildProtoShowers(const pandora::ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const;

    void FindShowerCores(const pandora::ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const;

    void FindShowerStartPositions(const pandora::ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const;

    void FindConnectionPathways(const pandora::ParticleFlowObject *const pShowerPfo, ProtoShowerVector &protoShowerVector) const;

    bool IsElectronPathway(const ProtoShower &protoShower);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_maxDistanceForConnection;
    float m_growingFitInitialLength;
    float m_macroSlidingFitWindow;
    float m_growingFitSegmentLength;
    float m_distanceToLine;
    int m_maxFittingHits;
    //private:
};

        class SortByDistanceToPoint
        {
        public:
            /**
             *  @brief  Constructor
             *
             *  @param  referencePoint the point relative to which constituent hits are ordered
             */
        SortByDistanceToPoint(const pandora::CartesianVector referencePoint) : m_referencePoint(referencePoint)
            {
            }

            /**
             *  @brief  Sort constituent hits by their position relative to a referencePoint
             *
             *  @param  lhs first constituent hit
             *  @param  rhs second constituent hit
             *
             *  @return  whether lhs hit is closer to the referencePoint than the rhs hit
             */
            bool operator()(const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs);

        private:
            const pandora::CartesianVector m_referencePoint; ///< The point relative to which constituent hits are ordered
        };

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------



} // namespace lar_content

#endif // #ifndef LAR_SHOWER_START_REFINEMENT_BASE_TOOL_H
