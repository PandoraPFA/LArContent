/**
 *  @file   larpandoracontent/LArShowerRefinement/ProtoShowerMatchingTool.cc
 *
 *  @brief  Implementation of the ProtoShower matching tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/ProtoShowerMatchingTool.h"

using namespace pandora;

namespace lar_content
{

ProtoShowerMatchingTool::ProtoShowerMatchingTool() :
    m_maxXSeparation(5.f),
    m_maxSeparation(5.f),
    m_maxAngularDeviation(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProtoShowerMatchingTool::Run(pandora::Algorithm *const pAlgorithm, const ElectronProtoShowerVector &protoShowerVectorU, 
    const ElectronProtoShowerVector &protoShowerVectorV, const ElectronProtoShowerVector &protoShowerVectorW, ProtoShowerMatchVector &protoShowerMatchVector)
{
    IntVector usedProtoShowersU, usedProtoShowersV, usedProtoShowersW; 

    for (unsigned int uIndex = 0; uIndex < protoShowerVectorU.size(); ++uIndex)
    {
        const ElectronProtoShower &protoShowerU(protoShowerVectorU.at(uIndex));

        for (unsigned int vIndex = 0; vIndex < protoShowerVectorV.size(); ++vIndex)
        {
            const ElectronProtoShower &protoShowerV(protoShowerVectorV.at(vIndex));

            for (unsigned int wIndex = 0; wIndex < protoShowerVectorW.size(); ++wIndex)
            {
                const ElectronProtoShower &protoShowerW(protoShowerVectorW.at(wIndex));

                if (std::find(usedProtoShowersU.begin(), usedProtoShowersU.end(), uIndex) != usedProtoShowersU.end())
                    continue;

                if (std::find(usedProtoShowersV.begin(), usedProtoShowersV.end(), vIndex) != usedProtoShowersV.end())
                    continue;             

                if (std::find(usedProtoShowersW.begin(), usedProtoShowersW.end(), wIndex) != usedProtoShowersW.end())
                    continue;

                Consistency consistency(Consistency::POSITION);

                if (!this->ArePathwaysConsistent(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, consistency))
                    continue;

                usedProtoShowersU.push_back(uIndex);
                usedProtoShowersV.push_back(vIndex);
                usedProtoShowersW.push_back(wIndex);

                std::cout << "uIndex: " << uIndex << std::endl;
                std::cout << "vIndex: " << vIndex << std::endl;
                std::cout << "wIndex: " << wIndex << std::endl;
                std::cout << "consistency: " << consistency << std::endl;


                protoShowerMatchVector.push_back(ProtoShowerMatch(protoShowerU, protoShowerV, protoShowerW, consistency));
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProtoShowerMatchingTool::ArePathwaysConsistent(pandora::Algorithm *const pAlgorithm, const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, 
    Consistency &consistency)
{
    if (LArConnectionPathwayHelper::AreShowerStartsConsistent(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, m_maxXSeparation, m_maxSeparation))
    {
        consistency = Consistency::POSITION;
    }
    else if (LArConnectionPathwayHelper::AreDirectionsConsistent(pAlgorithm, protoShowerU, protoShowerV, protoShowerW, m_maxAngularDeviation))
    {
        consistency = Consistency::DIRECTION;
    }
    else
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProtoShowerMatchingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxXSeparation", m_maxXSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxSeparation", m_maxSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxAngularDeviation", m_maxAngularDeviation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
