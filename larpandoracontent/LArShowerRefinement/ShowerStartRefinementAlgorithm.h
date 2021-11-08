/**
 *  @file   larpandoracontent/LArShowerRefinement/ShowerStartRefinementAlgorithm.h
 *
 *  @brief  Header file for the shower start refinement algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_START_REFINEMENT_ALGORITHM_H
#define LAR_SHOWER_START_REFINEMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

namespace lar_content
{

class ShowerStartRefinementBaseTool;

class ShowerStartRefinementAlgorithm : public pandora::Algorithm
{
public:
    ShowerStartRefinementAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_pfoListNames;
    std::string m_neutrinoVertexListName;

    typedef std::vector<ShowerStartRefinementBaseTool *> ShowerStartRefinementToolVector;
    ShowerStartRefinementToolVector m_algorithmToolVector;
    //private:
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------



} // namespace lar_content

#endif // #ifndef LAR_SHOWER_START_REFINEMENT_BASE_TOOL_H
