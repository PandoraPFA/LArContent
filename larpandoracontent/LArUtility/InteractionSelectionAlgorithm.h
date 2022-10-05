/**
 *  @file   larpandoracontent/LArUtility/InteractionSelectionAlgorithm.h
 *
 *  @brief  Header file for the post processing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_INTERACTION_SELECTION_ALGORITHM_H
#define LAR_INTERACTION_SELECTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  InteractionSelectionAlgorithm class
 */
class InteractionSelectionAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    InteractionSelectionAlgorithm();

private:
    pandora::StatusCode Run();

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::IntVector m_interactionIds;
};

} // namespace lar_content

#endif // #ifndef LAR_INTERACTION_SELECTION_ALGORITHM_H
