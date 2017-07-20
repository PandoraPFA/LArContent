/**
 *  @file   larpandoracontent/LArUtility/ParentAlgorithm.h
 *
 *  @brief  Header file for the parent algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PARENT_ALGORITHM_H
#define LAR_PARENT_ALGORITHM_H 1

#include "larpandoracontent/LArUtility/ParentSlicingBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ParentAlgorithm class
 */
class ParentAlgorithm : public ParentSlicingBaseAlgorithm
{
private:
    pandora::StatusCode Run();
    void FastReconstruction() const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_content

#endif // #ifndef LAR_PARENT_ALGORITHM_H
