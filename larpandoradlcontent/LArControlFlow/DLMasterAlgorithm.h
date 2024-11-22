/**
 *  @file   larpandoradlcontent/LArControlFlow/DLMasterAlgorithm.h
 *
 *  @brief  Header file for the master algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_MASTER_ALGORITHM_H
#define LAR_DL_MASTER_ALGORITHM_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  MasterAlgorithm class
 */
class DLMasterAlgorithm : public lar_content::MasterAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DLMasterAlgorithm() = default;

private:
    pandora::StatusCode Run();

protected:
    pandora::StatusCode RegisterCustomContent(const pandora::Pandora *const pPandora) const;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_MASTER_ALGORITHM_H
