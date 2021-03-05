/**
 *  @file   larpandoradlcontent/LArMonitoring/DlHitValidationAlgorithm.h
 *
 *  @brief  Header file for the deep learning track shower id validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_HIT_VALIDATION_ALGORITHM_H
#define LAR_DL_HIT_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DlHitValidationlgorithm class
 */
class DlHitValidationAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DlHitValidationAlgorithm();

    virtual ~DlHitValidationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector     m_caloHitListNames;    ///< Name of input calo hit list
    int                       m_confusionU[2][2];    ///< Confusion matrix for the U view
    int                       m_confusionV[2][2];    ///< Confusion matrix for the V view
    int                       m_confusionW[2][2];    ///< Confusion matrix for the W view
};

} // namespace lar_dl_content

#endif // LAR_DL_HIT_VALIDATION_ALGORITHM_H
