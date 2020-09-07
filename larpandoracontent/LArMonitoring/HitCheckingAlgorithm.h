/**
 *  @file   larpandoracontent/LArMonitoring/HitCheckingAlgorithm.h
 *
 *  @brief  Header file for the hit checking algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_HIT_CHECKING_ALGORITHM_H
#define LAR_HIT_CHECKING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  HitCheckingAlgorithm class
 */
class HitCheckingAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    HitCheckingAlgorithm();

    ~HitCheckingAlgorithm();
    
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    std::string m_caloHitListName;
    bool m_writeToTree;
    std::string m_treeName;
    std::string m_fileName;
    int m_eventNumber;
};

} // namespace lar_content

#endif // LAR_HIT_CHECKING_ALGORITHM_H
