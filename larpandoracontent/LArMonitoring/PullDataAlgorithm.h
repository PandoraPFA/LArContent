/**
 *  @file   larpandoracontent/LArMonitoring/PullDataAlgorithm.h
 *
 *  @brief Header file for the performance assessment algorithm
 *
 * $Log: $
 */

#ifndef LAR_PULL_DATA_ALGORITHM_H
#define LAR_PULL_DATA_ALGORITHM_H

#include "Pandora/Algorithm.h"


namespace lar_content {

  class PullDataAlgorithm : public pandora::Algorithm {

  public:

    PullDataAlgorithm();
    ~PullDataAlgorithm();

  private:

    pandora::StatusCode Run();
    void GetLArSoftAngles(const pandora::CartesianVector &vector, float &theta0XZ, float &theta0YZ);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    int m_eventNumber;

    std::string m_treeName;
    std::string m_fileName;

  };



} //namespace lar_content


#endif 
