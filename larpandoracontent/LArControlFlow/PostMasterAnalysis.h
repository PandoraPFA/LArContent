#ifndef POSTMASTERANALYSIS_H
#define POSTMASTERANALYSIS_H 1

//#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/ExternallyConfiguredAlgorithm.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

//#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"

#include "Pandora/Algorithm.h"
#include <string>
#include <vector>
#include "Objects/CaloHit.h"


namespace lar_content
{

  class PostMasterAnalysis : public pandora::Algorithm
  {
  public:
    PostMasterAnalysis();
    virtual ~PostMasterAnalysis();

   
  private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    std::string                     m_pfoListName;
    
    
  };



}
#endif //POSTMASTERANALYSIS_H
