#ifndef UG_PROJECT_ALGORITHM_H
#define UG_PROJECT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  UGProjectAlgorithm class
 */
class UGProjectAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    UGProjectAlgorithm();

    ~UGProjectAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    // Member variables here
    bool m_writeToTree;
    int m_eventId;
    std::string m_treeName;
    std::string m_fileName;
    std::string m_mcParticleListName;
    std::string m_caloHitListName;
    std::string m_inputPfoListName;
};

} // namespace lar_content

#endif // #ifndef UG_PROJECT_ALGORITHM_H
