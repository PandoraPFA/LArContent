/**
 *  @file   larpandoracontent/LArUtility/ParentBaseAlgorithm.h
 *
 *  @brief  Header file for the parent base algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PARENT_BASE_ALGORITHM_H
#define LAR_PARENT_BASE_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  ParentBaseAlgorithm class
 */
class ParentBaseAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Destructor
     */
    virtual ~ParentBaseAlgorithm();

    typedef std::vector<pandora::HitType> HitTypeList;
    typedef std::map<pandora::HitType, std::string> HitTypeToNameMap;

protected:
    pandora::StatusCode Initialize();

    /**
     *  @brief  Run provided algorithm
     *
     *  @param  algorithmName
     */
    void RunAlgorithm(const std::string &algorithmName) const;

    /**
     *  @brief  Run algorithms in provided list
     *
     *  @param  algorithmNames
     */
    void RunAlgorithms(const pandora::StringVector &algorithmNames) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    HitTypeList                 m_hitTypeList;                      ///< The hit type list
    HitTypeToNameMap            m_caloHitListNames;                 ///< The hit type to calo hit list name map
    HitTypeToNameMap            m_clusterListNames;                 ///< The hit type to cluster list name map

    std::string                 m_caloHitListNameU;                 ///< The name of the u input calo hit list
    std::string                 m_caloHitListNameV;                 ///< The name of the v input calo hit list
    std::string                 m_caloHitListNameW;                 ///< The name of the w input calo hit list

    std::string                 m_clusterListNameU;                 ///< The name of the u working cluster list
    std::string                 m_clusterListNameV;                 ///< The name of the v working cluster list
    std::string                 m_clusterListNameW;                 ///< The name of the w working cluster list
};

} // namespace lar_content

#endif // #ifndef LAR_PARENT_BASE_ALGORITHM_H
