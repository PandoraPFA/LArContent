/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoMopUp/ShowerHierarchyMopUpAlgorithm.h
 *
 *  @brief  Header file for the shower hierarchy mop up algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_SHOWER_HIERARCHY_MOP_UP_ALGORITHM_H
#define LAR_SHOWER_HIERARCHY_MOP_UP_ALGORITHM_H 1

#include "larpandoracontent/LArUtility/PfoMopUpBaseAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  ShowerHierarchyMopUpAlgorithm class
 */
class ShowerHierarchyMopUpAlgorithm : public PfoMopUpBaseAlgorithm
{
private:
    pandora::StatusCode Run();

    /**
     *  @brief  Starting with provided leading pfos, find all shower pfos that themselves have daughter pfos
     *
     *  @param  pLeadingPfoList the list of leading pfos
     *  @param  parentShowerPfos to receive the list of parent shower pfos
     */
    void FindParentShowerPfos(const pandora::PfoList *const pLeadingPfoList, pandora::PfoList &parentShowerPfos) const;

    /**
     *  @brief  Starting with provided pfo, find all downstream shower pfos that themselves have daughter pfos
     *
     *  @param  pPfo the address of a pfo
     *  @param  parentShowerPfos to receive the list of parent shower pfos
     */
    void FindParentShowerPfos(const pandora::Pfo *const pLeadiPfo, pandora::PfoList &parentShowerPfos) const;

    /**
     *  @brief  For each parent shower pfo, merge all downstream pfos back into the parent shower
     *
     *  @param  parentShowerPfos the list of parent shower pfos
     */
    void PerformPfoMerges(const pandora::PfoList &parentShowerPfos) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_leadingPfoListName;  ///< The input leading pfo list name (e.g. list of neutrino or testbeam pfos)
};

} // namespace lar_content

#endif // #ifndef LAR_SHOWER_HIERARCHY_MOP_UP_ALGORITHM_H
