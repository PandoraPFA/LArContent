/**
 *  @file   ExampleContent/include/ExampleAlgorithms/ExampleDirectionAlgorithm.h
 * 
 *  @brief  Header file for the access lists algorithm class.
 * 
 *  $Log: $
 */
#ifndef EXAMPLE_DIRECTION_ALGORITHM_H
#define EXAMPLE_DIRECTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "larpandoracontent/LArControlFlow/TrackDirectionTool.h"

namespace lar_content
{

/**
 *  @brief  ExampleDirectionAlgorithm class
 */
class ExampleDirectionAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Factory class for instantiating algorithm
     */
    class Factory : public pandora::AlgorithmFactory
    {
    public:
        pandora::Algorithm *CreateAlgorithm() const;
    };

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string             m_clusterListName;         ///< The requested calo hit list name
    std::string             m_pfoListName;         ///< The requested calo hit list name
    TrackDirectionTool      *m_pTrackDirectionTool;
    typedef std::vector<TrackDirectionTool*> TrackDirectionToolVector;
    TrackDirectionToolVector           m_trackDirectionToolVector;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *ExampleDirectionAlgorithm::Factory::CreateAlgorithm() const
{
    return new ExampleDirectionAlgorithm();
}

} // namespace lar_content

#endif // #ifndef EXAMPLE_DIRECTION_ALGORITHM_H

