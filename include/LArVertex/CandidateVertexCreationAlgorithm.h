/**
 *  @file   LArContent/include/LArVertex/CandidateVertexCreationAlgorithm.h
 * 
 *  @brief  Header file for the candidate vertex creation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H
#define LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  CandidateVertexCreationAlgorithm::Algorithm class
 */
class CandidateVertexCreationAlgorithm : public pandora::Algorithm
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

    /**
     *  @brief  Default constructor
     */
    CandidateVertexCreationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string     m_outputVertexListName;         ///< The name under which to save the output vertex list
    bool            m_replaceCurrentVertexList;     ///< Whether to replace the current vertex list with the output list
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *CandidateVertexCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new CandidateVertexCreationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_CANDIDATE_VERTEX_CREATION_ALGORITHM_H
