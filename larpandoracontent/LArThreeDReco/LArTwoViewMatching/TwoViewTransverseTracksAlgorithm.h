/**
 *  @file   larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewTransverseTracksAlgorithm.h
 *
 *  @brief  Header file for the two view transverse tracks algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TWO_VIEW_TRANSVERSE_TRACKS_ALGORITHM_H
#define LAR_TWO_VIEW_TRANSVERSE_TRACKS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArThreeDReco/LArThreeDBase/TwoViewTrackMatchingAlgorithm.h"

namespace lar_content
{

class TransverseMatrixTool;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TwoViewTransverseTracksAlgorithm class
 */
class TwoViewTransverseTracksAlgorithm : public TwoViewTrackMatchingAlgorithm<float>
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewTransverseTracksAlgorithm();

private:
    void CalculateOverlapResult(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2);
    void ExamineOverlapContainer();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::vector<TransverseMatrixTool*> MatrixToolVector;
    MatrixToolVector            m_algorithmToolVector;      ///< The algorithm tool vector

    unsigned int                m_nMaxMatrixToolRepeats;    ///< The maximum number of repeat loops over matrix tools
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TransverseMatrixTool class
 */
class TransverseMatrixTool : public pandora::AlgorithmTool
{
public:
    typedef TwoViewTransverseTracksAlgorithm::MatrixType MatrixType;
    typedef std::vector<MatrixType::ElementList::const_iterator> IteratorList;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  overlapMatrix the overlap matrix
     *
     *  @return whether changes have been made by the tool
     */
    virtual bool Run(TwoViewTransverseTracksAlgorithm *const pAlgorithm, MatrixType &overlapMatrix) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_TWO_VIEW_TRANSVERSE_TRACKS_ALGORITHM_H
