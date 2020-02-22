/**
 *  @file   larpandoracontent/DLVertexCreationAlgorithm.h
 * 
 *  @brief  Header file for the DLVertexCreation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_DLVERTEXCREATION_ALGORITHM_H
#define LAR_DLVERTEXCREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <torch/script.h>

namespace lar_content
{

/**
 *  @brief  DLVertexCreationAlgorithm class
 */
class DLVertexCreationAlgorithm : public pandora::Algorithm
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
      * @brief Default constructor
      */
    DLVertexCreationAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~DLVertexCreationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode Initialize();

    /**
     *  @brief  Get the DL Vertex position for chosen view
     *
     *  @param  pClusterList, address of the input cluster list for the specific view
     *  @param  view, string recording what view the clusters are of
     *  @param  positionInput, DL vertex position already calculated for chosen view or (0.f, 0.f, 0.f)
     *  @param  lenVecIndex, index for m_lenVec
     *  @param  vertReconCount, count of the reconstructed 2D DL vertices
     *
     *  @return The position of the DL vertex for chosen view
     */
    pandora::CartesianVector GetDLVertexForView(const pandora::ClusterList *pClusterList, const std::string &view, 
                             const pandora::CartesianVector &positionInput, const int &lenVecIndex, int &vertReconCount) const;

    /**
     *  @brief  Use Deep Learning on input image to get vertex pixel position for chosen view
     *
     *  @param  out2dVec, the image
     *  @param  view, string recording what view the image represents
     *  @param  lenVecIndex, index for m_lenVec
     *
     *  @return The pixel position of the DL vertex for chosen view
     */
    pandora::CartesianVector DeepLearning(const std::vector<std::vector<double>> &out2dVec, const std::string &view,
                             const int &lenVecIndex) const;

    /**
     *  @brief  Create the training files.
     *
     *  @param  out2dVec, the image
     *  @param  view, string recording what view the image represents
     *  @param  minx, the minimum x coordinate represented by the image
     *  @param  nstepx, the length of a pixel in the x direction
     *  @param  minz, the minimum z coordinate represented by the image
     *  @param  nstepz, the length of a pixel in the z direction
     */
    void CreateTrainingFiles(const std::vector<std::vector<double>> &out2dVec, const std::string &view,
         const float &minx, const float &nstepx, const float &minz, const float &nstepz) const;

    // Member variables here
    std::string             m_outputVertexListName;           ///< The name under which to save the output vertex list
    pandora::StringVector   m_inputClusterListNames;          ///< The list of cluster list names
    std::string             m_filePathEnvironmentVariable;    ///< The environment variable providing paths to model files
    std::string             m_modelFileNamePrefix;            ///< The model file name prefix
    unsigned int            m_numClusterCaloHitsPar;          ///< The number of cluster calo hits parameter 
    unsigned int            m_npixels;                        ///< The number of pixels, N, in the N*N square image 
    std::vector<float>      m_lenVec;                         ///< Vector of lengths of the images
    bool                    m_trainingSetMode;                ///< Whether to train
    int                     m_trainingLenVecIndex;            ///< Index for m_lenVec to use for output training images
    float                   m_lenBuffer;                      ///< Length of the buffer for scaled images
    std::vector<std::shared_ptr<torch::jit::script::Module>> m_pModule; ///< Vector of Pointers to Torch models
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *DLVertexCreationAlgorithm::Factory::CreateAlgorithm() const
{
    return new DLVertexCreationAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_DLVERTEXCREATION_ALGORITHM_H
