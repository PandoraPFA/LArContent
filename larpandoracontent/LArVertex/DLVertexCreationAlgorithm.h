/**
 *  @file   larpandoracontent/LArVertex/DLVertexCreationAlgorithm.h
 * 
 *  @brief  Header file for the DLVertexCreation algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_DLVERTEXCREATION_ALGORITHM_H
#define LAR_DLVERTEXCREATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <memory>

namespace torch
{
    namespace jit
    {
        namespace script
        {
            struct Module;
        }
    }
}

namespace lar_content
{

/**
 *  @brief  DLVertexCreationAlgorithm class
 */
class DLVertexCreationAlgorithm : public pandora::Algorithm
{
public:
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

    typedef std::vector<std::vector<double>> TwoDImage;
    typedef std::vector<double> DoubleVector;

    /**
     *  @brief  Get the DL Vertex position for chosen view
     *
     *  @param  pClusterList, address of the input cluster list for the specific view
     *  @param  view, pandora::HitType recording what view the clusters are of
     *  @param  positionInput, DL vertex position already calculated for chosen view or (0.f, 0.f, 0.f)
     *  @param  imgLenVecIndex, index for m_imgLenVec
     *  @param  vertReconCount, count of the reconstructed 2D DL vertices
     *  @param  ssBuf, buffer for training file writing
     *
     *  @return The position of the DL vertex for chosen view
     */
    pandora::CartesianVector GetDLVertexForView(const pandora::ClusterList *const pClusterList, const pandora::HitType &view, 
        const pandora::CartesianVector &positionInput, const int imgLenVecIndex, unsigned int &vertReconCount,
        std::stringstream ssBuf[6]) const;

    /**
     *  @brief  Use Deep Learning on input image to get vertex pixel position for chosen view
     *
     *  @param  out2dVec, the image
     *  @param  view, pandora::HitType recording what view the image represents
     *  @param  imgLenVecIndex, index for m_imgLenVec
     *
     *  @return The pixel position of the DL vertex for chosen view
     */
    pandora::CartesianVector DeepLearning(const TwoDImage &out2dVec, const pandora::HitType &view,
        const int imgLenVecIndex) const;

    /**
     *  @brief  Create the training files
     *
     *  @param  out2dVec, the image
     *  @param  view, pandora::HitType recording what view the image represents
     *  @param  minx, the minimum x coordinate represented by the image
     *  @param  nstepx, the length of a pixel in the x direction
     *  @param  minz, the minimum z coordinate represented by the image
     *  @param  nstepz, the length of a pixel in the z direction
     *  @param  ssBuf, buffer for training file writing
     *
     *  @return A flag recording if the end of the function was reached
     */
    int CreateTrainingFiles(const TwoDImage &out2dVec, const pandora::HitType &view,
        const float minx, const float nstepx, const float minz, const float nstepz, std::stringstream ssBuf[6]) const;

    /**
     *  @brief  Write the training files
     *
     *  @param  ssBuf, buffer for training file writing
     */
    void WriteTrainingFiles(std::stringstream ssBuf[6]) const;

    /**
     *  @brief  Check if a 3D position is in the detector 
     *
     *  @param  position3D, the 3D position to check
     *
     *  @return A flag recording if the 3D position is in the detector 
     */
    bool DetectorCheck(const pandora::CartesianVector &position3D) const;

    /**
     *  @brief  Check if the event view passes conditions
     *
     *  @param  pClusterList, address of the input cluster list for the specific view
     *
     *  @return A flag recording if the event view passed the conditions
     */
    bool EventViewCheck(const pandora::ClusterList *const pClusterList) const;

    /**
     *  @brief  Check if the event passes training conditions
     *
     *  @return A flag recording if the event passed the training conditions
     */
    bool TrainEventCheck() const;

    // Member variables here
    std::string             m_outputVertexListName;           ///< The name under which to save the output vertex list
    pandora::StringVector   m_inputClusterListNames;          ///< The list of cluster list names
    std::string             m_filePathEnvironmentVariable;    ///< The environment variable providing paths to model files
    std::string             m_modelFileNamePrefix;            ///< The model file name prefix
    unsigned int            m_numClusterCaloHitsPar;          ///< The number of cluster calo hits parameter 
    unsigned int            m_npixels;                        ///< The number of pixels, N, in the N*N square image 
    pandora::FloatVector    m_imgLenVec;                      ///< Vector of lengths of the images
    bool                    m_trainingSetMode;                ///< Whether to train
    int                     m_trainingImgLenVecIndex;         ///< Index for m_imgLenVec to use for output training images
    float                   m_lenBuffer;                      ///< Length of the buffer for scaled images
    unsigned int            m_numViews;                       ///< Number of available 2D views
    std::string             m_trainingDataFileName;           ///< The name of the training data output file
    std::string             m_trainingLabelsFileName;         ///< The name of the training labels output file
    float                   m_vertexXCorrection;              ///< The correction to the x-coordinate of the MC vertex position
    double                  m_hitWidthZ;                      ///< The width of the hits in the z-direction
    std::vector<std::shared_ptr<torch::jit::script::Module>> m_vecPModule; ///< Vector of Pointers to Torch models
};

} // namespace lar_content

#endif // #ifndef LAR_DLVERTEXCREATION_ALGORITHM_H
