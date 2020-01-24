/**
 *  @file   larpandoracontent/MyDlVtxAlgorithm.h
 * 
 *  @brief  Header file for the mydlvtx algorithm class.
 * 
 *  $Log: $
 */
#ifndef LAR_MYDLVTX_ALGORITHM_H
#define LAR_MYDLVTX_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <torch/script.h>

namespace lar_content
{

/**
 *  @brief  MyDlVtxAlgorithm class
 */
class MyDlVtxAlgorithm : public pandora::Algorithm
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
    MyDlVtxAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~MyDlVtxAlgorithm();

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
     *  @param  length, square length of the image
     *
     *  @return The position of the DL vertex for chosen view
     */
    pandora::CartesianVector GetDlVtxForView(const pandora::ClusterList *pClusterList, const std::string &view, 
                             const pandora::CartesianVector &positionInput, const double &length, int &vertReconCount) const;

    /**
     *  @brief  Use Deep Learning on input image to get vertex pixel position for chosen view
     *
     *  @param  out2darr, the image
     *  @param  view, string recording what view the image represents
     *  @param  length, square length of the image
     *
     *  @return The pixel position of the DL vertex for chosen view
     */
    pandora::CartesianVector DeepLearning(const std::vector<std::vector<double>> &out2darr, const std::string &view, const double &length) const;

    // Member variables here
    std::string             m_outputVertexListName;           ///< The name under which to save the output vertex list
    pandora::StringVector   m_inputClusterListNames;          ///< The list of cluster list names
    std::string             m_filePathEnvironmentVariable;    ///< The environment variable providing paths to model files
    std::string             m_modelFileNamePrefix;            ///< The model file name prefix
    std::vector<std::shared_ptr<torch::jit::script::Module>> m_pModule; ///< Vector of Pointers to Torch models.

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::Algorithm *MyDlVtxAlgorithm::Factory::CreateAlgorithm() const
{
    return new MyDlVtxAlgorithm();
}

} // namespace lar_content

#endif // #ifndef LAR_MYDLVTX_ALGORITHM_H
