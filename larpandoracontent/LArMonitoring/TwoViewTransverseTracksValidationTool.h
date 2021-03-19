/**
 *  @file   larpandoracontent/LArMonitoring/TwoViewTransverseTracksValidationTool.h
 *
 *  @brief  Header file for the two view transverse tracks validation tool
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_TRANSVERSE_TRACKS_VALIDATION_TOOL_H
#define TWO_VIEW_TRANSVERSE_TRACKS_VALIDATION_TOOL_H 1

#include "PandoraMonitoringApi.h"
#include "Objects/Cluster.h"
#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArDiscreteProbabilityVector.h"
#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"

#include <any>
#include <memory>
#include <type_traits>

namespace lar_content
{

/**
 *  @brief  TwoViewTransverseTracksValidationTool class
 */
class TwoViewTransverseTracksValidationTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Default constructor
     */
    TwoViewTransverseTracksValidationTool();

    /**
     *  @brief  Default destructor
     */
    ~TwoViewTransverseTracksValidationTool();

    /**
     *  @brief  Run the validation tool
     *
     *  @param  pCluster1 the first cluster
     *  @param  pCluster2 the second cluster
     *  @param  discreteProbabilityVector1 the resampled charge profile for the first cluster
     *  @param  discreteProbabilityVector2 the resampled charge profile for the second cluster
     *  @param  overlapResult the result containing the matching metrics
     *
     *  @return processing success
     */
    bool Run(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2, const DiscreteProbabilityVector &discreteProbabilityVector1,
        const DiscreteProbabilityVector &discreteProbabilityVector2, const TwoViewTransverseOverlapResult &overlapResult);

private:
    /**
     *  @brief  AbstractDatum class
     */
    class AbstractDatum
    {
    public:
        /**
         *  @brief  Add data to the data structure.  Only functional when the underlying structure is a vector
         *
         *  @param  T the data to be added to the data structure
         */
        template <typename T>
        void emplace_back(const T &t);

    private:
        /**
         *  @brief  Underlying pure virtual implementation of emplace_back.   
         *
         *  @param  value the data to be added to the data structure
         */
        virtual void emplace_back_impl(const std::any &value) = 0;
    };

    /**
     *  @brief  TreeDatum class
     */
    template <typename T>
    class TreeDatum final : public AbstractDatum
    {
    public:
        /**
         *  @brief  Constructor
         *  @param  t the data to be saved to the class
         *  @param  pTool the tool that requested the data be saved
         *  @param  datumName the label for the data
         */
        TreeDatum(const T &t, const TwoViewTransverseTracksValidationTool *const pTool, const std::string datumName);

        /**
         *  @brief  Overriding implementation of emplace_back.   
         *
         *  @param  value the data to be added to the data structure
         */
        void emplace_back_impl(const std::any &value) override;

    private:
        T m_Datum;      ///< the datum that is destined for a TTree
    };

    /**
     *  @brief  is_std_vector failing type trait struct
     */
    template <typename T>
    struct is_std_vector final : std::false_type
    {
    };

    /**
     *  @brief  is_std_vector true type trait struct
     */
    template <typename... T>
    struct is_std_vector<std::vector<T...>> final : std::true_type
    {
    };

    typedef std::map<std::string, std::unique_ptr<AbstractDatum>> TreeDataBox;

    /**
     *  @brief  Store some data in a box and register it with the pandora monitoring api for writing to a TTree 
     *
     *  @param  t the data to be saved
     *  @param  datumName the data label
     *  @param  treeDataBox the box where the data is stored before TTree writing
     */
    template <typename T>
    void StoreAndRegisterDatum(const T &t, const std::string datumName, TreeDataBox &treeDataBox);

    /**
     *  @brief  Register some data with the pandora monitoring api for writing to a TTree 
     *
     *  @param  t the data to be saved
     *  @param  datumName the data label
     */
    template <typename T>
    void RegisterTreeDatum(T &t, const std::string datumName) const;

    /**
     *  @brief  Convenience function for collecting and registering truth information with the pandora monitoring api for writing to a TTree 
     *
     *  @param  pCluster the cluster whose truth information needs to be stored
     *  @param  treeDataBox the box where the truth information will be stored
     *  @param  clusterName a cluster label, needed for TTree branch naming
     */
    void CollectTruthInformation(const pandora::Cluster *const pCluster, TreeDataBox &treeDataBox, std::string clusterName);

    /**
     *  @brief  Convenience function for collecting and registering charge profile information with the pandora monitoring api for writing to a TTree 
     *
     *  @param  discreteProbabilityVector the charge profile whose information needs to be stored 
     *  @param  treeDataBox the box where the charge profile information will be stored
     *  @param  clusterName a cluster label, needed for TTree branch naming
     */
    void CollectChargeProfileInformation(const DiscreteProbabilityVector &discreteProbabilityVector, TreeDataBox &treeDataBox, std::string clusterName);

    /**
     *  @brief  Check whether two clusters were created by the same true particle 
     *
     * @param  pCluster1 the first cluster
     * @param  pCluster2 the second cluster
     */
    bool IsSameMCParticle(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_treeName;       ///< name of output tree
    std::string m_outputFileName; ///< name of output file
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TwoViewTransverseTracksValidationTool::AbstractDatum::emplace_back(const T &t)
{
    emplace_back_impl(std::any(t));
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
TwoViewTransverseTracksValidationTool::TreeDatum<T>::TreeDatum(
    const T &t, const TwoViewTransverseTracksValidationTool *const pTool, const std::string datumName) :
    m_Datum(t)
{
    pTool->RegisterTreeDatum(m_Datum, datumName);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TwoViewTransverseTracksValidationTool::TreeDatum<T>::emplace_back_impl(const std::any &value)
{
    if constexpr (is_std_vector<T>::value)
    {
        m_Datum.emplace_back((std::any_cast<typename T::value_type>(value)));
    }
    else
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_ALLOWED);
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TwoViewTransverseTracksValidationTool::StoreAndRegisterDatum(const T &t, const std::string datumName, TreeDataBox &treeDataBox)
{
    auto [iterator, inserted] = treeDataBox.try_emplace(datumName, std::make_unique<TreeDatum<T>>(t, this, datumName));
    if (!inserted)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
void TwoViewTransverseTracksValidationTool::RegisterTreeDatum(T &t, const std::string datumName) const
{
    if constexpr (is_std_vector<T>::value)
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), datumName.c_str(), &t));
    }
    else
    {
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), datumName.c_str(), t));
    }
    return;
}

} // namespace lar_content

#endif // #ifndef TWO_VIEW_TRANSVERSE_TRACKS_VALIDATION_TOOL_H
