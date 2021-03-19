/**
 *  @file   larpandoracontent/LArMonitoring/TwoViewTransverseTracksValidationTool.h
 *
 *  @brief  Header file for the two view transverse tracks validation tool
 *
 *  $Log: $
 */
#ifndef TWO_VIEW_TRANSVERSE_TRACKS_VALIDATION_TOOL_H
#define TWO_VIEW_TRANSVERSE_TRACKS_VALIDATION_TOOL_H 1

#include <any>
#include <memory>
#include <type_traits>

#include "PandoraMonitoringApi.h"

#include "Objects/Cluster.h"

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArObjects/LArDiscreteProbabilityVector.h"
#include "larpandoracontent/LArObjects/LArTrackTwoViewOverlapResult.h"
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

    ~TwoViewTransverseTracksValidationTool();

    bool Run(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        const DiscreteProbabilityVector &discreteProbabilityVector1,
        const DiscreteProbabilityVector &discreteProbabilityVector2,
        const TwoViewTransverseOverlapResult &overlapResult);

private:
    class AbstractDatum
    {
        public:
            template <typename T>
            void emplace_back(const T &t);
        private:
            virtual void emplace_back_impl(const std::any &value) = 0;
    };

    template <typename T> 
    class TreeDatum : public AbstractDatum
    {
        public:
            TreeDatum(const T &t, const TwoViewTransverseTracksValidationTool *const pTool, const std::string datumName);

            void emplace_back_impl(const std::any &value) override;

        private:
            T m_Datum;
    };

    template<typename U> struct is_std_vector final : std::false_type {};
    template<typename... U> struct is_std_vector<std::vector<U...> > final : std::true_type {};

    typedef std::map<std::string, std::unique_ptr<AbstractDatum> > TreeDataBox;

    template <typename T>
    void StoreAndRegisterDatum(const T &t, const std::string datumName, TreeDataBox &treeDataBox);

    template <typename T>
    void RegisterTreeDatum(T &t, const std::string datumName) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_treeName;            ///< name of output tree
    std::string m_outputFileName;            ///< name of output file
};

template <typename T>
void TwoViewTransverseTracksValidationTool::AbstractDatum::emplace_back(const T &t)
{
    emplace_back_impl(std::any(t));
}

template <typename T>
TwoViewTransverseTracksValidationTool::TreeDatum<T>::TreeDatum(const T &t, const TwoViewTransverseTracksValidationTool *const pTool, 
    const std::string datumName) :
    m_Datum(t)
{
    pTool->RegisterTreeDatum(m_Datum, datumName);
}

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

template <typename T>
void TwoViewTransverseTracksValidationTool::StoreAndRegisterDatum(const T &t, const std::string datumName, TreeDataBox &treeDataBox)
{
    auto [iterator, inserted] = treeDataBox.try_emplace(datumName, (new TreeDatum<T>(t, this, datumName)));
    if (!inserted)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
    return;
}


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
}

} // namespace lar_content

#endif // #ifndef TWO_VIEW_TRANSVERSE_TRACKS_VALIDATION_TOOL_H
