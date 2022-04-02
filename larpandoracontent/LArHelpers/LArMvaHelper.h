/**
 *  @file   larpandoracontent/include/LArHelpers/LArMvaHelper.h
 *
 *  @brief  Header file for the lar mva helper class.
 *
 *  $Log: $
 */
#ifndef LAR_MVA_HELPER_H
#define LAR_MVA_HELPER_H 1

#include "larpandoracontent/LArObjects/LArMvaInterface.h"

#include "Api/PandoraContentApi.h"

#include "Helpers/XmlHelper.h"

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmTool.h"
#include "Pandora/StatusCodes.h"
// BH - commenting now but may be necessary (pandora::StringVector)
#include "Pandora/PandoraInternal.h"
#include "Pandora/StatusCodes.h"

#include <chrono>
#include <ctime>
#include <fstream>

namespace lar_content
{

/**
 *  @brief  MvaFeatureTool class template
 */
template <typename... Ts>
class MvaFeatureTool : public pandora::AlgorithmTool
{
public:
    typedef std::vector<MvaFeatureTool<Ts...> *> FeatureToolVector;
    typedef std::map<std::string, MvaFeatureTool<Ts...> *> FeatureToolMap;

    /**
     *  @brief  Default constructor.
     */
    MvaFeatureTool() = default;

    /**
     *  @brief  Run the algorithm tool
     *
     *  @param  featureVector the vector of features to append
     *  @param  args arguments to pass to the tool
     */
    virtual void Run(MvaTypes::MvaFeatureVector &featureVector, Ts... args) = 0;
    virtual void Run(const std::string featureToolName, MvaTypes::MvaFeatureMap &featureMap, Ts... args){ return; };
};

template <typename... Ts>
using MvaFeatureToolVector = std::vector<MvaFeatureTool<Ts...> *>;

// ?
template <typename... Ts>
using MvaFeatureToolMap = std::map<std::string, MvaFeatureTool<Ts...> *>;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  LArMvaHelper class
 */
class LArMvaHelper
{
public:
    typedef MvaTypes::MvaFeature MvaFeature;
    typedef MvaTypes::MvaFeatureVector MvaFeatureVector;
    typedef std::map<std::string, double> MvaFeatureMap;

    typedef pandora::StringVector StringVector;
    typedef MvaTypes::MvaFeatureMap MvaFeatureMap;
    typedef std::map<std::string, pandora::AlgorithmTool *> AlgorithmToolMap; // idea would be to put this in PandoraInternal.h at some point in PandoraSDK

    /**
     *  @brief  Produce a training example with the given features and result - meant for use with vectors primarily
     *
     *  @param  trainingOutputFile the file to which to append the example
     *  @param  featureLists the lists of features
     *
     *  @return success
     */
    template <typename... TLISTS>
    static pandora::StatusCode ProduceTrainingExample(const std::string &trainingOutputFile, const bool result, TLISTS &&... featureLists);

    /**
     *  @brief  Produce a training example with the given features and result - meant for use with maps of features
     *
     *  @param  trainingOutputFile the file to which to append the example
     *  @param  featureMap the map of features
     *  @param  featureToolOrder the vector of strings corresponding to ordered list of keys
     *  @param  featureLists the lists of features - in this mode forced to be not more than parameter in the pack
     */
    template <typename TLIST>
    static pandora::StatusCode ProduceTrainingExample(const std::string &trainingOutputFile, const bool result, const StringVector &featureToolOrder, TLIST && featureList);

    /**
     *  @brief  Use the trained classifier to predict the boolean class of an example
     *
     *  @param  classifier the classifier
     *  @param  featureLists the lists of features
     *
     *  @return the predicted boolean class of the example
     */
    template <typename... TLISTS>
    static bool Classify(const MvaInterface &classifier, TLISTS &&... featureLists);

    // BH - for now let's make this NOT templated. You need to give it an already concatenated map...
    // TODO: make more widely useful via additional templating...
    /**
     *  @brief  Use the trained classifier to predict the boolean class of an example
     *
     *  @param  classifier the classifier
     *  @param  featureMap the map of features
     *  @param  featureToolOrder the vector of strings corresponding to ordered list of keys
     *
     *  @return success
     */
    template <typename... TLISTS>
    static bool Classify(const MvaInterface &classifier, const MvaFeatureMap &featureMap, const StringVector &featureToolOrder);

    /**
     *  @brief  Use the trained classifer to calculate the classification score of an example (>0 means boolean class true)
     *
     *  @param  classifier the classifier
     *  @param  featureLists the lists of features
     *
     *  @return the classification score
     */
    template <typename... TLISTS>
    static double CalculateClassificationScore(const MvaInterface &classifier, TLISTS &&... featureLists);

    /**
     *  @brief  Use the trained mva to calculate a classification probability for an example
     *
     *  @param  classifier the classifier
     *  @param  featureLists the lists of features
     *
     *  @return the classification probability
     */
    template <typename... TLISTS>
    static double CalculateProbability(const MvaInterface &classifier, TLISTS &&... featureLists);

    // BH - for now let's make this NOT templated. You need to give it an already concatenated map...
    // TODO: make more widely useful via additional templating...
    /**
     *  @brief  Use the trained mva to calculate a classification probability for an example
     *
     *  @param  classifier the classifier
     *  @param  featureMap the map of features
     *  @param  featureToolOrder the vector of strings corresponding to ordered list of keys
     *
     *  @return the classification probability
     */
    template <typename... TLISTS>
    static double CalculateProbability(const MvaInterface &classifier, const MvaFeatureMap &featureMap, const StringVector &featureToolOrder );

    /**
     *  @brief  Calculate the features in a given feature tool vector
     *
     *  @param  featureToolVector the feature tool vector
     *  @param  args arguments to pass to the tool
     *
     *  @return the vector of features
     */
    template <typename... Ts, typename... TARGS>
    static MvaFeatureVector CalculateFeatures(const MvaFeatureToolVector<Ts...> &featureToolVector, TARGS &&... args);

    /**
     *  @brief  Calculate the features in a given feature tool map, and fill a MvaFeatureMap and MvaFeatureVector
     *
     *  @param  featureToolMap the feature tool map
     *  @param  featureToolOrder vector of strings of the ordered keys
     *  @param  args arguments to pass to the tool
     *
     *  @return the map of features
     */
    template <typename... Ts, typename... TARGS>
    static MvaFeatureMap CalculateFeatures(const MvaFeatureToolMap<Ts...> &featureToolMap, const StringVector &featureToolOrder, TARGS &&... args);

    /**
     *  @brief  Calculate the features of a given derived feature tool type in a feature tool vector
     *
     *  @param  featureToolVector the feature tool vector
     *  @param  args arguments to pass to the tool
     *
     *  @return the vector of features
     */
    template <typename T, typename... Ts, typename... TARGS>
    static MvaFeatureVector CalculateFeaturesOfType(const MvaFeatureToolVector<Ts...> &featureToolVector, TARGS &&... args);

    /**
     *  @brief  Add a feature tool to a vector of feature tools
     *
     *  @param  pFeatureTool the feature tool
     *  @param  featureToolVector the vector to append
     *
     *  @return success
     */
    template <typename... Ts>
    static pandora::StatusCode AddFeatureToolToVector(pandora::AlgorithmTool *const pFeatureTool, MvaFeatureToolVector<Ts...> &featureToolVector);

    /**
     *  @brief  Add a feature tool to a map of feature tools
     * 
     *  @param  pFeatureTool the feature tool
     *  @param  pFeatureToolName the name of the feature tool
     *  @param  featureToolMap the map to append
     * 
     *  @return success
     */
    template <typename... Ts>
    static pandora::StatusCode AddFeatureToolToMap(pandora::AlgorithmTool *const pFeatureTool, std::string pFeatureToolName, MvaFeatureToolMap<Ts...> &featureToolMap);

    /**
     *  @brief  Process a list of algorithms tools in an xml file, using a map. Idea is for this to go to XmlHelper in PandoraSDK eventually as an overload to ProcessAlgorithmToolList
     * 
     *  @param  algorithm the parent algorithm calling this function
     *  @param  xmlHandle the relevant xml handle
     *  @param  listName the name of the algorithm tool list
     *  @param  algorithmToolMap to receive the vector of addresses of the algorithm tool instances, but also keep the name
     */
    static pandora::StatusCode ProcessAlgorithmToolListToMap(const pandora::Algorithm &algorithm, const pandora::TiXmlHandle &xmlHandle, const std::string &listName,
							     StringVector &algorithToolNameVector, AlgorithmToolMap &algorithmToolMap);

private:
    /**
     *  @brief  Get a timestamp string for this point in time
     *
     *  @return a timestamp string
     */
    static std::string GetTimestampString();

    /**
     *  @brief  Recursively write the features of the given lists to file
     *
     *  @param  outfile the std::ofstream object to use
     *  @param  delimiter the delimiter string
     *  @param  featureList a list of features to write
     *  @param  featureLists optional further lists of features to write
     *
     *  @return success
     */
    template <typename TLIST, typename... TLISTS>
    static pandora::StatusCode WriteFeaturesToFile(std::ofstream &outfile, const std::string &delimiter, TLIST &&featureList, TLISTS &&... featureLists);

    /**
     *  @brief  Recursively write the features of the given lists to file (terminating method)
     *
     *  @return success
     */
    static pandora::StatusCode WriteFeaturesToFile(std::ofstream &, const std::string &);

    /**
     *  @brief  Write the features of the given list to file (implementation method)
     *
     *  @param  outfile the std::ofstream object to use
     *  @param  delimiter the delimiter string
     *  @param  featureList a list of features to write
     *
     *  @return success
     */
    template <typename TLIST>
    static pandora::StatusCode WriteFeaturesToFileImpl(std::ofstream &outfile, const std::string &delimiter, TLIST &&featureList);

    /**
     *  @brief  Recursively concatenate vectors of features
     *
     *  @param  featureList a list of features
     *  @param  featureLists optional further lists of features
     *
     *  @return the concatenated vector of features
     */
    template <typename TLIST, typename... TLISTS>
    static MvaFeatureVector ConcatenateFeatureLists(TLIST &&featureList, TLISTS &&... featureLists);

    /**
     *  @brief  Recursively concatenate vectors of features (terminating method)
     */
    static MvaFeatureVector ConcatenateFeatureLists();
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename... TLISTS>
pandora::StatusCode LArMvaHelper::ProduceTrainingExample(const std::string &trainingOutputFile, const bool result, TLISTS &&... featureLists)
{
    std::ofstream outfile;
    outfile.open(trainingOutputFile, std::ios_base::app); // always append to the output file

    if (!outfile.is_open())
    {
        std::cout << "LArMvaHelper: could not open file for training examples at " << trainingOutputFile << std::endl;
        return pandora::STATUS_CODE_FAILURE;
    }

    std::string delimiter(",");
    outfile << GetTimestampString() << delimiter;

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, WriteFeaturesToFile(outfile, delimiter, featureLists...));
    outfile << static_cast<int>(result) << '\n';

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TLIST>
pandora::StatusCode LArMvaHelper::ProduceTrainingExample(const std::string &trainingOutputFile, const bool result,
							 const LArMvaHelper::StringVector &featureToolOrder, TLIST && featureList)
{
    static_assert(std::is_same<typename std::decay<TLIST>::type, LArMvaHelper::MvaFeatureMap>::value,
		  "LArMvaHelper: Could not Produce Training Example as the list was not a map of MvaFeatures, yet this structure of parameters is being used.");

    // Make a feature vector from the map and calculate the features
    LArMvaHelper::MvaFeatureVector featureVector;

    for ( auto const& pFeatureToolName : featureToolOrder ) {
        if ( featureList.find( pFeatureToolName ) == featureList.end() )
	    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
	featureVector.push_back( featureList.at( pFeatureToolName ) );
    }

    return ProduceTrainingExample( trainingOutputFile, result, featureVector );
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename... TLISTS>
bool LArMvaHelper::Classify(const MvaInterface &classifier, TLISTS &&... featureLists)
{
    return classifier.Classify(ConcatenateFeatureLists(std::forward<TLISTS>(featureLists)...));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename... TLISTS>
bool LArMvaHelper::Classify(const MvaInterface &classifier, const LArMvaHelper::MvaFeatureMap &featureMap, const LArMvaHelper::StringVector &featureToolOrder )
{
    // Make a feature vector from the map and calculate the features
    LArMvaHelper::MvaFeatureVector featureVector;

    for ( auto const& pFeatureToolName : featureToolOrder ) {
        if ( featureMap.find( pFeatureToolName ) == featureMap.end() )
	    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
	featureVector.push_back( featureMap.at( pFeatureToolName ) );
    }

    return Classify(classifier, featureVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename... TLISTS>
double LArMvaHelper::CalculateClassificationScore(const MvaInterface &classifier, TLISTS &&... featureLists)
{
    return classifier.CalculateClassificationScore(ConcatenateFeatureLists(std::forward<TLISTS>(featureLists)...));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename... TLISTS>
double LArMvaHelper::CalculateProbability(const MvaInterface &classifier, TLISTS &&... featureLists)
{
    return classifier.CalculateProbability(ConcatenateFeatureLists(std::forward<TLISTS>(featureLists)...));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename... TLISTS>
double LArMvaHelper::CalculateProbability(const MvaInterface &classifier, const LArMvaHelper::MvaFeatureMap &featureMap, const LArMvaHelper::StringVector &featureToolOrder )
{
    // Make a feature vector from the map and calculate the features
    LArMvaHelper::MvaFeatureVector featureVector;

    for ( auto const& pFeatureToolName : featureToolOrder ) {
        if ( featureMap.find( pFeatureToolName ) == featureMap.end() )
	    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
	featureVector.push_back( featureMap.at( pFeatureToolName ) );
    }

    return CalculateProbability(classifier, featureVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename... Ts, typename... TARGS>
LArMvaHelper::MvaFeatureVector LArMvaHelper::CalculateFeatures(const MvaFeatureToolVector<Ts...> &featureToolVector, TARGS &&... args)
{
    LArMvaHelper::MvaFeatureVector featureVector;

    for (MvaFeatureTool<Ts...> *const pFeatureTool : featureToolVector)
        pFeatureTool->Run(featureVector, std::forward<TARGS>(args)...);

    return featureVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename... Ts, typename... TARGS>
LArMvaHelper::MvaFeatureMap LArMvaHelper::CalculateFeatures(const MvaFeatureToolMap<Ts...> &featureToolMap, const LArMvaHelper::StringVector &featureToolOrder, TARGS &&... args)
{
  LArMvaHelper::MvaFeatureMap featureMap;

    for ( auto const& pFeatureToolName : featureToolOrder ) {
      if ( featureToolMap.find( pFeatureToolName ) == featureToolMap.end() )
	throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
      featureToolMap.at(pFeatureToolName)->Run(pFeatureToolName, featureMap, std::forward<TARGS>(args)...);
    }

    return featureMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T, typename... Ts, typename... TARGS>
LArMvaHelper::MvaFeatureVector LArMvaHelper::CalculateFeaturesOfType(const MvaFeatureToolVector<Ts...> &featureToolVector, TARGS &&... args)
{
    using TD = typename std::decay<T>::type;
    LArMvaHelper::MvaFeatureVector featureVector;

    for (MvaFeatureTool<Ts...> *const pFeatureTool : featureToolVector)
    {
        if (TD *const pCastFeatureTool = dynamic_cast<TD *const>(pFeatureTool))
            pCastFeatureTool->Run(featureVector, std::forward<TARGS>(args)...);
    }

    return featureVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename... Ts>
pandora::StatusCode LArMvaHelper::AddFeatureToolToVector(pandora::AlgorithmTool *const pFeatureTool, MvaFeatureToolVector<Ts...> &featureToolVector)
{
    if (MvaFeatureTool<Ts...> *const pCastFeatureTool = dynamic_cast<MvaFeatureTool<Ts...> *const>(pFeatureTool))
    {
        featureToolVector.push_back(pCastFeatureTool);
        return pandora::STATUS_CODE_SUCCESS;
    }

    return pandora::STATUS_CODE_FAILURE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename... Ts>
pandora::StatusCode LArMvaHelper::AddFeatureToolToMap(pandora::AlgorithmTool *const pFeatureTool, std::string pFeatureToolName, MvaFeatureToolMap<Ts...> &featureToolMap)
{
    if (MvaFeatureTool<Ts...> *const pCastFeatureTool = dynamic_cast<MvaFeatureTool<Ts...> *const>(pFeatureTool))
    {
        featureToolMap[ pFeatureToolName ] = pCastFeatureTool;
        return pandora::STATUS_CODE_SUCCESS;
    }

    return pandora::STATUS_CODE_FAILURE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

/*
pandora::StatusCode LArMvaHelper::ProcessAlgorithmToolListToMap(const pandora::Algorithm &algorithm, const pandora::TiXmlHandle &xmlHandle, const std::string &listName,
								AlgorithmToolMap &algorithmToolMap)
{
    if ("algorithm" != xmlHandle.ToNode()->ValueStr())
        return pandora::STATUS_CODE_NOT_ALLOWED;

    const pandora::TiXmlHandle algorithmListHandle = pandora::TiXmlHandle(xmlHandle.FirstChild(listName).Element());

    std::cout << "ALG NAMES: ";

    for (pandora::TiXmlElement *pXmlElement = algorithmListHandle.FirstChild("tool").Element(); nullptr != pXmlElement;
	 pXmlElement = pXmlElement->NextSiblingElement("tool"))
    {
        pandora::AlgorithmTool *pAlgorithmTool(nullptr);
	PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateAlgorithmTool(algorithm, pXmlElement, pAlgorithmTool));
	std::cout << pXmlElement->Attribute("type") << " ";
	algorithmToolMap[ pXmlElement->Attribute("type") ] = pAlgorithmTool;
    }

    std::cout << std::endl;

    return pandora::STATUS_CODE_SUCCESS;
}
*/

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::string LArMvaHelper::GetTimestampString()
{
    std::time_t timestampNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    struct tm *pTimeInfo(NULL);
    char buffer[80];

    pTimeInfo = localtime(&timestampNow);
    strftime(buffer, 80, "%x_%X", pTimeInfo);

    std::string timeString(buffer);

    if (!timeString.empty() && timeString.back() == '\n') // last char is always a newline
        timeString.pop_back();

    return timeString;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TLIST, typename... TLISTS>
inline pandora::StatusCode LArMvaHelper::WriteFeaturesToFile(
    std::ofstream &outfile, const std::string &delimiter, TLIST &&featureList, TLISTS &&... featureLists)
{
    static_assert(std::is_same<typename std::decay<TLIST>::type, LArMvaHelper::MvaFeatureVector>::value,
        "LArMvaHelper: Could not write training set example because a passed parameter was not a vector of MvaFeatures");

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, WriteFeaturesToFileImpl(outfile, delimiter, featureList));
    return WriteFeaturesToFile(outfile, delimiter, featureLists...);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode LArMvaHelper::WriteFeaturesToFile(std::ofstream &, const std::string &)
{
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TLIST>
pandora::StatusCode LArMvaHelper::WriteFeaturesToFileImpl(std::ofstream &outfile, const std::string &delimiter, TLIST &&featureList)
{
    for (const MvaFeature &feature : featureList)
        outfile << feature.Get() << delimiter;

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TLIST, typename... TLISTS>
LArMvaHelper::MvaFeatureVector LArMvaHelper::ConcatenateFeatureLists(TLIST &&featureList, TLISTS &&... featureLists)
{
    static_assert(std::is_same<typename std::decay<TLIST>::type, LArMvaHelper::MvaFeatureVector>::value,
        "LArMvaHelper: Could not concatenate feature lists because one or more lists was not a vector of MvaFeatures");

    LArMvaHelper::MvaFeatureVector featureVector;

    for (const MvaFeature &feature : featureList)
        featureVector.push_back(feature);

    LArMvaHelper::MvaFeatureVector newFeatureVector = ConcatenateFeatureLists(std::forward<TLISTS>(featureLists)...);
    featureVector.insert(featureVector.end(), newFeatureVector.begin(), newFeatureVector.end());

    return featureVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArMvaHelper::MvaFeatureVector LArMvaHelper::ConcatenateFeatureLists()
{
    return LArMvaHelper::MvaFeatureVector();
}

} // namespace lar_content

#endif // #ifndef LAR_MVA_HELPER_H
