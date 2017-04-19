/**
 *  @file   larcontent/include/LArHelpers/LArSVMHelper.h
 *
 *  @brief  Header file for the lar SVM helper class.
 *
 *  $Log: $
 */
#ifndef LAR_SVM_HELPER_H
#define LAR_SVM_HELPER_H 1

#include "larpandoracontent/LArObjects/LArSupportVectorMachine.h"

#include <fstream>
#include <chrono>
#include <ctime>

namespace lar_content
{

/**
 *  @brief  SVMHelper class
 */
class SVMHelper
{
public:
    /**
     *  @brief  Produce a training example with the given features and result
     *
     *  @param  trainingOutputFile the file to which to append the example
     *  @param  featureLists the lists of features
     *
     *  @return success
     */
    template <typename ...TLISTS>
    static pandora::StatusCode ProduceTrainingExample(const std::string &trainingOutputFile, const bool result, TLISTS &&... featureLists);

    /**
     *  @brief  Use the trained SVM to predict the boolean class of an example
     *
     *  @param  sVMachine the support vector machine
     *  @param  featureLists the lists of features
     *
     *  @return the predicted boolean class of the example
     */
    template <typename ...TLISTS>
    static bool Classify(const SupportVectorMachine &sVMachine, TLISTS &&... featureLists);

    /**
     *  @brief  Use the trained SVM to calculate the classification score of an example (>0 means boolean class true)
     *
     *  @param  sVMachine the support vector machine
     *  @param  featureLists the lists of features
     *
     *  @return the classification score
     */
    template <typename ...TLISTS>
    static double CalculateClassificationScore(const SupportVectorMachine &sVMachine, TLISTS &&... featureLists);

    /**
     *  @brief  Calculate the features in a given feature tool vector
     *
     *  @param  featureToolVector the feature tool vector
     *  @param  args arguments to pass to the tool
     *
     *  @return the vector of features
     */
    template <typename ...Ts, typename ...TARGS>
    static SupportVectorMachine::DoubleVector CalculateFeatures(const SVMFeatureToolVector<Ts...> &featureToolVector, TARGS &&... args);

    /**
     *  @brief  Calculate the features of a given derived feature tool type in a feature tool vector
     *
     *  @param  featureToolVector the feature tool vector
     *  @param  args arguments to pass to the tool
     *
     *  @return the vector of features
     */
    template <typename T, typename ...Ts, typename ...TARGS>
    static SupportVectorMachine::DoubleVector CalculateFeaturesOfType(const SVMFeatureToolVector<Ts...> &featureToolVector, TARGS &&... args);

    /**
     *  @brief  Add a feature tool to a vector of feature tools
     *
     *  @param  pFeatureTool the feature tool
     *  @param  featureToolVector the vector to append
     *
     *  @return success
     */
    template <typename ...Ts>
    static pandora::StatusCode AddFeatureToolToVector(pandora::AlgorithmTool *const pFeatureTool, SVMFeatureToolVector<Ts...> &featureToolVector);

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
    template <typename TLIST, typename ...TLISTS>
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
    template <typename TLIST, typename ...TLISTS>
    static SupportVectorMachine::DoubleVector ConcatenateFeatureLists(TLIST &&featureList, TLISTS &&... featureLists);

    /**
     *  @brief  Recursively concatenate vectors of features (terminating method)
     */
    static SupportVectorMachine::DoubleVector ConcatenateFeatureLists();
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename ...TLISTS>
pandora::StatusCode SVMHelper::ProduceTrainingExample(const std::string &trainingOutputFile, const bool result, TLISTS &&... featureLists)
{
    std::ofstream outfile;
    outfile.open(trainingOutputFile, std::ios_base::app); // always append to the output file

    if (!outfile.is_open())
    {
        std::cout << "SVMHelper: could not open file for training examples at " << trainingOutputFile << std::endl;
        return pandora::STATUS_CODE_FAILURE;
    }

    std::string delimiter(",");
    outfile << GetTimestampString() << delimiter;

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, WriteFeaturesToFile(outfile, delimiter, featureLists...));
    outfile << static_cast<int>(result) << '\n';

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename ...TLISTS>
bool SVMHelper::Classify(const SupportVectorMachine &sVMachine, TLISTS &&... featureLists)
{
    return sVMachine.Classify(ConcatenateFeatureLists(std::forward<TLISTS>(featureLists)...));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename ...TLISTS>
double SVMHelper::CalculateClassificationScore(const SupportVectorMachine &sVMachine, TLISTS &&... featureLists)
{
    return sVMachine.CalculateClassificationScore(ConcatenateFeatureLists(std::forward<TLISTS>(featureLists)...));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename ...Ts, typename ...TARGS>
SupportVectorMachine::DoubleVector SVMHelper::CalculateFeatures(const SVMFeatureToolVector<Ts...> &featureToolVector, TARGS &&... args)
{
    SupportVectorMachine::DoubleVector featureVector;

    for (SVMFeatureTool<Ts...> *const pFeatureTool : featureToolVector)
        pFeatureTool->Run(featureVector, std::forward<TARGS>(args)...);

    return featureVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T, typename ...Ts, typename ...TARGS>
SupportVectorMachine::DoubleVector SVMHelper::CalculateFeaturesOfType(const SVMFeatureToolVector<Ts...> &featureToolVector, TARGS &&... args)
{
    using TD = typename std::decay<T>::type;
    SupportVectorMachine::DoubleVector featureVector;

    for (SVMFeatureTool<Ts...> *const pFeatureTool : featureToolVector)
    {
        if (TD *const pCastFeatureTool = dynamic_cast<TD *const>(pFeatureTool))
            pCastFeatureTool->Run(featureVector, std::forward<TARGS>(args)...);
    }

    return featureVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename ...Ts>
pandora::StatusCode SVMHelper::AddFeatureToolToVector(pandora::AlgorithmTool *const pFeatureTool, SVMFeatureToolVector<Ts...> &featureToolVector)
{
    if (SVMFeatureTool<Ts...> *const pCastFeatureTool = dynamic_cast<SVMFeatureTool<Ts...> *const>(pFeatureTool))
    {
        featureToolVector.push_back(pCastFeatureTool);
        return pandora::STATUS_CODE_SUCCESS;
    }

    return pandora::STATUS_CODE_FAILURE;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::string SVMHelper::GetTimestampString()
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

template <typename TLIST, typename ...TLISTS>
inline pandora::StatusCode SVMHelper::WriteFeaturesToFile(std::ofstream &outfile, const std::string &delimiter, TLIST &&featureList, TLISTS &&... featureLists)
{
    static_assert(std::is_same<typename std::decay<TLIST>::type, SupportVectorMachine::DoubleVector>::value, 
        "SVMHelper: Could not write training set example because a passed parameter was not a vector of doubles");
    
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, WriteFeaturesToFileImpl(outfile, delimiter, featureList));
    return WriteFeaturesToFile(outfile, delimiter, featureLists...);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode SVMHelper::WriteFeaturesToFile(std::ofstream &, const std::string &)
{
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TLIST>
pandora::StatusCode SVMHelper::WriteFeaturesToFileImpl(std::ofstream &outfile, const std::string &delimiter, TLIST &&featureList)
{
    for (const double feature : featureList)
        outfile << feature << delimiter;

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TLIST, typename ...TLISTS>
SupportVectorMachine::DoubleVector SVMHelper::ConcatenateFeatureLists(TLIST &&featureList, TLISTS &&... featureLists)
{
    static_assert(std::is_same<typename std::decay<TLIST>::type, SupportVectorMachine::DoubleVector>::value,
        "SVMHelper: Could not concatenate feature lists because one or more lists was not a vector of doubles");
    
    SupportVectorMachine::DoubleVector featureVector;

    for (const double feature : featureList)
        featureVector.push_back(feature);

    SupportVectorMachine::DoubleVector newFeatureVector = ConcatenateFeatureLists(std::forward<TLISTS>(featureLists)...);
    featureVector.insert(featureVector.end(), newFeatureVector.begin(), newFeatureVector.end());
    
    return featureVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline SupportVectorMachine::DoubleVector SVMHelper::ConcatenateFeatureLists()
{
    return SupportVectorMachine::DoubleVector();
}

} // namespace lar_content

#endif // #ifndef LAR_SVM_HELPER_H