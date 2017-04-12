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
private:
    /**
     *  @brief  Get a timestamp string for this point in time
     *
     *  @return a timestamp string
     */
    static std::string GetTimestampString();

    /**
     *  @brief  Add all the feature tools of a given return type T (int, float, double) to a given feature tool map
     *
     *  @param  pFeatureTool address of the feature tool
     *  @param  featureToolMap the feature tool map to append
     *
     *  @return whether the feature tool was added with this T
     */
    template <typename T, typename TALG, typename ...Ts>
    static bool AddFeatureToolToMapImpl(pandora::AlgorithmTool *const pFeatureTool, SVMFeatureToolMap<TALG, Ts...> &featureToolMap);
    
    /**
     *  @brief  Calculate the vector of features from a set of algorithm tools matching type T. 
     * 
     *  @param  featureToolMap the feature tool map
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  args arguments to pass to the algorithm tool
     * 
     *  @return the vector of features
     */
    template <typename T, typename TALG, typename ...TARGS>
    static typename T::FeatureList CalculateFeaturesImpl(std::true_type, const typename T::FeatureToolMap &featureToolMap, 
        const TALG *const pAlgorithm, TARGS &&... args);
    
    /**
     *  @brief  Calculate the vector of features from a set of algorithm tools with return type T and base type SVMFeatureToolBase<TALG, Ts...>.
     * 
     *  @param  featureToolMap the feature tool map
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  args arguments to pass to the algorithm tool
     * 
     *  @return the vector of features
     */
    template <typename T, typename TALG, typename ...Ts, typename ...TARGS>
    static typename SVMFeatureTool<typename std::decay<T>::type, TALG, Ts...>::FeatureList 
    CalculateFeaturesImpl(std::false_type, const SVMFeatureToolMap<TALG, Ts...> &featureToolMap, const TALG *const pAlgorithm, TARGS &&... args);
      
    /**
     *  @brief  Calculate the vector of features from a set of algorithm tools matching type T.
     * 
     *  @param  featureToolMap the feature tool map
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  args arguments to pass to the algorithm tool
     * 
     *  @return the vector of features
     */
    template <typename T, typename TALG, typename ...TARGS>
    static typename T::FeatureList GetVectorOfFeatures(const typename T::FeatureToolMap &featureToolMap, const TALG *const pAlgorithm, TARGS &&... args);

    /**
     *  @brief  Recursively write the features of the given lists to file
     *
     *  @param  outFile the std::ofstream object to use
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
     *  @param  outFile the std::ofstream object to use
     *  @param  delimiter the delimiter string
     *  @param  featureList a list of features to write
     *
     *  @return success
     */
    template <typename TLIST>
    static pandora::StatusCode WriteFeaturesToFileImpl(std::ofstream &outfile, const std::string &delimiter, TLIST &&featureList);
    
    /**
     *  @brief  Create an compile-time constant-sized array of doubles from a set of lists of features (implementation method)
     *
     *  @param  featureList a list of features to write
     *  @param  featureLists optional further lists of features to write
     *
     *  @return an array of doubles
     */
    template <std::size_t NFEATURES, typename TLIST, typename ...TLISTS>
    static typename SupportVectorMachine<NFEATURES>::FeatureArray GetFeatureArrayImpl(TLIST &&featureList, TLISTS &&... featureLists);

    /**
     *  @brief  Recursively concatenate sets of feature lists into an array
     *
     *  @param  featureArray the array to create
     *  @param  featureCounter the feature counter
     *  @param  featureList a list of features to write
     *  @param  featureLists optional further lists of features to write
     */
    template <std::size_t NFEATURES, typename TLIST, typename ...TLISTS>
    static void ConcatenateFeatureLists(typename SupportVectorMachine<NFEATURES>::FeatureArray &featureArray, std::size_t &featureCounter,
        TLIST &&featureList, TLISTS &&... featureLists);

    /**
     *  @brief  Recursively concatenate sets of feature lists into an array (terminating method)
     */
    template <std::size_t NFEATURES>
    static void ConcatenateFeatureLists(typename SupportVectorMachine<NFEATURES>::FeatureArray &, std::size_t &);
    
public:
    /**
     *  @brief  Produce a training example with the given features and result
     *
     *  @param  trainingOutputFile the file to which to append the example
     * 
     *  @param  featureList a list of int, float or double features
     *  @param  featureLists optional further lists of int, float or double features
     *
     *  @return  success
     */
    template <typename TLIST, typename ...TLISTS>
    static pandora::StatusCode ProduceTrainingExample(const std::string &trainingOutputFile, const bool result, TLIST &&featureList, TLISTS &&... featureLists);
    
    /**
     *  @brief  Use the trained SVM to predict the boolean class of an example
     *
     *  @param  sVMachine the support vector machine
     *  @param  featureList the list of features
     *  @param  featureLists optional further lists of features
     *
     *  @return  The predicted boolean class of the example
     */
    template <std::size_t NFEATURES, typename TLIST, typename ...TLISTS>
    static bool Classify(const SupportVectorMachine<NFEATURES> &sVMachine, TLIST &&featureList, TLISTS &&... featureLists);
    
    /**
     *  @brief  Use the trained SVM to calculate the classification score of an example (>0 means boolean class true)
     *
     *  @param  sVMachine the support vector machine
     *  @param  featureList the list of features
     *  @param  featureLists optional further lists of features
     *
     *  @return  the classification score
     */
    template <std::size_t NFEATURES, typename TLIST, typename ...TLISTS>
    static double CalculateClassificationScore(const SupportVectorMachine<NFEATURES> &sVMachine, TLIST &&featureList, TLISTS &&... featureLists);

    /**
     *  @brief  Use a feature tool map to calculate the features of a given return type (int, float, double, bool) or derived FeatureTool type T for a set of arguments
     *
     *  @param  featureToolMap the feature tool map
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  args arguments to pass to the tool
     *
     *  @return  the list of features
     */
    template <typename T, typename TALG, typename TALGARG, typename ...Ts, typename ...TARGS>
    static auto CalculateFeatures(const SVMFeatureToolMap<TALG, Ts...> &featureToolMap, const TALGARG *const pAlgorithm, TARGS &&... args) 
        -> decltype(CalculateFeaturesImpl<T>(typename std::is_base_of<pandora::AlgorithmTool, T>::type(), featureToolMap, pAlgorithm, std::forward<TARGS>(args)...)); // this line can disappear in C++14
    
    /**
     *  @brief  Use a feature tool map to calculate the feature from a derived FeatureTool T
     *
     *  @param  featureToolMap the feature tool map
     *  @param  pAlgorithm address of the calling algorithm
     *  @param  args arguments to pass to the tool
     *
     *  @return  the list of features
     */
    template <typename T, typename TALG, typename TALGARG, typename ...Ts, typename ...TARGS, typename std::enable_if<std::is_base_of<pandora::AlgorithmTool, T>{}, int>::type = 0>
    static typename T::ReturnType CalculateFeature(const SVMFeatureToolMap<TALG, Ts...> &featureToolMap, const TALGARG *const pAlgorithm, TARGS &&... args);
    
    /**
     *  @brief  Add all the feature tools of a given return type T (int, float, double) to a given feature tool map
     *
     *  @param  pFeatureTool address of the feature tool
     *  @param  featureToolMap the feature tool map to append
     *
     *  @return whether the feature tool was added with this T
     */
    template <typename TALG, typename ...Ts>
    static pandora::StatusCode AddFeatureToolToMap(pandora::AlgorithmTool *const pFeatureTool, SVMFeatureToolMap<TALG, Ts...> &featureToolMap);
    
    /**
     *  @brief  Create an compile-time constant-sized array of doubles from a set of lists of features
     *
     *  @param  featureList a list of features to write
     *  @param  featureLists optional further lists of features to write
     *
     *  @return an array of doubles
     */
    template <std::size_t NFEATURES, typename TLIST, typename ...TLISTS>
    static typename SupportVectorMachine<NFEATURES>::FeatureArray GetFeatureArray(TLIST &&featureList, TLISTS &&... featureLists);
};

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TLIST, typename ...TLISTS>
pandora::StatusCode SVMHelper::ProduceTrainingExample(const std::string &trainingOutputFile, const bool result, TLIST &&featureList, TLISTS &&... featureLists)
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

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, WriteFeaturesToFile(outfile, delimiter, featureList, featureLists...));
    outfile << static_cast<int>(result) << '\n';
  
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES, typename TLIST, typename ...TLISTS>
bool SVMHelper::Classify(const SupportVectorMachine<NFEATURES> &sVMachine, TLIST &&featureList, TLISTS &&... featureLists)
{    
    return sVMachine.Classify(GetFeatureArray<NFEATURES>(std::forward<TLIST>(featureList), std::forward<TLISTS>(featureLists)...));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES, typename TLIST, typename ...TLISTS>
double SVMHelper::CalculateClassificationScore(const SupportVectorMachine<NFEATURES> &sVMachine, TLIST &&featureList, TLISTS &&... featureLists)
{    
    return sVMachine.CalculateClassificationScore(GetFeatureArray<NFEATURES>(std::forward<TLIST>(featureList), std::forward<TLISTS>(featureLists)...));
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T, typename TALG, typename TALGARG, typename ...Ts, typename ...TARGS>
auto SVMHelper::CalculateFeatures(const SVMFeatureToolMap<TALG, Ts...> &featureToolMap, const TALGARG *const pAlgorithm, TARGS &&... args) 
    -> decltype(CalculateFeaturesImpl<T>(typename std::is_base_of<pandora::AlgorithmTool, T>::type(), featureToolMap, pAlgorithm, std::forward<TARGS>(args)...)) // this line can disappear in C++14
{
    return CalculateFeaturesImpl<T>(typename std::is_base_of<pandora::AlgorithmTool, T>::type(), featureToolMap, pAlgorithm, std::forward<TARGS>(args)...);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T, typename TALG, typename TALGARG, typename ...Ts, typename ...TARGS, typename std::enable_if<std::is_base_of<pandora::AlgorithmTool, T>{}, int>::type>
typename T::ReturnType SVMHelper::CalculateFeature(const SVMFeatureToolMap<TALG, Ts...> &featureToolMap, const TALGARG *const pAlgorithm, TARGS &&... args)
{
    using TD = typename std::decay<T>::type;
    
    auto iterPair = featureToolMap.equal_range(typeid(typename TD::ReturnType));    
    for (auto iter = iterPair.first; iter != iterPair.second; ++iter)
    {
        if (TD *const pFeatureTool = dynamic_cast<TD *const>(iter->second))
           return pFeatureTool->Run(pAlgorithm, std::forward<TARGS>(args)...);
    }
           
    std::cout << "SVMHelper: could not find the feature tool in the map" << std::endl;
    throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TALG, typename ...Ts>
pandora::StatusCode SVMHelper::AddFeatureToolToMap(pandora::AlgorithmTool *const pFeatureTool, SVMFeatureToolMap<TALG, Ts...> &featureToolMap)
{    
    if (AddFeatureToolToMapImpl<double>(pFeatureTool, featureToolMap) ||
        AddFeatureToolToMapImpl<float>(pFeatureTool, featureToolMap) ||
        AddFeatureToolToMapImpl<int>(pFeatureTool, featureToolMap) ||
        AddFeatureToolToMapImpl<bool>(pFeatureTool, featureToolMap))
    {
        return pandora::STATUS_CODE_SUCCESS;
    }
    
    std::cout << "SVMHelper: Could not add FeatureTool to map (return type must be double, float, int or bool)" << std::endl;
    return pandora::STATUS_CODE_INVALID_PARAMETER;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline std::string SVMHelper::GetTimestampString()
{
    std::time_t timestampNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    
    struct tm * timeinfo;
    char buffer [80];

    timeinfo = localtime (&timestampNow);
    strftime(buffer, 80, "%x_%X", timeinfo);
    
    std::string timeString(buffer);
    
    if (!timeString.empty() && timeString.back() == '\n') // last char is always a newline
        timeString.pop_back();
        
    return timeString;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES, typename TLIST, typename ...TLISTS>
inline typename SupportVectorMachine<NFEATURES>::FeatureArray
SVMHelper::GetFeatureArray(TLIST &&featureList, TLISTS &&... featureLists)
{    
    return GetFeatureArrayImpl<NFEATURES>(std::forward<TLIST>(featureList), std::forward<TLISTS>(featureLists)...);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T, typename TALG, typename ...TARGS>
typename T::FeatureList 
SVMHelper::CalculateFeaturesImpl(std::true_type, const typename T::FeatureToolMap &featureToolMap, const TALG *const pAlgorithm, TARGS &&... args)
{
    return GetVectorOfFeatures<T>(featureToolMap, pAlgorithm, std::forward<TARGS>(args)...);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T, typename TALG, typename ...Ts>
bool SVMHelper::AddFeatureToolToMapImpl(pandora::AlgorithmTool *const pFeatureTool, SVMFeatureToolMap<TALG, Ts...> &featureToolMap)
{
    using TD = typename std::decay<T>::type;
    if (SVMFeatureTool<TD, TALG, Ts...> *const pCastFeatureTool = dynamic_cast<SVMFeatureTool<TD, TALG, Ts...> *const>(pFeatureTool))
    {
        featureToolMap.emplace(typeid(TD), pCastFeatureTool);
        return true;
    }
    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T, typename TALG, typename ...Ts, typename ...TARGS>
typename SVMFeatureTool<typename std::decay<T>::type, TALG, Ts...>::FeatureList 
SVMHelper::CalculateFeaturesImpl(std::false_type, const SVMFeatureToolMap<TALG, Ts...> &featureToolMap, const TALG *const pAlgorithm, TARGS &&... args)
{
    using TD = typename std::decay<T>::type;
    static_assert((std::is_same<TD, double>::value || std::is_same<TD, float>::value || std::is_same<TD, int>::value || std::is_same<TD, bool>::value), 
              "SVMHelper: Features may only be double, float, int or bool");
              
    return GetVectorOfFeatures<SVMFeatureTool<T, TALG, Ts...>>(featureToolMap, pAlgorithm, std::forward<TARGS>(args)...);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T, typename TALG, typename ...TARGS>
typename T::FeatureList 
SVMHelper::GetVectorOfFeatures(const typename T::FeatureToolMap &featureToolMap, const TALG *const pAlgorithm, TARGS &&... args)
{
    using TD = typename std::decay<T>::type;
    
    typename TD::FeatureList features;
    auto iterPair = featureToolMap.equal_range(typeid(typename TD::ReturnType));
    
    for (auto iter = iterPair.first; iter != iterPair.second; ++iter)
    {
        if (TD *const pFeatureTool = dynamic_cast<TD *>(iter->second))
           features.push_back(pFeatureTool->Run(pAlgorithm, std::forward<TARGS>(args)...));
    }
            
    return features;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename TLIST, typename ...TLISTS>
inline pandora::StatusCode SVMHelper::WriteFeaturesToFile(std::ofstream &outfile, const std::string &delimiter, TLIST &&featureList, TLISTS &&... featureLists)
{
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
    for (const auto feature : featureList)
        outfile << feature << delimiter;
    
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES, typename TLIST, typename ...TLISTS>
typename SupportVectorMachine<NFEATURES>::FeatureArray SVMHelper::GetFeatureArrayImpl(TLIST &&featureList, TLISTS &&... featureLists)
{    
    typename SupportVectorMachine<NFEATURES>::FeatureArray featureArray;
    std::size_t featureCounter(0);
    ConcatenateFeatureLists<NFEATURES>(featureArray, featureCounter, featureList, featureLists...);

    if (featureCounter != NFEATURES)
    {
        std::cout << "SVMHelper: could not create feature array because the number of features (" << featureCounter << ") did not match "
                  << "the expected number of features (" << NFEATURES << ")" << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
    }

    return featureArray;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES, typename TLIST, typename ...TLISTS>
void SVMHelper::ConcatenateFeatureLists(typename SupportVectorMachine<NFEATURES>::FeatureArray &featureArray, std::size_t &featureCounter, 
                                        TLIST &&featureList, TLISTS &&... featureLists)
{
    if ((featureCounter + featureList.size()) > featureArray.size())
    {
        std::cout << "SVMHelper: could not create feature array because there were more than the expected " << featureArray.size() 
                  << " features" << std::endl;
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
    }
    
    for (const auto feature : featureList)
        featureArray[featureCounter++] = static_cast<double>(feature);
        
    ConcatenateFeatureLists<NFEATURES>(featureArray, featureCounter, std::forward<TLISTS>(featureLists)...);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <std::size_t NFEATURES>
void SVMHelper::ConcatenateFeatureLists(typename SupportVectorMachine<NFEATURES>::FeatureArray &, std::size_t &)
{
    return;
}

} // namespace lar_content

#endif // #ifndef LAR_SVM_HELPER_H