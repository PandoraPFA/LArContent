/**
 *  @file   larpandoracontent/LArUtility/RANSAC/RANSAC.h
 *
 *  @brief  Header and implementation for the RANSAC.
 *          Single header file since since that is easier for something so templated.
 *          Code is from https://github.com/drsrinathsridhar/GRANSAC
 *          Original license file can be found in larpandoracontent/LArUtility/RANSAC/LICENSE.
 *
 *  $Log: $
 */

#ifndef LAR_RANSAC_ALGO_TEMPLATED_H
#define LAR_RANSAC_ALGO_TEMPLATED_H 1

#include <random>
#include <vector>

#include "larpandoracontent/LArUtility/RANSAC/AbstractModel.h"
#include "larpandoracontent/LArHelpers/LArRandomHelper.h"

namespace lar_content
{

template <class T, int t_numParams>
class RANSAC
{

public:

    /**
     *  @brief  Default constructor.
     *
     *  @param threshold Threshold to be used for model to point distance.
     *  @param numIterations How many iterations to run at most.
     */
    RANSAC(double threshold, int numIterations = 1000) :
        m_numIterations(numIterations),
        m_threshold(threshold)
    {
        this->Reset();
    };

    /**
     *  @brief  Reset all the data.
     */
    void Reset() { m_data.clear(); };

    /**
     *  @brief Destructor
     */
    virtual ~RANSAC() {};

    /**
     *  @brief Get the best model, where the best model has the most inliers.
     */
    std::shared_ptr<T> GetBestModel() { return m_bestModel; };

    /**
     *  @brief Get the parameters that are inliers for the best model.
     */
    ParameterVector& GetBestInliers() { return m_bestInliers; };

    /**
     *  @brief Get the second best model, where the second best model has the second most unique inliers, compared to the best model.
     *         This prevents the second best model being a model that is right next to the best model, and is instead somewhat more
     *         unique.
     */
    std::shared_ptr<T> GetSecondBestModel() { return m_secondBestModel; };

    /**
     *  @brief Get the parameters that are inliers for the second best model.
     */
    ParameterVector& GetSecondBestInliers() { return m_secondBestInliers; };

//------------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Given a vector of data, get the best and second best model for it.
     *
     *  @param  data The data to use for this RANSAC run.
     *
     *  @return bool If the models were populated.
     */
    bool Estimate(const ParameterVector &data)
    {
        if (data.size() <= t_numParams)
            return false;

        m_data = data;
        this->GenerateSamples();

        std::vector<double> inlierFrac(m_numIterations, 0.0);
        std::vector<ParameterVector> inliers(m_numIterations);
        std::vector<std::shared_ptr<T>> sampledModels(m_numIterations);

        RANSAC::CheckModel(inlierFrac, inliers, sampledModels);
        double bestModelScore(-1);

        for (unsigned int i = 0; i < sampledModels.size(); ++i)
        {
            if (inlierFrac[i] == 0.0)
                continue;

            if (inlierFrac[i] > bestModelScore)
            {
                bestModelScore = inlierFrac[i];
                m_bestModel = sampledModels[i];
                m_bestInliers = inliers[i];
            }
        }

        double secondModelScore(-1);

        for (unsigned int i = 0; i < sampledModels.size(); ++i)
        {
            if (inlierFrac[i] == 0.0)
                continue;

            if (inlierFrac[i] > (secondModelScore * 0.9))
            {
                ParameterVector diff;
                this->CompareToBestModel(inliers[i], diff);

                if (diff.size() > m_secondUniqueParamCount)
                {
                    secondModelScore = inlierFrac[i];
                    m_secondUniqueParamCount = diff.size();
                    m_secondBestModel = sampledModels[i];
                    m_secondBestInliers = inliers[i];
                }
            }
        }

        this->Reset();

        return true;
    }

//------------------------------------------------------------------------------------------------------------------------------------------

private:
    ParameterVector m_data;                    ///< The data for the RANSAC model.
    std::vector<ParameterVector> m_samples;    ///< Samples to use for model generation.

    std::shared_ptr<T> m_secondBestModel;      ///< Second best model, with most unique parameters compared to first.
    std::shared_ptr<T> m_bestModel;            ///< Pointer to best model, valid only after Estimate() is called.

    ParameterVector m_bestInliers;             ///< The parameters in the best model.
    ParameterVector m_secondBestInliers;       ///< The parameters in the second model.

    unsigned int m_secondUniqueParamCount = 0; ///< A count of how many unique parameters are in the second model.
    unsigned int m_numIterations;              ///< Number of RANSAC iterations.
    double m_threshold;                        ///< Threshold for model consensus.

//------------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Compare a candidate model to the current best model, to decide if it is sufficiently unique compared to the current best
     *  model, and the existing most-unique model.
     *
     *  @param  candidateInliers  The current inlying parameters to consider.
     *  @param  diff A parameter vector to fill with unique parameters.
     */
    void CompareToBestModel(ParameterVector &candidateInliers, ParameterVector &diff)
    {
        if (candidateInliers.size() < m_secondUniqueParamCount)
            return;

        for (unsigned int i = 0; i < candidateInliers.size(); ++i)
        {
            const auto p1 = candidateInliers[i];
            const unsigned int maxDiffSize(diff.size() + (candidateInliers.size() - i));

            // ATTN: Early return.
            if (maxDiffSize < m_secondUniqueParamCount)
                return;

            // INFO: If current parameter is in the best model, skip it.
            if (m_bestModel->ComputeDistanceMeasure(p1) < m_threshold)
                continue;

            diff.push_back(p1);
        }
    };

//------------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Generate a vector of samples to generate models from. These samples are randomly sampled from the full input data.
     */
    void GenerateSamples()
    {
        std::mt19937 eng(m_data.size());
        m_samples.clear();

        for (unsigned int i = 0; i < m_numIterations; ++i)
        {
            ParameterVector currentParameters(t_numParams);

            for (unsigned int j = 0; j < t_numParams; ++j)
                currentParameters[j] = m_data[GetIntsInRange(0, m_data.size() - 1, eng)];

            m_samples.push_back(currentParameters);
        }
    };

//------------------------------------------------------------------------------------------------------------------------------------------

    /**
     *  @brief  Generate and check a model based on some samples.
     *
     *  @param  inlierFrac Vector to store the percentage of inlying parameters for each model.
     *  @param  inliers The vector of inlying parameters for each model.
     *  @param  sampledModels Vector of the actual models that were generated and evaluated.
     */
    void CheckModel(std::vector<double> &inlierFrac, std::vector<ParameterVector> &inliers,
            std::vector<std::shared_ptr<T>> &sampledModels)
    {
        for (unsigned int i = 0; i < m_numIterations; ++i)
        {
            // Evaluate the current model, so that its performance can be checked later.
            const std::shared_ptr<T> randomModel = std::make_shared<T>(m_samples[i]);
            const std::pair<double, ParameterVector> evalPair = randomModel->Evaluate(m_data, m_threshold);

            // Push back into history.
            inliers[i] = evalPair.second;
            sampledModels[i] = randomModel;
            inlierFrac[i] = evalPair.first;

            // If the model contained every data point, stop.
            if (evalPair.first == m_data.size())
                break;
        }
    };
};

} // namespace lar_content
#endif // LAR_RANSAC_ALGO_TEMPLATED_H
