/**
 *  @file   larpandoracontent/LArUtility/KalmanFilter.h
 *
 *  @brief  Kalman filter for predicting 3D trajectories.
 *
 *  $Log: $
 */
#ifndef LAR_KALMAN_FILTER_H
#define LAR_KALMAN_FILTER_H 1

#include "Eigen/Dense"

namespace lar_content
{

/**
 *  @brief  KalmanFilter3D class
 */
class KalmanFilter3D
{
public:
    /**
     *  @brief  Default constructor
     *
     *  @param  dt Time step
     *  @param  processVariance Process variance
     *  @param  measurementVariance Measurement variance
     *  @param  x Initial position
     */
    KalmanFilter3D(double dt, double processVariance, double measurementVariance, const Eigen::VectorXd &x);

    /**
     *  @brief  Destructor
     */
    ~KalmanFilter3D();

    /**
     * @brief  Uses the state transition matrix to predict the next state
     */
    void Predict();

    /**
     * @brief  Uses the measurement matrix to update the state estimate
     *
     * @param  z The position to use for updating the state estimate
     */
    void Update(const Eigen::VectorXd &z);

    /**
     * @brief  Get the current state vector
     *
     * @return The current state vector
     */
    const Eigen::VectorXd &GetState() const;

    /**
     * @brief  Get the current position vector
     *
     * @return The current postion vector
     */
    const Eigen::VectorXd GetPosition() const;

    /**
     * @brief  Get the current direction vector
     *
     * @return The current direction vector
     */
    const Eigen::VectorXd GetDirection() const;

    /**
     * @brief  Get the temporary state vector
     *
     * @return The temporary state vector
     */
    const Eigen::VectorXd &GetTemporaryState() const;

private:
    double m_dt;             ///< Time step
    Eigen::VectorXd m_x;     ///< State vector (tracks position and 'velocity')
    Eigen::MatrixXd m_P;     ///< Covariance matrix
    Eigen::MatrixXd m_F;     ///< State transition matrix
    Eigen::MatrixXd m_H;     ///< Measurement matrix
    Eigen::MatrixXd m_R;     ///< Measurement covariance matrix
    Eigen::MatrixXd m_Q;     ///< Process covariance matrix
    Eigen::VectorXd m_xTemp; ///< Temporary state vector
    Eigen::MatrixXd m_PTemp; ///< Temporary covariance matrix
};

/**
 *  @brief  KalmanFilter2D class
 */
class KalmanFilter2D
{
public:
    /**
     *  @brief  Default constructor
     *
     *  @param  dt Time step
     *  @param  processVariance Process variance
     *  @param  measurementVariance Measurement variance
     *  @param  x Initial position
     */
    KalmanFilter2D(double dt, double processVariance, double measurementVariance, const Eigen::VectorXd &x);

    /**
     *  @brief  Destructor
     */
    ~KalmanFilter2D();

    /**
     * @brief  Uses the state transition matrix to predict the next state
     */
    void Predict();

    /**
     * @brief  Uses the measurement matrix to update the state estimate
     *
     * @param  z The position to use for updating the state estimate
     */
    void Update(const Eigen::VectorXd &z);

    /**
     * @brief  Get the current state vector
     *
     * @return The current state vector
     */
    const Eigen::VectorXd &GetState() const;

    /**
     * @brief  Get the current position vector
     *
     * @return The current postion vector
     */
    const Eigen::VectorXd GetPosition() const;

    /**
     * @brief  Get the current direction vector
     *
     * @return The current direction vector
     */
    const Eigen::VectorXd GetDirection() const;

    /**
     * @brief  Get the temporary state vector
     *
     * @return The temporary state vector
     */
    const Eigen::VectorXd &GetTemporaryState() const;

private:
    double m_dt;             ///< Time step
    Eigen::VectorXd m_x;     ///< State vector (tracks position and 'velocity')
    Eigen::MatrixXd m_P;     ///< Covariance matrix
    Eigen::MatrixXd m_F;     ///< State transition matrix
    Eigen::MatrixXd m_H;     ///< Measurement matrix
    Eigen::MatrixXd m_R;     ///< Measurement covariance matrix
    Eigen::MatrixXd m_Q;     ///< Process covariance matrix
    Eigen::VectorXd m_xTemp; ///< Temporary state vector
    Eigen::MatrixXd m_PTemp; ///< Temporary covariance matrix
};

} // namespace lar_content

#endif // #ifndef LAR_KALMAN_FILTER_H
