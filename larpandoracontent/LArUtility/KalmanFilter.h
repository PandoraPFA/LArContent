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
 *  @brief  KalmanFilter class
 */
template <int DIM>
class KalmanFilter
{
    static_assert(DIM > 0, "KalmanFilter dimensionality must be positive");

public:
    using StateVector = Eigen::Matrix<double, 2 * DIM, 1>;
    using MeasurementVector = Eigen::Matrix<double, DIM, 1>;
    using StateMatrix = Eigen::Matrix<double, 2 * DIM, 2 * DIM>;
    using MeasurementMatrix = Eigen::Matrix<double, DIM, 2 * DIM>;
    using CovarianceMatrix = Eigen::Matrix<double, DIM, DIM>;
    using PositionVector = Eigen::Matrix<double, DIM, 1>;
    /**
     *  @brief  Default constructor
     *
     *  @param  dt Time step
     *  @param  processVariance Process variance
     *  @param  measurementVariance Measurement variance
     *  @param  x Initial position
     */
    KalmanFilter(double dt, double processVariance, double measurementVariance, const PositionVector &x, const double initialCovariance = 1.0);

    ~KalmanFilter() = default;

    KalmanFilter(const KalmanFilter &) = default;
    KalmanFilter &operator=(const KalmanFilter &) = default;
    KalmanFilter(KalmanFilter &&) = default;
    KalmanFilter &operator=(KalmanFilter &&) = default;

    /**
     * @brief  Uses the state transition matrix to predict the next state
     */
    void Predict();

    /**
     * @brief  Uses the measurement matrix to update the state estimate
     *
     * @param  z The position to use for updating the state estimate
     */
    void Update(const MeasurementVector &z);

    /**
     * @brief  Get the current state vector
     *
     * @return The current state vector
     */
    const StateVector &GetState() const;

    /**
     * @brief  Get the current position vector
     *
     * @return The current postion vector
     */
    const PositionVector GetPosition() const;

    /**
     * @brief  Get the current direction vector
     *
     * @return The current direction vector
     */
    const PositionVector GetDirection() const;

    /**
     * @brief  Get the temporary state vector
     *
     * @return The temporary state vector
     */
    const StateVector &GetTemporaryState() const;

private:
    double m_dt;            ///< Time step
    int m_stateSize;        ///< Size of the state vector
    int m_measurementSize;  ///  Size of the measurement vector
    StateVector m_x;        ///< State vector (tracks position and 'velocity')
    StateMatrix m_P;        ///< Covariance matrix
    StateMatrix m_F;        ///< State transition matrix
    MeasurementMatrix m_H;  ///< Measurement matrix
    CovarianceMatrix m_R;   ///< Measurement covariance matrix
    StateMatrix m_Q;        ///< Process covariance matrix
    StateMatrix m_identity; ///< Identity matrix
    StateVector m_xTemp;    ///< Temporary state vector
    StateMatrix m_PTemp;    ///< Temporary covariance matrix
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <int DIM>
KalmanFilter<DIM>::KalmanFilter(double dt, double processVariance, double measurementVariance, const KalmanFilter<DIM>::PositionVector &x,
    const double initialCovariance) :
    m_dt{dt},
    m_stateSize{2 * DIM},
    m_measurementSize{DIM}
{
    m_x = StateVector::Zero();
    m_x.head(DIM) = x;
    m_P = StateMatrix::Identity() * initialCovariance;
    m_F = StateMatrix::Identity();
    m_F.topRightCorner(DIM, DIM) = Eigen::Matrix<double, DIM, DIM>::Identity() * m_dt;
    m_H = MeasurementMatrix::Zero();
    m_H.leftCols(DIM) = Eigen::Matrix<double, DIM, DIM>::Identity();
    Eigen::Matrix<double, 2 * DIM, DIM> G = Eigen::Matrix<double, 2 * DIM, DIM>::Zero();
    G.topRows(DIM) = Eigen::Matrix<double, DIM, DIM>::Identity() * (0.5 * m_dt * m_dt);
    G.bottomRows(DIM) = Eigen::Matrix<double, DIM, DIM>::Identity() * m_dt;
    m_Q = G * G.transpose() * processVariance;
    m_R = CovarianceMatrix::Identity() * measurementVariance;

    m_xTemp = m_x;
    m_PTemp = m_P;
    m_identity = StateMatrix::Identity();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <int DIM>
void KalmanFilter<DIM>::Predict()
{
    m_xTemp = m_F * m_x;
    m_PTemp = m_F * m_P * m_F.transpose() + m_Q;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <int DIM>
void KalmanFilter<DIM>::Update(const KalmanFilter<DIM>::MeasurementVector &z)
{
    m_x = m_xTemp;
    m_P = m_PTemp;
    MeasurementVector y = z - m_H * m_x;
    auto Ht = m_H.transpose();
    auto S = m_H * m_P * Ht + m_R;
    auto K = m_P * Ht * S.inverse();
    m_x += K * y;
    m_P = (m_identity - K * m_H) * m_P;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <int DIM>
const typename KalmanFilter<DIM>::StateVector &KalmanFilter<DIM>::GetState() const
{
    return m_x;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <int DIM>
const typename KalmanFilter<DIM>::PositionVector KalmanFilter<DIM>::GetPosition() const
{
    return m_x.head(DIM);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <int DIM>
const typename KalmanFilter<DIM>::PositionVector KalmanFilter<DIM>::GetDirection() const
{
    return m_x.tail(DIM).normalized();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <int DIM>
const typename KalmanFilter<DIM>::StateVector &KalmanFilter<DIM>::GetTemporaryState() const
{
    return m_xTemp;
}

using KalmanFilter2D = KalmanFilter<2>;
using KalmanFilter3D = KalmanFilter<3>;

} // namespace lar_content

#endif // #ifndef LAR_KALMAN_FILTER_H
