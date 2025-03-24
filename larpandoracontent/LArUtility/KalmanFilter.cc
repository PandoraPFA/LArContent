/**
 *  @file   larpandoracontent/LArUtility/KalmanFilter.cc
 *
 *  @brief  Kalman filter for predicting 3D trajectories.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArUtility/KalmanFilter.h"

using namespace Eigen;

namespace lar_content
{

KalmanFilter3D::KalmanFilter3D(double dt, double processVariance, double measurementVariance, const VectorXd &x) :
    m_dt{dt}
{
    m_x = VectorXd::Zero(6);
    m_x.head(3) = x;
    m_P = MatrixXd::Identity(6, 6);
    m_F = MatrixXd::Identity(6, 6);
    for (int i = 0; i < 3; ++i)
    {
        m_F(i, i + 3) = dt;
    }
    m_H = MatrixXd::Zero(3, 6);
    m_H.block<3,3>(0,0) = MatrixXd::Identity(3,3);
    m_Q = MatrixXd::Zero(6, 6);
    for (int i = 0; i < 3; ++i)
    {
        m_Q(i, i) = (0.25 * dt * dt * dt * dt) * processVariance;
        m_Q(i, i + 3) = (0.5 * dt * dt * dt) * processVariance;
        m_Q(i + 3, i) = (0.5 * dt * dt * dt) * processVariance;
        m_Q(i + 3, i + 3) = (dt * dt) * processVariance;
    }
    m_R = MatrixXd::Identity(3, 3) * measurementVariance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

KalmanFilter3D::~KalmanFilter3D()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanFilter3D::Predict()
{
    m_xTemp = m_F * m_x;
    m_PTemp = m_F * m_P * m_F.transpose() + m_Q;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanFilter3D::Update(const VectorXd &z)
{
    m_x = m_xTemp;
    m_P = m_PTemp;
    VectorXd y = z - m_H * m_x;
    MatrixXd S = m_H * m_P * m_H.transpose() + m_R;
    MatrixXd K = m_P * m_H.transpose() * S.inverse();
    m_x += K * y;
    m_P = (MatrixXd::Identity(6, 6) - K * m_H) * m_P;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VectorXd &KalmanFilter3D::GetState() const
{
    return m_x;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VectorXd KalmanFilter3D::GetPosition() const
{
    return m_x.head(3);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VectorXd KalmanFilter3D::GetDirection() const
{
    return m_x.tail(3).normalized();
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VectorXd &KalmanFilter3D::GetTemporaryState() const
{
    return m_xTemp;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

KalmanFilter2D::KalmanFilter2D(double dt, double processVariance, double measurementVariance, const VectorXd &x) :
    m_dt{dt}
{
    m_x = VectorXd::Zero(4);
    m_x.head(2) = x;
    m_P = MatrixXd::Identity(4, 4);
    m_F = MatrixXd::Identity(4, 4);
    for (int i = 0; i < 2; ++i)
    {
        m_F(i, i + 2) = dt;
    }
    m_H = MatrixXd::Zero(2, 4);
    m_H.block<2,2>(0,0) = MatrixXd::Identity(3,3);
    m_Q = MatrixXd::Zero(4, 4);
    for (int i = 0; i < 2; ++i)
    {
        m_Q(i, i) = (0.25 * dt * dt * dt * dt) * processVariance;
        m_Q(i, i + 2) = (0.5 * dt * dt * dt) * processVariance;
        m_Q(i + 2, i) = (0.5 * dt * dt * dt) * processVariance;
        m_Q(i + 2, i + 2) = (dt * dt) * processVariance;
    }
    m_R = MatrixXd::Identity(2, 2) * measurementVariance;
}

//------------------------------------------------------------------------------------------------------------------------------------------

KalmanFilter2D::~KalmanFilter2D()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanFilter2D::Predict()
{
    m_xTemp = m_F * m_x;
    m_PTemp = m_F * m_P * m_F.transpose() + m_Q;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanFilter2D::Update(const VectorXd &z)
{
    m_x = m_xTemp;
    m_P = m_PTemp;
    VectorXd y{z - m_H * m_x};
    MatrixXd S{m_H * m_P * m_H.transpose() + m_R};
    MatrixXd K{m_P * m_H.transpose() * S.inverse()};
    m_x += K * y;
    m_P = (MatrixXd::Identity(4, 4) - K * m_H) * m_P;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VectorXd &KalmanFilter2D::GetState() const
{
    return m_x;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VectorXd KalmanFilter2D::GetPosition() const
{
    return m_x.head(2);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VectorXd KalmanFilter2D::GetDirection() const
{
    return m_x.tail(2).normalized();
}

//------------------------------------------------------------------------------------------------------------------------------------------

const VectorXd &KalmanFilter2D::GetTemporaryState() const
{
    return m_xTemp;
}

} // namespace lar_content
