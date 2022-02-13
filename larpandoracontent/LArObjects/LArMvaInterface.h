/**
 *  @file   larpandoracontent/LArObjects/LArMvaInterface.h
 *
 *  @brief  Header file for the lar multivariate analysis interface class.
 *
 *  $Log: $
 */
#ifndef LAR_MVA_INTERFACE_H
#define LAR_MVA_INTERFACE_H 1

#include "Pandora/PandoraInputTypes.h"

#include <vector>

namespace lar_content
{

/**
 *  @brief  MvaTypes class
 */
class MvaTypes
{
public:
    /**
     *   @brief  InitializedDouble class used to define mva features
     */
    class InitializedDouble
    {
    public:
        /**
         *  @brief  Default constructor.
         */
        InitializedDouble();

        /**
         *  @brief  Constructor.
         *
         *  @param  number to hold in class
         */
        InitializedDouble(const double number);

        /**
         *  @brief  Copy constructor
         *
         *  @param  rhs the initialized double to copy
         */
        InitializedDouble(const InitializedDouble &rhs);

        /**
         *  @brief  Assignment operator
         *
         *  @param  number the double to assign
         */
        InitializedDouble &operator=(const double number);

        /**
         *  @brief  Assignment operator
         *
         *  @param  rhs the initialized double to assign
         */
        InitializedDouble &operator=(const InitializedDouble rhs);

        /**
         *  @brief  Get number held in class
         *
         *  @return number held in class
         */
        double Get() const;

        /**
         *  @brief  Check number has been initialized
         *
         *  @return whether number has been initialized
         */
        bool IsInitialized() const;

    private:
        double m_number;      ///< Number held by class
        bool m_isInitialized; ///< Whether the number has been initialized
    };

    typedef InitializedDouble MvaFeature;
    typedef std::vector<MvaFeature> MvaFeatureVector;
    typedef std::map<std::string, MvaFeature> MvaFeatureMap;
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  MvaInterface class
 */
class MvaInterface
{
public:
    /**
     *  @brief  Classify the set of input features based on the trained model
     *
     *  @param  features the input features
     *
     *  @return the classification
     */
    virtual bool Classify(const MvaTypes::MvaFeatureVector &features) const = 0;

    /**
     *  @brief  Calculate the classification score for a set of input features, based on the trained model
     *
     *  @param  features the input features
     *
     *  @return the classification score
     */
    virtual double CalculateClassificationScore(const MvaTypes::MvaFeatureVector &features) const = 0;

    /**
     *  @brief  Calculate the classification probability for a set of input features, based on the trained model
     *
     *  @param  features the input features
     *
     *  @return the classification probability
     */
    virtual double CalculateProbability(const MvaTypes::MvaFeatureVector &features) const = 0;

    /**
     *  @brief  Destructor
     */
    virtual ~MvaInterface() = default;
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

inline MvaTypes::InitializedDouble::InitializedDouble() : m_number(0.), m_isInitialized(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline MvaTypes::InitializedDouble::InitializedDouble(const double number) : m_number(number), m_isInitialized(true)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline MvaTypes::InitializedDouble::InitializedDouble(const InitializedDouble &rhs) :
    m_number(rhs.m_number),
    m_isInitialized(rhs.m_isInitialized)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline MvaTypes::InitializedDouble &MvaTypes::InitializedDouble::operator=(const double number)
{
    m_number = number;
    m_isInitialized = true;

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline MvaTypes::InitializedDouble &MvaTypes::InitializedDouble::operator=(const InitializedDouble rhs)
{
    if (this != &rhs)
    {
        m_number = rhs.m_number;
        m_isInitialized = rhs.m_isInitialized;
    }

    return *this;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline double MvaTypes::InitializedDouble::Get() const
{
    if (!m_isInitialized)
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return m_number;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool MvaTypes::InitializedDouble::IsInitialized() const
{
    return m_isInitialized;
}

} // namespace lar_content

#endif // #ifndef LAR_MVA_INTERFACE_H
