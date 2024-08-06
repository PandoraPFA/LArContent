/**
 *  @file   larpandoracontent/LArUtility/KDTreeLinkerToolsT.h
 *
 *  @brief  Header file for the kd tree linker tools template class
 *
 *  $Log: $
 */
#ifndef LAR_KD_TREE_LINKER_TOOLS_TEMPLATED_H
#define LAR_KD_TREE_LINKER_TOOLS_TEMPLATED_H 1

#include "Api/PandoraContentApi.h"

#include "Objects/CaloHit.h"
#include "Objects/CartesianVector.h"

#include "Pandora/PandoraInternal.h"

#include <array>
#include <vector>

namespace pandora
{
class Algorithm;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_content
{

/**
 *  @brief  Box structure used to define 2D field. It's used in KDTree building step to divide the detector space (ECAL, HCAL...) and
 *          in searching step to create a bounding box around the demanded point (Track collision point, PS projection...).
 */
template <unsigned DIM>
class KDTreeBoxT
{
public:
    /**
     *  @brief  Default constructor
     */
    KDTreeBoxT();

    /**
     *  @brief  Constructor
     *
     *  @param  dimargs
     */
    template <typename... Ts>
    KDTreeBoxT(Ts... dimargs);

    std::array<float, DIM> dimmin; ///<
    std::array<float, DIM> dimmax; ///<
};

typedef KDTreeBoxT<2> KDTreeBox;
typedef KDTreeBoxT<3> KDTreeCube;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Data stored in each KDTree node. The dim1/dim2 fields are usually the duplication of some PFRecHit values
 *          (eta/phi or x/y). But in some situations, phi field is shifted by +-2.Pi
 */
template <typename DATA, unsigned DIM>
class KDTreeNodeInfoT
{
public:
    /**
     *  @brief  Default constructor
     */
    KDTreeNodeInfoT();

    /**
     *  @brief  Constructor
     *
     *  @param  d
     *  @param  dimargs
     */
    template <typename... Ts>
    KDTreeNodeInfoT(const DATA &d, Ts... dimargs);

    DATA data;                   ///<
    std::array<float, DIM> dims; ///<
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  KDTree node
 */
template <typename DATA, unsigned DIM>
class KDTreeNodeT
{
public:
    /**
     *  @brief  Default constructor
     */
    KDTreeNodeT();

    /**
     *  @brief  setAttributs
     *
     *  @param  regionBox
     *  @param  infoToStore
     */
    void setAttributs(const KDTreeBoxT<DIM> &regionBox, const KDTreeNodeInfoT<DATA, DIM> &infoToStore);

    /**
     *  @brief  setAttributs
     *
     *  @param  regionBox
     */
    void setAttributs(const KDTreeBoxT<DIM> &regionBox);

    KDTreeNodeInfoT<DATA, DIM> info; ///< Data
    KDTreeNodeT<DATA, DIM> *left;    ///< Left son
    KDTreeNodeT<DATA, DIM> *right;   ///< Right son
    KDTreeBoxT<DIM> region;          ///< Region bounding box.
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  kdtree_type_adaptor
 */
template <typename T>
class kdtree_type_adaptor
{
public:
    /**
     *  @brief  position
     *
     *  @param  t
     *
     *  @return position
     */
    static const pandora::CartesianVector &position(const T *const t);
};

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  minmax
 *
 *  @param  a
 *  @param  b
 *
 *  @return minmax
 */
std::pair<float, float> minmax(const float a, const float b);

/**
 *  @brief  fill_and_bound_2d_kd_tree
 *
 *  @param  points
 *  @param  nodes
 *
 *  @return KDTreeCube
 */
template <typename T>
KDTreeBox fill_and_bound_2d_kd_tree(const MANAGED_CONTAINER<const T *> &points, std::vector<KDTreeNodeInfoT<const T *, 2>> &nodes);

/**
 *  @brief  fill_and_bound_3d_kd_tree
 *
 *  @param  points
 *  @param  nodes
 *
 *  @return KDTreeCube
 */
template <typename T>
KDTreeCube fill_and_bound_3d_kd_tree(const MANAGED_CONTAINER<const T *> &points, std::vector<KDTreeNodeInfoT<const T *, 3>> &nodes);

/**
 *  @brief  build_2d_kd_search_region
 *
 *  @param  point
 *  @param  x_span
 *  @param  z_span
 *
 *  @return KDTreeBox
 */
KDTreeBox build_2d_kd_search_region(const pandora::CaloHit *const point, const float x_span, const float z_span);

/**
 *  @brief  build_2d_kd_search_region
 *
 *  @param  pos
 *  @param  x_span
 *  @param  z_span
 *
 *  @return KDTreeBox
 */
KDTreeBox build_2d_kd_search_region(const pandora::CartesianVector &pos, const float x_span, const float z_span);

/**
 *  @brief  build_3d_kd_search_region
 *
 *  @param  point
 *  @param  x_span
 *  @param  y_span
 *  @param  z_span
 *
 *  @return KDTreeCube
 */
KDTreeCube build_3d_kd_search_region(const pandora::CaloHit *const point, const float x_span, const float y_span, const float z_span);

/**
 *  @brief  build_3d_kd_search_region
 *
 *  @param  pos
 *  @param  x_span
 *  @param  y_span
 *  @param  z_span
 *
 *  @return KDTreeCube
 */
KDTreeCube build_3d_kd_search_region(const pandora::CartesianVector &pos, const float x_span, const float y_span, const float z_span);

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <unsigned DIM>
inline KDTreeBoxT<DIM>::KDTreeBoxT()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <unsigned DIM>
template <typename... Ts>
inline KDTreeBoxT<DIM>::KDTreeBoxT(Ts... dimargs)
{
    static_assert(sizeof...(dimargs) == 2 * DIM, "Constructor requires 2*DIM args");
    std::vector<float> dims = {dimargs...};

    for (unsigned i = 0; i < DIM; ++i)
    {
        dimmin[i] = dims[2 * i];
        dimmax[i] = dims[2 * i + 1];
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline KDTreeNodeInfoT<DATA, DIM>::KDTreeNodeInfoT() :
    data()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
template <typename... Ts>
inline KDTreeNodeInfoT<DATA, DIM>::KDTreeNodeInfoT(const DATA &d, Ts... dimargs) :
    data(d),
    dims{{dimargs...}}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline KDTreeNodeT<DATA, DIM>::KDTreeNodeT() :
    left(nullptr),
    right(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline void KDTreeNodeT<DATA, DIM>::setAttributs(const KDTreeBoxT<DIM> &regionBox, const KDTreeNodeInfoT<DATA, DIM> &infoToStore)
{
    info = infoToStore;
    region = regionBox;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline void KDTreeNodeT<DATA, DIM>::setAttributs(const KDTreeBoxT<DIM> &regionBox)
{
    region = regionBox;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
inline const pandora::CartesianVector &kdtree_type_adaptor<T>::position(const T *const t)
{
    return t->GetPosition();
}

template <>
inline const pandora::CartesianVector &kdtree_type_adaptor<const pandora::CaloHit>::position(const pandora::CaloHit *const t)
{
    return t->GetPositionVector();
}

template <>
inline const pandora::CartesianVector &kdtree_type_adaptor<const pandora::CartesianVector>::position(const pandora::CartesianVector *const t)
{
    return *t;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
KDTreeBox fill_and_bound_2d_kd_tree(const MANAGED_CONTAINER<const T *> &points, std::vector<KDTreeNodeInfoT<const T *, 2>> &nodes)
{
    std::array<float, 2> minpos{{0.f, 0.f}}, maxpos{{0.f, 0.f}};

    unsigned i = 0;

    for (const T *const point : points)
    {
        const pandora::CartesianVector &pos = kdtree_type_adaptor<const T>::position(point);
        nodes.emplace_back(point, pos.GetX(), pos.GetZ());

        if (0 == i)
        {
            minpos[0] = pos.GetX();
            minpos[1] = pos.GetZ();
            maxpos[0] = pos.GetX();
            maxpos[1] = pos.GetZ();
        }
        else
        {
            minpos[0] = std::min(pos.GetX(), minpos[0]);
            minpos[1] = std::min(pos.GetZ(), minpos[1]);
            maxpos[0] = std::max(pos.GetX(), maxpos[0]);
            maxpos[1] = std::max(pos.GetZ(), maxpos[1]);
        }

        ++i;
    }

    return KDTreeBox(minpos[0], maxpos[0], minpos[1], maxpos[1]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename T>
KDTreeCube fill_and_bound_3d_kd_tree(const MANAGED_CONTAINER<const T *> &points, std::vector<KDTreeNodeInfoT<const T *, 3>> &nodes)
{
    std::array<float, 3> minpos{{0.f, 0.f, 0.f}}, maxpos{{0.f, 0.f, 0.f}};

    unsigned i = 0;

    for (const T *const point : points)
    {
        const pandora::CartesianVector &pos = kdtree_type_adaptor<const T>::position(point);
        nodes.emplace_back(point, pos.GetX(), pos.GetY(), pos.GetZ());

        if (0 == i)
        {
            minpos[0] = pos.GetX();
            minpos[1] = pos.GetY();
            minpos[2] = pos.GetZ();
            maxpos[0] = pos.GetX();
            maxpos[1] = pos.GetY();
            maxpos[2] = pos.GetZ();
        }
        else
        {
            minpos[0] = std::min(pos.GetX(), minpos[0]);
            minpos[1] = std::min(pos.GetY(), minpos[1]);
            minpos[2] = std::min(pos.GetZ(), minpos[2]);
            maxpos[0] = std::max(pos.GetX(), maxpos[0]);
            maxpos[1] = std::max(pos.GetY(), maxpos[1]);
            maxpos[2] = std::max(pos.GetZ(), maxpos[2]);
        }

        ++i;
    }

    return KDTreeCube(minpos[0], maxpos[0], minpos[1], maxpos[1], minpos[2], maxpos[2]);
}

} // namespace lar_content

#endif // LAR_KD_TREE_LINKER_TOOLS_TEMPLATED_H
