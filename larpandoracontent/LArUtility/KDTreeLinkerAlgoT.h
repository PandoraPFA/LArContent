/**
 *  @file   larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h
 *
 *  @brief  Header file for the kd tree linker algo template class
 *
 *  $Log: $
 */
#ifndef LAR_KD_TREE_LINKER_ALGO_TEMPLATED_H
#define LAR_KD_TREE_LINKER_ALGO_TEMPLATED_H

#include "KDTreeLinkerToolsT.h"

#include <vector>

namespace lar_content
{

/**
 *  @brief  Class that implements the KDTree partition of 2D space and a closest point search algorithm
 */
template <typename DATA, unsigned DIM = 2>
class KDTreeLinkerAlgo
{
public:
    /**
     *  @brief  Default constructor
     */
    KDTreeLinkerAlgo();

    /**
     *  @brief  Destructor calls clear
     */
    ~KDTreeLinkerAlgo();

    /**
     *  @brief  Build the KD tree from the "eltList" in the space define by "region"
     *
     *  @param  eltList
     *  @param  region
     */
    void build(std::vector<KDTreeNodeInfoT<DATA, DIM>> &eltList, const KDTreeBoxT<DIM> &region);

    /**
     *  @brief  Search in the KDTree for all points that would be contained in the given searchbox
     *          The founded points are stored in resRecHitList
     *
     *  @param  searchBox
     *  @param  resRecHitList
     */
    void search(const KDTreeBoxT<DIM> &searchBox, std::vector<KDTreeNodeInfoT<DATA, DIM>> &resRecHitList);

    /**
     *  @brief  findNearestNeighbour
     *
     *  @param  point
     *  @param  result
     *  @param  distance
     */
    void findNearestNeighbour(const KDTreeNodeInfoT<DATA, DIM> &point, const KDTreeNodeInfoT<DATA, DIM> *&result, float &distance);

    /**
     *  @brief  Whether the tree is empty
     *
     *  @return boolean
     */
    bool empty();

    /**
     *  @brief  Return the number of nodes + leaves in the tree (nElements should be (size() +1) / 2)
     *
     *  @return the number of nodes + leaves in the tree
     */
    int size();

    /**
     *  @brief  Clear all allocated structures
     */
    void clear();

private:
    /**
     *  @brief  Get the next node from the node pool
     *
     *  @return the next node from the node pool
     */
    KDTreeNodeT<DATA, DIM> *getNextNode();

    /**
     *  @brief  Fast median search with Wirth algorithm in eltList between low and high indexes.
     *
     *  @param  low
     *  @param  high
     *  @param  treeDepth
     */
    int medianSearch(int low, int high, int treeDepth);

    /**
     *  @brief  Recursive kdtree builder. Is called by build()
     *
     *  @param  low
     *  @param  high
     *  @param  depth
     *  @param  region
     */
    KDTreeNodeT<DATA, DIM> *recBuild(int low, int high, int depth, const KDTreeBoxT<DIM> &region);

    /**
     *  @brief  Recursive kdtree search. Is called by search()
     *
     *  @param  current
     *  @param  trackBox
     */
    void recSearch(const KDTreeNodeT<DATA, DIM> *current, const KDTreeBoxT<DIM> &trackBox);

    /**
     *  @brief  Recursive nearest neighbour search. Is called by findNearestNeighbour()
     *
     *  @param  depth
     *  @param  current
     *  @param  point
     *  @param  best_match
     *  @param  best_dist
     */
    void recNearestNeighbour(unsigned depth, const KDTreeNodeT<DATA, DIM> *current, const KDTreeNodeInfoT<DATA, DIM> &point,
        const KDTreeNodeT<DATA, DIM> *&best_match, float &best_dist);

    /**
     *  @brief  Add all elements of an subtree to the closest elements. Used during the recSearch().
     *
     *  @param  current
     */
    void addSubtree(const KDTreeNodeT<DATA, DIM> *current);

    /**
     *  @brief  dist2
     *
     *  @param  a
     *  @param  b
     *
     *  @return dist2
     */
    float dist2(const KDTreeNodeInfoT<DATA, DIM> &a, const KDTreeNodeInfoT<DATA, DIM> &b) const;

    /**
     *  @brief  Frees the KDTree.
     */
    void clearTree();

    KDTreeNodeT<DATA, DIM> *root_;     ///< The KDTree root
    KDTreeNodeT<DATA, DIM> *nodePool_; ///< Node pool allows us to do just 1 call to new for each tree building
    int nodePoolSize_;                 ///< The node pool size
    int nodePoolPos_;                  ///< The node pool position

    std::vector<KDTreeNodeInfoT<DATA, DIM>> *closestNeighbour; ///< The closest neighbour
    std::vector<KDTreeNodeInfoT<DATA, DIM>> *initialEltList;   ///< The initial element list
};

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline KDTreeLinkerAlgo<DATA, DIM>::KDTreeLinkerAlgo() :
    root_(nullptr), nodePool_(nullptr), nodePoolSize_(-1), nodePoolPos_(-1), closestNeighbour(nullptr), initialEltList(nullptr)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline KDTreeLinkerAlgo<DATA, DIM>::~KDTreeLinkerAlgo()
{
    this->clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline void KDTreeLinkerAlgo<DATA, DIM>::build(std::vector<KDTreeNodeInfoT<DATA, DIM>> &eltList, const KDTreeBoxT<DIM> &region)
{
    if (eltList.size())
    {
        initialEltList = &eltList;
        const size_t mysize = initialEltList->size();

        nodePoolSize_ = mysize * 2 - 1;
        nodePool_ = new KDTreeNodeT<DATA, DIM>[nodePoolSize_];

        // Here we build the KDTree
        root_ = this->recBuild(0, mysize, 0, region);
        initialEltList = nullptr;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline int KDTreeLinkerAlgo<DATA, DIM>::medianSearch(int low, int high, int treeDepth)
{
    // We should have at least 1 element to calculate the median...
    //assert(low < high);

    const int nbrElts = high - low;
    int median = nbrElts / 2 - (1 - 1 * (nbrElts & 1));
    median += low;

    int l = low;
    int m = high - 1;

    while (l < m)
    {
        KDTreeNodeInfoT<DATA, DIM> elt = (*initialEltList)[median];
        int i = l;
        int j = m;

        do
        {
            // The even depth is associated to dim1 dimension, the odd one to dim2 dimension
            const unsigned thedim = treeDepth % DIM;
            while ((*initialEltList)[i].dims[thedim] < elt.dims[thedim])
                ++i;
            while ((*initialEltList)[j].dims[thedim] > elt.dims[thedim])
                --j;

            if (i <= j)
            {
                std::swap((*initialEltList)[i], (*initialEltList)[j]);
                i++;
                j--;
            }
        } while (i <= j);

        if (j < median)
            l = i;
        if (i > median)
            m = j;
    }

    return median;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline void KDTreeLinkerAlgo<DATA, DIM>::search(const KDTreeBoxT<DIM> &trackBox, std::vector<KDTreeNodeInfoT<DATA, DIM>> &recHits)
{
    if (root_)
    {
        closestNeighbour = &recHits;
        this->recSearch(root_, trackBox);
        closestNeighbour = nullptr;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline void KDTreeLinkerAlgo<DATA, DIM>::recSearch(const KDTreeNodeT<DATA, DIM> *current, const KDTreeBoxT<DIM> &trackBox)
{
    // By construction, current can't be null
    //assert(current != 0);
    // By Construction, a node can't have just 1 son.
    //assert (!(((current->left == 0) && (current->right != 0)) || ((current->left != 0) && (current->right == 0))));

    if ((current->left == nullptr) && (current->right == nullptr))
    {
        // Leaf case
        // If point inside the rectangle/area
        bool isInside = true;

        for (unsigned i = 0; i < DIM; ++i)
        {
            const auto thedim = current->info.dims[i];
            isInside = isInside && thedim >= trackBox.dimmin[i] && thedim <= trackBox.dimmax[i];
        }

        if (isInside)
            closestNeighbour->push_back(current->info);
    }
    else
    {
        // Node case
        // If region( v->left ) is fully contained in the rectangle
        bool isFullyContained = true;
        bool hasIntersection = true;

        for (unsigned i = 0; i < DIM; ++i)
        {
            const auto regionmin = current->left->region.dimmin[i];
            const auto regionmax = current->left->region.dimmax[i];
            isFullyContained = isFullyContained && (regionmin >= trackBox.dimmin[i] && regionmax <= trackBox.dimmax[i]);
            hasIntersection = hasIntersection && (regionmin < trackBox.dimmax[i] && regionmax > trackBox.dimmin[i]);
        }

        if (isFullyContained)
        {
            this->addSubtree(current->left);
        }
        else if (hasIntersection)
        {
            this->recSearch(current->left, trackBox);
        }

        //if region( v->right ) is fully contained in the rectangle
        isFullyContained = true;
        hasIntersection = true;

        for (unsigned i = 0; i < DIM; ++i)
        {
            const auto regionmin = current->right->region.dimmin[i];
            const auto regionmax = current->right->region.dimmax[i];
            isFullyContained = isFullyContained && (regionmin >= trackBox.dimmin[i] && regionmax <= trackBox.dimmax[i]);
            hasIntersection = hasIntersection && (regionmin < trackBox.dimmax[i] && regionmax > trackBox.dimmin[i]);
        }

        if (isFullyContained)
        {
            this->addSubtree(current->right);
        }
        else if (hasIntersection)
        {
            this->recSearch(current->right, trackBox);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline void KDTreeLinkerAlgo<DATA, DIM>::findNearestNeighbour(
    const KDTreeNodeInfoT<DATA, DIM> &point, const KDTreeNodeInfoT<DATA, DIM> *&result, float &distance)
{
    if (nullptr != result || distance != std::numeric_limits<float>::max())
    {
        result = nullptr;
        distance = std::numeric_limits<float>::max();
    }

    if (root_)
    {
        const KDTreeNodeT<DATA, DIM> *best_match = nullptr;
        this->recNearestNeighbour(0, root_, point, best_match, distance);

        if (distance != std::numeric_limits<float>::max())
        {
            result = &(best_match->info);
            distance = std::sqrt(distance);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline void KDTreeLinkerAlgo<DATA, DIM>::recNearestNeighbour(unsigned int depth, const KDTreeNodeT<DATA, DIM> *current,
    const KDTreeNodeInfoT<DATA, DIM> &point, const KDTreeNodeT<DATA, DIM> *&best_match, float &best_dist)
{
    const unsigned int current_dim = depth % DIM;

    if (current->left == nullptr && current->right == nullptr)
    {
        best_match = current;
        best_dist = this->dist2(point, best_match->info);
        return;
    }
    else
    {
        const float dist_to_axis = point.dims[current_dim] - current->info.dims[current_dim];

        if (dist_to_axis < 0.f)
        {
            this->recNearestNeighbour(depth + 1, current->left, point, best_match, best_dist);
        }
        else
        {
            this->recNearestNeighbour(depth + 1, current->right, point, best_match, best_dist);
        }

        // If we're here we're returned so best_dist is filled. Compare to this node and see if it's a better match. If it is, update result
        const float dist_current = this->dist2(point, current->info);

        if (dist_current < best_dist)
        {
            best_dist = dist_current;
            best_match = current;
        }

        // Now we see if the radius to best crosses the splitting axis
        if (best_dist > dist_to_axis * dist_to_axis)
        {
            // if it does we traverse the other side of the axis to check for a new best
            const KDTreeNodeT<DATA, DIM> *check_best = best_match;
            float check_dist = best_dist;

            if (dist_to_axis < 0.f)
            {
                this->recNearestNeighbour(depth + 1, current->right, point, check_best, check_dist);
            }
            else
            {
                this->recNearestNeighbour(depth + 1, current->left, point, check_best, check_dist);
            }

            if (check_dist < best_dist)
            {
                best_dist = check_dist;
                best_match = check_best;
            }
        }
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline void KDTreeLinkerAlgo<DATA, DIM>::addSubtree(const KDTreeNodeT<DATA, DIM> *current)
{
    // By construction, current can't be null
    //assert(current != 0);

    if ((current->left == nullptr) && (current->right == nullptr))
    {
        // Leaf case
        closestNeighbour->push_back(current->info);
    }
    else
    {
        // Node case
        this->addSubtree(current->left);
        this->addSubtree(current->right);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline float KDTreeLinkerAlgo<DATA, DIM>::dist2(const KDTreeNodeInfoT<DATA, DIM> &a, const KDTreeNodeInfoT<DATA, DIM> &b) const
{
    double d = 0.;

    for (unsigned i = 0; i < DIM; ++i)
    {
        const double diff = a.dims[i] - b.dims[i];
        d += diff * diff;
    }

    return (float)d;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline void KDTreeLinkerAlgo<DATA, DIM>::clearTree()
{
    delete[] nodePool_;
    nodePool_ = nullptr;
    root_ = nullptr;
    nodePoolSize_ = -1;
    nodePoolPos_ = -1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline bool KDTreeLinkerAlgo<DATA, DIM>::empty()
{
    return (nodePoolPos_ == -1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline int KDTreeLinkerAlgo<DATA, DIM>::size()
{
    return (nodePoolPos_ + 1);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline void KDTreeLinkerAlgo<DATA, DIM>::clear()
{
    if (root_)
        this->clearTree();
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline KDTreeNodeT<DATA, DIM> *KDTreeLinkerAlgo<DATA, DIM>::getNextNode()
{
    ++nodePoolPos_;

    // The tree size is exactly 2 * nbrElts - 1 and this is the total allocated memory.
    // If we have used more than that....there is a big problem.
    //assert(nodePoolPos_ < nodePoolSize_);

    return &(nodePool_[nodePoolPos_]);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <typename DATA, unsigned DIM>
inline KDTreeNodeT<DATA, DIM> *KDTreeLinkerAlgo<DATA, DIM>::recBuild(int low, int high, int depth, const KDTreeBoxT<DIM> &region)
{
    const int portionSize = high - low;

    // By construction, portionSize > 0 can't happen.
    //assert(portionSize > 0);

    if (portionSize == 1)
    {
        // Leaf case
        KDTreeNodeT<DATA, DIM> *leaf = this->getNextNode();
        leaf->setAttributs(region, (*initialEltList)[low]);
        return leaf;
    }
    else
    {
        // The even depth is associated to dim1 dimension, the odd one to dim2 dimension
        int medianId = this->medianSearch(low, high, depth);

        // We create the node
        KDTreeNodeT<DATA, DIM> *node = this->getNextNode();
        node->setAttributs(region);
        node->info = (*initialEltList)[medianId];

        // Here we split into 2 halfplanes the current plane
        KDTreeBoxT<DIM> leftRegion = region;
        KDTreeBoxT<DIM> rightRegion = region;

        const unsigned thedim = depth % DIM;
        auto medianVal = (*initialEltList)[medianId].dims[thedim];
        leftRegion.dimmax[thedim] = medianVal;
        rightRegion.dimmin[thedim] = medianVal;

        ++depth;
        ++medianId;

        // We recursively build the son nodes
        node->left = this->recBuild(low, medianId, depth, leftRegion);
        node->right = this->recBuild(medianId, high, depth, rightRegion);
        return node;
    }
}

} // namespace lar_content

#endif // LAR_KD_TREE_LINKER_ALGO_TEMPLATED_H
