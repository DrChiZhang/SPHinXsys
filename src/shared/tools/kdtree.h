/* -------------------------------------------------------------------------*
 *                                SPHinXsys                                 *
 * -------------------------------------------------------------------------*
 * Copyright (c) 2024, The AUTHORS.                                         *
 * Authors: Chi Zhang                                                       *
 * Contributors:                                                            *
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at                                                   *
 *          http://www.apache.org/licenses/LICENSE-2.0.                     *
 * Unless required by applicable law or agreed to in writing, software      *
 * distributed under the License is distributed on an "AS IS" BASIS,        *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
 * See the License for the specific language governing permissions and      *
 * limitations under the License.                                           *
 * ------------------------------------------------------------------------*/

/**
 * @file 	  kdtree.h
 * @brief 	kd tree for neighbor searching. 
 * @author	Chi ZHang Email: zhangchi0118@gmail.com
 */

#pragma once

#include "data_type.h"
#include <nanoflann.hpp>

namespace nanoflann 
{
    struct metric_L2;
    template <class MatrixType, int DIM, class Distance, bool row_major>
    struct KDTreeEigenMatrixAdaptor;
}

namespace SPH
{
    class PolygonMesh;

    class KDTree
    {
    public:
        /** Default Constructor. */
        KDTree(){};
        /** 
         * Parameterized Constructor.
         * Paraeters
         * [in] data Provides set of data points for KDTree construction.
        */
        KDTree(const Eigen::MatrixXd &data) { setMatrixData(data); };
        /** 
         * Parameterized Constructor.
         * Paraeters
         * [in] geometry Provides geometry from which KDTree is constructed.
        */
        KDTree(const PolygonMesh &polymesh) { setGeometry( polymesh ); };
        /** Default Destructor. */
        ~KDTree(){};

    public:
        /** 
         * Sets the data for the KDTree from a eigen matrix.
         * Paraeters
         * [in] data Input Eigen matrix data. 
        */
        bool setMatrixData( const Eigen::MatrixXd &data );

        /** 
         * Sets the data for the KDTree from geometry, viz. Polygon mesh. 
         * Paraeters
         * [in] poly_mesh Input Polygon mesh. 
        */
        bool setGeometry( const PolygonMesh &polymesh );
        
        /**
         * Find the "num_closest" nearest neighbors to the query_point, here as Vecd.
         * Their indices and distances are stored in the provided pointers to
         * array/vector.
         *
         * \note If L2 norms are used, all returned distances are actually squared
         *       distances.
         *
         * \note Only the first `N` entries in `out_indices` and `out_distances`
         *       will be valid. Return is less than `num_closest` only if the
         *       number of elements in the tree is less than `num_closest`.
         * 
         * Returns
         * Number `N` of valid points in the result set.
         */
        int searchKNN(const Vecd &query 
                    , int knn
                    , std::vector<int> &indices
                    , std::vector<Real> &distance2) const;

        /**
         * Find all the neighbors to query_point, here as Vecd, within a maximum
         * radius. The output is given as a vector of pairs, of which the first
         * element is a point index and the second the corresponding distance.
         * Previous contents of IndicesDists are cleared.
         *
         *  If searchParams.sorted==true, the output list is sorted by ascending
         * distances.
         *
         *  For a better performance, it is advisable to do a .reserve() on the
         * vector if you have any wild guess about the number of expected matches.
         *
         *
         * \note If L2 norms are used, search radius and all returned distances
         *       are actually squared distances.
         * 
         * Returns
         * The number of points within the given radius (i.e. indices.size() or dists.size() )
         */
        int searchRadius(const Vecd &query
                        , Real radius
                        , std::vector<int> &indices
                        , std::vector<Real> &distance2) const;

        /**
         * This is optimized code for heavily repeated search.
         * \note It is also the recommended setting for search Other flann::Index::radiusSearch() implementations lose performance due
         *  to memory allocation/deallocation.
         * 
         * Returns
         * Number `N` of valid points in the result set.
         */
        int searchHybrid(const Vecd &query
                        , double radius
                        , int max_nn
                        , std::vector<int> &indices
                        , std::vector<double> &distance2) const;

        private:
        /**
         * Sets the KDTree data from the data provided by the other methods.
         * Parameters
         * [in] data Input data, e.g., Eigen matrix or Polygon mesh. 
        */
        inline bool setRawData(const Eigen::Map<const Eigen::MatrixXd> &data);
        
        protected:
        using Adaptor = nanoflann::KDTreeEigenMatrixAdaptor<
                                                            Eigen::Map<const Eigen::MatrixXd>
                                                            , -1
                                                            , nanoflann::metric_L2
                                                            , false>;

        std::vector<Real> data_;
        std::unique_ptr<Eigen::Map<const Eigen::MatrixXd>> data_ptr_;
        std::unique_ptr<Adaptor> adaptor_ptr_;
        size_t dimension_ = 0;
        size_t dataset_size_ = 0;
    };

    inline bool KDTree::setRawData(const Eigen::Map<const Eigen::MatrixXd> &data)
    {
        dimension_ = data.rows();
        dataset_size_ = data.cols();
        if (dimension_ == 0 || dataset_size_ == 0) 
        {
            std::cout << " KDTree::setRawData Failed due to no data." << std::endl;
            return false;
        }

        data_.resize( dataset_size_ * dimension_ );
        memcpy( data_.data(), data.data(), dataset_size_ * dimension_ * sizeof(Real) );
        data_ptr_.reset( new Eigen::Map<const Eigen::MatrixXd>(data) );
        adaptor_ptr_.reset( new Adaptor(dimension_, std::cref(*data_ptr_), 15) );
        adaptor_ptr_->index_->buildIndex();
        return true;
    };

}
