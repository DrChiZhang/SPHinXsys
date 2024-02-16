#include "kdtree.h"
#include "polygon_mesh.h"

//=====================================================================================================//
namespace SPH
{
    //=================================================================================================//
    int KDTree::searchKNN(
        const Vecd &query 
        , int knn
        , std::vector<int> &indices
        , std::vector<Real> &distance2) const
    {
        if ( data_.empty() || dataset_size_ <= 0 || size_t(query.rows()) != dimension_ || knn < 0 ) 
            return -1;

        indices.resize(knn);
        distance2.resize(knn);
        std::vector<Eigen::Index> indices_eigen(knn);
        int k = adaptor_ptr_->index_->knnSearch( query.data(), knn, indices_eigen.data(), distance2.data() );
        indices.resize(k);
        distance2.resize(k);
        std::copy_n( indices_eigen.begin(), k, indices.begin() );
        return k;
    }
    //=================================================================================================//
    int KDTree::searchRadius(
        const Vecd &query
        , Real radius
        , std::vector<int> &indices
        , std::vector<Real> &distance2) const 
    {
        if ( data_.empty() || dataset_size_ <= 0 || size_t(query.rows()) != dimension_ ) 
            return -1;

        std::vector<nanoflann::ResultItem<Eigen::Index, double>> indices_dists;
        int k = adaptor_ptr_->index_->radiusSearch( query.data(), radius * radius, indices_dists, nanoflann::SearchParameters(-1, 0.0) );
        indices.resize(k);
        distance2.resize(k);
        for ( int i = 0; i < k; ++i ) 
        {
            indices[i] = indices_dists[i].first;
            distance2[i] = indices_dists[i].second;
        }

        return k;
    }
    //=================================================================================================//
    int KDTree::searchHybrid(
        const Vecd &query
        , double radius
        , int max_nn
        , std::vector<int> &indices
        , std::vector<double> &distance2) const 
    {
        if ( data_.empty() || dataset_size_ <= 0 || size_t(query.rows()) != dimension_ || max_nn < 0 ) 
            return -1;

        distance2.resize(max_nn);
        std::vector<Eigen::Index> indices_eigen(max_nn);
        int k = adaptor_ptr_->index_->knnSearch( query.data(), max_nn, indices_eigen.data(), distance2.data() );
        k = std::distance( distance2.begin(), 
                               std::lower_bound( distance2.begin(), distance2.begin() + k, radius * radius ) 
                             );
        indices.resize(k);
        distance2.resize(k);
        std::copy_n(indices_eigen.begin(), k, indices.begin());

        return k;
    }
	//=================================================================================================//
    bool KDTree::setMatrixData(const Eigen::MatrixXd &data)
    {
        return setRawData( Eigen::Map<const Eigen::MatrixXd>( 
            data.data()
            , data.rows()
            , data.cols() 
        ) );
    }
    //=================================================================================================//
    bool KDTree::setGeometry(const PolygonMesh &polymesh)
    {
        return setRawData(Eigen::Map<const Eigen::MatrixXd>(
            (const Real *)((const PolygonMesh &)polymesh).vertices_.data()
            , 3
            , ((const PolygonMesh &)polymesh).vertices_.size()
        ));
    }
    //=================================================================================================//
}
