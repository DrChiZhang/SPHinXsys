#include "polygon_mesh.h"

//=====================================================================================================//
namespace SPH
{
	//=================================================================================================//
    bool PolygonMesh::loadFile( const std::string& pathname ) 
    {
        bool is_loaded = false;
        fs::path filePath(pathname);

        if ( !filePath.has_extension() )
        {
                std::cout << "PolygonMesh::loadFile():  " << pathname << " has no extension " << std::endl;
                exit( EXIT_FAILURE );
        } else
        {
            if ( filePath.extension() == ".obj" )
            {
                std::cout << " Loading a obj file ..." << std::endl;
                loadObjFile(pathname);
                is_loaded = true;
            }else if( filePath.extension() == ".vtp" )
            {
                std::cout << " Loading a vtp file ..." << std::endl;
                loadVtpFile(pathname);    
                is_loaded = true;           
            }else if( filePath.extension() == ".stl" )
            {
                std::cout << " Loading a stl file ..." << std::endl;
                loadStlFile(pathname);  
                is_loaded = true;
            }else {
                std::cout << "PolygonMesh::loadFile() Unrecognized file extension: " <<   pathname << std::endl;
                exit( EXIT_FAILURE );
            }
        }
        
        bool is_manifold = isVertexManifold();
        if( is_loaded  && is_manifold )
        {
            computeTriangleNormals();
            buildKDTree();
        }

        return is_loaded && is_manifold;
    }
    //=================================================================================================//
    void PolygonMesh::loadObjFile( 	const std::string &pathname )
    {
        std::ifstream ifs(pathname);

        if ( !ifs.good() )
		{
			std::cout << "ObjParser::loadObjFile() failed to open file: " <<  pathname << std::endl;
			exit( EXIT_FAILURE );
		}

        loadObjFile(ifs);
        ifs.close();
    } 	
    //=================================================================================================//
    void PolygonMesh::loadVtpFile( 	const std::string &pathname )
    {
        std::ifstream ifs(pathname);

        if ( !ifs.good() )
		{
				std::cout << "ObjParser::loadObjFile() failed to open file: " <<  pathname << std::endl;
				exit( EXIT_FAILURE );
		}

        ifs.close();
    }
    //=================================================================================================//
    void PolygonMesh::loadStlFile( 	const std::string &pathname )
        {
        std::ifstream ifs(pathname);

        if ( !ifs.good() )
		{
			std::cout << "ObjParser::loadObjFile() failed to open file: " <<  pathname << std::endl;
			exit( EXIT_FAILURE );
		}
        
        ifs.close();
    }
    //=================================================================================================//
    bool PolygonMesh::writeFile( const std::string& pathname ) 
    {
         bool is_written = false;
        fs::path filePath(pathname);

        if ( !filePath.has_extension() )
        {
            std::cout << "PolygonMesh::writeFile():  " << pathname << " has no extension " << std::endl;
            exit( EXIT_FAILURE );
        } else
        {
            if ( fs::exists( pathname ) )
			{
				fs::remove( pathname );
			}

            if ( filePath.extension() == ".obj" )
            {
                std::cout << " Writing a obj file ..." << std::endl;
                writeObjFile(pathname);
                is_written = true;
            }else if( filePath.extension() == ".vtp" )
            {
                std::cout << " Writing a vtp file ..." << std::endl;
                writeVtpFile(pathname);    
                is_written = true;           
            }else if( filePath.extension() == ".stl" )
            {
                std::cout << " Writing a stl file ..." << std::endl;
                writeStlFile(pathname);  
                is_written = true;
            }else {
                std::cout << "PolygonMesh::loadFile() Unrecognized file extension: " <<   pathname << std::endl;
                exit( EXIT_FAILURE );
            }
        }

        return is_written;
    }
    //=================================================================================================//
    void PolygonMesh::writeObjFile( 	const std::string &pathname )
    {
        std::ofstream ofs( pathname );

        int num_vertices = NumVertices();
        for (int i = 0; i != num_vertices; i++)
            ofs << "v " << vertices_[i][0] << " " << vertices_[i][1] << " " <<  vertices_[i][2] << "\n";
  
        for (int i = 0; i != vertex_textures_.size(); i++)
            ofs << "vt " << vertex_textures_[i][0] << " " << vertex_textures_[i][1] << "\n";

        for (int i = 0; i != vertex_normals_.size(); i++)
            ofs << "vn " << vertex_normals_[i][0] << " " << vertex_normals_[i][1] << " " <<  vertex_normals_[i][2] << "\n";
            

        int num_triangle = Numtriangles();
        for (int i = 0; i != num_triangle; i++)
        {
            if( has_vetex_texture_ && has_vertex_normal_ )
            {
                ofs << "f " << triangles_[i][0] + 1 << "/" << triangle_texture_id_[i][0] + 1 << "/" << triangle_normal_id_[i][0] + 1 
                    << " "  << triangles_[i][1] + 1 << "/" << triangle_texture_id_[i][1] + 1 << "/" << triangle_normal_id_[i][1] + 1
                    << " "  << triangles_[i][2] + 1 << "/" << triangle_texture_id_[i][2] + 1 << "/" << triangle_normal_id_[i][2] + 1
                    << "\n"; 
            }
            else if ( has_vetex_texture_ && !has_vertex_normal_ )
            {
                ofs << "f " << triangles_[i][0] + 1 << "/" << triangle_texture_id_[i][0] + 1
                    << " "  << triangles_[i][1] + 1 << "/" << triangle_texture_id_[i][1] + 1
                    << " "  << triangles_[i][2] + 1 << "/" << triangle_texture_id_[i][2] + 1
                    << "\n"; 
            }
            else if ( has_vertex_normal_ && !has_vetex_texture_ )
            {
                ofs << "f " << triangles_[i][0] + 1 << "//" << triangle_normal_id_[i][0] + 1
                    << " "  << triangles_[i][1] + 1 << "//" << triangle_normal_id_[i][1] + 1
                    << " "  << triangles_[i][2] + 1 << "//" << triangle_normal_id_[i][2] + 1
                    << "\n"; 
            }
            else 
            {
                ofs << "f " << triangles_[i][0] + 1  
                    << " "  << triangles_[i][1] + 1
                    << " "  << triangles_[i][2] + 1
                    << "\n"; 
            }
        }

        ofs.close();
    }
    //=================================================================================================//	
    void PolygonMesh::writeVtpFile( 	const std::string &pathname )
    {
        std::ofstream ofs( pathname );
        // TO DO
        ofs.close();
    }
    //=================================================================================================//
    void PolygonMesh::writeStlFile( 	const std::string &pathname )
    {
        std::ofstream ofs( pathname );
        // TO DO
        ofs.close();
    }
    //=================================================================================================//
    void PolygonMesh::normalizeNormals() 
    {
        for ( size_t i = 0; i < triangle_normals_.size(); i++ ) 
        {
            triangle_normals_[i].normalize();
            if ( std::isnan( triangle_normals_[i](0) ) ) 
            {
                triangle_normals_[i] = Vecd::Identity();
            }
        }
    }
    //=================================================================================================//
    void PolygonMesh::computeTriangleNormals( bool normalized ) 
    {
        triangle_normals_.resize( Numtriangles() );

        for ( size_t i = 0; i < Numtriangles(); i++ ) 
        {
            auto &triangle = triangles_[i];
            Vecd v01 = vertices_[triangle(1)] - vertices_[triangle(0)];
            Vecd v02 = vertices_[triangle(2)] - vertices_[triangle(0)];
            triangle_normals_[i] = v01.cross(v02);
        }

        if ( normalized ) 
        {
            normalizeNormals();
        }
    }
    //=================================================================================================//
    void PolygonMesh::computeVertexNormals( bool normalized ) 
    {
        computeTriangleNormals( normalized );

        if( !has_vertex_normal_ )
        {
            vertex_normals_.resize( vertices_.size(), Vecd::Zero() );
            triangle_normal_id_.resize( triangles_.size(), Veci::Zero() );

            for ( size_t i = 0; i < triangles_.size(); i++ ) 
            {
                auto &triangle = triangles_[i];
                vertex_normals_[triangle(0)] += triangle_normals_[i];
                vertex_normals_[triangle(1)] += triangle_normals_[i];
                vertex_normals_[triangle(2)] += triangle_normals_[i];
                triangle_normal_id_[i] = Veci( triangle(0), triangle(1), triangle(2) );
            }

            has_vertex_normal_ = true;            
        }
    }
    //=================================================================================================//
    void PolygonMesh::buildKDTree() 
    { 
        kdtree_ = kdtree_ptr_keeper_.createPtr<KDTree>( *this );
    }
    //=================================================================================================//
    std::vector<int> PolygonMesh::findTriangles(const int & vidx)
    {
        std::vector<int> triangles;
        for ( int i = 0; i < triangles_.size(); i++ )  
        {
            const auto &triangle = triangles_[i];
            if ( triangle(0) == vidx || triangle(1) == vidx  || triangle(2) == vidx )
                triangles.push_back( i );
        }

        return triangles;
    }
    //=================================================================================================//
    Vecd PolygonMesh::findNearestPoint( 	const Vecd &position, Vecd &normal, bool &inside )
    {
        Real distance = MaxRealNumber;
        Vecd nearestpoint = Vecd::Zero();
        int face; 

        int knn = 1;
        std::vector<int> indices;
        std::vector<Real> distances;

        int result = kdtree_->searchKNN(position, knn, indices, distances);

        std::vector<int> triangles;
        for (int i = 0; i < result; i++)
        {
            std::vector<int> triangles = findTriangles( indices[i] );
            for ( int j = 0; j < triangles.size(); j++ )
            {
                int tidx = triangles[j];
                const auto &triangle = triangles_[tidx];
                Vecd p_find = findNearestPointToTriangle( position, tidx );
                Real dist = ( p_find - position ).squaredNorm();

                if ( dist < distance ) 
                {
                    nearestpoint = p_find;
                    distance = dist;
                    face = tidx;
                }
            }
        }

        Vecd delta = position-nearestpoint;
        inside = ( delta.dot( triangle_normals_[face] ) < 0 );
        
        return nearestpoint;
    }
    //=================================================================================================//
    bool PolygonMesh::isInside( 	const Vecd &position )
    {
        Vecd norm = Vecd::Zero();
        bool is_inside = false;
        Vecd nearestpoint = findNearestPoint( position, norm, is_inside);
 
        return is_inside;
    }
    //=================================================================================================//
    std::shared_ptr<PolygonMesh> PolygonMesh::createSphere(
        Real radius 
        , int resolution )
    {
        auto mesh = std::make_shared<PolygonMesh>();

        if ( radius <= 0 ) {
            std::cout << "PolygonMesh::createSphere(): radius <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }
        if ( resolution <= 0 ) {
            std::cout << "PolygonMesh::createSphere(): resolution <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }

        mesh->vertices_.resize( 2 * resolution * (resolution - 1) + 2 );

        mesh->vertices_[0] = Vecd( 0.0, 0.0,  radius );
        mesh->vertices_[1] = Vecd( 0.0, 0.0, -radius );

        Real step = M_PI / (Real)resolution;
        for ( int i = 1; i < resolution; i++ ) 
        {
            Real alpha = step * i;
            int base = 2 + 2 * resolution * ( i - 1 );

            for ( int j = 0; j < 2 * resolution; j++ ) 
            {
                Real theta = step * j;
                mesh->vertices_[base + j] = radius *
                    Vecd( sin(alpha) * cos(theta), sin(alpha) * sin(theta), cos(alpha) );
            }
        }

        // Triangles for poles.
        for ( int j = 0; j < 2 * resolution; j++ ) 
        {
            int j1 = ( j + 1 ) % ( 2 * resolution );
            int base = 2;
            mesh->triangles_.push_back( Veci(0, base + j, base + j1) );
            base = 2 + 2 * resolution * ( resolution - 2 );
            mesh->triangles_.push_back( Veci(1, base + j1, base + j) );
        }

        // Triangles for non-polar region.
        for (int i = 1; i < resolution - 1; i++) 
        {
            int base1 = 2 + 2 * resolution * ( i - 1 );
            int base2 = base1 + 2 * resolution;
            for ( int j = 0; j < 2 * resolution; j++ ) 
            {
                int j1 = ( j + 1 ) % ( 2 * resolution );
                mesh->triangles_.push_back( Veci(base2 + j, base1 + j1, base1 + j) );
                mesh->triangles_.push_back( Veci(base2 + j, base2 + j1, base1 + j1));
            }
        }

        bool is_manifold = mesh->isVertexManifold();
        if( is_manifold )
        {
            mesh->computeVertexNormals();
            mesh->buildKDTree();
        }

        return mesh;
    }
    //=================================================================================================//
    std::shared_ptr<PolygonMesh> PolygonMesh::createBox(
        Real width 
        , Real height
        , Real depth
        , bool map_texture_to_each_triangle )
    {
        auto mesh = std::make_shared<PolygonMesh>();  

        if ( width <= 0 ) {
            std::cout << "PolygonMesh::createBox():  width <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }
        if ( height <= 0 ) {
            std::cout << "PolygonMesh::createBox():  height <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }
        if ( depth <= 0 ) {
            std::cout << "PolygonMesh::createBox():  depth <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }

        // Vertices.
        mesh->vertices_.resize(8);
        mesh->vertices_[0] = Vecd( 0.0, 0.0, 0.0 );
        mesh->vertices_[1] = Vecd( width, 0.0, 0.0 );
        mesh->vertices_[2] = Vecd( 0.0, 0.0, depth );
        mesh->vertices_[3] = Vecd( width, 0.0, depth );
        mesh->vertices_[4] = Vecd( 0.0, height, 0.0 );
        mesh->vertices_[5] = Vecd( width, height, 0.0 );
        mesh->vertices_[6] = Vecd( 0.0, height, depth );
        mesh->vertices_[7] = Vecd( width, height, depth );

        // Triangles.
        mesh->triangles_ = {{4, 7, 5}, {4, 6, 7}, {0, 2, 4}, {2, 6, 4},
                            {0, 1, 2}, {1, 3, 2}, {1, 5, 7}, {1, 7, 3},
                            {2, 3, 7}, {2, 7, 6}, {0, 4, 1}, {1, 4, 5}};

        bool is_manifold = mesh->isVertexManifold();
        if( is_manifold )
        {
            mesh->computeVertexNormals();
            mesh->buildKDTree();
        }

        return mesh;
    }
    //=================================================================================================//
    std::shared_ptr<PolygonMesh> PolygonMesh::createCylinder(
        Real radius 
        , Real height 
        , int resolution 
        , int split )
    {
        auto mesh = std::make_shared<PolygonMesh>();

        if ( radius <= 0 ) {
            std::cout << "PolygonMesh::createSphere():  radius <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }
        if ( height <= 0 ) {
            std::cout << "PolygonMesh::createSphere():  height <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }
        if ( resolution <= 0 ) {
            std::cout << "PolygonMesh::createSphere():  resolution <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }
        if ( split <= 0 ) {
            std::cout << "PolygonMesh::createSphere():  split <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }

        mesh->vertices_.resize( resolution * (split + 1) + 2 );
        mesh->vertices_[0] = Vecd( 0.0, 0.0, height * 0.5 );
        mesh->vertices_[1] = Vecd( 0.0, 0.0, -height * 0.5 );
        Real step = M_PI * 2.0 / (Real)resolution;
        Real h_step = height / (Real)split;
        for ( int i = 0; i <= split; i++ ) 
        {
            for ( int j = 0; j < resolution; j++ ) 
            {
                Real theta = step * j;
                mesh->vertices_[2 + resolution * i + j] =
                            Vecd( cos(theta) * radius, 
                                  sin(theta) * radius,
                                  height * 0.5 - h_step * i );
            }
        }

        // Triangles for top and bottom face.
        for ( int j = 0; j < resolution; j++ ) 
        {
            int j1 = (j + 1) % resolution;
            int base = 2;
            mesh->triangles_.push_back( Veci(0, base + j, base + j1) );
            base = 2 + resolution * split;
            mesh->triangles_.push_back( Veci(1, base + j1, base + j) );
        }

        // Triangles for cylindrical surface.
        for ( int i = 0; i < split; i++ ) 
        {
            int base1 = 2 + resolution * i;
            int base2 = base1 + resolution;
            for ( int j = 0; j < resolution; j++ ) 
            {
                int j1 = (j + 1) % resolution;
                mesh->triangles_.push_back( Veci(base2 + j, base1 + j1, base1 + j) );
                mesh->triangles_.push_back( Veci(base2 + j, base2 + j1, base1 + j1) );
            }
        }

        bool is_manifold = mesh->isVertexManifold();
        if( is_manifold )
        {
            mesh->computeVertexNormals();
            mesh->buildKDTree();
        }

        return mesh;
    }
    //=================================================================================================//
    std::shared_ptr<PolygonMesh> PolygonMesh::createCone(
        Real radius 
        , Real height 
        , int resolution 
        , int split )
    {
        auto mesh = std::make_shared<PolygonMesh>();

        if ( radius <= 0 ) {
            std::cout << "PolygonMesh::createCone():  radius <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }
        if ( height <= 0 ) {
            std::cout << "PolygonMesh::createCone():  height <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }
        if ( resolution <= 0 ) {
            std::cout << "PolygonMesh::createCone():  resolution <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }
        if ( split <= 0 ) {
            std::cout << "PolygonMesh::createCone():  split <= 0 " << std::endl;
            exit( EXIT_FAILURE );
        }

        mesh->vertices_.resize(resolution * split + 2);

        mesh->vertices_[0] = Vecd( 0.0, 0.0, 0.0 );
        mesh->vertices_[1] = Vecd( 0.0, 0.0, height );

        Real step = M_PI * 2.0 / (Real)resolution;
        Real h_step = height / (Real)split;
        Real r_step = radius / (Real)split;

        std::unordered_map<int64_t, std::pair<Real, Real>> map_vertices_to_uv;
        for ( int i = 0; i < split; i++ ) 
        {
            int base = 2 + resolution * i;
            Real r = r_step * ( split - i );
            for ( int j = 0; j < resolution; j++ ) 
            {
                Real theta = step * j;
                mesh->vertices_[base + j] = Vecd( cos(theta) * r, sin(theta) * r, h_step * i );
            }
        }

        for ( int j = 0; j < resolution; j++ ) 
        {
            int j1 = ( j + 1 ) % resolution;
            // Triangles for bottom surface.
            int base = 2;
            mesh->triangles_.push_back( Veci(0, base + j1, base + j) );

            // Triangles for top segment of conical surface.
            base = 2 + resolution * ( split - 1 );
            mesh->triangles_.push_back( Veci(1, base + j, base + j1) );
        }

        // Triangles for conical surface other than top-segment.
        for ( int i = 0; i < split - 1; i++ ) 
        {
            int base1 = 2 + resolution * i;
            int base2 = base1 + resolution;
            for ( int j = 0; j < resolution; j++ ) 
            {
                int j1 = ( j + 1 ) % resolution;
                mesh->triangles_.push_back( Veci(base2 + j1, base1 + j, base1 + j1) );
                mesh->triangles_.push_back( Veci(base2 + j1, base2 + j, base1 + j) );
            }
        }

        bool is_manifold = mesh->isVertexManifold();
        if( is_manifold )
        {
            mesh->computeVertexNormals();
            mesh->buildKDTree();
        }

        return mesh;
    };
    //=================================================================================================//
}
