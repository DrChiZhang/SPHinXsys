/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	Polygon_mesh.h
 * @brief 	Polygon mesh, deal with OBJ file parser with writer, loader and parser.
 * @author	Chi Zhang
 */

#ifndef POLYGON_MESH_H
#define POLYGON_MESH_H

#include "base_data_type.h"
#include "large_data_containers.h"

#include <xml_parser.h>
#include "kdtree.h"

#include <iostream>
#include <string>
#include <cstdio>
#include <charconv>
#include <sstream>
#include <assert.h> 
#include <unordered_set>
#include <queue>

#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

namespace SPH
{
    inline std::vector<std::string> tokenize(const std::string& str, const std::string& delimiters = " ")
    {
        std::vector<std::string> tokens; 

        std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
        std::string::size_type pos = str.find_first_of(delimiters, lastPos);

        while (std::string::npos != pos || std::string::npos != lastPos)
        {
            tokens.push_back(str.substr(lastPos, pos - lastPos));
            lastPos = str.find_first_not_of(delimiters, pos);
            pos = str.find_first_of(delimiters, lastPos);
        }

        return tokens;
    }
    /**
	 * @class   PolygonMesh 
	 * @note    The PolygonMesh contains a class of format paser, e.g. stl, obj and  vtp. 
	 */
	class PolygonMesh 
	{
    private:
        UniquePtrKeeper<KDTree> kdtree_ptr_keeper_;
        KDTree* kdtree_;
	protected:
        StdLargeVec<Vec2d> vertex_textures_;
        StdLargeVec<Vec3d> vertex_normals_;

        StdLargeVec<Array3i> triangles_;
        StdLargeVec<Array3i> triangle_texture_id_;
        StdLargeVec<Array3i> triangle_normal_id_;
        StdLargeVec<Vec3d> triangle_normals_;

        bool has_vetex_texture_ = false;
        bool has_vertex_normal_  = false;

    public:
        StdLargeVec<Vec3d> vertices_;

	public:
		/** Default Constructor.  */
		PolygonMesh(){};
        /** Default Destructor. */
		virtual ~PolygonMesh()
        {
            vertex_textures_.clear();
            vertex_normals_.clear();

            triangles_.clear();
            triangle_texture_id_.clear();
            triangle_normal_id_.clear();
            triangle_normals_.clear();
        };

        /** 
         * Get the number of vertices parsered in the given file. 
         * Returns
         * int Number of vertices. 
        */
        size_t NumVertices() const{ return vertices_.size(); }

        /** 
         * Get the number of triangles parsered in the given file. 
         * Returns
         * int Number of triangles. 
        */
        size_t Numtriangles() const { return triangles_.size();}

        /** 
         * Get the specifed vertices parsered in the given file. 
         * Parameters
         * [in]	    vertex_idex	Vertex idex. 
         * Returns
         * Vec3d Vertex postion.  
        */
        Vec3d getVertexPosition(int vertex_idex) const
        {
            assert(0 <= vertex_idex && vertex_idex < NumVertices());
            return vertices_[vertex_idex];
        }

        /** 
         * Scale the  Polygon mesh with given paramer. 
         * Parameters
         * [in]	    scale	Scale. 
        */
        void scaleMesh( const Real scale) 
        {
            size_t num_vertices = NumVertices();
            for (size_t i = 0; i != num_vertices; i++)
                vertices_[i] *= scale;
        }

        /** 
         * Translate the  Polygon mesh with given vector. 
         * Parameters
         * [in]	    translation	Translation parameter.  
        */
        void translateMesh( const Vec3d translation ) 
        {
            size_t num_vertices = NumVertices();
            for (size_t i = 0; i != num_vertices; i++)
                vertices_[i] += translation;
        }

        /** 
         * Transform the  Polygon mesh with given matrix. 
         * Parameters
         * [in]	    transform	Transform matrix. 
        */
        void transformMesh( const Transform &transform )
        {
            // int num_vertices = NumVertices();
            // for (int i = 0; i != num_vertices; i++)
            //     vertices_[i] = transform * vertices_[i];         
        } 	

        /** 
         * Given the path of the a file, check the extension and load it if the extension is valid. 
         * Parameters
         * [in]	    pathname	Path to the file. 
        */
        bool loadFile( const std::string& pathname );

        /** 
         * Given the path of the Obj file, read the data and store it for future usage. 
         * Parameters
         * [in]	    pathname	Path to the Obj file. 
        */
        void loadObjFile( 	const std::string &pathname ); 	

        /** 
         * Given the path of the Obj file, read the data and store it for future usage. 
         * Parameters
         * [in]	    filestream	ifstream of the importing path file. 
        */
        inline void loadObjFile( std::ifstream& filestream );	

        /** 
         * Given the path of the Obj file, read the data and store it for future usage. 
         * Parameters
         * [in]	    pathname	Path to the Obj file. 
        */
        void loadVtpFile( 	const std::string &pathname ); 

        /** 
         * Given the path of the Obj file, read the data and store it for future usage. 
         * Parameters
         * [in]	    pathname	Path to the Obj file. 
        */
        void loadStlFile( 	const std::string &pathname );

        /** 
         * Given the path of the a file, check the extension and write it in specified format. 
         * Parameters
         * [in]	    pathname	Path to the writing file. 
        */
        bool writeFile( const std::string& pathname );

        /** 
         * Given the path of the file, write the data in obj format. 
         * Parameters
         * [in]	    pathname	Path to the file. 
        */
        void writeObjFile( 	const std::string &pathname ); 		

        /** 
         * Given the path of the file, write the data in vtp format. 
         * Parameters
         * [in]	    pathname	Path to the file. 
         */
        void writeVtpFile( 	const std::string &pathname ); 

        /** 
         * Given the path of the file, write the data in stl format. 
         * Parameters
         * [in]	    pathname	Path to the file. 
         */
        void writeStlFile( 	const std::string &pathname );
         /** 
         * Function to Normalize both triangle normals and vertex normals to length 1. 
         */
        void normalizeNormals();

        /** 
         * Function to compute triangle normals. 
         * Parameters
         * [in]	    normalized	Normalize the normals. 
         */
        void computeTriangleNormals( bool normalized = true );

        /** 
         * Function to compute vertex normals. 
         * Parameters
         * [in]	    normalized	Normalize the normals. 
         */
        void computeVertexNormals( bool normalized = true );

        /** 
         * Function that checks if all vertices in the triangle mesh are manifold.
         * \note A vertex is manifold if its star is edge‐manifold and edge‐connected.
         * (Two or more faces connected only by a vertex and not by an edge.)
         * 
         * Returns
         * Ture or not . 
         */
        bool isVertexManifold() { return getNonManifoldVertices().empty(); }

        /** 
         * Function that returns a list of non-manifold vertex indices.
         * \note  A vertex is manifold if its star is edge‐manifold and edge‐connected.
         *  (Two or more faces connected only by a vertex and not by an edge.)
         * 
         * Returns
         * Vectore of non manifold vertices.
         */
        inline std::vector<int> getNonManifoldVertices();

        /** 
         * Function for building a KD tree for nearest neighbor search.
         * 
         */
        void buildKDTree();

        /** 
         * Function that returns a list of triangles that connected to the given vertex. 
         * 
         * Parameters
         * [in] vidx The Given index of vertex
         * Returns
         * Vectore of triangles.
         */
        std::vector<int> findTriangles(const int & vidx);

        /** 
         * Given a point, find the nearest point on the surtriangle of this object presented by the Obj file.
         * If multiple points on the surtriangle are equally close to the given point, this may return any of them.
         * Parameters
         * [in]	    position	The point in question.
         * [out]	inside	    On exit, this is set to true if the given point is inside this object, false otherwise.
         * [out]	normal	    On exit, this contains the surtriangle normal at the returned point.
         * Returns
         * The point on the surtriangle of the object which is closest to the given point. 
         */
        Vec3d findNearestPoint( const Vec3d &position, Vec3d &normal, bool &inside ); 	

        /** 
         * Given a point, check the given point is inside or outside the object presented by the Obj file.
         * Parameters
         * [in]	    position	The point in question.
         * Returns
         * TRUE or FALSE, inside or not. 
         */
        bool isInside( 	const Vec3d &position );

        /** 
         * Given a point, calculate the distance between the point in space and a triangle.
         * \note This algorithm is based on a description by David Eberly found at 
         *  http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf.
         * Parameters
         * [in]	    position	The point in question.
         * [in]     triangle    The triangle in question.
         * 
         * Returns
         * The nearest point on the face. 
         */
        inline Vec3d findNearestPointToTriangle( const Vec3d &position, int face); 

        /** 
         * Gcreate a sphere mesh. 
         * Parameters
         * [in]	    radius Defines radius of the sphere.
         * [in]     resolution Defines the resolution of the sphere.
         * Return 
         * PolygonMesh* SmartPointer to the generated mesh. 
         */
        std::shared_ptr<PolygonMesh> createSphere(
              Real radius = 1.0
            , int resolution = 20);

        /** 
         * Gcreate a box mesh. 
         * Parameters
         * [in]	    width Defines x-directional length.
         * [in]     height Defines y-directional length.
         * [in]	    depth Defines z-directional length.
         * [in]     texture_to_triangle If true, maps the entire texture image
         *          to each triangle. If false, sets the default uv map to the mesh.
         * Return 
         * PolygonMesh* SmartPointer to the generated mesh. 
         */
        std::shared_ptr<PolygonMesh> createBox(
              Real width = 1.0
            , Real height = 1.0
            , Real depth = 1.0
            , bool map_texture_to_each_triangle = false);

        /** 
         * Gcreate a Cone mesh. 
         * Parameters
         * [in]	    radius defines the radius of the cylinder.
         * [in]     height defines the height of the cylinder.
         * [in]	    resolution defines that the circle will be split into resolution segments.
         * [in]     split defines that the height will be split into split segments.
         * Return 
         * PolygonMesh* SmartPointer to the generated mesh. 
         */
        std::shared_ptr<PolygonMesh> createCylinder(
              Real radius = 1.0
            , Real height = 2.0
            , int resolution = 20
            , int split = 4);

        /** 
         * Gcreate a Cone mesh. 
         * Parameters
         * [in]	    radius Defines the radius of the cone.
         * [in]     height Defines the height of the cone.
         * [in]	    resolution defines that the circle will be split into resolution segments.
         * [in]     split defines that the height will be split into split segments.
         * Return 
         * PolygonMesh* SmartPointer to the generated mesh. 
         */
        std::shared_ptr<PolygonMesh> createCone(
              Real radius = 1.0
            , Real height = 2.0
            , int resolution = 20
            , int split = 1);
    
    };

    /** 
     * Given the path of the Obj file, read the data and store it for future usage. 
     * Parameters
     * [in]	    filestream	ifstream of the importing file. 
     */
    inline void PolygonMesh::loadObjFile( std::ifstream& filestream ) 
    {
        std::string line;
        std::vector<std::string> sub_buffer;

        while( std::getline( filestream, line ) )
        {
            std::stringstream str_stream(line);
			std::string type_str;
			str_stream >> type_str;

			if ( type_str == "v" )
			{
                std::istringstream v_ss( line.substr( line.find("v") + 1 ) );
                Real x, y, z;
                v_ss >> x; v_ss >> y; v_ss >> z;
                vertices_.push_back( Vec3d( x, y, z ) );
			} 
            else if( type_str == "vt" )
            {
                std::istringstream vt_ss( line.substr( line.find("vt") + 2 ) );
                Real t_x, t_y;
                vt_ss >> t_x; vt_ss >> t_y;
                vertex_textures_.push_back( Vec2d( t_x, t_y ) );       
                has_vetex_texture_ = true;  
            } 
            else if( type_str == "vn" )
            {
                std::istringstream vn_ss( line.substr( line.find("vn") + 2 ) );
                Real n_x, n_y, n_z;
                vn_ss >> n_x; vn_ss >> n_y; vn_ss >> n_z;
                vertex_normals_.push_back( Vec3d( n_x, n_y, n_z ) );       
                has_vertex_normal_ = true;        
            } 
            else if( type_str == "f" )
            {
                Array3i triangle = Array3i::Zero();
                Array3i triangle_normal = Array3i::Zero();
                Array3i triangle_texture = Array3i::Zero();
                std::string parse_str = line.substr( line.find("f") + 1 );

                if( has_vetex_texture_ && has_vertex_normal_ )
                {
                    sub_buffer.clear();
                    sub_buffer = tokenize( parse_str );
					for ( int i = 0; i < 3; ++i )
					{
						std::vector<std::string> buffer = tokenize( sub_buffer[i], "/" );

						triangle[i] = stoi( buffer[0] ) - 1;
						triangle_texture[i] = stoi( buffer[1] ) - 1;
						triangle_normal[i] = stoi( buffer[2] ) - 1;
					}                  
                } 
                else if ( has_vetex_texture_ && !has_vertex_normal_ )
                {
                    sub_buffer.clear();
                    sub_buffer = tokenize( parse_str );
					for ( int i = 0; i < 3; ++i )
					{
						std::vector<std::string> buffer = tokenize( sub_buffer[i], "/" );

						triangle[i] = stoi( buffer[0] ) - 1;
						triangle_texture[i] = stoi( buffer[1] ) - 1;
					} 
                }      
                else if ( has_vertex_normal_ &&  !has_vetex_texture_)
                {
                    sub_buffer.clear();
                    sub_buffer = tokenize( parse_str );
					for (int i = 0; i < 3; ++i)
					{
						std::vector<std::string> buffer = tokenize( sub_buffer[i], "/" );

						triangle[i] = stoi( buffer[0] ) - 1;
						triangle_normal[i] = stoi( buffer[1] ) - 1;
					} 
                }  
                else
                {
                    sub_buffer.clear();
                    sub_buffer = tokenize( parse_str );
					for (int i = 0; i < 3; ++i)
					{
						triangle[i] = stoi( sub_buffer[i] ) - 1;
					} 
                }

                triangles_.push_back( triangle ); 
                if( has_vetex_texture_ )
                    triangle_texture_id_.push_back( triangle_texture );
                if( has_vertex_normal_ )   
                    triangle_normal_id_.push_back( triangle_normal );     
            }
        }
    };
    /** 
     * Function that returns a list of non-manifold vertex indices.
     * \note  A vertex is manifold if its star is edge‐manifold and edge‐connected.
     *  (Two or more faces connected only by a vertex and not by an edge.)
     * 
     * Returns
     * Vectore of non manifold vertices.
     */

    inline std::vector<int> PolygonMesh::getNonManifoldVertices()
    {
        std::vector<std::unordered_set<int>> vertex_to_triangles( vertices_.size() );
        for ( size_t i = 0; i < triangles_.size(); ++i ) 
        {
            const auto &tria = triangles_[i];

            vertex_to_triangles[tria(0)].emplace( int(i) );
            vertex_to_triangles[tria(1)].emplace( int(i) );
            vertex_to_triangles[tria(2)].emplace( int(i) );
        }

        std::vector<int> non_manifold_verts;
        for (int vidx = 0; vidx < int(vertices_.size()); ++vidx) 
        {
            const auto &triangles = vertex_to_triangles[vidx];
            if ( triangles.size() == 0 ) 
            {
                continue;
            }

            // collect edges and vertices
            std::unordered_map<int, std::unordered_set<int>> edges;
            for (int i : triangles) 
            {
                const auto &triangle = triangles_[i];
                if ( triangle(0) != vidx && triangle(1) != vidx ) 
                {
                    edges[triangle(0)].emplace(triangle(1));
                    edges[triangle(1)].emplace(triangle(0));
                } else if ( triangle(0) != vidx && triangle(2) != vidx ) 
                {
                    edges[triangle(0)].emplace(triangle(2));
                    edges[triangle(2)].emplace(triangle(0));
                } else if ( triangle(1) != vidx && triangle(2) != vidx ) 
                {
                    edges[triangle(1)].emplace(triangle(2));
                    edges[triangle(2)].emplace(triangle(1));
                }
            }

            // test if vertices are connected
            std::queue<int> next;
            std::unordered_set<int> visited;
            next.push(edges.begin()->first);
            visited.emplace(edges.begin()->first);
            while ( !next.empty() ) 
            {
                int vert = next.front();
                next.pop();

                for (auto nb : edges[vert]) {
                    if (visited.count(nb) == 0) {
                        visited.emplace(nb);
                        next.emplace(nb);
                    }
                }
            }
            if ( visited.size() != edges.size() ) 
            {
                non_manifold_verts.push_back( vidx );
            }
        }

        return non_manifold_verts;
    };

    /** 
     * Given a point, calculate the distance between the point in space and a triangle.
     * \note This algorithm is based on a description by David Eberly found at 
     *  http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf.
     * Parameters
     * [in]	    position	The point in question.
     * [in]     triangle    The triangle in question.
     * 
     * Returns
     * The nearest point on the face. 
     */
    inline Vec3d PolygonMesh::findNearestPointToTriangle( const Vec3d &position, int face)
    {
        const auto &triangle = triangles_[face];
        const Vec3d& vert1 = vertices_[triangle(0)];
        const Vec3d& vert2 = vertices_[triangle(1)];
        const Vec3d& vert3 = vertices_[triangle(2)];

        const Vec3d e0 = vert2 - vert1;
        const Vec3d e1 = vert3 - vert1;
        const Vec3d delta = vert1 - position;
        
        const Real a = e0.squaredNorm();
        const Real b = e0.dot( e1 );
        const Real c = e1.squaredNorm();
        const Real d = e0.dot( delta );
        const Real e = e1.dot( delta );
        //const Real f = delta.squaredNorm();
        const Real det = a * c - b * b;

        Real s = b * e - c * d;
        Real t = b * d - a * e;

        if ( s+t <= det ) {
            if ( s < 0 ) {
                if ( t < 0 ) {
                    // Region 4

                    if ( d < 0 ) {
                        s = ( -d >= a ? 1 : -d / a );
                        t = 0;
                    }
                    else {
                        s = 0;
                        t = ( e >= 0 ? 0 : ( -e >= c ? 1 : -e / c ));
                    }
                }
                else {
                    // Region 3

                    s = 0;
                    t = ( e >= 0 ? 0 : ( -e >= c ? 1 : -e / c ));
                }
            }
            else if ( t < 0 ) {
                // Region 5

                s = ( d >= 0 ? 0 : ( -d >= a ? 1 : -d / a ));
                t = 0;
            }
            else {
                // Region 0

                const Real invDet = Real(1)/det;
                s *= invDet;
                t *= invDet;
            }
        }
        else {
            if ( s < 0 ) {
                // Region 2

                Real temp0 = b + d;
                Real temp1 = c + e;
                if ( temp1 > temp0 ) {
                    Real numer = temp1 - temp0;
                    Real denom = a - 2 * b + c;
                    s = ( numer >= denom ? 1 : numer / denom );
                    t = 1 - s;
                }
                else {
                    s = 0;
                    t = (temp1 <= 0 ? 1 : (e >= 0 ? 0 : -e / c));
                }
            }
            else if ( t < 0 ) {
                // Region 6

                Real temp0 = b + e;
                Real temp1 = a + d;
                if ( temp1 > temp0 ) {
                    Real numer = temp1 - temp0;
                    Real denom = a - 2 * b + c;
                    t = ( numer >= denom ? 1 : numer / denom );
                    s = 1 - t;
                }
                else {
                    s = ( temp1 <= 0 ? 1 : ( e >= 0 ? 0 : -d / a ));
                    t = 0;
                }
            }
            else {
                // Region 1

                const Real numer = c + e - b - d;
                if ( numer <= 0 )
                    s = 0;
                else {
                    const Real denom = a - 2 * b + c;
                    s = ( numer >= denom ? 1 : numer / denom );
                }
                t = 1 - s;
            }
        }
        
        return vert1 + s * e0 + t * e1;
    };
}
#endif //POLYGON_MESH_H
