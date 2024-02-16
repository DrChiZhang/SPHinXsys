#include <gtest/gtest.h>
#include "sphinxsys.h"
#include "polygon_mesh.h"

using namespace SPH;

TEST(test_poly_mesh, test_loader)
{
	PolygonMesh poly_mesh;
	bool is_loaded = poly_mesh.loadFile( "input/bunny.obj" );
	int num_face = poly_mesh.NumVertices();

	poly_mesh.scaleMesh( 0.5 );
	bool is_written = poly_mesh.writeFile( "output/bunny-scaled.obj" );

	ASSERT_TRUE( is_loaded );
	ASSERT_TRUE( is_written );
	EXPECT_EQ( 5550, num_face );
}

TEST(test_poly_mesh, test_box)
{
	PolygonMesh poly_mesh;
	auto box_mesh = poly_mesh.createBox( );
	bool is_written = box_mesh->writeFile( "output/box.obj" );

	ASSERT_TRUE( is_written );
}

TEST(test_poly_mesh, test_sphere)
{
	PolygonMesh poly_mesh;
	auto sphere_mesh = poly_mesh.createSphere( );
	bool is_written = sphere_mesh->writeFile( "output/sphere.obj" );

	ASSERT_TRUE( is_written );
}

TEST(test_poly_mesh, test_cylinder)
{
	PolygonMesh poly_mesh;
	auto cylinder_mesh = poly_mesh.createCylinder( );
	bool is_written = cylinder_mesh->writeFile( "output/cylinder.obj" );

	ASSERT_TRUE( is_written );
}

TEST(test_poly_mesh, test_cone)
{
	PolygonMesh poly_mesh;
	auto cone_mesh = poly_mesh.createCone( );
	bool is_written = cone_mesh->writeFile( "output/cone.obj" );

	ASSERT_TRUE( is_written );
}

TEST(test_poly_mesh, test_find_nearest_point)
{
	PolygonMesh poly_mesh;
	auto sphere_mesh = poly_mesh.createSphere( );

	Vecd point = Vecd(2.0, 0.0, 0.0);
	bool inside = true;
	Vecd normal = Vecd::Zero(); 

	Vecd nearestpoint = sphere_mesh->findNearestPoint( point, normal, inside );
	
	Real dist = ( nearestpoint - point ).squaredNorm();
	
	EXPECT_EQ( 1.0, dist );
	ASSERT_TRUE( !inside );
}

TEST(test_poly_mesh, test_find_nearest_point_bunny)
{
	PolygonMesh poly_mesh;
	bool is_loaded = poly_mesh.loadFile( "input/bunny.obj" );

	Vecd point = Vecd(-3.38, -4.5, 28.0);
	bool inside = false;
	Vecd normal = Vecd::Zero(); 

	Vecd nearestpoint = poly_mesh.findNearestPoint( point, normal, inside );
	
	Real dist = ( nearestpoint - point ).squaredNorm();
	bool is_written = poly_mesh.writeFile( "output/bunny.obj" );

	EXPECT_EQ( 1.0, dist );
	ASSERT_TRUE( is_loaded );
	ASSERT_TRUE( inside );
	ASSERT_TRUE( is_written );
}

int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
