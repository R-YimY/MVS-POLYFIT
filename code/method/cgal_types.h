#ifndef _CGAL_TYPES_H_
#define _CGAL_TYPES_H_

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/intersections.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <CGAL/compute_average_spacing.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/property_map.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include<CGAL/Aff_transformation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h> 
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>v

#include "../math/math_types.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT Float;
typedef Kernel::FT FT;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Line_3 Line_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Iso_cuboid_3 BoundBox;
typedef CGAL::Polygon_2<Kernel> Polygon_2;

typedef CGAL::Polygon_with_holes_2<Kernel> Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2> Pwh_list_2;


//邻域搜索定义
typedef CGAL::Search_traits_3<Kernel> SearchTraits_3;
typedef CGAL::Orthogonal_k_neighbor_search<SearchTraits_3> Neighbor_search;
typedef Neighbor_search::Tree KD_Tree;


//形状检测定义
typedef std::pair<Point_3, Vector_3> Point_with_normal;
typedef std::vector<Point_with_normal> Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;


typedef std::array<unsigned char, 3> Color_;
typedef std::tuple<Point_3, Vector_3, Color_, double> Point_NCI;
typedef CGAL::Nth_of_tuple_property_map<0, Point_NCI> PointMap;
typedef CGAL::Nth_of_tuple_property_map<1, Point_NCI> NormalMap;
typedef CGAL::Nth_of_tuple_property_map<2, Point_NCI> ColorMap;
typedef CGAL::Nth_of_tuple_property_map<3, Point_NCI> ImageMap;

typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel> CDT2;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT2> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT2, Criteria> Meshing_engine;


inline Point_2 to_cgal_point(const vec2& p)
{
	return Point_2(p.x, p.y);
}
inline Vector_2 to_cgal_vector(const vec2& v)
{
	return Vector_2(v.x, v.y);
}

inline Point_3 to_cgal_point(const vec3& p)
{
	return Point_3(p.x, p.y, p.z);
}
inline Vector_3 to_cgal_vector(const vec3& v)
{
	return Vector_3(v.x, v.y, v.z);
}

inline vec2 to_my_point(const Point_2& p)
{
	return vec2(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
}
inline vec2 to_my_vector(const Vector_2& v)
{
	return vec2(CGAL::to_double(v.x()), CGAL::to_double(v.y()));
}

inline vec3 to_my_point(const Point_3& p)
{
	return vec3(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
}
inline vec3 to_my_vector(const Vector_3& v)
{
	return vec3(CGAL::to_double(v.x()), CGAL::to_double(v.y()), CGAL::to_double(v.z()));
}

inline Plane_3 to_cgal_plane(const Plane3d& plane)
{
	return Plane_3(to_cgal_point(plane.point()), to_cgal_vector(plane.normal()));
}

// the coordinate system is defined by (orig, base1, base2)
inline Point_2 convert_to_2d(const Point_3& orig, const Vector_3& base1, const Vector_3& base2, const Point_3& p)
{
	Vector_3 vec = p - orig;
	Float x = vec * base1;
	Float y = vec * base2;
	return Point_2(x, y);
}

inline Point_3 convert_to_3d(const Point_3& orig, const Vector_3& base1, const Vector_3& base2, const Point_2& p)
{
	return orig + base1 * p.x() + base2 * p.y();
}



//比较运算符
template <typename T>
bool jin(const T & a, const T & x, const T & b) {
	return (a <= x && x <= b);
}

template <typename T>
bool jinr(const T & a, const T & x, const T & b) {
	return (a < x && x < b);
}

template <typename T>
const T & jmin(const T & a, const T & b) {
	return (a < b ? a : b);
}

template <typename T>
const T & jmax(const T & a, const T & b) {
	return (a > b ? a : b);
}

template <typename T>
const T & jclamp(const T & a, const T & x, const T & b) {
	return (x < a ? a : (x > b ? b : x));
}



#endif