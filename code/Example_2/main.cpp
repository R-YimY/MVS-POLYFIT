
/*
平面提取
*/

#include "../model/point_set.h"
#include "../model/map.h"
#include "../model/map_io.h"
#include "../model/point_set_io.h"
#include"../method/PlaneDetection.h"
#include "../basic/file_utils.h"
#include<CGAL/IO/read_ply_points.h>

PointSet * dataread(std::string filepath) {

	std::ifstream in(filepath);
	std::vector<Point_NCI> points;

	if (CGAL::IO::read_PLY_with_properties(in, std::back_inserter(points),
		CGAL::make_ply_point_reader(PointMap()), std::make_tuple(ColorMap(),
			CGAL::Construct_array(),
			CGAL::IO::PLY_property<unsigned char>("red"),
			CGAL::IO::PLY_property<unsigned char>("green"),
			CGAL::IO::PLY_property<unsigned char>("blue")),
		CGAL::IO::make_ply_normal_reader(NormalMap()),
		std::make_pair(ImageMap(), CGAL::PLY_property<double>("Image"))
	)) {
		std::cout << "pointcloud: " << points.size() << std::endl;
	}

	PointSet * pts = new PointSet;
	std::vector<vec3>&point = pts->points();
	std::vector<vec3>&normal = pts->normals();
	std::vector<vec3>&color = pts->colors();
	std::vector<double>&images = pts->images();

	for (int i(0); i < points.size(); i++) {
		point.push_back(to_my_point(get<0>(points[i])));
		normal.push_back(to_my_vector(get<1>(points[i])));
		vec3 co = vec3((int)get<2>(points[i]).at(0), (int)get<2>(points[i]).at(1), (int)get<2>(points[i]).at(2));
		color.push_back(co);
		images.push_back(get<3>(points[i]));
	}

	std::cout << "pts information: \n"
		<< "point size: " << point.size()
		<< "\n normal size: " << normal.size()
		<< "\n color size: " << color.size()
		<< "\n image size: " << images.size() << std::endl;

	return pts;
}




int main(int argc, char **argv) {
  
	int stop_iterations = 3;                    //停止迭代的次数
	double pd_normal_deviation = 0.9;			//法线阈值
	int pd_nn = 10;                             //K邻域值
	int pd_sigma(200);

	double pd_c = 0.6;
	double pd_s = 0.6;

	bool pd_out_vg = false;
	bool pd_out_alpha_shape = false;
	bool pd_out_convex_hull = false;
	bool save_refine = false;
	const std::string path_point_cloud = "E:/MvsDataBuild/database/Cloud32.ply";
	std::cout << "processing " << path_point_cloud << " ..." << std::endl;
	ShapeDetector* CS;

	PointSet *pset = dataread(path_point_cloud);
	CS = new ShapeDetector(pset);

	int num_0 = CS->get_points_number();

	double pd_epsilon_0 = 0.005 * CS->get_bbox_diagonal();
	double pd_epsilon = (round(1000 * pd_epsilon_0)) / 1000.0;

	std::cout << "shape_min_number: " << pd_sigma << std::endl;
	std::cout << "epsilon: " << pd_epsilon << std::endl;

	CS->set_detection_parameters(pd_sigma, pd_nn, pd_normal_deviation);
	CS->set_max_steps(stop_iterations);
	CS->set_lambda_r(pd_s);
	CS->set_lambda_c(pd_c);
	CS->set_epsilon(pd_epsilon);
	CS->refine_plane(save_refine);

	CS->detect_shapes();
	CS->set_primitives_simple();
	CS->planar_shape_detection_hybrid();
	CS->set_primitives_simple();
	//CS->save_alpha_shapes();

	PointSet *pset1 = new PointSet;
	std::vector<VertexGroup::Ptr> roofs;
	CS->Get_VertexGroup(pset1, roofs);

    return EXIT_FAILURE;
}