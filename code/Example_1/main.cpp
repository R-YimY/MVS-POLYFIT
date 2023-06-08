/*
test 1:
	引导滤波算法实验


*/

#include "../model/point_set.h"
#include "../model/map.h"
#include "../model/map_io.h"
#include "../model/point_set_io.h"
#include "../method/reconstruction.h"
#include"../method/Image_Guided_filtering.h"
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
		std::make_pair(ImageMap(), CGAL::PLY_property<double>("scalar_Image"))
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

void data_output(PointSet*pts, std::string filepath) {

	std::cout << "saving filtered pointcloud......" << std::endl;
	std::ofstream stm(filepath, std::ios::out);

	std::vector<vec3>&points = pts->points();
	std::vector<vec3>&normals = pts->normals();
	std::vector<vec3>&colors = pts->colors();
	std::vector<double>&images = pts->images();

	int num = pts->points().size();
	stm << "ply" << std::endl;
	stm << "format ascii 1.0" << std::endl;
	stm << "element vertex " << num << std::endl;
	stm << "property float x" << std::endl;
	stm << "property float y" << std::endl;
	stm << "property float z" << std::endl;
	stm << "property uchar red" << std::endl;
	stm << "property uchar green" << std::endl;
	stm << "property uchar blue" << std::endl;
	stm << "property float nx" << std::endl;
	stm << "property float ny" << std::endl;
	stm << "property float nz" << std::endl;
	stm << "property float scalar_Image" << std::endl;
	stm << "end_header" << std::endl;

	for (int i(0); i < num; i++) {
		stm << points[i].x << " " << points[i].y << " " << points[i].z << " " <<
			colors[i].x << " " << colors[i].y << " " << colors[i].z << " " <<
			normals[i].x << " " << normals[i].y << " " << normals[i].z << " " <<
			images[i] << std::endl;
	}


}




int main(int argc, char **argv) {
	std::string file = "E:/MvsDataBuild/model/model_5-2.ply";
	const std::string output = "G:/City3D-main/data/MVS/filtered.ply";

	PointSet *pset = dataread(file);


	Box3d box = pset->bbox();
	double epsilon = box.radius() * 0.001f;
	int knn = 6;

	std::cout << "epsilon:" << epsilon << std::endl;
	ImageGuidedFilter filter;
	PointSet* pointfilter = filter.ImageFilter(pset, epsilon, knn);

	data_output(pointfilter, output);

}