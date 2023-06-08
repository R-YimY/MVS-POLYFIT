/*
初始数据最好进行数据尺度转换，点间距不能过小
*/

#ifndef _IMAGE_GUIDED_FILTERING_H_
#define _IMAGE_GUIDED_FILTERING_H_

#include"../model/point_set.h"
#include"../model/vertex_group.h"
#include"../method/cgal_types.h"
#include<Eigen/Dense>
#include<Eigen/Core>

class PointSet;
class ImageGuidedFilter {
public:
	PointSet* ImageFilter(PointSet*pset, double epsilon, int KNN);
private:
	void DataReady(PointSet *pset,double _epsilon,int KNN); //数据格式转换
	double CurrentDisVar(std::vector<double> distance_neighbors); //计算当前点的距离方差
	double CurrentImgVar(std::vector<int> PointIndex);//计算当前点的影像方差

	double NeighborDisVar(std::vector<double> distance);//计算邻域点的距离方差
	double NeighborImgVar(std::vector<int> imgs);//计算邻域点的影像方差
private:
	int RadiuK;
	double epsilon;
	Pwn_vector points;//原始点云数据
	std::vector<vec3>ori_points;
	std::vector<vec3>ori_normals;
	std::vector<double>ori_images;
	std::vector<vec3>ori_colors;
	int ImgMAX;
	int ImgMIN;
};

#endif


