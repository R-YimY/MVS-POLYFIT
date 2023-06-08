/*
��ʼ������ý������ݳ߶�ת�������಻�ܹ�С
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
	void DataReady(PointSet *pset,double _epsilon,int KNN); //���ݸ�ʽת��
	double CurrentDisVar(std::vector<double> distance_neighbors); //���㵱ǰ��ľ��뷽��
	double CurrentImgVar(std::vector<int> PointIndex);//���㵱ǰ���Ӱ�񷽲�

	double NeighborDisVar(std::vector<double> distance);//���������ľ��뷽��
	double NeighborImgVar(std::vector<int> imgs);//����������Ӱ�񷽲�
private:
	int RadiuK;
	double epsilon;
	Pwn_vector points;//ԭʼ��������
	std::vector<vec3>ori_points;
	std::vector<vec3>ori_normals;
	std::vector<double>ori_images;
	std::vector<vec3>ori_colors;
	int ImgMAX;
	int ImgMIN;
};

#endif


