#include"../method/Image_Guided_filtering.h"
#include<map>


void ImageGuidedFilter::DataReady(PointSet *pset, double _epsilon, int KNN) {
	if (pset->points().size() == 0) {
		std::cout << "Unaviliable data for GuidedFiltering" << std::endl;
	}
	else {
		ori_points = pset->points();
		ori_normals = pset->normals();
		ori_images = pset->images();
		ori_colors = pset->colors();
		epsilon = _epsilon;
		RadiuK = KNN;

	}


	if (ori_points.size() != 0) {
		ImgMIN = FLT_MAX; ImgMAX = -FLT_MAX;
		for (int i(0); i < ori_points.size(); i++) {
			Point_3 pt(ori_points[i][0], ori_points[i][1], ori_points[i][2]);
			Vector_3 nr(ori_normals[i][0], ori_normals[i][1], ori_normals[i][2]);
			points.push_back(std::make_pair(pt, nr));
			if (ori_images[i] > ImgMAX) ImgMAX = ori_images[i];
			if (ori_images[i] < ImgMIN) ImgMIN = ori_images[i];
		}
	}

	//std::cout << "Total Points Number is: " << points.size() << std::endl;
}


//数据滤波处理
PointSet* ImageGuidedFilter::ImageFilter(PointSet * pset, double epsilon, int KNN)
{
	//step 1: initialization
	DataReady(pset,epsilon,KNN); //数据初始化
	PointSet *pointcloud = new PointSet;
	std::vector<vec3>&after_points = pointcloud->points();
	after_points.resize(points.size());

	std::map<Point_3, int>map_indice_point;
	std::list<Point_3>list_points;

	for (int i(0); i < points.size(); i++) {
		const Point_3 & pt = points[i].first;
		map_indice_point[pt] = i;
		list_points.push_back(pt);
	}
	KD_Tree tree(list_points.begin(), list_points.end()); //建立检索树
	std::cout << "Computing K nearest neighbors ";
	int ten_percent_of_points = points.size() / 10;


	std::cout << " Filtering: " ;

	//step 2 : Filtering every point
	for (int i(0); i < points.size(); i++) 
	{
		if (i % ten_percent_of_points == 0) std::cout << ". " ;//显示处理进程

		Point_3 query = points[i].first;
		Neighbor_search search(tree, query, RadiuK + 1);

		std::vector<int> index_of_neighbors;
		index_of_neighbors.resize((int)RadiuK);
		std::vector<double> distance_neighbors;			
		distance_neighbors.resize((int)RadiuK);

		for (Neighbor_search::iterator it = search.begin(); it != search.end(); it++) {
			std::map<Point_3, int>::iterator iter = map_indice_point.begin();
			iter = map_indice_point.find(it->first);
			if (iter != map_indice_point.end() && iter->second != i) {
				index_of_neighbors.push_back(iter->second);
				double d = sqrt((query - iter->first).squared_length());
				distance_neighbors.push_back(d);
			}
		}

		//邻域点的矩阵计算
		Eigen::MatrixXd neighbors_as_matrix(3, index_of_neighbors.size());
		for (std::size_t j(0); j < index_of_neighbors.size(); j++) {
			neighbors_as_matrix(0, j) = points[index_of_neighbors[j]].first.x();
			neighbors_as_matrix(1, j) = points[index_of_neighbors[j]].first.y();
			neighbors_as_matrix(2, j) = points[index_of_neighbors[j]].first.z();
		}

		Eigen::Vector3d mean;
		mean = neighbors_as_matrix.rowwise().mean();//计算质心
		neighbors_as_matrix.transposeInPlace();//矩阵转置
		Eigen::MatrixXd centered = neighbors_as_matrix.rowwise() - neighbors_as_matrix.colwise().mean();
		Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(neighbors_as_matrix.rows() - 1);


		//方差计算
		//当前点
		double cur_dis_var = CurrentDisVar(distance_neighbors);
		double cur_img_var = CurrentImgVar(index_of_neighbors);
		//周围点
		
		std::vector<double> all_distance;
		std::vector<int>all_image;
		for (std::size_t k(0); k < index_of_neighbors.size(); k++) {
			Point_3 neighbor_query = points[index_of_neighbors[k]].first;
			Neighbor_search search2(tree, neighbor_query, RadiuK + 1);
			for (Neighbor_search::iterator nei_it = search2.begin(); nei_it != search2.end(); nei_it++) 
			{
				std::map<Point_3, int>::iterator nei_iter = map_indice_point.begin();
				nei_iter = map_indice_point.find(nei_it->first);
				if (nei_iter != map_indice_point.end() && nei_iter->second != k) {
					double nei_d = sqrt((neighbor_query - nei_iter->first).squared_length());
					int nei_img = ori_images[index_of_neighbors[k]] - ori_images[nei_iter->second];
					all_distance.push_back(nei_d);
					all_image.push_back(nei_img);
				}
			}
		}
		
		double nei_dis_var = NeighborDisVar(all_distance);
		double nei_img_var=NeighborImgVar(all_image);
		

		//权重计算
		double weight = (cur_dis_var + cur_img_var) / (nei_dis_var + nei_img_var);
		weight = weight / (double(index_of_neighbors.size()));

		Eigen::MatrixXd e = (cov + (epsilon / weight) * Eigen::MatrixXd::Identity(3, 3));
		e = e.inverse();//求逆 
		Eigen::MatrixXd A = cov * e; //系数a
		Eigen::MatrixXd b = mean - A * mean; //系数b

		Eigen::Vector3d SearchPointType;
		SearchPointType[0] = query.x();
		SearchPointType[1] = query.y();
		SearchPointType[2] = query.z();
		  
		SearchPointType = A * SearchPointType + b;

		vec3 pt(SearchPointType[0], SearchPointType[1], SearchPointType[2]);//最后平滑的点
		after_points[i]=(pt);

	}

	pointcloud->normals() = ori_normals;
	pointcloud->images() = ori_images;
	pointcloud->colors() = ori_colors;
	std::cout << " Done: " << pointcloud->points().size() << std::endl;

	return pointcloud;

}

double ImageGuidedFilter::CurrentDisVar(std::vector<double> distance_neighbors)
{
	double avg(0);
	for (int i(0); i < distance_neighbors.size(); i++) {
		avg += distance_neighbors[i];
	}
	avg = avg / distance_neighbors.size();

	double Variance(0);
	for (int i(0); i < distance_neighbors.size(); i++) {
		Variance += std::pow((distance_neighbors[i] - avg), 2);
	}
	Variance = Variance / distance_neighbors.size();
	return Variance;
}

double ImageGuidedFilter::CurrentImgVar(std::vector<int> PointIndex) {
	
	double avg(0);
	for (std::size_t i(0); i < PointIndex.size(); i++) {
		avg += ori_images[PointIndex[i]];
	}
	avg = avg / double(PointIndex.size());

	double Variance(0);
	for (std::size_t i(0); i < PointIndex.size(); i++) {
		Variance += std::pow((ori_images[PointIndex[i]] - avg), 2);
	}
	Variance = Variance / double(PointIndex.size());
	return  Variance;
}


double ImageGuidedFilter::NeighborDisVar(std::vector<double> distance) {

	double avg(0);
	for (int i(0); i < distance.size(); i++) {
		avg += distance[i];
	}
	avg = avg / (double)distance.size();

	double variance(0);
	for (int i(0); i < distance.size(); i++) {
		variance += std::pow((avg - distance[i]), 2);
	}
	variance = variance / (double)distance.size();
	return variance;
}

double  ImageGuidedFilter::NeighborImgVar(std::vector<int> imgs) {

	double avg(0);
	for (int i(0); i < imgs.size(); i++) {
		avg += imgs[i];
	}
	avg = avg / (double)imgs.size();

	double variance(0);
	for (int i(0); i < imgs.size(); i++) {
		variance += std::pow((avg - imgs[i]), 2);
	}
	variance = variance / (double)imgs.size();
	return variance;

}