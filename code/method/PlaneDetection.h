#ifndef _POINTS_REGION_GROWING_H_
#define _POINTS_REGION_GROWING_H_

#include"../method/cgal_types.h"
#include"../method/alpha_shape.h"
#include"../model/point_set.h"
#include"../model/vertex_group.h"
#include"../model/map_attributes.h"
#include"../model/map.h"
#include"../model/map_circulators.h"
#include"../model/kdtree_search.h"
#include"../basic/file_utils.h"

#include<vector>
#include<string>
#include<random>
class PointSet;
class VertexGroup;

class ShapeDetector {

	//操作符重载
	struct Add_Comparator {
		bool operator()(std::pair<int, double> p1, std::pair<int, double> p2) {
			return p1.second > p2.second;
		}
	};
	struct Add_Comparator_normal {
		bool operator()(std::pair<int, double> p1, std::pair<int, double> p2) {
			return p1.second < p2.second;
		}
	};

	struct Remove_Comparator {
		bool operator()(std::pair<int, double> p1, std::pair<int, double> p2) {
			return p1.second < p2.second;
		}
	};
	struct Remove_Comparator_normal {
		bool operator()(std::pair<int, double> p1, std::pair<int, double> p2) {
			return p1.second > p2.second;
		}
	};

	struct Weight_Comparator_with_energy {
		bool operator()(std::pair<std::pair<int, std::vector<int>>, std::vector<double>> p1, std::pair<std::pair<int, std::vector<int>>, std::vector<double>> p2) {
			return p1.second[0] > p2.second[0];
		}
	};


public:
	//外部接口函数

	ShapeDetector(const std::string &filepath);
	ShapeDetector(PointSet*pset);
	~ShapeDetector();


	void set_detection_parameters(int _rg_min_points, int _knn, double _rg_normal_threshold);
	void detect_shapes();
	void planar_shape_detection_hybrid();
	void set_primitives_simple();

	bool read_ply();
	void load_PointSet(PointSet* pset);

	void Get_VertexGroup(PointSet*pset, std::vector<VertexGroup::Ptr> &results);

public:
	//参数设置
	void set_lambda_r(double db) { lambda_r = db; }
	void set_lambda_c(double db) { lambda_c = db; }
	void set_weight_m(int wm) { weight_mode = wm; }
	void set_epsilon(double _rg_epsilon) { epsilon = _rg_epsilon; }
	void set_max_steps(int si) { stop_iterations = si; }
	void refine_plane(bool re) { refine_planes = re; }
	double get_bbox_diagonal() { return bbox_diagonal; }
	void show_result(double t);
	double get_points_number() { return points.size(); }
	//保存结果
	void save_convex_hull();
	void save_alpha_shapes();
	void save_vg();
	void save_curvature_conf();
	void save_horsduff_ditance();


protected:
	//中间函数
	void set_extrema();
	void detect_planes();
	void compute_average_spacing_and_k_nearest_neighbors();
	void do_region_growing();
	void get_coverage_and_mean_error();
	std::vector<Plane_3>  refine_coplane();//共面检测及合并




private:
	//曲率计算
	void compute_curvature();
	void seed_sort();

	void merge(const std::vector<int> a, const std::vector<int >b);
	//基础函数	


	bool inliers_arent_aligned(const std::vector<int> & inds);
	void get_last_information();
	void test_connected_primitives();
	void get_distance_diviation();
	void get_good_points();
	void get_bad_points();
	double add_changed_error(int id, std::vector<int> point_ids);
	double energy_changed_second(double dis, double numb);
	double remove_changed_error(int id, std::vector<int> point_ids);
	double merge_distance_changed_with_epsilon(int i, int j, double & dif, std::vector<int> & move_ids);
	double energy_changed_merge(double dis, double numb, int nmoves);
	bool update_good_points(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	bool update_bad_points(int id_shape, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	double energy_changed(double dis, double numb);
	bool test_if_connected(int i, int j);
	bool separate_two_out(int id, std::vector<int> & max_list, std::vector<int> & min_list, double & dif);


	void local_operators();
	bool convert_form_merge_l2(int i, int j, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	bool convert_form_split_l2(int i, std::vector<std::vector<int>>& max_list_vector, std::vector<std::vector<int>>& min_list_vector, std::pair<std::pair<int, std::vector<int>>, std::vector<double>> & one_element);
	bool convert_form_exclude_l2(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);
	bool convert_form_insert_l2(int i, std::pair<std::pair<int, std::vector<int>>, std::vector<double>>& one_element);

	void get_distance_diviation_show_merge_info(double t);
	void transfer_operator_normal_for_hybrid();
	void get_distance_diviation_show_normal_info(double t);
	void get_coverage_and_mean_error_pure();





	//初始化参数
private:
	std::vector<Point_NCI> points;

private:
	int number_iterations;//迭代次数
	int min_points;//最少点数
	int knn;
	double lambda_c;//完整项系数
	double lambda_r;//
	int stop_iterations;//停止迭代次数
	double normal_threshold;
	double epsilon;
	int weight_mode;
	int non_coplanar_planes;//非共面平面数量
	int number_of_insert_exclude; //设置的外点数量


	int t_l;
	double mean_error;
	double all_error;
	double ori_all_error;
	int t_m;
	int t_merge;
	int t_split;
	int t_insert;
	int t_exlude;
	int all_t_transfer;
	int all_t_merge;
	int all_t_split;
	int all_t_insert;
	int all_t_exlude;

	double mean_normal_diviation;
	double mean_distance_diaviation;
	double interval_all;
	double all_distance_diaviation;
	double all_normal_diaviation;
	std::size_t number_of_assigned_points;
	double size_current_primitives;
	double mean_distance_current;
	double mean_normal_current;
	double ori_mean_normal_diviation;
	double number_inlier_before_opers;
	double ori_coverage;
	double ori_mean_error;
	int ori_inliers_number;
	double mean_distance_diaviation_before_transfer;
	double mean_normal_diviation_befor_transfer;


	int ori_primitives_number;
	int last_primitives_number;
	double last_normal_deviation;
	double last_coverage;
	double last_mean_error;

	std::vector<int> region_type;
	std::vector<int> region_type_m;
	std::vector<int> points_if_added;
	std::vector<int> points_changed;
	std::vector<int> bad_points;
	std::vector<int> if_insert_candidate;
	std::vector<std::pair<int, int>> good_points;
	std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> good_points_shape;
	std::vector<std::pair<std::pair<int, std::vector<int>>, std::vector<double>>> bad_points_shape;



	int old_size_current_primitives;
	double old_mean_distance_diaviation;
	double old_mean_normal_diaviation;
	double old_coverage;

protected:

	std::string path_point_cloud;
	std::string base_path_point_cloud;
	std::string path_point_cloud_basename;

	bool if_constraint;
	bool spacing_is_known;
	bool should_compute_knn;
	bool should_compute_neighborhood;
	bool refine_planes;
private:
	double x_min, x_max, y_min, y_max, z_min, z_max, img_max, img_min;
	double average_spacing;
	double bbox_diagonal;
	double coverage;
	int primitives_number;

	std::vector<Plane_3> planes_0;
	std::vector<Plane_3> planes_1;
	std::vector<Plane_3> planes_2;
	std::vector<Plane_3> planes_3;
	std::vector<Plane_3> planes_centroids;
	std::vector<CGAL::Color> planes_to_colors;

	std::vector<std::vector<int> > spherical_neighborhood_former;
	std::vector<std::vector<int> > spherical_neighborhood;//环形邻域
	std::vector<double>spherical_curvature; //种子点可信度值

	std::vector<std::vector<int> > planes_to_inliers;//外层为平面数，内层为平面内的点索引
	std::vector<int> inliers_to_planes; //表示在有效平面内的点
	std::vector<std::vector<int>> primitive_connection;

	std::vector<std::vector<Point_3> > alpha_shapes_pts;
	std::vector<CGAL::Color> alpha_shapes_colors;

	std::vector<std::vector<Point_3> > convex_hulls_pts;
	std::vector<CGAL::Color> convex_hulls_colors;

	std::vector<std::pair<int, int> > convex_hulls_seeds;


};
#endif