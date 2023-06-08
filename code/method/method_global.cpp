#include "method_global.h"


namespace Method
{

	double lambda_data_fitting = 0.35;
	double lambda_model_image = 0.15;
	double lambda_model_complexity = 0.5;

    double number_region_growing=100;
    double image_resolution =0.05;
	double ground_height = -10.0;
	double coincident_threshold = 1e-7;

	//________________ names for various quality measures ____________________

	std::string facet_attrib_supporting_vertex_group = "facet_supporting_vertex_group";
	std::string facet_attrib_supporting_point_num = "facet_supporting_point_num";
	std::string facet_attrib_facet_area = "facet_area";
	std::string facet_attrib_covered_area = "facet_covered_area";

	double dist_thresh = 0.005;
	double bitmap_reso = 0.05;
	double normal_thresh = 0.95;

	double Image_max = 15.0;
	double Image_min = 3.0;



}
