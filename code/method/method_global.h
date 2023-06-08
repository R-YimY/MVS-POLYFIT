#ifndef _METHOD_GLOBAL_H_
#define _METHOD_GLOBAL_H_



#include "../math/linear_program_solver.h"

#include <string>


namespace Method {

	extern double lambda_data_fitting;
	extern double lambda_model_image;
	extern double lambda_model_complexity;

    extern double number_region_growing;
    extern double image_resolution;
	extern double ground_height;
	extern double coincident_threshold;

	//________________ names for various quality measures ____________________

	extern std::string facet_attrib_supporting_vertex_group;
	extern std::string facet_attrib_supporting_point_num;
	extern std::string facet_attrib_facet_area;
	extern std::string facet_attrib_covered_area;


	extern double dist_thresh;
	extern double bitmap_reso;
	extern double normal_thresh;


	extern double Image_max;
	extern double Image_min;


}


#endif