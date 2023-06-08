#pragma once

#include "../model/vertex_group.h"
#include "../model/map.h"
#include "../model/point_set.h"
#include "../math/linear_program_solver.h"


#include <vector>

class PolyFitInfo;

class Reconstruction
{
 public:
	Reconstruction() {}
	~Reconstruction() {}

	void Initialize(PointSet* pset_, Map* foot_print_, std::vector<VertexGroup::Ptr> Group_);
	void Extract_Building_planes();
	bool Reconstruct(Map *result, LinearProgramSolver::SolverName solver_name,bool update_display=false);


 private:

	void Get_PlaneGroup( VertexGroup *building);

	PointSet * create_roof_point_set(const std::vector<VertexGroup::Ptr> &segments, VertexGroup *building);
	PointSet * create_projected_point_set( const PointSet *roof);


	std::vector<vec3> compute_line_segment(PointSet *seg_pset, PointSet *roof_pset, Map::Facet *footprint);
	std::vector<std::vector<int>> detect_height_jump(PointSet *pset, Map::Facet *foot_print, double min_height, double max_height);
	std::vector<std::vector<int>>compute_height_field(PointSet *pset, Map::Facet *footprint);

	Map *reconstruct_single_building(PointSet *roof_pset,
		std::vector<vec3> line_segments,
		Map::Facet *footprint,
		LinearProgramSolver::SolverName solver_name);

	void extrude_boundary_to_ground(Map *model, const Plane3d &ground, PolyFitInfo *polyfit_info);

	Map *generate_polygon(PointSet *pSet, double footprint_height, double density);
private:

	bool is_simple_polygon(const Map::Facet *f);



private:
	PointSet* pset;
	Map* foot_print;
	std::vector<VertexGroup::Ptr> PlaneGroup;

	double image_x_min;
	double image_y_min;
	double image_delta_res;

	int width = 400, height = 400;

};

