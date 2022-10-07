#ifndef EXTRACT_FROM_ODB_H
#define EXTRACT_FROM_ODB_H

#include <stdio.h>
#include<vector>
#include "Input_Reader.h"

//Include for Abaqus
#include <odb_API.h>

const double Zero = 1e-10;

//Data structure for a cuboid
struct cuboid
{
	//Cuboid's minimum corrdinates
	double min_x, min_y, min_z;
    //Cuboid's length, width and height
    double lx, ly, lz;
    //Cuboid's maximum corrdinates
    double max_x, max_y, max_z;
};

class Extract_From_ODB
{
public:

	//Data members

	//Constructor
	Extract_From_ODB() {};

	//Function members
	int Extract_data_from_odb(const Input& Init);
	int Is_step_in_odb(const odb_Odb& odb, const string& step_name);
	int Get_number_of_cnts_and_points(int& N, vector<int>& n_points);
	int Extract_cnt_data_for_frame(const int& frame, const int& n_cnts, const vector<int>& n_points, odb_Assembly& root_assy, odb_FieldOutput& current_fieldU);
	int Extract_incremental_cnt_data_for_frame(const int& frame, const int& n_cnts, const vector<int>& n_points, odb_Assembly& root_assy, odb_FieldOutput& previous_fieldU, odb_FieldOutput& current_fieldU);
	string Get_cnt_set_name(const int& cnt_i);
	int Get_number_of_gnps_and_vertices_inside(int& N, vector<vector<int> >& vertices_in);
	int Get_sample_cuboid(cuboid& sample);
	int Calculate_vertex_coordinates(const double vertex[], const double& Cx, const double& Cy, const double& Cz, const double& cosy, const double& siny, const double& cosz, const double& sinz, double new_vertex[]);
	bool Is_vertex_inside_cuboid(const double vertex[], const cuboid& sample);
	int Extract_gnp_data_for_frame(const int& frame, const int& n_gnps, const vector<vector<int> >& vertices_in, odb_Assembly& root_assy, odb_FieldOutput& current_fieldU);
	string Get_gnp_set_name(const int& gnp_i, const int& vertex);
	vector<int> Vector_for_all_vertices_inside();
	int Extract_matrix_data_for_frame(const int& frame, odb_Assembly& root_assy, odb_FieldOutput& current_fieldU);
	int Get_displacement_for_node_set(const string& set_name, odb_Assembly& root_assy, odb_FieldOutput& current_fieldU, double disp[]);
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
