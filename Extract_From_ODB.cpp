
#include "Extract_From_ODB.h"

//Main function to extract data from ODB
int Extract_From_ODB::Extract_data_from_odb(const Input& Init)
{
    time_t itt0, itt1;
    itt0 = time(NULL);
    //Initialize Abaqus C++ API
    odb_initializeAPI();

    //Open Abaqus database using the name/path intidated in the input file
    //Use a C string, since that seems to be equivalent to (or able to be cast as) an odb_String
    odb_Odb& odb = openOdb(Init.extract_para.odb_file.c_str());
    //cout << "open odb=" << Init.extract_para.odb_file.c_str() << endl;

    //Access the root assembly
    odb_Assembly& root_assy = odb.rootAssembly();
    //cout << "root_assy" << endl;

    //Make sure the step indicated in the input file is in the odb file
    if (!Is_step_in_odb(odb, Init.extract_para.step_name))
    {
        cout << "Error in Nanoparticle_resistor_network_from_odb when calling Is_step_in_odb." << endl;
        return 0;
    }

    //Get all frames from the step and save them in a (pointer) variable
    odb_SequenceFrame& allFramesInStep = odb.steps()[Init.extract_para.step_name.c_str()].frames();
    //Get the number of frames in the database
    int n_frames = allFramesInStep.size();
    //Check the first frame indicated in the input parameters is valid
    if (Init.extract_para.first_frame >= n_frames)
    {
        cout << "Error in Extract_data_from_odb: First frame number input is greater than the maximum frame number." << endl;
        cout << "First frame number input is " << Init.extract_para.first_frame << ", maximum frame number is " << n_frames - 1 << " (there are " << n_frames << " frames)." << endl;
        return 0;
    }
    cout << endl << "There are " << n_frames << " frames in the Abaqus database." << endl;
    itt1 = time(NULL);
    cout << "Initialization of Abaqus API and variables time: " << (int)(itt1 - itt0) << " secs." << endl;

    //Data needed in case there are CNTs
    int n_cnts = 0;
    vector<int> n_points;

    //Booleans to identify the fillers in the sample
    bool cnts_present = Init.extract_para.filler_type == "CNTs_only" || Init.extract_para.filler_type == "GNP_CNT_mix";
    bool gnps_present = Init.extract_para.filler_type == "GNPs_only" || Init.extract_para.filler_type == "GNP_CNT_mix";

    //Get the CNT data if there are CNTs
    if (cnts_present)
    {
        time_t it0, it1;
        it0 = time(NULL);

        //Get CNT data
        cout << endl << "Reading CNT data from csv files ..." << endl;
        if (!Get_number_of_cnts_and_points(n_cnts, n_points))
        {
            cout << "Error in Extract_data_from_odb when calling Get_number_of_cnts_and_points." << endl;
            return 0;
        }
        it1 = time(NULL);
        cout << "Read CNT data from csv files in " << (int)(it1 - it0) << " secs." << endl << endl;
    }

    //Data needed in case there are GNPs
    int n_gnps = 0;
    vector<vector<int> > vertices_in;

    //Get the GNP data if there are GNPs
    if (gnps_present)
    {
        time_t it0, it1;
        it0 = time(NULL);

        if (!cnts_present)
            cout << endl;

        //Get GNP data
        cout << "Reading GNP data from csv files ..." << endl;
        if (!Get_number_of_gnps_and_vertices_inside(n_gnps, vertices_in))
        {
            cout << "Error in Extract_data_from_odb when calling Get_number_of_gnps_and_vertices_inside." << endl;
            return 0;
        }
        it1 = time(NULL);
        cout << "Read GNP data from csv files in " << (int)(it1 - it0) << " secs." << endl << endl;
    }

    //Iterate overt the frames, ignoring the first one which has no deformation
    for (int i = Init.extract_para.first_frame; i < n_frames; i++)
    {
        cout << "============================================================================" << endl;
        cout << "============================================================================" << endl;
        cout << "Frame " << i << endl;
        time_t it0, it1;
        it0 = time(NULL);

        //Access displacement field ("U") in the current frame
        //hout << "current_fieldU" << endl;
        odb_FieldOutput& current_fieldU = allFramesInStep[i].fieldOutputs()["U"];
        //Access displacement field ("U") in the previous frame
        //hout << "previous_fieldU" << endl;
        //odb_FieldOutput& previous_fieldU = allFramesInStep[i-1].fieldOutputs()["U"];

        //Extract the CNT data if needed
        if (cnts_present)
        {
            time_t cnt0, cnt1;
            cnt0 = time(NULL);
            cout << "Extracting CNT data..." << endl;
            if (!Extract_cnt_data_for_frame(i, n_cnts, n_points, root_assy, current_fieldU))
            {
                cout << "Error in Extract_data_from_odb when calling Extract_cnt_data_for_frame." << endl;
                return 0;
            }
            cnt1 = time(NULL);
            cout << "Extracted CNT data in " << (int)(cnt1 - cnt0) << " secs." << endl;
        }

        //Extract the GNP data if needed
        if (gnps_present)
        {
            time_t gnp0, gnp1;
            gnp0 = time(NULL);
            cout << "Extracting GNP data..." << endl;
            if (!Extract_gnp_data_for_frame(i, n_gnps, vertices_in, root_assy, current_fieldU))
            {
                cout << "Error in Extract_data_from_odb when calling Extract_gnp_data_for_frame." << endl;
                return 0;
            }
            gnp1 = time(NULL);
            cout << "Extracted GNP data in " << (int)(gnp1 - gnp0) << " secs." << endl;
        }

        //Extract the matrix data
        if (!Extract_matrix_data_for_frame(i, root_assy, current_fieldU))
        {
            cout << "Error in Extract_data_from_odb when calling Extract_matrix_data_for_frame." << endl;
            return 0;
        }

        it1 = time(NULL);
        cout << "Frame " << i << " time: " << (int)(it1 - it0) << " secs." << endl;
    }

    //Close Abaqus database
    odb.close();

    //Finalize usage of Abaqus C++ API
    odb_finalizeAPI();

	return 1;
}
//This function determines if the step indicated in the input file is in the odb file
int Extract_From_ODB::Is_step_in_odb(const odb_Odb& odb, const string& step_name)
{
    //Check if there is at least one step
    if (odb.steps().size() < 1)
    {
        cout << "Error in Is_step_in_odb: There are no steps in the odb file." << endl;
        return 0;
    }

    //Check that the step exists in the odb file
    odb_StepRepositoryIT stepIter(odb.steps());
    //cout << "from input file:" << step_name << endl;
    for (stepIter.first(); !stepIter.isDone(); stepIter.next())
    {
        if (step_name == stepIter.currentKey().CStr())
            return 1;
        //cout << stepIter.currentKey().CStr() << endl;
    }
    //cout << endl;

    cout << "Error in Is_step_in_odb: Step name indicated in input file was not found in odb file." << endl;
    cout << "Step name in input file is: " << step_name << "." << endl;
    cout << "This is a list of all steps in odb file:" << endl;
    for (stepIter.first(); !stepIter.isDone(); stepIter.next())
    {
        cout << stepIter.currentKey().CStr() << endl;
    }

    return 0;
}
//This functiongets the number of CNTs and the number of points per CNT from the 
//structure csv file
int Extract_From_ODB::Get_number_of_cnts_and_points(int& N, vector<int>& n_points)
{
    //Open the structure file
    ifstream struc_file;
    struc_file.open("cnt_struct.csv");
    if (!struc_file) {
        cout << "Error in Get_number_of_cnts_and_points: Failed to open CNT structure file cnt_struct.csv." << endl;
        return 0;
    }
    //String to read lines from file
    string line;

    //Integers to store the number of CNTs
    int n_cnts, n_cnts_check;

    //Read the first line line form the file and store it in a string stream
    getline(struc_file, line);
    stringstream ss_n_cnts(line);

    //Read the number of CNTs from the string stream and ignore the comma
    ss_n_cnts >> n_cnts; ss_n_cnts.ignore();
    ss_n_cnts >> n_cnts_check;

    //Check that both numbers in the first line of the structure file are the same
    if (n_cnts != n_cnts_check) {
        cout << "Error in Get_number_of_cnts_and_points: The first line of file cnt_struct.csv indicates different number of CNTs: " << n_cnts << " and " << n_cnts_check << "." << endl;
        return 0;
    }

    //Set the number of CNTs in the simulation
    N = n_cnts;

    //Set the size of the vector with the number of points
    n_points.assign(n_cnts, -1);

    //Iterate over the number of CNTs
    for (int i = 0; i < n_cnts; i++)
    {
        //Check if end-of-file has been reached
        if (struc_file.eof())
        {
            cout << "Error in Get_number_of_cnts_and_points. The end-of-file of cnt_struct.csv has been reached before reading all CNT data." << endl;
            return 0;
        }

        //Read a line form the file and store it in a string stream
        getline(struc_file, line);
        stringstream ss(line);

        //Read the number of CNT points in CNT i from the string stream and
        //save it in the vector of number of points
        ss >> n_points[i];
    }

    //Close structure file
    struc_file.close();

    return 1;
}
//This function extracts the CNT data for a given frame and saves it in a file
int Extract_From_ODB::Extract_cnt_data_for_frame(const int& frame, const int& n_cnts, const vector<int>& n_points, odb_Assembly& root_assy, odb_FieldOutput& current_fieldU)
{
    //Filename for current frame
    string filename_points = "CNTs_disp_F" + to_string(frame) + ".dat";

    //Open binary file for incremental displacements of CNT nodes in current frame
    ofstream otec_points(filename_points.c_str(), ios::binary | ios::out);

    //Get the size of a double
    streamsize double_size = sizeof(double);

    //Variable to store the total number of CNT points in the network
    long int Np = 0;

    //Iterate over the number of CNTs
    for (int i = 0; i < n_cnts; i++)
    {
        //Variable to store the number of points in CNT i
        int np = n_points[i];

        //Get the set name for CNT i
        string set_name = Get_cnt_set_name(i + 1);

        //Variable to store displacement object
        odb_FieldOutput cnt_disp;

        //Sometimes abaqus does not generate the set, so catch that error if this happens
        try {
            //Access set from root assembly
            odb_Set& cnt_set = root_assy.nodeSets()[set_name.c_str()];

            //Get the displacement objects of the set
            cnt_disp = current_fieldU.getSubset(cnt_set);
            //hout << "current_fieldU" << endl;
        }
        catch (...) {
            cout << "Error in Extract_cnt_data_for_frame." << endl;
            cout << "Error while accessing set " << set_name.c_str() << "." << endl;
            return 0;
        }
        //hout << "nodeSets" << endl;

        //Get the values of the displacement object
        const odb_SequenceFieldValue& vals = cnt_disp.values();
        //hout << "cnt_disp" << endl;
        //Output size of values
        //cout << "values.size=" << vals.size() <<" structure[CNT="<<i<<"].size()="<< structure[i].size() << endl;

        //Check there is the same number of points in the CNT and in the set
        if (np != vals.size())
        {
            cout << "Error in Apply_displacements_to_cnts. The number of points in CNT " << i << " (" << np << " points) is different from the number of points in set " << set_name << " (" << vals.size() << " points)." << endl;
            return 0;
        }

        //Iterate over the values, i.e., the points in the current CNT
        //Note that nodes are actually in reverse order given the way
        //CNTs are generated (and meshed) in Abaqus
        //This this loop goes in reverse order
        for (int j = vals.size() - 1; j >= 0; j--)
        {
            //Get values at node j
            const odb_FieldValue val = vals[j];
            //Output node label
            //cout << "  Node: " << val.nodeLabel() << endl;

            //Get the data of the displacements
            int numComp = 0; //This integer is needed to call data() in the line below
            const float* const data = val.data(numComp);

            //Displacements are in data, where:
            //data[0] corresponds to displacement in x
            //data[1] corresponds to displacement in y
            //data[2] corresponds to displacement in z
            // 
            //Note data[i] is the total displacement along diraction i for a given frame, 
            //not the displacement with respect to the previous frame
            double disp_x = (double)data[0];
            double disp_y = (double)data[1];
            double disp_z = (double)data[2];

            //Save displacements in the binary file
            otec_points.write((char*)&disp_x, double_size);
            otec_points.write((char*)&disp_y, double_size);
            otec_points.write((char*)&disp_z, double_size);
        }
    }

    //Close binary file
    otec_points.close();

    return 1;
}
//This function extracts the CNT data for a given frame and saves it in a file
int Extract_From_ODB::Extract_incremental_cnt_data_for_frame(const int& frame, const int& n_cnts, const vector<int>& n_points, odb_Assembly& root_assy, odb_FieldOutput& previous_fieldU, odb_FieldOutput& current_fieldU)
{
    //Filename for current frame
    string filename_points = "CNTs_inc_disp_F" + to_string(frame) + ".dat";

    //Open binary file for incremental displacements of CNT nodes in current frame
    ofstream otec_points(filename_points.c_str(), ios::binary | ios::out);

    //Get the size of a double
    streamsize double_size = sizeof(double);

    //Iterate over the number of CNTs
    for (int i = 0; i < n_cnts; i++)
    {
        //Variable to store the number of points in CNT i
        int np = n_points[i];

        //Get the set name for CNT i
        string set_name = Get_cnt_set_name(i + 1);

        //Variables to store displacement objects
        odb_FieldOutput cnt_disp, cnt_disp_prev;

        //Sometimes abaqus does not generate the set, so catch that error if this happens
        try {
            //Access set from root assembly
            odb_Set& cnt_set = root_assy.nodeSets()[set_name.c_str()];

            //Get the displacement objects of the set
            cnt_disp = current_fieldU.getSubset(cnt_set);
            //hout << "current_fieldU" << endl;
            cnt_disp_prev = previous_fieldU.getSubset(cnt_set);
            //hout << "previous_fieldU" << endl;
        }
        catch (...) {
            cout << "Error in Extract_cnt_data_for_frame." << endl;
            cout << "Error while accessing set " << set_name.c_str() << "." << endl;
            return 0;
        }
        //hout << "nodeSets" << endl;

        //Get the values of the displacement object
        const odb_SequenceFieldValue& vals = cnt_disp.values();
        //hout << "cnt_disp" << endl;
        const odb_SequenceFieldValue& vals_prev = cnt_disp_prev.values();
        //hout << "cnt_disp_prev" << endl;
        //Output size of values
        //cout << "values.size=" << vals.size() <<" structure[CNT="<<i<<"].size()="<< structure[i].size() << endl;

        //Check there is the same number of points in the CNT and in the set
        if (np != vals.size())
        {
            cout << "Error in Apply_displacements_to_cnts. The number of points in CNT " << i << " (" << np << " points) is different from the number of points in set " << set_name << " (" << vals.size() << " points)." << endl;
            return 0;
        }

        //Iterate over the values, i.e., the points in the current CNT
        //Note that nodes are actually in reverse order given the way
        //CNTs are generated (and meshed) in Abaqus
        //This this loop goes in reverse order
        for (int j = vals.size() - 1; j >= 0; j--)
        {
            //Get values at node j
            const odb_FieldValue val = vals[j];
            //Get previous values
            const odb_FieldValue val_prev = vals_prev[j];
            //Output node label
            //cout << "  Node: " << val.nodeLabel() << endl;
            
            //Get the data of the displacements
            int numComp = 0; //This integer is needed to call data() in the line below
            const float* const data = val.data(numComp);
            //Get the data of the previous displacements
            int numComp_prev = 0; //This integer is needed to call data() in the line below
            const float* const data_prev = val_prev.data(numComp_prev);

            //Displacements are in data, where:
            //data[0] corresponds to displacement in x
            //data[1] corresponds to displacement in y
            //data[2] corresponds to displacement in z
            // 
            //Note data[i] is the total displacement along diraction i for a given frame, 
            //not the displacement with respect to the previous frame
            //Thus, calculate the increment of displacement with respect to the previous frame
            //Otherwise I would need the initial position of the points for each frame
            double disp_x = (double)data[0] - (double)data_prev[0];
            double disp_y = (double)data[1] - (double)data_prev[1];
            double disp_z = (double)data[2] - (double)data_prev[2];

            //Save displacements in the binary file
            otec_points.write((char*)&disp_x, double_size);
            otec_points.write((char*)&disp_y, double_size);
            otec_points.write((char*)&disp_z, double_size);
        }
    }

    //Close binary file
    otec_points.close();

    return 1;
}
//This function generates the set name of CNT i, which follows the convention: CNT-i-NODES
string Extract_From_ODB::Get_cnt_set_name(const int& cnt_i)
{
    return ("CNT-" + to_string(cnt_i) + "-NODES");
}
//This function gets the number of GNPs in the sample 
int Extract_From_ODB::Get_number_of_gnps_and_vertices_inside(int& N, vector<vector<int> >& vertices_in)
{
    //Get the sample cuboid
    cuboid sample;
    if (!Get_sample_cuboid(sample))
    {
        cout << "Error in Get_number_of_gnps_and_vertices_inside while calling Get_sample_cuboid" << endl;
        return 0;
    }

    //Open the file with the GNP geometric data
    ifstream gnp_file;
    gnp_file.open("gnp_data.csv");
    if (!gnp_file) {
        cout << "Error in Get_number_of_gnps_and_vertices_inside: Failed to open file with GNP geometric data gnp_data.csv." << endl;
        return 0;
    }

    //String to read lines from file
    string line;

    //Variable to count the GNPs
    int n_gnps = 0;

    //Read the file
    while (getline(gnp_file, line))
    {
        //cout << "GS=" << n_gnps << endl;
        //Read a line from the file and it in a string stream
        stringstream ss(line);

        //Variables to store the GNP geometry
        double lx, ly, t;

        //Variables to store the angles
        double angle_y, angle_z;

        //Variable to store the GNP center
        double Cx, Cy, Cz;

        //Read the values while ignoring the commas
        ss >> lx; ss.ignore();
        ss >> ly;  ss.ignore();
        ss >> t; ss.ignore();
        ss >> angle_y; ss.ignore();
        ss >> angle_z; ss.ignore();
        ss >> Cx; ss.ignore();
        ss >> Cy; ss.ignore();
        ss >> Cz;

        //Calculate sine and cosines of angles
        double cosy = cos(angle_y);
        double siny = sin(angle_y);
        double cosz = cos(angle_z);
        double sinz = sin(angle_z);

        //Calculate half values of length and thickenss
        double lx_half = 0.5 * lx;
        double ly_half = 0.5 * ly;
        double t_half = 0.5 * t;

        //Calculate the vertex coodinates with the GNP centered at the origin
        double GS_vertices[][3] = { {lx_half, ly_half, t_half},
            {-lx_half, ly_half, t_half},
            {-lx_half, -ly_half, t_half},
            {lx_half, -ly_half, t_half},
            {lx_half, ly_half, -t_half},
            {-lx_half, ly_half, -t_half},
            {-lx_half, -ly_half, -t_half},
            {lx_half, -ly_half, -t_half} };

        //Vector to store the vertices inside the sample
        vector<int> v_inside;

        //Iterate over the number of vetices
        for (int i = 0; i < 8; i++)
        {
            //Calculate coordinates of vertex i
            double new_vertex[] = {0.0, 0.0, 0.0};
            if (!Calculate_vertex_coordinates(GS_vertices[i], Cx, Cy, Cz, cosy, siny, cosz, sinz, new_vertex))
            {
                cout << "Error in Get_number_of_gnps_and_vertices_inside when calling Calculate_vertex_coordinates." << endl;
                return 0;
            }
            //cout << "v=" << i << ": " << new_vertex[0] << ", " << new_vertex[1] << ", " << new_vertex[2] << endl;

            //Check if vertex i is inside the sample
            if (Is_vertex_inside_cuboid(new_vertex, sample)) 
            {
                //Vertex i is inside the sample, so add it to the vector of
                //vertices inside the sample
                v_inside.push_back(i);
                //cout << "    v=" << i << " inside" << endl;
            }
        }

        //Check if there were vertices outside the sample
        if (v_inside.size() == 8)
        {
            //All vertices are inside the sample, so add an empty vector
            vertices_in.push_back(vector<int>());
        }
        else
        {
            //Not all vertices are inside the sample,so add the vector v_inside
            vertices_in.push_back(v_inside);
        }

        //Increase the number of GNPs
        n_gnps++;
    }

    //Set the number of GNPs
    N = n_gnps;

    return 1;
}
//This function reads the sample geomtry from a csv file (sample_geom.csv)
int Extract_From_ODB::Get_sample_cuboid(cuboid& sample)
{
    //Open the file with the sample geometry data
    ifstream sample_file;
    sample_file.open("sample_geom.csv");
    if (!sample_file) {
        cout << "Error in Get_sample_cuboid: Failed to open file with sample geometry data sample_geom.csv." << endl;
        return 0;
    }

    //String to read lines from file
    string line;

    //Read the first line form the file and store it in a string stream
    getline(sample_file, line);
    stringstream ss_point(line);

    //Read the lower left corner of the sample while ignoring the commas in between
    ss_point >> sample.min_x;
    ss_point.ignore();
    ss_point >> sample.min_y;
    ss_point.ignore();
    ss_point >> sample.min_z;

    //Read the second line form the file and store it in a string stream
    getline(sample_file, line);
    stringstream ss_size(line);

    //Read the dimensions of the sample along each direction while ignoring the commas in between
    ss_size >> sample.lx;
    ss_size.ignore();
    ss_size >> sample.ly;
    ss_size.ignore();
    ss_size >> sample.lz;

    //Calculate the coordinates of the sample's boundaries opposite to those given by the coordinates of origin
    sample.max_x = sample.min_x + sample.lx;
    sample.max_y = sample.min_y + sample.ly;
    sample.max_z = sample.min_z + sample.lz;

    //Close file
    sample_file.close();

    return 1;
}
//This function calculates the coordinates of the GNP vertex in the final position
int Extract_From_ODB::Calculate_vertex_coordinates(const double vertex[], const double& Cx, const double& Cy, const double& Cz, const double& cosy, const double& siny, const double& cosz, const double& sinz, double new_vertex[])
{
    //Calculate the coordinates of new_vertex after being rotated
    double tmp1 = cosy * vertex[0];
    double tmp2 = siny * vertex[2];
    new_vertex[0] = cosz * tmp1 - sinz * vertex[1] + cosz * tmp2;
    new_vertex[1] = sinz * tmp1 + cosz * vertex[1] + sinz * tmp2;
    new_vertex[2] = cosy * vertex[2] - siny * vertex[0];

    //Apply tanstaltion 
    new_vertex[0] = new_vertex[0] + Cx;
    new_vertex[1] = new_vertex[1] + Cy;
    new_vertex[2] = new_vertex[2] + Cz;

    return 1;
}
//Thin function checks if a vertex is inside the sample
bool Extract_From_ODB::Is_vertex_inside_cuboid(const double vertex[], const cuboid& sample)
{
    if (vertex[0] - sample.min_x > Zero &&
        vertex[0] - sample.max_x < Zero &&
        vertex[1] - sample.min_y > Zero &&
        vertex[1] - sample.max_y < Zero &&
        vertex[2] - sample.min_z > Zero &&
        vertex[2] - sample.max_z < Zero)
    {
        return true;
    }

    return false;
}
//This function extracts the GNP data
int Extract_From_ODB::Extract_gnp_data_for_frame(const int& frame, const int& n_gnps, const vector<vector<int> >& vertices_in, odb_Assembly& root_assy, odb_FieldOutput& current_fieldU)
{
    //Filename for current frame
    string filename_points = "GNPs_disp_F" + to_string(frame) + ".dat";

    //Open binary file for incremental displacements of CNT nodes in current frame
    ofstream otec_points(filename_points.c_str(), ios::binary | ios::out);

    //Get the size of a double
    streamsize double_size = sizeof(double);

    //Create vector with numbers from 0 to 8
    vector<int> all_vertices = Vector_for_all_vertices_inside();

    //Iterate over the number of GNPs
    for (int i = 0; i < n_gnps; i++)
    {
        //Get the number of vertices of the GNP that are inside the sample
        int n_v = (vertices_in[i].empty()) ? 8: (int)vertices_in[i].size();
        //cout << "GNP=" << i << " n_v=" << n_v << endl;

        //Get the vector with the vertices inside the sample
        vector<int> v_inside = (vertices_in[i].empty()) ? all_vertices : vertices_in[i];

        //Iterate over the vertices of the GNP
        for (int j = 0; j < n_v; j++)
        {
            //Get the vertex number
            int v = v_inside[j];

            //Get the set name for the node v in GNP i
            string set_name = Get_gnp_set_name(i, v);
            //cout << "v=" << v << " set_name=" << set_name << endl;

            //Variable to get the displacement objects of the set
            odb_FieldOutput set_disp;

            //Sometimes abaqus does not generate the set, so catch that error if this happens
            try {
                //Access set from root assembly
                odb_Set& gnp_set = root_assy.nodeSets()[set_name.c_str()];

                //Get the displacement objects of the set
                set_disp = current_fieldU.getSubset(gnp_set);
                //hout << "current_fieldU" << endl;
            }
            catch (...) {
                cout << "Error in Extract_gnp_data_for_frame." << endl;
                cout << "Error while accessing set " << set_name.c_str() << "." << endl;
                return 0;
            }
            //hout << "nodeSets" << endl;

            //Get the sequence of values of the displacement object for matrix1
            const odb_SequenceFieldValue& vals = set_disp.values();

            //Check the size of vals is 1 (since it contains one node)
            if (vals.size() != 1)
            {
                cout << "Error in Extract_gnp_data_for_frame. The number of node set " << set_name << " is not 1. Size is " << vals.size() << endl;
                return 0;
            }

            //Get the acutal values
            const odb_FieldValue val = vals[0];
            //Output node label
            //cout << "  Node: " << val.nodeLabel() << endl;
            //Get the data of the displacements
            int numComp = 0; //This integer is needed to call data() in the line below
            const float* const data = val.data(numComp);
            //hout << "disp=" << data[0] << " " << data[1] << " " << data[2] << endl;

            //Get the displacements as a double
            double disp_vx = (double)data[0];
            double disp_vy = (double)data[1];
            double disp_vz = (double)data[2];

            //Save displacements in the binary file
            otec_points.write((char*)&disp_vx, double_size);
            otec_points.write((char*)&disp_vy, double_size);
            otec_points.write((char*)&disp_vz, double_size);
        }
    }

    //Close file
    otec_points.close();

    return 1;
}
//This function generates a vector with the numbers from 1 to 8
vector<int> Extract_From_ODB::Vector_for_all_vertices_inside()
{
    vector<int> tmp(8, 0);

    //Set values 1 to 7
    tmp[1] = 1;
    tmp[2] = 2;
    tmp[3] = 3;
    tmp[4] = 4;
    tmp[5] = 5;
    tmp[6] = 6;
    tmp[7] = 7;

    return tmp;
}
//This function generates the set name of GNP gnp_i, that contains the node vertex
//Set naming follows the convention: GS-i-_N-vertex
string Extract_From_ODB::Get_gnp_set_name(const int& gnp_i, const int& vertex)
{
    return ("GS-" + to_string(gnp_i) + "_N-" + to_string(vertex));
}
//This function extracts the datat for the matrix sets
int Extract_From_ODB::Extract_matrix_data_for_frame(const int& frame, odb_Assembly& root_assy, odb_FieldOutput& current_fieldU)
{
    //Filename for current frame
    string filename_points = "Matrix_disp_F" + to_string(frame) + ".dat";

    //Open binary file for incremental displacements of CNT nodes in current frame
    ofstream otec_points(filename_points.c_str(), ios::binary | ios::out);

    //Get the size of a double
    streamsize double_size = sizeof(double);

    //Names of the sets for the matrix corners (it is hard coded in the python scritp too)
    string set0 = "MATRIX0";
    string set1 = "MATRIX1";

    //Variable to store displacements
    double disp[] = {0.0, 0.0, 0.0};

    //Get the displacements for set0
    if (!Get_displacement_for_node_set(set0, root_assy, current_fieldU, disp))
    {
        cout << "Erro in Extract_matrix_data_for_frame when calling Get_displacement_for_node_set (set0)." << endl;
        return 0;
    }

    //Save displacements to file
    otec_points.write((char*)&disp[0], double_size);
    otec_points.write((char*)&disp[1], double_size);
    otec_points.write((char*)&disp[2], double_size);

    //Get the displacements for set1
    if (!Get_displacement_for_node_set(set1, root_assy, current_fieldU, disp))
    {
        cout << "Erro in Extract_matrix_data_for_frame when calling Get_displacement_for_node_set (set1)." << endl;
        return 0;
    }

    //Save displacements to file
    otec_points.write((char*)&disp[0], double_size);
    otec_points.write((char*)&disp[1], double_size);
    otec_points.write((char*)&disp[2], double_size);

    //Close the file
    otec_points.close();

    return 1;
}
//This function gets the displacements for a given set node that contains only one node
int Extract_From_ODB::Get_displacement_for_node_set(const string& set_name, odb_Assembly& root_assy, odb_FieldOutput& current_fieldU, double disp[])
{

    //Get the displacement objects of the set
    odb_FieldOutput set_disp;

    //Access set from root assembly
    try
    {
        odb_Set& sub_set = root_assy.nodeSets()[set_name.c_str()];

        //Get the displacement objects of the set
        set_disp = current_fieldU.getSubset(sub_set);
    }
    catch (...)
    {
        cout << "Error while accessing set " << set_name << " form Abaqus odb file." << endl;
        return 0;
    }

    //Get the sequence of values of the displacement object for matrix1
    const odb_SequenceFieldValue& vals = set_disp.values();

    //Check the size of vals is 1 (since it contains one node)
    if (vals.size() != 1)
    {
        cout << "Error in Get_displacement_for_node_set. The number of nodes in set " << set_name << " is not 1. Size is " << vals.size() << endl;
        return 0;
    }

    //Get the acutal values
    const odb_FieldValue val = vals[0];
    //Output node label
    //cout << "  Node: " << val.nodeLabel() << endl;
    //Get the data of the displacements
    int numComp = 0; //This integer is needed to call data() in the line below
    const float* const data = val.data(numComp);
    //hout << "disp=" << data[0] << " " << data[1] << " " << data[2] << endl;

    //Use the values stored in data1 as the components of the point disp
    disp[0] = (double)data[0];
    disp[1] = (double)data[1];
    disp[2] = (double)data[2];

    return 1;
}
