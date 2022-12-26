
#include "Input_Reader.h"

//---------------------------------------------------------------------------
//Read data
int Input::Read_input_file(ifstream& infile)
{
	//Read parameters from input file
	if (!Read_extraction_parameters(extract_para, infile))
	{
		cout << "Error in Read_input_file when calling Read_extraction_parameters." << endl;
		return 0;
	}

	return 1;
}
//Read parameters for data extraction
int Input::Read_extraction_parameters(Extraction_para& extract_para, ifstream& infile)
{
	//Read line for particle type
	istringstream istr1(Get_Line(infile));
	istr1 >> extract_para.filler_type;

	//Check a valid filler type was input
	if (extract_para.filler_type != "CNTs_only" &&
		extract_para.filler_type != "GNPs_only" &&
		extract_para.filler_type != "GNP_CNT_mix")
	{
		cout << "Error in Read_extraction_parameters: Invalid filler type, valid option are: CNTs_only, GNPs_only or GNP_CNT_mix." << endl;
		return 0;
	}
	//cout << "extract_para.filler_type=" << extract_para.filler_type << endl;

	//Read line for the path to obd file or name if it is in the same folder as the executable
	istringstream istr2(Get_Line(infile));
	istr2 >> extract_para.odb_file;

	//Check that the odb file exists
	ifstream file(extract_para.odb_file.c_str());
	if (!file.good())
	{
		//File does not exists
		cout << "Error in Read_extraction_parameters: odb file was not found." << endl;
		cout << "Please double check odb filename or path." << endl;
		cout << "Input odb filename/path: " << extract_para.odb_file << endl;
		return 0;
	}
	//cout << "extract_para.odb_file=" << extract_para.odb_file << endl;

	//Read line for the step name
	istringstream istr_step(Get_Line(infile));
	istr_step >> extract_para.step_name;
	//cout << "extract_para.step_name=" << extract_para.step_name << endl;

	//Read line for the first frame
	istringstream istr_first_frame(Get_Line(infile));
	istr_first_frame >> extract_para.first_frame;

	//Check the first frame is a valid frame
	if (extract_para.first_frame < 0)
	{
		cout << "Error in Read_extraction_parameters: Invalid first frame number. " << endl;
		cout << "First frame number should be greater or equal to zero. Input was: " << extract_para.first_frame << "." << endl;
		return 0;
	}

	//Check if first frame is frame 0
	if (extract_para.first_frame == 0)
	{
		//Set first frame as frame 1 and send a message
		extract_para.first_frame = 1;
		cout << "Warning: First frame is frame 0. Frame 0 will be ignored and first frame will be set as 1." << endl;
	}

	return 1;
}
//Read the input data in a whole line (to skip over the comment line starting with a '%')
string Input::Get_Line(ifstream& infile)const
{
	string s;
	//Read the input data in a whole line
	getline(infile, s);
	//to skip over the comment line starting with a '%'
	while (!infile.eof() && s.substr(0, 1) == "%")
		getline(infile, s);
	return s;
}
//===========================================================================
