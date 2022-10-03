
#include "Input_Reader.h"

//---------------------------------------------------------------------------
//Read data
int Input::Read_input_file(ifstream& infile)
{
	cout << "Reading input file..." << endl;

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
