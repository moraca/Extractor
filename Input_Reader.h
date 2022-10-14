#ifndef INPUT_READER_H
#define INPUT_READER_H

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

//Parameters to extract data from odb file
struct Extraction_para 
{
    //Type of fillers or fillers in the simulation
    string filler_type;
    //Path to obd file or name if it is in the same folder as the executable
    string odb_file;
    //Name of the simulation step in Abaqus
    string step_name;
    //Starting frame to extract data from
    int first_frame;
};


class Input
{
public:

    //Data members
    Extraction_para extract_para;

    //Constructor
    Input() {};

    //Member functions
    int Read_input_file(ifstream& infile);
    int Read_extraction_parameters(Extraction_para& extract_para, ifstream& infile);
    string Get_Line(ifstream& infile)const;


};
//---------------------------------------------------------------------------
#endif
//===========================================================================
