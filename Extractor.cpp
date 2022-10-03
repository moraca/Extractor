// Extractor.cpp : Este archivo contiene la función "main". La ejecución del programa comienza y termina ahí.
//

#include <iostream>

#include "Input_Reader.h"
#include "Extract_From_ODB.h"

using namespace std;

int main(int argc, char** argv)
{
    //Read input file name into in_file
    string in_file;
    if (argc > 1)
    {
        //Set the input file name
        in_file = argv[1];
    }
    else
    {
        //If no input file name give, use the default one
        in_file = "in.txt";
        cout << "The input file name is:  " << in_file << endl;
    }

    //Open the input file
    ifstream infile;
    infile.open(in_file.c_str());
    if (!infile) { 
        cout << "Failed to open input file: " << in_file << endl;  
        return 0;
    }

    //Time markers for simulation
    time_t it_begin, it_end;
    it_begin = time(NULL);

    //----------------------------------------------------------------------
    //Input file reader
    cout << "======================================================" << endl;
    cout << "Reading input file......" << endl;
    Input Init;
    //Read the input file
    if (!Init.Read_input_file(infile)) 
    {
        cout << "Error when reading input file." << endl;
        return 0;
    }
    it_end = time(NULL);
    cout << "Input file read in " << (int)(it_end - it_begin) << " secs." << endl;

    //Close input file
    infile.close();

    //Extract data from odb file and save into file
    Extract_From_ODB Odb_extractor;
    it_begin = time(NULL);
    if (!Odb_extractor.Extract_data_from_odb(Init))
    {
        cout << "Error when calling Extract_data_from_odb." << endl;
        return 0;
    }
    it_end = time(NULL);
    cout << "======================================================" << endl;
    cout << "Data extracted in " << (int)(it_end - it_begin) << " secs." << endl;

    return 1;
}
