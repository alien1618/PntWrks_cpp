#include "../pntwrks.h"

//----------------------------------------
// FUNCTIONS FOR PRINTING POINTS TO FILE
//----------------------------------------

void POINTSET::printTXT(vector<double> U, string filename)
{
    cout << "Printing solution to file..." << endl;
    string direc;
    direc = output_dir+"/" + filename +"_txt/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    fname = output_dir+"/" + filename +"_txt/" + filename+".txt";
    ofstream outfile(fname.c_str());
    for (int i = 1; i <= TotalPoints; i++)
    {outfile << points[i].x << "\t" << points[i].y<<  "\t"  << points[i].z << "\t" << U[i] << endl;}
    outfile.close();
}
void POINTSET::printTXT(vector<double> U, string filename, int t)
{
    cout << "Printing solution to file..." << endl;
    string direc;
    direc = output_dir+"/" + filename +"_txt/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    stringstream tt;
    tt << t;
    if (t>=1 && t <= 9)
    {fname = direc + filename+ "00"+tt.str()+".txt";}
    else if (t >= 10 && t <= 99)
    {fname = direc + filename+ "0"+tt.str()+".txt";}
    else
    {fname = direc + filename+ tt.str()+".txt";}
    ofstream outfile(fname.c_str());
    for (int i = 1; i <= TotalPoints; i++)
    {outfile << points[i].x << "\t" << points[i].y<<  "\t"  << points[i].z << "\t" << U[i] << endl;}
    outfile.close();
}
void POINTSET::printTXT(vector<int> U, string filename, int t)
{
    cout << "Printing solution to file..." << endl;
    string direc;
    direc = output_dir+"/" + filename +"_txt/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    stringstream tt;
    tt << t;
    if (t>=1 && t <= 9)
    {fname = direc + filename+ "00"+tt.str()+".txt";}
    else if (t >= 10 && t <= 99)
    {fname = direc + filename+ "0"+tt.str()+".txt";}
    else
    {fname = direc + filename+ tt.str()+".txt";}
    ofstream outfile(fname.c_str());
    for (int i = 1; i <= TotalPoints; i++)
    {outfile << points[i].x << "\t" << points[i].y<<  "\t"  << points[i].z << "\t" << U[i] << endl;}
    outfile.close();
}
void POINTSET::printVTK(vector<double> U, string variable, int t)
{
    cout << "Printing solution to file..." << endl;
    string direc;
    direc = output_dir+"/" + variable +"_vtk/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    stringstream tt;
    tt << t;

    if (t>=1 && t <= 9)
    {fname = direc + variable + "00"+tt.str()+".vtk";}
    else if (t >= 10 && t <= 99)
    {fname = direc + variable + "0"+tt.str()+".vtk";}
    else
    {fname = direc + variable + tt.str()+".vtk";}

    ofstream outfile(fname.c_str());
    outfile << "# vtk DataFile Version 2.0" << endl;
    outfile << "Unstructured Grid Example" << endl;
    outfile << "ASCII" << endl;
    outfile << "DATASET UNSTRUCTURED_GRID" << endl;
    outfile << "POINTS " << TotalPoints << " float" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile << points[i].x << "\t" << points[i].y <<  "\t"  << points[i].z << endl;}

    outfile << endl << "POINT_DATA " << TotalPoints << endl;
    outfile << "SCALARS data float 1"<< endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile << U[i] << endl;}
    outfile.close();
}
void POINTSET::printVTK(vector<int> U, string variable, int t)
{
    cout << "Printing solution to file..." << endl;
    string direc;
    direc = output_dir+"/" + variable +"_vtk/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    stringstream tt;
    tt << t;

    if (t>=1 && t <= 9)
    {fname = direc + variable+ "00"+tt.str()+".vtk";}
    else if (t >= 10 && t <= 99)
    {fname = direc + variable+ "0"+tt.str()+".vtk";}
    else
    {fname = direc + variable+ tt.str()+".vtk";}

    ofstream outfile(fname.c_str());
    outfile << "# vtk DataFile Version 2.0" << endl;
    outfile << "Unstructured Grid Example" << endl;
    outfile << "ASCII" << endl;
    outfile << "DATASET UNSTRUCTURED_GRID" << endl;
    outfile << "POINTS " << TotalPoints << " float" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile << points[i].x << "\t" << points[i].y <<  "\t"  << points[i].z << endl;}

    outfile << endl << "POINT_DATA " << TotalPoints << endl;
    outfile << "SCALARS data float 1"<< endl;
    outfile << "LOOKUP_TABLE default" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile << U[i] << endl;}
    outfile.close();
}
