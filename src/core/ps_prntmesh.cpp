#include "../pntwrks.h"

using namespace std;

//-------------------------------------------------------------------------------
// FUNCTIONS FOR PRINTING MESH IN VTK FORMAT
//-------------------------------------------------------------------------------
void POINTSET::printMeshVTK(vector<double> U, string filename, int t)
{
    cout << "Printing mesh data to " << filename << "..." << endl;
    string direc;
    direc = output_dir+"/" + filename +"_vtk/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    stringstream tt;
    tt << t;

    if (t>=1 && t <= 9)
    {fname =  direc + filename+ "00"+tt.str()+".vtk";}
    else if (t >= 10 && t <= 99)
    {fname =  direc + filename+ "0"+tt.str()+".vtk";}
    else
    {fname =  direc + filename+ tt.str()+".vtk";}
    int sum = 0;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        int num=0;
        if (SurfaceShape == 1)
        {num = 4;}
        else if (SurfaceShape == 2)
        {num = 5;}
        else
        {cout << "ERROR2: unrecognized element type" << endl;}
        sum = sum + num;
    }
    ofstream outfile1(fname.c_str());
    outfile1 << "# vtk DataFile Version 1.0" << endl;
    outfile1 << "Cube example" << endl;
    outfile1 << "ASCII" << endl;
    outfile1 << endl;
    outfile1 <<"DATASET POLYDATA" << endl;
    outfile1 << "POINTS " << TotalPoints << " double" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << points[i].x << "\t" << points[i].y <<  "\t"  << points[i].z << endl;}
    outfile1 << endl;
    outfile1 << "POLYGONS " << TotalSurfaces << "\t" << sum << endl;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        outfile1 << SurfaceNds/ElementOrder << "\t";
        for (int j = 1; j <= SurfaceNds/ElementOrder; j++)
        {outfile1 << (surfaces[i][j]-1)<< "\t";}
        outfile1 << endl;
    }
    outfile1 << endl;
    outfile1 <<"POINT_DATA " << TotalPoints << endl;
    outfile1 << "SCALARS myscalars double"<< endl;
    outfile1 << "LOOKUP_TABLE custom_table" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << U[i] << endl;}
    outfile1.close();
}
void POINTSET::printMeshVTK(vector<int> U, string filename, int t)
{
    cout << "Printing mesh data to " << filename << "..." << endl;
    string direc;
    direc = output_dir+"/" + filename +"_vtk/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    stringstream tt;
    tt << t;

    if (t>=1 && t <= 9)
    {fname =  direc + filename+ "00"+tt.str()+".vtk";}
    else if (t >= 10 && t <= 99)
    {fname =  direc + filename+ "0"+tt.str()+".vtk";}
    else
    {fname =  direc + filename+ tt.str()+".vtk";}
    int sum = 0;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        int num=0;
        if (SurfaceShape == 1)
        {num = 4;}
        else if (SurfaceShape == 2)
        {num = 5;}
        else
        {cout << "ERROR2: unrecognized element type" << endl;}
        sum = sum + num;
    }
    ofstream outfile1(fname.c_str());
    outfile1 << "# vtk DataFile Version 1.0" << endl;
    outfile1 << "Cube example" << endl;
    outfile1 << "ASCII" << endl;
    outfile1 << endl;
    outfile1 <<"DATASET POLYDATA" << endl;
    outfile1 << "POINTS " << TotalPoints << " double" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << points[i].x << "\t" << points[i].y <<  "\t"  << points[i].z << endl;}
    outfile1 << endl;
    outfile1 << "POLYGONS " << TotalSurfaces << "\t" << sum << endl;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        outfile1 << SurfaceNds/ElementOrder << "\t";
        for (int j = 1; j <= SurfaceNds/ElementOrder; j++)
        {outfile1 << (surfaces[i][j]-1)<< "\t";}
        outfile1 << endl;
    }
    outfile1 << endl;
    outfile1 <<"POINT_DATA " << TotalPoints << endl;
    outfile1 << "SCALARS myscalars double"<< endl;
    outfile1 << "LOOKUP_TABLE custom_table" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << U[i] << endl;}
    outfile1.close();
}
void POINTSET::printMeshVTK(vector<vector<double> > u, RESOLUTION g, string filename, int t)
{
    string direc;
    direc = output_dir+"/" + filename +"_vtk/";
    mkdir(direc.c_str(), 0777);
    
    if (g.x*g.y*g.z != TotalPoints)
    {cout << "ERROR in POINTSET::printMeshVTK. Grid density must be equal to mesh total points." << endl; exit(0);}
    else
    {
        int m = 0;
        vector<double> U(TotalPoints+1);
        for (int j = 1; j <= g.y; j++)
        {
            for (int i = 1; i <= g.x; i++)
            {
                m++;
                U[m] = u[i][j];
            }
        }

        cout << "Printing " << filename << " mesh data to file..." << endl;
        string fname;
        stringstream tt;
        tt << t;

        if (t>=1 && t <= 9)
        {fname =  direc + filename+ "00"+tt.str()+".vtk";}
        else if (t >= 10 && t <= 99)
        {fname =  direc + filename+ "0"+tt.str()+".vtk";}
        else
        {fname =  direc + filename+ tt.str()+".vtk";}
        int sum = 0;
        for (int i = 1; i <= TotalSurfaces; i++)
        {
            int num=0;
            if (SurfaceShape == 1)
            {num = 4;}
            else if (SurfaceShape == 2)
            {num = 5;}
            else
            {cout << "ERROR2: unrecognized element type" << endl;}
            sum = sum + num;
        }
        ofstream outfile1(fname.c_str());
        outfile1 << "# vtk DataFile Version 1.0" << endl;
        outfile1 << "Cube example" << endl;
        outfile1 << "ASCII" << endl;
        outfile1 << endl;
        outfile1 <<"DATASET POLYDATA" << endl;
        outfile1 << "POINTS " << TotalPoints << " double" << endl;
        for (int i = 1; i <= TotalPoints; i++)
        {outfile1 << points[i].x << "\t" << points[i].y <<  "\t"  << points[i].z << endl;}
        outfile1 << endl;
        outfile1 << "POLYGONS " << TotalSurfaces << "\t" << sum << endl;
        for (int i = 1; i <= TotalSurfaces; i++)
        {
            outfile1 << SurfaceNds/ElementOrder << "\t";
            for (int j = 1; j <= SurfaceNds/ElementOrder; j++)
            {outfile1 << (surfaces[i][j]-1)<< "\t";}
            outfile1 << endl;
        }
        outfile1 << endl;
        outfile1 <<"POINT_DATA " << TotalPoints << endl;
        outfile1 << "SCALARS myscalars double"<< endl;
        outfile1 << "LOOKUP_TABLE custom_table" << endl;
        for (int i = 1; i <= TotalPoints; i++)
        {outfile1 << U[i] << endl;}
        outfile1.close();
    }
}
void POINTSET::printVectorsVTK(vector<double> U, vector<double> Vx, vector<double> Vy, vector<double> Vz, string filename, int t)
{
    cout << "Printing mesh data to " << filename << "..." << endl;
    string direc;
    direc = output_dir+"/" + filename +"_vtk/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    stringstream tt;
    tt << t;

    if (t>=1 && t <= 9)
    {fname = direc + filename+ "00"+tt.str()+".vtk";}
    else if (t >= 10 && t <= 99)
    {fname = direc + filename+ "0"+tt.str()+".vtk";}
    else
    {fname = direc + filename+ tt.str()+".vtk";}
    
    int sum = 0;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        int num=0;
        if (SurfaceShape == 1)
        {num = 4;}
        else if (SurfaceShape == 2)
        {num = 5;}
        else
        {cout << "ERROR2: unrecognized element type. SurfaceShape = " << SurfaceShape << endl;}
        sum = sum + num;
    }
    ofstream outfile1(fname.c_str());
    outfile1 << "# vtk DataFile Version 1.0" << endl;
    outfile1 << "Cube example" << endl;
    outfile1 << "ASCII" << endl;
    outfile1 << endl;
    outfile1 <<"DATASET POLYDATA" << endl;
    outfile1 << "POINTS " << TotalPoints << " double" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << points[i].x << "\t" << points[i].y <<  "\t"  << points[i].z << endl;}
    outfile1 << endl;
    outfile1 << "POLYGONS " << TotalSurfaces << "\t" << sum << endl;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        outfile1 << SurfaceNds/ElementOrder << "\t";
        for (int j = 1; j <= SurfaceNds/ElementOrder; j++)
        {outfile1 << (surfaces[i][j]-1)<< "\t";}
        outfile1 << endl;
    }
    outfile1 << endl;
    outfile1 <<"POINT_DATA " << TotalPoints << endl;
    outfile1 << "SCALARS myscalars double"<< endl;
    outfile1 << "LOOKUP_TABLE custom_table" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << U[i] << endl;}
    outfile1 << endl;
    outfile1 << "VECTORS vectors double" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << Vx[i] << "\t" << Vy[i] <<"\t" << Vz[i] << endl;}
    outfile1.close();
}

void POINTSET::printUnstructuredGridVTK(vector<double> U, string filename, int t)
{
    cout << "Printing structured grid data to " << filename << "..." << endl;
    string direc;
    direc = output_dir+"/" + filename +"_vtk/";
    mkdir(direc.c_str(), 0777);

    string fname;
    stringstream tt;
    tt << t;

    if (t>=1 && t <= 9)
    {fname = direc + filename+ "00"+tt.str()+".vtk";}
    else if (t >= 10 && t <= 99)
    {fname = direc + filename+ "0"+tt.str()+".vtk";}
    else
    {fname = direc + filename+ tt.str()+".vtk";}

    int sum=0;
    for (int i = 1; i <= TotalVolumes; i++)
    {
        int num=0;
        if (VolumeShape == 1)
        {num = 4;}
        else if (VolumeShape == 2)
        {num = 8;}
        else
        {cout << "ERROR2: unrecognized element type. VolumeShape = " << VolumeShape << endl;}
        sum = sum + num;
    }
    ofstream outfile1(fname.c_str());
    outfile1 << "# vtk DataFile Version 1.0" << endl;
    outfile1 << "Cube example" << endl;
    outfile1 << "ASCII" << endl;
    outfile1 << endl;
    outfile1 <<"DATASET UNSTRUCTURED_GRID" << endl;
    outfile1 << "POINTS " << TotalPoints << " double" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << points[i].x << "\t" << points[i].y <<  "\t"  << points[i].z << endl;}
    outfile1 << endl;
    outfile1 << "CELLS " << TotalVolumes << "\t" << (sum + TotalVolumes)<< endl;
    for (int i = 1; i <= TotalVolumes; i++)
    {
        outfile1 << VolumeNds<< "\t";
        for (int j = 1; j <= VolumeNds; j++)
        {outfile1 << (volumes[i][j]-1)<< "\t";}
        outfile1 << endl;
    }
    outfile1 << endl;
    outfile1 << "CELL_TYPES " << TotalVolumes << endl;
    for (int i = 1; i <= TotalVolumes; i++)
    {outfile1 << 12 << endl;}
    outfile1 << endl;
    outfile1 <<"POINT_DATA " << TotalPoints << endl;
    outfile1 << "SCALARS myscalars double"<< endl;
    outfile1 << "LOOKUP_TABLE default" << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << U[i] << endl;}
    outfile1.close();
}
void printStructuredGridVTK(vector<double> U, RESOLUTION divs, string filename, int t)
{
    cout << "Printing structured grid data to " << filename << "..." << endl;
    string direc;
    direc = output_dir+"/" + filename +"_vtk/";
    mkdir(direc.c_str(), 0777);

    string fname;
    stringstream tt;
    tt << t;
    if (t>=1 && t <= 9)
    {fname = direc + filename+ "00"+tt.str()+".vtk";}
    else if (t >= 10 && t <= 99)
    {fname = direc + filename+ "0"+tt.str()+".vtk";}
    else
    {fname = direc + filename+ tt.str()+".vtk";}
    ofstream outfile1(fname.c_str());
    outfile1 << "# vtk DataFile Version 1.0" << endl;
    outfile1 << "Cube example" << endl;
    outfile1 << "ASCII" << endl;
    outfile1 << endl;
    outfile1 <<"DATASET STRUCTURED_POINTS" << endl;
    outfile1 << "DIMENSIONS " << divs.x << " " << divs.y << " " << divs.z << endl;
    outfile1 << "ORIGIN " << "0.0 0.0 0.0" << endl;
    outfile1 << "SPACING " << "1.0 1.0 1.0" << endl;
    int points = divs.x*divs.y*divs.z;
    outfile1 << "POINT_DATA " << points << endl;
    outfile1 <<"SCALARS data double" << endl;
    outfile1 << "LOOKUP_TABLE default" << endl;
    int m = 0;
    for (int k = 1; k <= divs.z; k++)
    {
        for (int j = 1; j <= divs.y; j++)
        {
            for (int i = 1; i <= divs.x; i++)
            {
                m++;
                outfile1 << U[m] << "\t";
            }
        }
    }
    outfile1.close();
}
void printElements(vector<vector<int> > surfaces, int SurfaceNds, int TotalSurfaces, string filename)
{
    cout << "Printing mesh data to " << filename << "..." << endl;
    string direc = output_dir+"/gmtry/";
    mkdir(direc.c_str(), 0777);
    string fname;
    fname = direc + filename + ".txt";
    ofstream outfile1(fname.c_str());
    for (int i = 1; i <= TotalSurfaces; i++)
    {
        for (int j = 1; j <= SurfaceNds; j++)
        {outfile1 << (surfaces[i][j]-1)<< "\t";}
        outfile1 << endl;
    }
    outfile1.close();
}
void pltctrl(vector<POINT> points, int TotalPoints, int n_t, int printfrequency)
{
    string direc = output_dir+"/gmtry/";
    mkdir(direc.c_str(), 0777);

    string fname;
    double xmin = points[1].x;
    double xmax = points[1].x;
    double ymin = points[1].y;
    double ymax = points[1].y;
    for (int i = 1; i <= TotalPoints; i++)
    {
        if (points[i].x <= xmin)
        {xmin =  points[i].x;}
        if (points[i].x >= xmax)
        {xmax =  points[i].x;}
        if (points[i].y <= ymin)
        {ymin =  points[i].y;}
        if (points[i].y >= ymax)
        {ymax =  points[i].y;}
    }
    fname = output_dir+"/gmtry/"+"pltctrl.txt";
    ofstream outfile(fname.c_str());
    outfile << n_t << "\t" << printfrequency << "\t" << xmin << "\t" <<  xmax << "\t" << ymin << "\t" << ymax << endl;
    outfile.close();
}
void POINTSET::printMesh()
{
    string direc;
    direc = output_dir+"/" + "msh/";
    mkdir(direc.c_str(), 0777);
    
    string fname;
    fname = direc + "nds.txt";
    ofstream outfile1(fname.c_str());
    outfile1 << TotalPoints << "\t" << 3 << endl;
    for (int i = 1; i <= TotalPoints; i++)
    {outfile1 << points[i].x << "\t" << points[i].y <<  "\t"  << points[i].z << endl;}
        
        
    fname = direc + "surfs.txt";
    ofstream outfile2(fname.c_str());
    outfile2 << TotalSurfaces << "\t" << SurfaceNds << "\t" << SurfaceShape << "\t" << ElementOrder << endl;
    for (int i = 1; i <= TotalSurfaces; i++)
    {
		for (int j = 1; j <= SurfaceNds; j++)
		{outfile2 << surfaces[i][j] << "\t";}
		outfile2 << endl;
	}   
	
	    fname = direc + "vols.txt";
    ofstream outfile3(fname.c_str());
    outfile3 << TotalVolumes << "\t" << VolumeNds << "\t" << VolumeShape << "\t" << ElementOrder << endl;
    for (int i = 1; i <= TotalVolumes; i++)
    {
		for (int j = 1; j <= VolumeNds; j++)
		{outfile3 << volumes[i][j] << "\t";}
		outfile3 << endl;
	}   
	
    outfile1.close();
    outfile2.close();
    outfile3.close();
}
