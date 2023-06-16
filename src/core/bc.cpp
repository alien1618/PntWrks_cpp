#include "../pntwrks.h"

BOUNDARY_CONDITION::BOUNDARY_CONDITION() //constructor
{
    total = 0;
    points.resize(1);
	values.resize(1);
}
void BOUNDARY_CONDITION::assignDBC(POINTSET pointset, string salomeface_fname, double value)
{
    cout << "Assigning Surface DBCs..." << endl;
    if (salomeface_fname == "")
    {cout << "ERROR: Direchlet boundary conditions File Name is Invalid..." << endl;exit(0);}
    else
    {
        POINTSET boundary(salomeface_fname);
        if (boundary.TotalPoints == 0)
        {cout << "ERROR: Direchlet boundary conditions Points are Empty..." << endl;exit(0);}
        else
        {
            vector<int> pnt_num(boundary.TotalPoints+1);
            for (int i = 1; i <= boundary.TotalPoints; i++)
            {
                for (int j = 1; j <= pointset.TotalPoints; j++)
                {
                    double d = boundary.points[i].distance(pointset.points[j]);
                    if (d <= 1e-5)
                    {pnt_num[i] = j;}
                }
            }
            for (int i = 1; i <= boundary.TotalPoints; i++)
            {
				total++;
				points.push_back(1);
				values.push_back(1);
				points[total] = pnt_num[i];
				values[total] = value;

			}
            cout << "SUCCESS: Direchlet boundary conditions read successfully..." << endl;
        }
    }
}
void BOUNDARY_CONDITION::assignDBC(vector<POINT> pnts, int TotalPoints, string edge, double location, double value)
{
    cout << "Assigning Grid DBCs..." << endl;
    vector<int> ndnum(1);
    int m = 0;
    if (edge == "x")
    {
        for (int i = 1; i <= TotalPoints; i++)
        {
            if ((pnts[i].x >= location-0.01*abs(location)) & (pnts[i].x <= location+0.01*abs(location)))
            {
                m++;
                ndnum.push_back(1);
                ndnum[m] = i;
            }
        }
    }
    if (edge == "y")
    {
        for (int i = 1; i <= TotalPoints; i++)
        {
            if ((pnts[i].y >= location-0.01*abs(location)) & (pnts[i].y <= location+0.01*abs(location)))
            {
                m++;
                ndnum.push_back(1);
                ndnum[m] = i;
            }
        }
    }
    if (edge == "z")
    {
        for (int i = 1; i <= TotalPoints; i++)
        {
            if ((pnts[i].z >= location-0.01*abs(location)) & (pnts[i].z <= location+0.01*abs(location)))
            {
                m++;
                ndnum.push_back(1);
                ndnum[m] = i;
            }
        }
    }

    if (m == 0)
    {cout << "ERROR: Direchlet boundary conditions Points are Empty..." << endl;exit(0);}
    else
    {	
        for (int i = 1; i <= m; i++)
        {
			total++;
			points.push_back(1);
			values.push_back(1);
			points[total] = ndnum[i];
			values[total] = value;
		}
        cout << "SUCCESS: Direchlet boundary conditions assigned successfully..." << endl;
    }
}
BOUNDARY_CONDITION::~BOUNDARY_CONDITION()
{
    //destructor
}
BOUNDARY_CONDITIONS::BOUNDARY_CONDITIONS()
{
	Vx.total = 0;
	Vx.points.resize(1);
	Vx.values.resize(1);

	Vy.total = 0;
	Vy.points.resize(1);
	Vy.values.resize(1);


	Vz.total = 0;
	Vz.points.resize(1);
	Vz.values.resize(1);


	P.total = 0;
	P.points.resize(1);
	P.values.resize(1);


	U.total = 0;
	U.points.resize(1);
	U.values.resize(1);

	phi.total = 0;
	phi.points.resize(1);
	phi.values.resize(1);
}
BOUNDARY_CONDITIONS::~BOUNDARY_CONDITIONS()
{
    //destructor
}
vector<double> setVector(vector<double> U, BOUNDARY_CONDITION direchletBC)
{
	 for (int l = 1; l <= direchletBC.total; l++)
     {U[direchletBC.points[l]] = direchletBC.values[l];}
     return U;
}
vector<double> setVector(int TotalPoints, double val)
{
	vector<double> U(TotalPoints+1);
    for (int r = 1; r <= TotalPoints; r++)
    {U[r] = val;}
    return U;
}
