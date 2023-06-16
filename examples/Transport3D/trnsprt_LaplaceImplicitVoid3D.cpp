#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "3D laplace equation with implicit void" << endl;
    cout << "-----------------------------------------------------------------" << endl;

    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN size(1,1,0.1);
    RESOLUTION divs(30,30,2);
    POINTSET pointset(p, size, divs, "T4");

    cout << "Assigning implicit obstacles..."<< endl;
    POINT cntr(0.5, 0.5, 0);
    pointset.addCylinderZ(cntr,0.25);
    pointset.computeVOF(0.5*pointset.AvgPntSpacing);
	for (int i =1; i<=pointset.TotalPoints; i++)
	{pointset.phi[i] = 1-pointset.phi[i];}
    
    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.D[1] = 1;
    materials.D[2] = 0.0001;

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
	
    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITION bc;
    bc.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 1);
    bc.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);

    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "SP"; //critical for getting correct solution for implicit void
    settings.kernel.radius_ratio = 3;
    settings.kernel.approximation_method = "RBF";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "S4";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.dt = 0.00001;
    settings.nt = 20000;
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    solveTransportSteadyState(pointset, materials, bc, settings, U, Vx, Vy, Vz, Q);
    
    return 0;
}
