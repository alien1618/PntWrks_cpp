#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D flow past implicit obstacles" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(5,1,0.0);
    RESOLUTION divs(250,50,1);
    POINTSET pointset(p, s, divs, "T3");

    cout << "Assigning implicit obstacles..."<< endl;
    POINT cntr1(1, 0.5, -0.5);
    pointset.addCylinderZ(cntr1, 0.125);
    pointset.printMeshVTK(pointset.phi, "phi", 0);

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(1);
    materials.Nu[1] = 0.1;
    materials.ro[1] = 1000;

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
	vector<double> Fx = setVector(pointset.TotalPoints,0);
	vector<double> Fy = setVector(pointset.TotalPoints,0);
	vector<double> Fz = setVector(pointset.TotalPoints,0);
	
    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 10);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.P.assignDBC(pointset.points, pointset.TotalPoints, "x", 5, 0);

    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 2;
    settings.kernel.approximation_method = "RBF";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "LS";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.dt = 0.00001; // use for RBF, MLS, WLS, NSPH
    settings.nt = 600000;
    settings.P_factor = 1000;
    settings.AV_factor = 0.1;
    settings.prnt_freq = 1000;

    cout << "Running solver..."<< endl;
    solveNavierStokes1Phase(pointset, materials, bc, settings, U, Vx, Vy, Vz, Q, Fx, Fy, Fz);
    
    return 0;
}

