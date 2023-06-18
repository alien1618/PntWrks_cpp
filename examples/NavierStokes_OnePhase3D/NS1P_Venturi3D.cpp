#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "3D flow past cylinder" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINTSET pointset("gmtry/Venturi3D/mesh_1450.dat");

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
    double vel = 10;
    BOUNDARY_CONDITIONS bc;
    bc.Vx.assignDBC(pointset, "gmtry/Venturi3D/inlet_1450.dat", vel);
    bc.Vy.assignDBC(pointset, "gmtry/Venturi3D/inlet_1450.dat", 0);
    bc.Vz.assignDBC(pointset, "gmtry/Venturi3D/inlet_1450.dat", 0);
    bc.Vx.assignDBC(pointset, "gmtry/Venturi3D/walls_1450.dat", 0);
    bc.Vy.assignDBC(pointset, "gmtry/Venturi3D/walls_1450.dat", 0);
    bc.Vz.assignDBC(pointset, "gmtry/Venturi3D/walls_1450.dat", 0);
    bc.P.assignDBC(pointset, "gmtry/Venturi3D/outlet_1450.dat", 0);

    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 1.3;
    settings.kernel.approximation_method = "RBF";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "GE";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.dt = 0.00001;// use for RBF, MLS, WLS, NSPH
    settings.nt = 600000;
    settings.P_factor = 1000;
    settings.AV_factor = 0.1;
    settings.prnt_freq = 1000;
    pointset.AvgPntSpacing = 0.2; // essential to work

    cout << "Running solver..."<< endl;
    double h = 1;
    double mu = materials.ro[1]*materials.Nu[1];
    double Re = (materials.ro[1]*vel*h)/mu;
    cout << "Reynolds number = " << Re << endl;
    solveNavierStokes1Phase(pointset, materials, bc, settings, U, Vx, Vy, Vz, Q, Fx, Fy, Fz);
    return 0;
}

