#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D dam break" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(1,1,0);
    RESOLUTION divs(30,30,1);
    POINTSET pointset(p, s, divs, "Q4");

    cout << "Assigning phases..."<< endl;
    POINT o(-0.05,-0.05);
    DOMAIN s2(0.45, 0.45);
    pointset.addCube(o, s2);
    pointset.computeVOF(0.5*pointset.AvgPntSpacing);
    pointset.printMeshVTK(pointset.phi, "phi", 0);

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.Nu[1] = 1000*1.48e-5;
    materials.ro[1] = 1;
    materials.Nu[2] = 1000*1e-6;
    materials.ro[2] = 1000;

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
	vector<double> Fx = setVector(pointset.TotalPoints,0);
	vector<double> Fy = setVector(pointset.TotalPoints,-9.8);
	vector<double> Fz = setVector(pointset.TotalPoints,0);
	
    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.solver  = "AC";
    interf_prmtrs.PFM.mobility = 5;
    /*
    interf_prmtrs.sharpen_interface = true;
    interf_prmtrs.sharpen_interface_beta = 10*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_eps = 2*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_dt =  1000*0.00001*0.25;
    interf_prmtrs.sharpen_interface_nt =  5;
    */

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc; //must be all wall boundary conditions to get correct velocity profile
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);

    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 1.5;
    settings.kernel.approximation_method = "SPH";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "LS";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";

    cout << "Smoothing out the interface corners and converting to vof..."<< endl;
    settings.dt = 0.0001;
    settings.nt = 300;
    settings.prnt_freq = 100;
    pointset = solveCurvatureDrivenInterface(pointset, interf_prmtrs, settings);

    cout << "Running solver..."<< endl;
    settings.dt =  0.00001*0.25;
    settings.nt = 5000000;
    settings.P_factor = 1000;//1.0/settings.dt;
    settings.gravity.y = -9.8;
    settings.prnt_freq = 2000;
    solveNavierStokes2Phase(pointset, materials, interf_prmtrs, bc, settings, U, Vx, Vy, Vz, Q, Fx, Fy, Fz);
    return 0;
}

