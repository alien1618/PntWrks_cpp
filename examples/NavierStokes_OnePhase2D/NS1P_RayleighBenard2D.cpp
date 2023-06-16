#include "src.pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D Rayleigh-Benard convection" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(1,1,0);
    RESOLUTION divs(50,50,1);
    POINTSET pointset(p, s, divs, "T3");

    cout << "Assigning phase..."<< endl;
    POINT cntr(0.5,0);
    pointset.addSphere(cntr, 0.05);

    //POINT o3(.25,.25);
    //DOMAIN s3(0.5, .25);
    //pointset.addCube(o3, s3);

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.Nu[1] = 0.1;
    materials.ro[1] = 1;
    materials.D[1] = 1;
    materials.Nu[2] = 0.1;
    materials.ro[2] = 1;
    materials.D[2] = 1;

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
	vector<double> Fx = setVector(pointset.TotalPoints,0);
	vector<double> Fy = setVector(pointset.TotalPoints,0.001);
	vector<double> Fz = setVector(pointset.TotalPoints,0);
	
    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);

    bc.U.assignDBC(pointset.points, pointset.TotalPoints, "x", 1, 0);
    bc.U.assignDBC(pointset.points, pointset.TotalPoints, "x", 0, 0);
    bc.U.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);
    bc.U.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 1);

    //bc.P.assignDBC(pointset.points, pointset.TotalPoints, "y", 1, 0);

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
    settings.T_factor = 500; //use 5000 for 200X100
    settings.P_factor = 1000;
    settings.Thot = 1;
    settings.Tcold = 0;
    settings.AV_factor = 0.1;
    settings.Pr = 1.0;
    settings.Ra = 20000.0;
    settings.dt = 0.00005; //use 0.00001 for 200X100
    settings.nt = 1000000;
    settings.prnt_freq = 2000;
    settings.gravity.y = 0.001;

    cout << "Running solver..."<< endl;
    solveNavierStokes1Phase(pointset, materials, bc, settings, U, Vx, Vy, Vz, Q, Fx, Fy, Fz);
    return 0;
}
