#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D rising bubble" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(1,1,0);
    RESOLUTION divs(50,50,1);
    POINTSET pointset(p, s, divs, "Q4");

    cout << "Assigning phases..."<< endl;
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {pointset.phi[i] = -(pointset.points[i].y -0.4);}
    POINT cntr(0.5,0.20,0);
    pointset.addCylinderZ(cntr, 0.15);
    pointset.computeVOF(0.5*pointset.AvgPntSpacing);

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.Nu[1] = 1000*1e-6;
    materials.ro[1] = 1000;
    materials.Nu[2] = 1000*1.48e-5;
    materials.ro[2] = 1;

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
    interf_prmtrs.solver  = "VOF";
    interf_prmtrs.interp_at_bndry = true;
    interf_prmtrs.sharpen_interface = true;
    interf_prmtrs.sharpen_interface_beta = 10*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_eps = 2*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_dt =  1000*0.00001;
    interf_prmtrs.sharpen_interface_nt =  5;

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
    settings.dt = 0.00001;
    settings.nt = 5000000;
    settings.P_factor = 1.0/settings.dt;
    settings.AV_factor = 0.1;
    settings.gravity.y = -9.8;
    settings.prnt_freq = 2000;

    cout << "Running solver..."<< endl;
    solveNavierStokes2Phase(pointset, materials, interf_prmtrs, bc, settings, U, Vx, Vy, Vz, Q, Fx, Fy, Fz);
    return 0;
}

