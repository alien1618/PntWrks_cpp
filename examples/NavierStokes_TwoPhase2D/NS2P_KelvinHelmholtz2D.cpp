#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D Kelvin-Helmholtz instability" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0);
    DOMAIN d(1,0.5);
    RESOLUTION r(100,50);
    POINTSET pointset(p,d,r,"T3");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.Nu[1] = 1.48e-5;
    materials.ro[1] = 2;
    materials.Nu[2] = 1e-6;
    materials.ro[2] = 1;

    cout << "Assigning phases..."<< endl;
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {pointset.phi[i] = ((pointset.points[i].y - (0.01*cos(25*pointset.points[i].x+0.1)+0.0))-0.25);}
    pointset.computeVOF(pointset.AvgPntSpacing);

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,0);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	vector<double> Q = setVector(pointset.TotalPoints,0);
	vector<double> Fx = setVector(pointset.TotalPoints,0);
	vector<double> Fy = setVector(pointset.TotalPoints,-9.8);
	vector<double> Fz = setVector(pointset.TotalPoints,0);
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        if (pointset.phi[i] <= 0.5)
        {Vx[i] = 10;}
        else
        {Vx[i] = -10;} 
    }
    
    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.solver  = "VOF";
    interf_prmtrs.PFM.mobility = 100;
    interf_prmtrs.sharpen_interface = true;
    interf_prmtrs.sharpen_interface_beta = 10*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_eps = 2*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_dt =  50*0.00001;
    interf_prmtrs.sharpen_interface_nt =  5;

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, -10);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 0, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 0.5, 10);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 0.5, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 0.5, 0);

    bc.phi.total = 0;
    bc.phi.points.resize(1);
	bc.phi.values.resize(1);
	for (int i= 1; i<= pointset.TotalPoints; i++)
    {
        if (pointset.points[i].x <= 0.001 || pointset.points[i].x >= 0.999)
        {
            bc.phi.total++;
			bc.phi.points.push_back(1);
			bc.phi.values.push_back(1);
            bc.phi.points[bc.phi.total] = i;
            if (interf_prmtrs.solver == "CH2"||interf_prmtrs.solver == "AC2")
            {bc.phi.values[bc.phi.total] = 2*pointset.phi[i]-1;}
            else
            {bc.phi.values[bc.phi.total] = pointset.phi[i];}
        }
    }
    
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
    settings.kernel.upwind = false;
    settings.kernel.upwind_ratio = 0.75;
    settings.dt = 0.0002;
    settings.nt = 300000;
    settings.P_factor = 2000;//1/(settings.temporal.dt);
    settings.P_iterations = 1;
    settings.AV_factor = 0.1;
    settings.prnt_freq = 1000;
    settings.gravity.x = 0;
    settings.gravity.y = -9.8;
    settings.gravity.z = 0;

    cout << "Running solver..."<< endl;
    solveNavierStokes2Phase(pointset, materials, interf_prmtrs, bc, settings, U, Vx, Vy, Vz, Q, Fx, Fy, Fz);

    return 0;
}

