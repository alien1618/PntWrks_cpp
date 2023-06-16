#include "src/pntwrks.h"

int main ()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D dendritic solidification with flow using the Allen-Cahn" << endl;
    cout << "phase field model" << endl;
    cout << "-----------------------------------------------------------------" << endl;  

    cout << "Constructing pointset..."<< endl;
    POINT p(-3,-3,0);
    DOMAIN s(6,6,0);
    RESOLUTION divs(100,100,1);
    POINTSET pointset(p, s, divs, "Q4");

    cout << "Assigning phases..."<< endl;
    POINT cntr(0,0);
    pointset.addSphere(cntr, 0.05);
    pointset.computeVOF(pointset.AvgPntSpacing);
    pointset.printMeshVTK(pointset.phi,"phi", 0);
    
    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.D[1] = 1;
    materials.Nu[1] = 0.1;
    materials.ro[1] = 1;
    materials.D[2] = 1;
    materials.Nu[2] = 1;
    materials.ro[2] = 1;

    cout << "Initializing field variables..."<< endl;
    double vel = 3;
	vector<double> U = setVector(pointset.TotalPoints,-1);
	vector<double> Vx = setVector(pointset.TotalPoints,vel);
	vector<double> Vy = setVector(pointset.TotalPoints,0.0);
	vector<double> Vz = setVector(pointset.TotalPoints,0.0);
	for (int i = 1; i<= pointset.TotalPoints; i++)
	{
		if (pointset.phi[i] >= 0.5)
		{U[i] = 0;}
	}
	
    cout << "Assigning interface properties..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.thermal.T_eq = 0;
    interf_prmtrs.surface_tension.b = 4.0;
    interf_prmtrs.surface_tension.A = 0.018;        //the smaller it is the more dendrite splitting
    interf_prmtrs.surface_tension.d0 = 0.01;        //directly related to the interface thickness
    interf_prmtrs.surface_tension.theta = 0;//0.785;        //preferential growth angle

    interf_prmtrs.PFM.tau0 = 0.0003;
    interf_prmtrs.PFM.alpha = 0.9;         //smaller it is the smaller the solidification rate
    interf_prmtrs.PFM.gamma = 10.0;        //the higher it is the more non-constant the interface thickness. smaller it is the more the surface tension

    interf_prmtrs.thermal.fusion_latent_heat = 1.8;

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;
    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "x", -3, vel);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "x", -3, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "x", -3, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", 3, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", 3, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", 3, 0);

    bc.Vx.assignDBC(pointset.points, pointset.TotalPoints, "y", -3, 0);
    bc.Vy.assignDBC(pointset.points, pointset.TotalPoints, "y", -3, 0);
    bc.Vz.assignDBC(pointset.points, pointset.TotalPoints, "y", -3, 0);

    //bc.P.assignDBC(pointset.points, pointset.TotalPoints, "x", 3, 0);

    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.simulation_type = "thermal";
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 2;
    settings.kernel.approximation_method = "WLS";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "GE";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.P_factor = 1000;
    settings.dt = 0.00005;
    settings.nt = 5000;
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    solveAllenCahn(pointset, materials, interf_prmtrs, bc, settings, U, Vx, Vy, Vz);    return 0;
}
