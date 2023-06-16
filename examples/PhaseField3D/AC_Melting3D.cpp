#include "src/pntwrks.h"

int main ()
{  
    cout << "-----------------------------------------------------------------" << endl;
    cout << "3D melting using the Allen-Cahn phase field model" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(-3,-3);
    DOMAIN s(6,6);
    RESOLUTION divs(200,200);
    POINTSET pointset(p, s, divs, "T3");
    //POINTSET pointset("geometry/CircleR3.0/mesh_50000.dat");
    //POINTSET pointset("geometry/ClosedSpline_v3/mesh_50000.dat");

    cout << "Assigning phases..."<< endl;
    POINT cntr(-1.5, -1.5);
    DOMAIN s2(3,3);
    pointset.addCube(cntr, s2);
    pointset.computeVOF(pointset.AvgPntSpacing);

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.D[1] = 1;
    materials.D[2] = 1;

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,1);
	vector<double> V0 = setVector(pointset.TotalPoints,0);
	for (int i = 1; i<= pointset.TotalPoints; i++)
	{
		if (pointset.phi[i] >= 0.5)
		{U[i] = -1;}
	}
	
    cout << "Assigning interface properties..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.thermal.T_eq = 0;
    interf_prmtrs.surface_tension.b = 4.0;
    interf_prmtrs.surface_tension.A = 0.02;        //the smaller it is the more dendrite splitting
    interf_prmtrs.surface_tension.d0 = 0.01;        //directly related to the interface thickness
    interf_prmtrs.surface_tension.theta = 0.2;        //preferential growth angle

    interf_prmtrs.PFM.tau0 = 0.0003;
    interf_prmtrs.PFM.alpha = 0.9;         //smaller it is the smaller the solidification rate
    interf_prmtrs.PFM.gamma = 10.0;        //the higher it is the more non-constant the interface thickness. smaller it is the more the surface tension

    interf_prmtrs.thermal.fusion_latent_heat = 1.8;
    
    cout << "Assigning bounday conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;
    //bc.U.assign_SurfaceDBC("geometry/CircleR3.0/bndr_50000.dat", -1);
    //bc.U.assign_SurfaceDBC("geometry/ClosedSpline_v3/walls_50000.dat", -1);

    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.simulation_type = "thermal";
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 1.5;
    settings.kernel.approximation_method = "RBF";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "LS";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.dt = 0.0001;
    settings.nt = 10000;
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    solveAllenCahn(pointset, materials, interf_prmtrs, bc, settings, U, V0, V0, V0);
    return 0;
}
