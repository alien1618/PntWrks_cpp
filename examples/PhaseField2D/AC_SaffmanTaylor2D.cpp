#include "src/pntwrks.h"

int main ()
{  
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D simulation of the Saffman-Taylor instability using the" << endl;
    cout << "Allen-Cahn phase field model" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    
    cout << "Constructing pointset..."<< endl;
    POINT p(-3,-3,0);
    DOMAIN s(6,6,0);
    RESOLUTION divs(100,100,1);
    POINTSET pointset(p, s, divs, "T3");

    cout << "Assigning phases..."<< endl;
    POINT cntr(0, 0, 0);
    pointset.addSphere(cntr, 0.1);
    pointset.computeVOF(0.5*pointset.AvgPntSpacing);

    cout << "Assigning material properties..."<< endl;
    double P_inf = -1;
    double P_vap = 0;
    double P_eq = 0;
    double Nu_vap= 1;
    double Nu_liq= 1;

    MATERIALS materials(2);
    materials.D[1] = Nu_liq;
    materials.D[2] = Nu_vap;

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,P_inf);
	vector<double> Vx = setVector(pointset.TotalPoints,0);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	for (int i = 1; i<= pointset.TotalPoints; i++)
	{
		if (pointset.phi[i] >= 0.5)
		{U[i] = P_vap;}
	}
	
    cout << "Assigning interface properties..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.thermal.T_eq = P_eq;
    interf_prmtrs.surface_tension.b = 8.0;
    interf_prmtrs.surface_tension.A = 0.01;        //the smaller it is the more dendrite splitting
    interf_prmtrs.surface_tension.d0 = 0.01;        //directly related to the interface thickness
    interf_prmtrs.surface_tension.theta = 0.2;        //preferential growth angle

    interf_prmtrs.PFM.tau0 = 0.0003;
    interf_prmtrs.PFM.alpha = 0.9;         //smaller it is the smaller the solidification rate
    interf_prmtrs.PFM.gamma = 10.0;        //the higher it is the more non-constant the interface thickness. smaller it is the more the surface tension

    double cell_gap_size = 4.64;
    interf_prmtrs.thermal.fusion_latent_heat = pow(cell_gap_size,2)/12;

    cout << "Assigning bounday conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;

    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.simulation_type = "hele_shaw";
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 2.5;
    settings.kernel.approximation_method = "SPH";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "LS";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.dt = 0.00005;
    settings.nt = 10000;
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    solveAllenCahn(pointset, materials, interf_prmtrs, bc, settings, U, Vx, Vy, Vz);
    
    return 0;
}
