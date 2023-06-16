#include "src/pntwrks.h"

int main ()
{  
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D dendritic solidification using the Allen-Cahn phase field model" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    
    cout << "Constructing pointset..."<< endl;
    POINT p(-3,-3,0);
    DOMAIN s(6,6,0);
    RESOLUTION divs(225,225,1);
    POINTSET pointset(p, s, divs, "Q4");

    cout << "Assigning phases..."<< endl;
    POINT cntr(0, 0, 0);
    pointset.addSphere(cntr, 0.025);
    pointset.computeVOF(pointset.AvgPntSpacing);
    pointset.printTXT(pointset.phi, "phi", 0);

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.D[1] = 1;
    materials.D[2] = 1;

    cout << "Initializing field variables..."<< endl;
	vector<double> U = setVector(pointset.TotalPoints,-1);
	vector<double> Vx = setVector(pointset.TotalPoints,0.0);
	vector<double> Vy = setVector(pointset.TotalPoints,0.0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
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
    interf_prmtrs.surface_tension.theta = 0.0;        //preferential growth angle

    interf_prmtrs.PFM.tau0 = 0.0002;//0.0003;
    interf_prmtrs.PFM.alpha = 0.9;         //smaller it is the smaller the solidification rate
    interf_prmtrs.PFM.gamma = 9.0;        //the higher it is the more non-constant the interface thickness. smaller it is the more the surface tension

    interf_prmtrs.thermal.fusion_latent_heat = 1.8;

    cout << "Assigning bounday conditions..."<< endl;
    BOUNDARY_CONDITIONS bc;

    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.simulation_type = "thermal";
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 2;
    settings.kernel.approximation_method = "RBF";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "GE";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "GE";
    settings.dt = 0.00005;
    settings.nt = 3000;
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    clock_t start, end;
    start = clock();
    solveAllenCahn(pointset, materials, interf_prmtrs, bc, settings, U, Vx, Vy, Vz);
    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by program is : " << time_taken << setprecision(10) << " seconds" << endl;    return 0;
}
