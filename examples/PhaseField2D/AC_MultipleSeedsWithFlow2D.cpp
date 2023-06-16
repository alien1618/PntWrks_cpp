#include "src/pntwrks.h"

int main ()
{  
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D dendritic solidification of multiple seeds with flow using" << endl;
    cout << "using the Allen-Cahn phase field model" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(-3,-3);
    DOMAIN s(6,6);
    RESOLUTION divs(225,225);
    POINTSET pointset(p, s, divs, "Q4");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.D[1] = 1;
    materials.Nu[1] = 0.1;
    materials.ro[1] = 1;
    materials.D[2] = 1;
    materials.Nu[2] = 1;
    materials.ro[2] = 1;

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
    interf_prmtrs.seeds(4);
    interf_prmtrs.seed_theta[1] = 0;
    interf_prmtrs.seed_theta[2] = 0.5;
    interf_prmtrs.seed_theta[3] = 2;
    interf_prmtrs.seed_theta[4] = 1.5;

    cout << "Assigning phases..."<< endl;
    int seeds = 4;
    vector<vector<double> > phi(pointset.TotalPoints+1, vector<double> (seeds + 1));
    POINT cntr(-1.5, -1.5, 0);
    pointset.addSphere(cntr, 0.1);
    for(int i = 1; i <= pointset.TotalPoints; i++)
    {
        phi[i][1] = 1-0.5*(tanh(pointset.phi[i]/interf_prmtrs.surface_tension.d0)+1);
        pointset.phi[i] = 1;
    }
    POINT cntr1(1.5, -1.5, 0);
    pointset.addSphere(cntr1, 0.1);
    for(int i = 1; i <= pointset.TotalPoints; i++)
    {
        phi[i][2] = 1-0.5*(tanh(pointset.phi[i]/interf_prmtrs.surface_tension.d0)+1);
        pointset.phi[i] = 1;
    }
    POINT cntr2(-1.5, 1.5, 0);
    pointset.addSphere(cntr2, 0.1);
    for(int i = 1; i <= pointset.TotalPoints; i++)
    {
        phi[i][3] = 1-0.5*(tanh(pointset.phi[i]/interf_prmtrs.surface_tension.d0)+1);
        pointset.phi[i] = 1;
    }

    POINT cntr3(1.5, 1.5, 0);
    pointset.addSphere(cntr3, 0.1);
    for(int i = 1; i <= pointset.TotalPoints; i++)
    {
        phi[i][4] = 1-0.5*(tanh(pointset.phi[i]/interf_prmtrs.surface_tension.d0)+1);
        pointset.phi[i] = 1;
    }

    vector<double> phi_tot(pointset.TotalPoints+1);
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        for(int s = 1;  s <= seeds; s++)
        {
            if (phi[i][s] > 0)
            {phi_tot[i] = phi[i][s];}
        }
    }
    pointset.printMeshVTK(phi_tot, "phi", 0);

    cout << "Initializing field variables..."<< endl;
    double vel = 3;
	vector<double> U = setVector(pointset.TotalPoints,-1);
	vector<double> Vx = setVector(pointset.TotalPoints,vel);
	vector<double> Vy = setVector(pointset.TotalPoints,0);
	vector<double> Vz = setVector(pointset.TotalPoints,0);
	for (int i = 1; i<= pointset.TotalPoints; i++)
	{
		if (phi_tot[i] >= 0.5)
		{U[i] = 0;}
	}
	
    cout << "Assigning bounday conditions..."<< endl;
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

    cout << "Assigning simulation parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.simulation_type = "thermal";
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 1.5;
    settings.kernel.approximation_method = "WLS";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "GE";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.P_factor = 1000;
    settings.dt = 0.00005;
    settings.nt = 3000;
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    solveAllenCahnMultipleSeeds(pointset, materials, interf_prmtrs, bc, settings, phi, U, Vx, Vy, Vz);
    return 0;
}
