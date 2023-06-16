#include "src/pntwrks.h"

int main ()
{  
    cout << "-----------------------------------------------------------------" << endl;
    cout << "3D phase separation using the Cahn-Hilliard phase field model" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(100,100,2);
    RESOLUTION divs(100,100,2);
    POINTSET pointset(p, s, divs, "H8");

    cout << "Assigning phases..."<< endl;
    double noise = 0.02;
    double c0 = 0.4;
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {pointset.phi[i] = c0 + noise*(0.5-random0to1());}
    pointset.printMeshVTK(pointset.phi, "phi", 0);

    cout << "Setting interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.surface_tension.d0 = 1*pointset.AvgPntSpacing;
    interf_prmtrs.PFM.mobility = 1;

    cout << "Assigning solver settings..."<< endl;
    vector<double> V0 = setVector(pointset.TotalPoints, 0);
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 2;
    settings.kernel.approximation_method = "RBF";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "LS";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.dt = 0.005;
    settings.nt = 200000;
    settings.prnt_freq = 100;

    cout << "Assigning solver..."<< endl;
    solveCahnHilliard(pointset, interf_prmtrs, settings, V0, V0, V0, V0);
    return 0;
}
