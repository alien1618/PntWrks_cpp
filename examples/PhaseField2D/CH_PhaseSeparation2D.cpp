#include "src/pntwrks.h"

int main ()
{  
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D phase separation using the Cahn-Hilliard phase field model" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(64,64,0);
    RESOLUTION divs(64,64,1);
    POINTSET pointset(p, s, divs, "Q4");

    cout << "Assigning phases..."<< endl;
    double noise = 0.02;
    double c0 = 0.4;
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {pointset.phi[i] = c0 + noise*(0.5-random0to1());}
    pointset.printMeshVTK(pointset.phi, "phi", 0);
    pointset.printTXT(pointset.phi, "phi", 0);

    cout << "Setting interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.surface_tension.d0 = 0.5*pointset.AvgPntSpacing;
    interf_prmtrs.PFM.mobility = 1;

    cout << "Assigning solver settings..."<< endl;
    vector<double> V0(pointset.TotalPoints+1);
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
    settings.kernel.SPH = "GE";
    settings.dt = 0.01;
    settings.nt = 20000;
    settings.prnt_freq = 100;

    cout << "Assigning solver..."<< endl;
    clock_t start, end;
    start = clock();
    solveCahnHilliard(pointset, interf_prmtrs, settings, V0, V0, V0, V0);
    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by program is : " << time_taken << setprecision(10) << " seconds" << endl;
    return 0;
}
