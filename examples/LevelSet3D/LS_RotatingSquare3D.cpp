#include "src/pntwrks.h"

int main ()
{  
    cout << "-----------------------------------------------------------------" << endl;
    cout << "3D level set evolution for a rotating square" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(-1,-1,-1);
    DOMAIN s(2,2,2);
    RESOLUTION divs(20,20,20);
    POINTSET pointset(p, s, divs, "H8");

    cout << "Assigning phases..."<< endl;
    POINT cntr(-0.25, -0.75, -0.25);
    DOMAIN s2(0.6, 0.6, 0.6);
    pointset.addCube(cntr, s2);
    pointset.computeVOF(pointset.AvgPntSpacing);
    pointset.printUnstructuredGridVTK(pointset.phi, "phi", 0);

    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.solver = "AC2";

    cout << "Initializing field variables..."<< endl;
    auto[Vx,Vy,Vz] = assignRotatingAdvectiveVelocity(pointset.points, pointset.TotalPoints, 1, 1, 0);
    vector<double> Vn_ext = setVector(pointset.TotalPoints,0);

    cout << "Assigning solver parameters..."<< endl;
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
    settings.dt = 0.0005;
    settings.nt = 20000;
    settings.prnt_freq = 100;

    cout << "Assigning solver..."<< endl;
    solveInterfaceEvolution(pointset, interf_prmtrs, settings, Vn_ext, Vx, Vy, Vz);
    return 0;
}
