#include "src/pntwrks.h"

int main ()
{  
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D level set evolution under vortex flow " << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0);
    DOMAIN s(1,1);
    RESOLUTION divs(75,75);
    POINTSET pointset(p, s, divs, "Q4");
    pointset.AvgPntSpacing = 1.0/75;

    cout << "Assigning phases..."<< endl;
    POINT cntr(0.5,0.25);
    pointset.addCylinderZ(cntr, 0.2);
    pointset.computeVOF(pointset.AvgPntSpacing);
    pointset.printUnstructuredGridVTK(pointset.phi, "phi", 0);

    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.solver = "VOF";
    interf_prmtrs.sharpen_interface = false;
    interf_prmtrs.sharpen_interface_beta = 10*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_eps = 2*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_dt =  10*0.00001;
    interf_prmtrs.sharpen_interface_nt =  5;

    cout << "Initializing field variables..."<< endl;
    auto[Vx,Vy,Vz] = assignSingleVortexVelocity(pointset.points, pointset.TotalPoints);
    vector<double> Vn_ext = setVector(pointset.TotalPoints,0);

    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 2;
    settings.kernel.approximation_method = "SPH";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "S4";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "GE";
    settings.dt = 0.0002;
    settings.nt = 20000;
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    clock_t start, end;
	start = clock();
    solveInterfaceEvolution(pointset, interf_prmtrs, settings, Vn_ext, Vx, Vy, Vz);
    end = clock();
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << "Time taken by program is : " << time_taken << setprecision(10) << " seconds" << endl;

    return 0;
}
