#include "src/pntwrks.h"

int main ()
{  
    cout << "-----------------------------------------------------------------" << endl;
    cout << "3D level set evolution under shear flow " << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0,0);
    DOMAIN s(1,1,1);
    RESOLUTION divs(20,20,20);
    POINTSET pointset(p, s, divs, "H8");

    cout << "Assigning phases..."<< endl;
    POINT cntr(0.5, 0.5, 0.5);
    pointset.addSphere(cntr, 0.25);
    pointset.computeVOF(pointset.AvgPntSpacing);
    pointset.printUnstructuredGridVTK(pointset.phi, "phi", 0);

    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.solver = "CH";
    interf_prmtrs.sharpen_interface = true;
    interf_prmtrs.sharpen_interface_beta = 10*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_eps = 2*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_dt =  10*0.00001;
    interf_prmtrs.sharpen_interface_nt =  5;

    cout << "Initializing field variables..."<< endl;
    auto[Vx,Vy,Vz] = assignExtremeDeformationVelocity(pointset.points, pointset.TotalPoints);
    vector<double> Vn_ext = setVector(pointset.TotalPoints,0);

    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
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
    settings.dt = 0.00001;
    settings.nt = 20000;
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    solveInterfaceEvolution(pointset, interf_prmtrs, settings, Vn_ext, Vx, Vy, Vz);
    
    return 0;
}
