#include "src/pntwrks.h"

int main ()
{  
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D level set evolution under curvature driven flow " << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p(0,0);
    DOMAIN s(1,1);
    RESOLUTION divs(50,50);
    POINTSET pointset(p, s, divs, "Q4");

    cout << "Assigning phases..."<< endl;
    POINT cntr(0.3,0.3);
    DOMAIN s2(0.4,0.4);
    pointset.addCube(cntr, s2);
    pointset.computeVOF(2*pointset.AvgPntSpacing);

	cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interf_prmtrs;
    interf_prmtrs.sharpen_interface = false;
    interf_prmtrs.sharpen_interface_beta = 5*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_eps = 1*pointset.AvgPntSpacing;
    interf_prmtrs.sharpen_interface_dt =  10*0.00001;
    interf_prmtrs.sharpen_interface_nt =  5;
        
    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 1.5;
    settings.kernel.approximation_method = "MLS";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "LS";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";
    settings.dt = 0.001;
    settings.nt = 20000;
    settings.prnt_freq = 100;

    cout << "Running solver..."<< endl;
    solveCurvatureDrivenInterface(pointset, interf_prmtrs, settings);
    
    return 0;
}
