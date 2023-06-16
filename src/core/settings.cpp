#include "../pntwrks.h"

MATERIALS::MATERIALS()
{
    total  = 1;
    D.resize(2);
	Nu.resize(2);
	ro.resize(2);
    D[1] = 0;
    Nu[1] = 0;
    ro[1] = 0;
}

MATERIALS::MATERIALS(int num)
{
	total = num;
	D.resize(num+1);
	Nu.resize(num+1);
	ro.resize(num+1);
    for (int i = 1; i <= num; i++)
    {
        D[i] = 0;
        Nu[i] = 0;
        ro[i] = 0;
    }
}

MATERIALS::~MATERIALS()
{
    //destructor
}
SOLVER_SETTINGS::SOLVER_SETTINGS() //constructor
{
    //DEFAULT SOLVER PARAMETERS
    kernel.recompute_approximants = true;
    kernel.support_method = "BF";
    kernel.radius_ratio = 1.5;
    kernel.approximation_method = "RBF";
    kernel.approximation_order = 1;
    kernel.grad_order = 1;
    kernel.RBF = "GE";
    kernel.RBF_alpha = 0.1;
    kernel.WLS = "GE";
    kernel.MLS = "S4";
    kernel.SPH = "WC4";
    kernel.upwind = false;
    kernel.upwind_ratio = 1;
    lagrangian = false;
    T_factor = 500;
    P_factor = 1000;
    P_iterations = 1;
    AV_factor = 0.1;
    Pr = 1.0;
    Ra = 20000.0;
    gravity.x = 0;
    gravity.y = 0;
    gravity.z = 0;
    prnt_freq = 1;
}
