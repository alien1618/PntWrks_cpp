#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D reaction-diffusion" << endl;
    cout << "-----------------------------------------------------------------" << endl;

    cout << "Constructing pointset..."<< endl;
    POINT p(-1,-1,0);
    DOMAIN size(2,2,0);
    RESOLUTION divs(100, 100, 1);
    POINTSET pointset(p, size, divs, "T3");

    cout << "Assigning phases..."<< endl;
    POINT cntr(0,0,0);
    pointset.addSphere(cntr, 0.025);
    /*
    POINT cntr1(-0.5,-0.5);
    pointset.addSphere(cntr1, 0.01);
    POINT cntr2(0.5,-0.5);
    pointset.addSphere(cntr2, 0.01);
    POINT cntr3(0.5,0.5);
    pointset.addSphere(cntr3, 0.01);
    POINT cntr4(-0.5,0.5);
    pointset.addSphere(cntr4, 0.01);
    */
    pointset.printMeshVTK(pointset.phi,"phi", 0);
    
    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.D[1] = 0.2*(pointset.AvgPntSpacing*pointset.AvgPntSpacing);
    materials.D[2] = 0.1*(pointset.AvgPntSpacing*pointset.AvgPntSpacing);

	cout << "Initializing field variables..."<< endl;
	vector<double> U_A = setVector(pointset.TotalPoints,1);
	vector<double> U_B = setVector(pointset.TotalPoints,0);

    cout << "Assigning boundary conditions..."<< endl;
    BOUNDARY_CONDITION bc;

    cout << "Assigning solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 1.5;
    settings.kernel.approximation_method = "RBF";
    settings.kernel.approximation_order = 1;
    settings.kernel.WLS = "GE";
    settings.kernel.MLS = "S4";
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.SPH = "WC4";
    settings.dt = 0.1;
    settings.nt = 50000;
    settings.prnt_freq = 100;

    double feed=0;
    double kill=0;
    int parameters = 7;

    if (parameters == 1)
    {
        feed = 0.0545;
        kill = 0.062;
    }
    else if (parameters == 2)
    {
        feed = 0.055;
        kill = 0.062;
    }
    else if (parameters == 3)
    {
        feed = 0.018;
        kill = 0.051;
    }
    else if (parameters == 4)
    {
        feed = 0.026;
        kill = 0.053;
    }
    else if (parameters == 5)
    {
        feed = 0.055;
        kill = 0.063;
    }
    else if (parameters == 6)
    {
        feed = 0.037;
        kill = 0.064;
    }
    else if (parameters == 7)//solitons
    {
        feed = 0.03;
        kill = 0.062;
    }
    else if (parameters == 8) //pulsating solitons
    {
        feed = 0.025;
        kill = 0.06;
    }
    else if (parameters == 9) //worms
    {
        feed = 0.078;
        kill = 0.061;
    }
    else if (parameters == 10) //maze
    {
        feed = 0.029;
        kill = 0.057;
    }
    else if (parameters == 11) // holes
    {
        feed = 0.039;
        kill = 0.058;
    }
    else if (parameters == 12) //chaos
    {
        feed = 0.026;
        kill = 0.051;
    }
    else if (parameters == 13) //chaos and holes
    {
        feed = 0.034;
        kill = 0.056;
    }
    else if (parameters == 14) //moving spots
    {
        feed = 0.014;
        kill = 0.054;
    }
    else if (parameters == 14) //spots and loops
    {
        feed = 0.018;
        kill = 0.051;
    }
    else if (parameters == 15) //wave
    {
        feed = 0.014;
        kill = 0.045;
    }
    else if (parameters == 16) //u-skate world
    {
        feed = 0.062;
        kill = 0.061;
    }
    else if (parameters == 17)
    {
        feed = 0.037;
        kill = 0.060;
    }

    cout << "Running solver..."<< endl;
    solveReactionDiffusion(pointset, materials, bc, settings,  feed, kill, U_A, U_B);
}