#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D hour glass" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "Constructing pointset..."<< endl;
    POINT p1(0.1,1.3);
    DOMAIN d1(0.8,0.4);
    RESOLUTION r1(60,30);
    POINTSET pointset1(p1,d1,r1);

    POINT p2(0.07,0.07);
    DOMAIN d2(0.865,0.2);
    RESOLUTION r2(65,15);   
    POINTSET pointset2(p2,d2,r2);

    POINTSET pointset;
    pointset.points.resize(pointset1.TotalPoints+pointset2.TotalPoints+1);
    pointset.points.resize(pointset1.TotalPoints+pointset2.TotalPoints+1);

    int m = 1;
    pointset.dim = 2;
    pointset.AvgPntSpacing = pointset1.AvgPntSpacing;
    pointset.TotalSurfaces = 0;
    pointset.TotalPoints = pointset1.TotalPoints+pointset2.TotalPoints;
    vector<double> phi(pointset.TotalPoints+1);
    for (int i = 1; i <= pointset1.TotalPoints; i++)
    {
        pointset.points[i].x = pointset1.points[i].x;
        pointset.points[i].y = pointset1.points[i].y;
        phi[i] = 0;
    }
    for (int i = pointset1.TotalPoints+1; i <= pointset.TotalPoints; i++)
    {
        pointset.points[i].x = pointset2.points[m].x;
        pointset.points[i].y = pointset2.points[m].y;
        phi[i] = 1;
        m++;
    }

    cout << "Constructing boundaries..."<< endl;
    POINTSET bndry("geometry/SPH_HourGlass/mesh_500.dat");

    cout << "Constructing background grid..."<< endl;
    POINT p3(-0.25,-0.125);
    DOMAIN d3(1.5, 2);
    RESOLUTION r3(10, 15); 
    POINTSET bg_grid(p3,d3,r3, "Q4");

    cout << "Assigning material properties..."<< endl;
    MATERIALS materials(2);
    materials.Nu[1] = 0.01;//1e-6;
    materials.ro[1] = 100;
    materials.Nu[2] = 0.01;//1e-6;
    materials.ro[2] = 1000;

    cout << "Assigning interface parameters..."<< endl;
    INTERFACE_PRMTRS interface_parameters;
    interface_parameters.solver = "NONE";
    interface_parameters.surface_tension.surface_tension_co = 0.0;

    cout << "Assigning boundary conditions..."<< endl;
    int s = 1;
    vector<int> bndry_pnt_num(pointset.TotalPoints+bndry.TotalPoints+1);
    vector<double> values(pointset.TotalPoints+bndry.TotalPoints+1);
    for (int i = pointset.TotalPoints+1; i <= pointset.TotalPoints+bndry.TotalPoints; i++)
    {
        bndry_pnt_num[s] = i;
        s++;
    }

    BOUNDARY_CONDITIONS bc;
    bc.Vx.total = bndry.TotalPoints;
    bc.Vx.points = bndry_pnt_num;
    bc.Vx.values = values;
    
    bc.Vy.total = bndry.TotalPoints;
    bc.Vy.points = bndry_pnt_num;
    bc.Vy.values = values;
    
    bc.Vz.total = bndry.TotalPoints;
    bc.Vz.points = bndry_pnt_num;
    bc.Vz.values = values;
    
    cout << "Assigning solver settings..."<< endl;
    SOLVER_SETTINGS settings;
    settings.dt = 0.00005;
    settings.nt = 50000;
    settings.kernel.radius_ratio = 1.5;
    settings.kernel.support_method = "BF";
    settings.kernel.approximation_method = "SPH";
    settings.P_factor = 200;
    settings.P_iterations = 1;
    settings.AV_factor = 0.01;
    settings.lagrangian = true;
    settings.particle_shifting = true;
    settings.particle_shifting_factor = 1;
    settings.particle_shifting_surf = 1;
    settings.gravity.x = 0;
    settings.gravity.y = -9.8;
    settings.gravity.z = 0;
    settings.prnt_freq = 100;

    cout << "Running SPH Solver..."<< endl;
    solveNavierStokesLagrangian(pointset, bndry, bg_grid, materials, interface_parameters, bc, settings, phi);
    
    return 0;
}
