#include "../pntwrks.h"

// solver works with 0<=phi<=1 given as input

void solveInterfaceEvolution(POINTSET pointset, INTERFACE_PRMTRS inter_param, SOLVER_SETTINGS settings, vector<double> Vn_ext, vector<double> Vx, vector<double> Vy, vector<double> Vz)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..." << endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}

    int method = 1;
    if (inter_param.solver == "VOF")
    {method = 1;}
    else if(inter_param.solver == "AC")
    {method = 2;}
    else if(inter_param.solver == "AC2")
    {method = 3;}
    else if(inter_param.solver == "CH")
    {method = 4;}
    else if(inter_param.solver == "CH2")
    {method = 5;}
    else
    {
        cout << "inter_param.solver = " << inter_param.solver << endl;
        cout << "ERROR: inter_param.solver is undefined..." << endl;
        exit(0);
    }

    int c = 1;
    int iters = 5;
    double kernel_radius = settings.kernel.radius_ratio*pointset.AvgPntSpacing;
    double W = 0.5*pointset.AvgPntSpacing;
    double M = inter_param.PFM.mobility*W*W;
    double tot_vof = 0;
    if (inter_param.solver == "AC" || inter_param.solver == "CH" || inter_param.solver == "VOF")
    {
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {tot_vof = tot_vof+pointset.phi[i];}
    }
    if (inter_param.solver == "AC2" || inter_param.solver == "CH2" )
    {
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {pointset.phi[i] = 2*pointset.phi[i]-1;}
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {
            double vof0 = 0.5*(pointset.phi[i]+1);
            tot_vof = tot_vof+vof0;
        }
    }

    cout << "Printing field variables at time 0..." << endl;
    pointset.printMeshVTK(pointset.phi, "phi", 0);

    cout << "Computing kernel approximation functions..." << endl;
    vector<KERNEL> kernels = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);

    cout << "Computing the transient solution..." << endl;
    for (int t = 1; t <= settings.nt; t++)
    {
        switch(method)
        {
            case 1:
            {pointset.phi = computeVolumeOfFluid(kernels, pointset.phi, pointset.TotalPoints, Vx, Vy, Vz, Vn_ext, iters, settings.dt);}break;
            case 2://sun-beckerman
            {pointset.phi = computeAllenCahn(W, M, kernels, pointset.phi, pointset.dim, pointset.points, pointset.TotalPoints,  Vx, Vy, Vz, Vn_ext, iters, settings.dt);}break;
            case 3:
            {pointset.phi = computeAllenCahn2(W, M, kernels, pointset.phi, pointset.dim, pointset.points, pointset.TotalPoints,  Vx, Vy, Vz, Vn_ext, iters, settings.dt);}break;
            case 4:
            {pointset.phi = computeCahnHilliard(W, M, kernels, pointset.phi, pointset.dim, pointset.points, pointset.TotalPoints,  Vx, Vy, Vz, Vn_ext, iters, settings.dt);}break;
            case 5:
            {pointset.phi = computeCahnHilliard2(W, M, kernels, pointset.phi, pointset.dim, pointset.points, pointset.TotalPoints,  Vx, Vy, Vz, Vn_ext, iters, settings.dt);}break;
        }
        if (inter_param.sharpen_interface == true)
        {
            if(method == 1 || method == 2 || method == 4)
            {pointset.phi = computeSharpenInterface(true, inter_param.sharpen_interface_beta, inter_param.sharpen_interface_eps, pointset.phi, pointset.TotalPoints, kernels, pointset.AvgPntSpacing, inter_param.sharpen_interface_dt, inter_param.sharpen_interface_nt);}
            else
            {pointset.phi = computeSharpenInterface(false, inter_param.sharpen_interface_beta, inter_param.sharpen_interface_eps, pointset.phi, pointset.TotalPoints, kernels, pointset.AvgPntSpacing, inter_param.sharpen_interface_dt, inter_param.sharpen_interface_nt);}
        }
        if(method == 1 || method == 2 || method == 4)
        {
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {
                if (pointset.phi[i] < 0)
                {pointset.phi[i] = 0.0;}
                else if (pointset.phi[i] > 1)
                {pointset.phi[i] = 1.0;}
            }
        }
        else
        {
            for (int i = 1; i <= pointset.TotalPoints; i++)
            {
                if (pointset.phi[i] < -1.0)
                {pointset.phi[i] = -1.0;}
                else if (pointset.phi[i] > 1.0)
                {pointset.phi[i] = 1.0;}
            }
        }
        if (t/settings.prnt_freq == c)
        {
            c++;
            cout << "Computation of time step " << t << " of " << settings.nt  << " complete..." << endl;
            if (pointset.TotalSurfaces > 0 || pointset.TotalVolumes > 0)
            {
                if (pointset.dim == 2)
                {
                    pointset.printMeshVTK(pointset.phi, "phi", t);
                    pointset.printTXT(pointset.phi, "phi", t);
                }
                else if (pointset.dim == 3)
                {pointset.printUnstructuredGridVTK(pointset.phi, "phi", t);}
            }
            else
            {pointset.printTXT(pointset.phi, "phi", t);}
            double tot_vof_new = 0;
            if (method == 1 || method == 2 || method == 4)
            {
                for (int i = 1; i <= pointset.TotalPoints; i++)
                {tot_vof_new = tot_vof_new+pointset.phi[i];}
            }
            else
            {
                for (int i = 1; i <= pointset.TotalPoints; i++)
                {
                    double vof = (pointset.phi[i]+1)/2;
                    tot_vof_new = tot_vof_new+vof;
                }
            }
            double vl = 100*(tot_vof-tot_vof_new)/tot_vof;
            cout << "Accumulated volume loss % = " << vl << endl;
        }
    }
    
    cout << "SOLVER COMPLETE" << endl;
}
POINTSET solveCurvatureDrivenInterface(POINTSET pointset, INTERFACE_PRMTRS inter_param, SOLVER_SETTINGS settings)
{
    cout << "SOLVER START..." << endl;
    cout << "Initializing variables..." << endl;
    pltctrl(pointset.points, pointset.TotalPoints, settings.nt, settings.prnt_freq);
    if (pointset.TotalSurfaces>0)
    {printElements(pointset.surfaces, pointset.SurfaceNds, pointset.TotalSurfaces, "elements");}

    int c = 1;
    double b = 0.01;
    double kernel_radius = pointset.AvgPntSpacing*settings.kernel.radius_ratio;
    vector<double> Vx = setVector(pointset.TotalPoints,0);
    vector<double> Vy = setVector(pointset.TotalPoints,0);
    vector<double> Vz = setVector(pointset.TotalPoints,0);
    vector<double> curvature;

    cout << "Printing field variables at time 0..." << endl;
    pointset.printMeshVTK(pointset.phi, "phi", 0);

    cout << "Computing kernel approximation functions..." << endl;
    vector<KERNEL> kernels(pointset.TotalPoints+1);
    kernels = computeKernels(pointset.points, pointset.phi, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);
    
    cout << "Computing transient solution..."<< endl;
    for (int t = 1; t <= settings.nt; t++)
    {
        //cout << "Computing curvature..." << endl;
        curvature = computeCurvature(pointset.phi, pointset.TotalPoints, kernels, pointset.dim);

        //cout << "Solving the level set equation..." << endl;
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {curvature[i] = curvature[i]*b;}
        pointset.phi = computeVolumeOfFluid(kernels, pointset.phi, pointset.TotalPoints, Vx, Vy, Vz, curvature, 5, settings.dt);

        //cout << "Sharpening the diffuse interface..." << endl;
        if (inter_param.sharpen_interface == true)
        {pointset.phi = computeSharpenInterface(true, inter_param.sharpen_interface_beta, inter_param.sharpen_interface_eps, pointset.phi, pointset.TotalPoints, kernels, pointset.AvgPntSpacing, inter_param.sharpen_interface_dt, inter_param.sharpen_interface_nt);}
        for (int i = 1; i <= pointset.TotalPoints; i++)
        {
            if (pointset.phi[i] <= 0.0001)
            {pointset.phi[i] = 0;}
            if (pointset.phi[i] >= 0.999)
            {pointset.phi[i] = 1;}
        }

        if (t/settings.prnt_freq == c)
        {
            c++;
            cout << "Computation of time step " << t << " of " << settings.nt  << " complete..." << endl;
            switch(pointset.dim)
            {
                case 2:
                {
                    pointset.printMeshVTK(pointset.phi, "phi", t);
                    pointset.printMeshVTK(curvature, "curvature", t);
                }break;
                case 3:
                {pointset.printUnstructuredGridVTK(pointset.phi, "phi", t);}break;
            }
        }
    }
    cout << "SOLVER COMPLETE" << endl;
    return pointset;
}
