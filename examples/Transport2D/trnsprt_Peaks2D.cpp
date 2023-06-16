#include "src/pntwrks.h"

int main()
{
    cout << "-----------------------------------------------------------------" << endl;
    cout << "2D peaks function approximation" << endl;
    cout << "-----------------------------------------------------------------" << endl;

    cout << "Constructing pointset..."<< endl;
    POINT p(-3, -3);
    DOMAIN size(6,6);
    RESOLUTION divs(100, 100);
    POINTSET pointset(p, size, divs, "T3");

    cout << "Calculating exact solution of the peaks function and its derivatives..."<< endl;
    auto[f,dfdx,dfdy,dfdxx,dfdyy] = computePeaksFunction(pointset.points, pointset.TotalPoints);

    cout << "Setting solver parameters..."<< endl;
    SOLVER_SETTINGS settings;
    settings.kernel.recompute_approximants = true;
    settings.kernel.support_method = "BF";
    settings.kernel.radius_ratio = 2;
    settings.kernel.approximation_method = "RBF";
    settings.kernel.approximation_order = 1;
    settings.kernel.RBF = "GE";
    settings.kernel.RBF_alpha = 0.1;
    settings.kernel.WLS = "GE";
    settings.kernel.MLS = "S4";
    settings.kernel.SPH = "WC4";

    cout << "Computing meshfree interpolants..." << endl;
    vector<double> zero(pointset.TotalPoints+1);
    double kernel_radius = settings.kernel.radius_ratio*pointset.AvgPntSpacing;
    double t_start = omp_get_wtime();
    vector<KERNEL> kernels = computeKernels(pointset.points, zero, pointset.TotalPoints, pointset.dim, settings.kernel, kernel_radius);
    double t_end = omp_get_wtime() - t_start;
    cout << "Time elapsed to compute interpolants = " << t_end << endl;

    cout << "Approximating solution of peaks function derivatives..." << endl;
    vector<double> f_approx(pointset.TotalPoints+1);
    vector<double> dfdx_approx(pointset.TotalPoints+1);
    vector<double> dfdy_approx(pointset.TotalPoints+1);
    vector<double> dfdxx_approx(pointset.TotalPoints+1);
    vector<double> dfdyy_approx(pointset.TotalPoints+1);
    vector<double> dfdxx_approx_direct(pointset.TotalPoints+1);
    vector<double> dfdyy_approx_direct(pointset.TotalPoints+1);

    vector<double> f_error(pointset.TotalPoints+1);
    vector<double> dfdx_error(pointset.TotalPoints+1);
    vector<double> dfdy_error(pointset.TotalPoints+1);
    vector<double> dfdxx_error(pointset.TotalPoints+1);
    vector<double> dfdyy_error(pointset.TotalPoints+1);
    vector<double> dfdxx_error_direct(pointset.TotalPoints+1);
    vector<double> dfdyy_error_direct(pointset.TotalPoints+1);
    
    double f_max_error = 0;
    double fx_max_error = 0;
    double fy_max_error = 0;
    double sum_exact_f = 0;
    double sum_exact_fx = 0;
    double sum_exact_fy = 0;
    double sum_exact_fxx = 0;
    double sum_exact_fyy = 0;
    double sum_f = 0;
    double sum_fx = 0;
    double sum_fy = 0;
    double sum_fxx = 0;
    double sum_fyy = 0;
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        f_approx[i] = 0;
        dfdx_approx[i] = 0;
        dfdy_approx[i] = 0;
        dfdxx_approx_direct[i] = 0;
        dfdyy_approx_direct[i] = 0;
        for (int s = 1; s <= kernels[i].TotalNbrs; s++)
        {
            int nbrpnt = kernels[i].nbrs[s];
            f_approx[i] = f_approx[i]+ kernels[i].N[s]*f[nbrpnt];
            dfdx_approx[i] = dfdx_approx[i]+ kernels[i].dNdx[s]*(f[nbrpnt]-f[i]);
            dfdy_approx[i] = dfdy_approx[i]+ kernels[i].dNdy[s]*(f[nbrpnt]-f[i]);
            dfdxx_approx_direct[i] = dfdxx_approx_direct[i]+ kernels[i].dNdxx[s]*(f[nbrpnt]-f[i]);
            dfdyy_approx_direct[i] = dfdyy_approx_direct[i]+ kernels[i].dNdyy[s]*(f[nbrpnt]-f[i]);
        }
        f_error[i] = 100*abs(f_approx[i]-f[i])/f[i];
        dfdx_error[i] = 100*abs(dfdx_approx[i]-dfdx[i])/dfdx[i];
        dfdy_error[i] = 100*abs(dfdy_approx[i]-dfdy[i])/dfdy[i];
        dfdxx_error_direct[i] = 100*abs(dfdxx_approx_direct[i]-dfdxx[i])/dfdxx[i];
        dfdyy_error_direct[i] = 100*abs(dfdyy_approx_direct[i]-dfdyy[i])/dfdyy[i];

        double f_num = abs(f[i]-f_approx[i]);
        if (f_num >= f_max_error)
        {f_max_error = f_num;}
        double fx_num = abs(dfdx[i]-dfdx_approx[i]);
        if (fx_num >= fx_max_error)
        {fx_max_error = fx_num;}
        double fy_num = abs(dfdy[i]-dfdy_approx[i]);
        if (fy_num >= fy_max_error)
        {fy_max_error = fy_num;}

        sum_exact_f = sum_exact_f + pow(f[i],2);
        sum_exact_fx = sum_exact_fx + pow(dfdx[i],2);
        sum_exact_fy = sum_exact_fy + pow(dfdy[i],2);
        sum_exact_fxx = sum_exact_fxx + pow(dfdxx[i],2);
        sum_exact_fyy = sum_exact_fyy + pow(dfdyy[i],2);
        sum_f = sum_f + pow((f_approx[i]-f[i]),2);
        sum_fx = sum_fx + pow((dfdx_approx[i]-dfdx[i]),2);
        sum_fy = sum_fy + pow((dfdy_approx[i]-dfdy[i]),2);
        sum_fxx = sum_fxx + pow((dfdxx_approx_direct[i]-dfdxx[i]),2);
        sum_fyy = sum_fyy + pow((dfdyy_approx_direct[i]-dfdyy[i]),2);
    }
    double error_f = pow((sum_f/sum_exact_f),0.5);
    double error_fx = pow((sum_fx/sum_exact_fx),0.5);
    double error_fy = pow((sum_fy/sum_exact_fy),0.5);
    double error_fxx = pow((sum_fxx/sum_exact_fxx),0.5);
    double error_fyy = pow((sum_fyy/sum_exact_fyy),0.5);

    double sum_fxx2 = 0;
    double sum_fyy2 = 0;
    for (int i = 1; i <= pointset.TotalPoints; i++)
    {
        dfdxx_approx[i] = 0;
        dfdyy_approx[i] = 0;
        for (int s = 1; s <= kernels[i].TotalNbrs; s++)
        {
            dfdxx_approx[i] = dfdxx_approx[i]+ kernels[i].dNdx[s]*dfdx_approx[kernels[i].nbrs[s]];
            dfdyy_approx[i] = dfdyy_approx[i]+ kernels[i].dNdy[s]*dfdy_approx[kernels[i].nbrs[s]];
        }
        sum_fxx2 = sum_fxx2 + pow((dfdxx_approx[i]-dfdxx[i]),2);
        sum_fyy2 = sum_fyy2 + pow((dfdyy_approx[i]-dfdyy[i]),2);
    }
    double error_fxx2 = pow((sum_fxx2/sum_exact_fxx),0.5);
    double error_fyy2 = pow((sum_fyy2/sum_exact_fyy),0.5);

    cout << "Solution summary..." << endl;
    cout << "compute time (sec) = " << t_end << endl;
    cout << "F error = " << error_f << endl;
    cout << "dFdx error = " << error_fx << endl;
    cout << "dFdy error = " << error_fy << endl;
    cout << "dFdxx error direct = " << error_fxx << endl;
    cout << "dFdxx error indirect = " << error_fxx2 << endl;
    cout << "dFdyy error direct = " << error_fyy << endl;
    cout << "dFdyy error indirect = " << error_fyy2 << endl;

    cout << "Printing solutions..." << endl;
    pointset.printMeshVTK(f, "f_exact", 0);
    pointset.printMeshVTK(dfdx,"dfdx_exact", 0);
    pointset.printMeshVTK(dfdy, "dfdy_exact", 0);
    pointset.printMeshVTK(dfdxx, "dfdxx_exact", 0);
    pointset.printMeshVTK(dfdyy, "dfdyy_exact", 0);

    pointset.printMeshVTK(f_approx,"f_approx", 0);
    pointset.printMeshVTK(dfdx_approx,"dfdx_approx", 0);
    pointset.printMeshVTK(dfdy_approx,"dfdy_approx", 0);
    pointset.printMeshVTK(dfdxx_approx_direct,"dfdxx_approx_direct", 0);
    pointset.printMeshVTK(dfdyy_approx_direct,"dfdyy_approx_direct", 0);

    pointset.printMeshVTK(f_error,"f_error", 0);
    pointset.printMeshVTK(dfdx_error,"dfdx_error", 0);
    pointset.printMeshVTK(dfdy_error,"dfdy_error", 0);
    pointset.printMeshVTK(dfdxx_error,"dfdxx_error", 0);
    pointset.printMeshVTK(dfdyy_error,"dfdyy_error", 0);
    pointset.printMeshVTK(dfdxx_approx,"dfdxx_approx_indirect", 0);
    pointset.printMeshVTK(dfdyy_approx,"dfdyy_approx_indirect", 0);

    return 0;
}
