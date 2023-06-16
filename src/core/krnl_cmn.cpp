#include "../pntwrks.h"

vector<KERNEL> computeKernels(vector<POINT> points, vector<double> phi, int TotalPoints, int dim, KERNEL_PRMTRS kernel_prmtrs, double kernel_radius)
{
    int order = kernel_prmtrs.approximation_order;
    string interp = kernel_prmtrs.approximation_method;
    string supp = kernel_prmtrs.support_method;
    int RBF, WLS, MLS, SPH;
    int interpolation = 0;

    if (kernel_prmtrs.RBF == "IMQ")
    {RBF = 1;}
    if (kernel_prmtrs.RBF == "GE")
    {RBF = 2;}
    if (kernel_prmtrs.RBF == "TPS")
    {RBF = 3;}

    if (kernel_prmtrs.MLS == "S3")
    {MLS = 1;}
    if (kernel_prmtrs.MLS == "S4")
    {MLS = 2;}
    if (kernel_prmtrs.MLS == "S5")
    {MLS = 3;}
    if (kernel_prmtrs.MLS == "RegS4")
    {MLS = 4;}

    if (kernel_prmtrs.WLS == "LS")
    {WLS = 0;}
    if (kernel_prmtrs.WLS == "GE")
    {WLS = 1;}
    if (kernel_prmtrs.WLS == "RegS4")
    {WLS = 2;}
    if (kernel_prmtrs.WLS == "S3")
    {WLS = 3;}
    if (kernel_prmtrs.WLS == "S4")
    {WLS = 4;}
    if (kernel_prmtrs.WLS == "S5")
    {WLS = 5;}

    if (kernel_prmtrs.SPH == "GE")
    {SPH = 1;}
    if (kernel_prmtrs.SPH == "S3")
    {SPH = 2;}
    if (kernel_prmtrs.SPH == "S5")
    {SPH = 3;}
    if (kernel_prmtrs.SPH == "WC4")
    {SPH = 4;}
    if (kernel_prmtrs.SPH == "WC6")
    {SPH = 5;}

    if (interp == "MLS")
    {interpolation = 1;}
    else if (interp == "RBF")
    {interpolation = 2;}
    else if (interp == "SPH")
    {interpolation = 3;}
    else if (interp == "WLS")
    {interpolation = 4;}
    else
    {
        cout << "ERROR: interpolation type is undefined." << endl;
        cout << "interpolation type = " << interpolation << endl;
        exit (0);
    }
    double eps = 1.0e-6;
    double factor = 1;
    vector<KERNEL> kernels(TotalPoints+1);
    cout << "kernel_prmtrs.grad_order = " << kernel_prmtrs.grad_order << endl;
    if (kernel_prmtrs.recompute_approximants == true)
    {
        int cc = 1;
#pragma omp parallel for
        for (int i = 1; i <= TotalPoints; i++)
        {
            if (supp == "BF")
            {kernels[i].collectNbrsBruteForce(points[i], points, TotalPoints, kernel_radius);}
            else if(supp == "SP")
            {kernels[i].collectNbrsSamePhaseBruteForce(points[i], phi[i], points, phi, TotalPoints, kernel_radius);}

            if (interpolation == 1)
            {kernels[i].computeMLS(points[i], points, kernels[i].nbrs, kernels[i].TotalNbrs, MLS, order, kernel_radius);}
            else if (interpolation == 2)
            {kernels[i].computeRBF(points[i], points, kernels[i].nbrs, kernels[i].TotalNbrs, RBF, order, 0.1*kernel_radius, kernel_prmtrs.RBF_alpha);}
            else if (interpolation == 3)
            {kernels[i].computeSPH(points[i], points, kernels[i].nbrs, kernels[i].TotalNbrs, SPH, kernel_radius, dim);}
            else if (interpolation == 4)
            {kernels[i].computeWLS(points[i], points, kernels[i].nbrs, kernels[i].TotalNbrs, WLS, order, kernel_radius);}

    	    kernels[i].dNdxx.resize(kernels[i].TotalNbrs+1);
    	    kernels[i].dNdyy.resize(kernels[i].TotalNbrs+1);
    	    kernels[i].dNdzz.resize(kernels[i].TotalNbrs+1);
    	    kernels[i].nabla2.resize(kernels[i].TotalNbrs+1);
    	    kernels[i].dNdxy.resize(kernels[i].TotalNbrs+1);
    	    kernels[i].dNdxz.resize(kernels[i].TotalNbrs+1);
    	    kernels[i].dNdyz.resize(kernels[i].TotalNbrs+1);
    	    if (kernel_prmtrs.grad_order == 1)
            {
    	    	for (int j = 1; j <= kernels[i].TotalNbrs; j++)
                {
                   	int nbr = kernels[i].nbrs[j];
            	   	double d = points[nbr].distance(points[i])+eps;
    	           	double nijx = (points[nbr].x-points[i].x)/d;
    	           	double nijy = (points[nbr].y-points[i].y)/d;
                	double nijz = (points[nbr].z-points[i].z)/d;
                    kernels[i].dNdxx[j] = factor*(nijx/d)*kernels[i].dNdx[j];
                    kernels[i].dNdyy[j] = factor*(nijy/d)*kernels[i].dNdy[j];
                	kernels[i].dNdzz[j] = factor*(nijz/d)*kernels[i].dNdz[j];
    				kernels[i].nabla2[j] = kernels[i].dNdxx[j]+kernels[i].dNdyy[j]+kernels[i].dNdzz[j];
            	    kernels[i].dNdxy[j] = factor*(nijx/d)*kernels[i].dNdy[j];
    	            kernels[i].dNdxz[j] = factor*(nijx/d)*kernels[i].dNdz[j];
                  	kernels[i].dNdyz[j] = factor*(nijy/d)*kernels[i].dNdz[j];
    	    	}
    	    }
    	    else if (kernel_prmtrs.grad_order == 0)	    
    	    {
    			for (int j = 1; j <= kernels[i].TotalNbrs; j++)
               	{
                   	int nbr = kernels[i].nbrs[j];
            	   	double d = points[nbr].distance(points[i])+eps;
                   	double xj_xi = points[nbr].x-points[i].x;
                   	double yj_yi = points[nbr].y-points[i].y;
                   	double zj_zi = points[nbr].z-points[i].z;
                   	kernels[i].dNdx[i] = factor*(xj_xi)*kernels[i].N[j];
                   	kernels[i].dNdy[i] = factor*(yj_yi)*kernels[i].N[j];
                   	kernels[i].dNdz[i] = factor*(zj_zi)*kernels[i].N[j];
    				kernels[i].nabla2[j] = factor*kernels[i].N[j]/(d*d);
    	    	}
    	    }
    	    if (i/1000 == cc)
            {
                cc++;
                cout  << "Processed " << i << " of " << TotalPoints << " points..."  << endl;
            }
        }
        writeKernels(kernels, TotalPoints);
    }
    else{kernels = readKernels(TotalPoints);}
    return kernels;
}
void KERNEL::bases(int order, double x, double y, double z)
{
    if (order == 1)
    {monomials = 4;}
    else if (order == 2)
    {monomials = 10;}
    else if (order == 3)
    {monomials = 20;}
    else
    {cout << "ERROR: order for basis construction is invalid" << endl;
    cout << "order = " << order << endl;exit (0);}
    double e = 0;

    double bases[21] =    {e, 1, x, y, z, x*x, y*y, z*z, x*y, x*z, y*z, x*x*x, y*y*y, z*z*z, x*x*y, x*x*z, y*y*x, y*y*z, z*z*x, z*z*y, x*y*z};
    double bases_dx[21] = {e, 0, 1, 0, 0, 2*x, 0, 0 , y, z, 0, 3*x*x, 0, 0, 2*x*y, 2*x*z, y*y, 0, z*z, 0, y*z};
    double bases_dy[21] = {e, 0, 0, 1, 0, 0, 2*y, 0, x, 0, z, 0, 3*y*y, 0, x*x, 0, 2*y*x, 2*y*z, 0, z*z, x*z};
    double bases_dz[21] = {e, 0, 0, 0, 1, 0, 0, 2*z, 0, x, y, 0, 0, 3*z*z, 0, x*x, 0, y*y, 2*z*x, 2*z*y, x*y};
    double bases_dxx[21] ={e, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 6*x, 0, 0, 2*y, 2*z, 0, 0, 0, 0, 0};
    double bases_dyy[21] ={e, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 6*y, 0, 0, 0, 2*x, 2*z, 0, 0, 0};
    double bases_dzz[21] = {e, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 6*z, 0, 0, 0, 0, 2*x, 2*y, 0};
    double bases_dxy[21] = {e, 0, 0, 0, 0, 0, 0, 0 , 1, 0, 0, 0, 0, 0, 2*x, 0, 2*y, 0, 0, 0, z};
    double bases_dxz[21] = {e, 0, 0, 0, 0, 0, 0, 0 , 0, 1, 0, 0, 0, 0, 0, 2*x, 0, 0, 2*z, 0, y};
    double bases_dyz[21] = {e, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2*y, 0, 2*z, x};

    N.resize(monomials+1);
    dNdx.resize(monomials+1);
    dNdy.resize(monomials+1);
    dNdz.resize(monomials+1);
    dNdxx.resize(monomials+1);
    dNdyy.resize(monomials+1);
    dNdzz.resize(monomials+1);
    dNdxy.resize(monomials+1);
    dNdxz.resize(monomials+1);
    dNdyz.resize(monomials+1);

    for (int i = 1; i <= monomials; i++)
    {
        N[i] = bases[i];
        dNdx[i] = bases_dx[i];
        dNdy[i] = bases_dy[i];
        dNdz[i] = bases_dz[i];
        dNdxx[i] = bases_dxx[i];
        dNdyy[i] = bases_dyy[i];
        dNdzz[i] = bases_dzz[i];
        dNdxy[i] = bases_dxy[i];
        dNdxz[i] = bases_dxz[i];
        dNdyz[i] = bases_dyz[i];
    }
}
KERNEL::~KERNEL()
{
    //destructor
}
