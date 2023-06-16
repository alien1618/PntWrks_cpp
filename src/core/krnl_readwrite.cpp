#include "../pntwrks.h"

void writeKernels(vector<KERNEL> kernels, int TotalPoints)
{
    int c = 1;
    //-------------------------------------------------
    cout << "Computing gradients..." << endl;
    //-------------------------------------------------
    ofstream ofile_N("logs/N.txt");
    ofstream ofile_dNdx("logs/dNdx.txt");
    ofstream ofile_dNdy("logs/dNdy.txt");
    ofstream ofile_dNdz("logs/dNdz.txt");
    ofstream ofile_dNdxx("logs/dNdxx.txt");
    ofstream ofile_dNdyy("logs/dNdyy.txt");
    ofstream ofile_dNdzz("logs/dNdzz.txt");
    ofstream ofile_dNdxy("logs/dNdxy.txt");
    ofstream ofile_dNdxz("logs/dNdxz.txt");
    ofstream ofile_dNdyz("logs/dNdyz.txt");
    ofstream ofile_supp("logs/supp.txt");
    for (int i = 1; i <= TotalPoints; i++)
    {
        //WARNING: DON'T CHAGE THE INTERPOLATION METHOD HERE. IT WAS FOUND TO YIELD BEST 2nd ORDER INTERPOLATION RESULTS
        ofile_N << i << "\t" << kernels[i].TotalNbrs << "\t";
        ofile_dNdx << i << "\t" << kernels[i].TotalNbrs << "\t";
        ofile_dNdy << i << "\t" << kernels[i].TotalNbrs << "\t";
        ofile_dNdz << i << "\t" << kernels[i].TotalNbrs << "\t";
        ofile_dNdxx << i << "\t" << kernels[i].TotalNbrs << "\t";
        ofile_dNdyy << i << "\t" << kernels[i].TotalNbrs << "\t";
        ofile_dNdzz << i << "\t" << kernels[i].TotalNbrs << "\t";
        ofile_dNdxy << i << "\t" << kernels[i].TotalNbrs << "\t";
        ofile_dNdxz << i << "\t" << kernels[i].TotalNbrs << "\t";
        ofile_dNdyz << i << "\t" << kernels[i].TotalNbrs << "\t";
        ofile_supp << i << "\t" << kernels[i].TotalNbrs << "\t";
        for (int j = 1; j <= kernels[i].TotalNbrs; j++)
        {
            ofile_N << kernels[i].N[j] << "\t";
            ofile_dNdx << kernels[i].dNdx[j] << "\t";
            ofile_dNdy << kernels[i].dNdy[j] << "\t";
            ofile_dNdz << kernels[i].dNdz[j] << "\t";
            ofile_dNdxx << kernels[i].dNdxx[j] << "\t";
            ofile_dNdyy << kernels[i].dNdyy[j] << "\t";
            ofile_dNdzz << kernels[i].dNdzz[j] << "\t";
            ofile_dNdxy << kernels[i].dNdxy[j] << "\t";
            ofile_dNdxz << kernels[i].dNdxz[j] << "\t";
            ofile_dNdyz << kernels[i].dNdyz[j] << "\t";
            ofile_supp << kernels[i].nbrs[j] << "\t";
        }
        ofile_N << endl;
        ofile_dNdx << endl;
        ofile_dNdy << endl;
        ofile_dNdz << endl;
        ofile_dNdxx << endl;
        ofile_dNdyy << endl;
        ofile_dNdzz << endl;
        ofile_dNdxy << endl;
        ofile_dNdxz << endl;
        ofile_dNdyz << endl;
        ofile_supp << endl;
        if (i/1000 == c)
        {
            c++;
            cout  << "Wrote " << i << " of " << TotalPoints << endl;
        }
    }
}
vector<KERNEL> readKernels(int total_nds)
{
    //------------------------------
    // READ SIMULATION PARAMETERS
    //------------------------------
    int size = 1e8;
    vector<KERNEL> kernels(2);
    vector<double> N(3);
    vector<double> dNdx(3);
    vector<double> dNdy(3);
    vector<double> dNdz(3);
    vector<double> dNdxx(3);
    vector<double> dNdyy(3);
    vector<double> dNdzz(3);
    vector<double> dNdxy(3);
    vector<double> dNdxz(3);
    vector<double> dNdyz(3);
    vector<double> supp(3);

    ifstream fin_N("logs/N.txt");
    ifstream fin_dNdx("logs/dNdx.txt");
    ifstream fin_dNdy("logs/dNdy.txt");
    ifstream fin_dNdz("logs/dNdz.txt");
    ifstream fin_dNdxx("logs/dNdxx.txt");
    ifstream fin_dNdyy("logs/dNdyy.txt");
    ifstream fin_dNdzz("logs/dNdzz.txt");
    ifstream fin_dNdxy("logs/dNdxy.txt");
    ifstream fin_dNdxz("logs/dNdxz.txt");
    ifstream fin_dNdyz("logs/dNdyz.txt");
    ifstream fin_supp("logs/supp.txt");
    if (!fin_N)
    {cout << ("ERROR: UNABLE TO OPEN INPUT FILE N.txt \n");exit(0);}
    else if(!fin_dNdx)
    {cout << ("ERROR: UNABLE TO OPEN INPUT FILE dNdx.txt \n");exit(0);}
    else if(!fin_dNdy)
    {cout << ("ERROR: UNABLE TO OPEN INPUT FILE dNdy.txt \n");exit(0);}
    else if(!fin_dNdz)
    {cout << ("ERROR: UNABLE TO OPEN INPUT FILE dNdz.txt \n");exit(0);}
    else if(!fin_dNdxx)
    {cout << ("ERROR: UNABLE TO OPEN INPUT FILE dNdxx.txt \n");exit(0);}
    else if(!fin_dNdyy)
    {cout << ("ERROR: UNABLE TO OPEN INPUT FILE dNdyy.txt \n");exit(0);}
    else if(!fin_dNdzz)
    {cout << ("ERROR: UNABLE TO OPEN INPUT FILE dNdzz.txt \n");exit(0);}
    else if(!fin_dNdxy)
    {cout << ("ERROR: UNABLE TO OPEN INPUT FILE dNdxy.txt \n");exit(0);}
    else if(!fin_dNdxz)
    {cout << ("ERROR: UNABLE TO OPEN INPUT FILE dNdxz.txt \n");exit(0);}
    else if(!fin_dNdyz)
    {cout << ("ERROR: UNABLE TO OPEN INPUT FILE dNdyz.txt \n");exit(0);}
    else
    {
        int i = 0;
        cout << "Scanning N.txt..." << endl;
        while (!fin_N.eof() && (i < size))
        {
            ++i;
            N.push_back(1);
            fin_N >> N[i];
        }
        i = 0;
        cout << "Scanning dNdx.txt..." << endl;
        while (!fin_dNdx.eof() && (i < size))
        {
            ++i;
            dNdx.push_back(1);
            fin_dNdx >> dNdx[i];
        }
        i = 0;
        cout << "Scanning dNdy.txt..." << endl;
        while (!fin_dNdy.eof() && (i < size))
        {
            ++i;
            dNdy.push_back(1);
            fin_dNdy >> dNdy[i];
        }
        i = 0;
        cout << "Scanning dNdz.txt..." << endl;
        while (!fin_dNdz.eof() && (i < size))
        {
            ++i;
            dNdz.push_back(1);
            fin_dNdz >> dNdz[i];
        }
        i = 0;
        cout << "Scanning dNdxx.txt..." << endl;
        while (!fin_dNdxx.eof() && (i < size))
        {
            ++i;
            dNdxx.push_back(1);
            fin_dNdxx >> dNdxx[i];
        }
        i = 0;
        cout << "Scanning dNdyy.txt..." << endl;
        while (!fin_dNdyy.eof() && (i < size))
        {
            ++i;
            dNdyy.push_back(1);
            fin_dNdyy >> dNdyy[i];
        }
        i = 0;
        cout << "Scanning dNdzz.txt..." << endl;
        while (!fin_dNdzz.eof() && (i < size))
        {
            ++i;
            dNdzz.push_back(1);
            fin_dNdzz >> dNdzz[i];
        }
        i = 0;
        cout << "Scanning dNdxy.txt..." << endl;
        while (!fin_dNdxy.eof() && (i < size))
        {
            ++i;
            dNdxy.push_back(1);
            fin_dNdxy >> dNdxy[i];
        }
        i = 0;
        cout << "Scanning dNdxz.txt..." << endl;
        while (!fin_dNdxz.eof() && (i < size))
        {
            ++i;
            dNdxz.push_back(1);
            fin_dNdxz >> dNdxz[i];
        }
        i = 0;
        cout << "Scanning dNdyz.txt..." << endl;
        while (!fin_dNdyz.eof() && (i < size))
        {
            ++i;
            dNdyz.push_back(1);
            fin_dNdyz >> dNdyz[i];
        }
        i = 0;
        cout << "Scanning supp.txt..." << endl;
        while (!fin_supp.eof() && (i < size))
        {
            ++i;
            supp.push_back(1);
            fin_supp >> supp[i];
        }
        vector<KERNEL> g(total_nds+1);
        kernels = g;
        int l = 1;

        for (int p = 1; p <= total_nds; p++)
        {
            int total_supp_nds_per_nd = int(N[l+1]);
            //cout << "nd_number = " << nd_number << " total_supp_nds_per_nd " << total_supp_nds_per_nd << endl;
            kernels[p].N.resize(total_supp_nds_per_nd+2);
            kernels[p].dNdx.resize(total_supp_nds_per_nd+2);
            kernels[p].dNdy.resize(total_supp_nds_per_nd+2);
            kernels[p].dNdz.resize(total_supp_nds_per_nd+2);
            kernels[p].dNdxx.resize(total_supp_nds_per_nd+2);
            kernels[p].dNdyy.resize(total_supp_nds_per_nd+2);
            kernels[p].dNdzz.resize(total_supp_nds_per_nd+2);
            kernels[p].dNdxy.resize(total_supp_nds_per_nd+2);
            kernels[p].dNdxz.resize(total_supp_nds_per_nd+2);
            kernels[p].dNdyz.resize(total_supp_nds_per_nd+2);
            kernels[p].TotalNbrs = total_supp_nds_per_nd;
            kernels[p].nbrs.resize(total_supp_nds_per_nd+2);

            int m = 0;
            for(int j = l+2; j <= l+1+total_supp_nds_per_nd; j++)
            {
                m++;
                kernels[p].N[m] = N[j];
                kernels[p].dNdx[m] = dNdx[j];
                kernels[p].dNdy[m] = dNdy[j];
                kernels[p].dNdz[m] = dNdz[j];
                kernels[p].dNdxx[m] = dNdxx[j];
                kernels[p].dNdyy[m] = dNdyy[j];
                kernels[p].dNdzz[m] = dNdzz[j];
                kernels[p].dNdxy[m] = dNdxy[j];
                kernels[p].dNdxz[m] = dNdxz[j];
                kernels[p].dNdyz[m] = dNdyz[j];
                kernels[p].nbrs[m] = supp[j];
            }
            l = l+1+total_supp_nds_per_nd+1;
        }
    }
    fin_N.close();
    fin_dNdx.close();
    fin_dNdy.close();
    fin_dNdz.close();
    fin_dNdxx.close();
    fin_dNdyy.close();
    fin_dNdzz.close();
    fin_dNdxy.close();
    fin_dNdxz.close();
    fin_dNdyz.close();

    return kernels;
}
