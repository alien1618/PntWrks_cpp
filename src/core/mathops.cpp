#include "../pntwrks.h"


vector<double> solve(vector<vector<double> > K, vector<double> F, int TotalPoints)
{
    vector<vector<double> > inv_K(TotalPoints+1, vector<double>(TotalPoints+1));
    vector<double> X(TotalPoints+1);
    inv_K = inv(K, TotalPoints);
    double mat_mult = 0;
    for(int i = 1; i<=TotalPoints ; i++)
    {
        for(int j = 1; j<=TotalPoints ; j++)
        {mat_mult = mat_mult + (inv_K[i][j] * F[j]);}
        X[i] = mat_mult;
        mat_mult = 0;
    }
    return X;
}

vector<vector<double> > inv(vector<vector<double> > a, int size)
{
    vector<vector<double> > b(size+1, vector<double>(size+1));
    double s,t;
#pragma omp parallel for
    for (int j = 1; j <= size; j++)
    {
#pragma omp parallel for
        for(int i = 1; i <= size; i++)
        {b[i][j] = 0;}
    }
#pragma omp parallel for
    for (int i = 1; i<= size; i++)
    {b[i][i] = 1;}
    //The following code actually performs the matrix inversion
    for (int j = 1; j <= size; j++)
    {
        for (int i = j; i<= size; i++)
        {
            if (a[i][j] != 0)
            {
                for (int k = 1; k<= size; k++)
                {
                    s = a[j][k];
                    a[j][k] = a[i][k];
                    a[i][k] = s;
                    s = b[j][k];
                    b[j][k] = b[i][k];
                    b[i][k] = s;
                }
                t = 1/a[j][j];

                for (int k = 1; k <= size; k++)
                {
                    a[j][k] = t * a[j][k];
                    b[j][k] = t * b[j][k];
                }
                for (int L = 1; L <= size; L++)
                {
                    if (L != j)
                    {
                        t = -a[L][j];
                        for (int k = 1; k<= size; k++)
                        {
                            a[L][k] = a[L][k] + t * a[j][k];
                            b[L][k] = b[L][k] + t * b[j][k];
                        }
                    }
                }
            }
            break;
        }
    }
    a.clear();
    return b;
}

double checkConvergence(int iternum, vector<double> U_temp, vector<double> U_temp_old, int TotalPoints)
{
    double max_diff=0;
    if (iternum > 1)
    {
        vector<double> difference(TotalPoints+1);
        //cout << "Checking Convergence..." << endl;
#pragma omp parallel for
        for (int i = 1; i <= TotalPoints; i++)
        {difference[i] = abs(U_temp[i] - U_temp_old[i]);}
        max_diff = maximum<double>(difference, TotalPoints);
    }
    return max_diff;
}


//==============================================================================
// RANDOM NUMBER GENERATION
//==============================================================================
double random0to1()
{return rand()/ double(RAND_MAX);}

void randomSeed()
{srand(time(0));}

double maxnum(double a, double b)
{
    double answer=0;
    if (a > b)
    {answer = a;}
    else if (b > a)
    {answer = b;}
    else if (b == a)
    {answer = b;}
    return answer;
}
double minnum(double a, double b)
{
    double answer=0;
    if (a < b)
    {answer = a;}
    else if (b < a)
    {answer = b;}
    else if (b == a)
    {answer = b;}
    return answer;
}
double mod(double number, double divnum)
{
    double ans = number/divnum;
    double remainder = number - double(int(ans) * divnum);
    return remainder;
}

//----------------------------------------------------------------------
// FUNCTION TO COMPUTE UNIQUE VALUES OF AN ARRAY VECTOR
//----------------------------------------------------------------------
tuple<vector<int>, int> UniqueNums(vector<int> V, int z)
{

    int m, k;
    vector <int> answer;
    vector<int> index(z+1);
    int entry;

    for (int i = 1; i <= z; i++)
    {index[i] = 0;}

    for (int i = 1; i <= z; i++)
    {
        m = 0;
        entry = V[i];
        for (int j = 1; j <= z; j++)
        {
            if ((entry == V[j]) &  (index[j] == 0))
            {
                m = m + 1;
                index[j] = m;
            }
        }
    }

    k = 0;
    answer.resize(1);
    for (int i = 1; i <= z; i++)
    {
        if (index[i] == 1)
        {
            k = k + 1;
            answer.push_back(1);
            answer[k] = V[i];
        }
    }
    return make_tuple(answer,k);
}

