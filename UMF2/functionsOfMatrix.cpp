#include "functions.h"

std::function <double(double)> u0 = [](double x) { return x * x; };
vector <std::function <double(double)>> lambda
{
   [](double u) { return u + 1; },
   [](double u) { return u + 1; },
};

vector <std::function <double(double)>> diffLambda
{
   [](double u) { return 1; },
   [](double u) { return u + 1;  },
};

vector <std::function<double(double, double)>> F
{
   [](double x, double t) { return -(6 * x * x + 2); },
   [](double x, double t) { return x; }
};

vector <std::function<double(double)>> ug
{
   [](double t) { return 0; },
   [](double t) { return 16; }
};

vector <vector<double>> C = { { 2, 1 },
                              { 1, 2 } };
vector <vector<double>> G = { { 1, -1 },
                              { -1, 1 } };

void createPortrait(matrix &A, grid &g)
{
   auto &ggl = A.ggl, &ggu = A.ggu;
   auto &di = A.di, &b = A.b;
   const int xSize = g.x.size();

   di.resize(xSize, 0);
   b.resize(xSize, 0);
   ggl.resize(xSize);
   ggu.resize(xSize);
   for (int i = 0; i < xSize; i++)
   {
      ggl[i].resize(1);
      ggu[i].resize(1);
   }
}

void addElementInMatrix(matrix &A, double elem, int i, int j)
{
   auto &ggl = A.ggl, &ggu = A.ggu;
   auto &di = A.di;

   if (i == j)
      di[i] += elem;
   else
      if (i > j)
         ggl[i][0] += elem;
      else
         ggu[j][0] += elem;
}

void begApproximation(vector <vector<double>> &q, grid &g)
{
   auto &x = g.x, &t = g.t;
   const int xSize = x.size(), tSize = t.size();
   q.resize(tSize);
   for (int i = 0; i < tSize; i++)
      q[i].resize(xSize);

   for (int i = 0; i < xSize; i++)
      q[0][i] = u0(x[i]);
}

void calcGlobalMatrix(double tValue, double deltaT, matrix &A, grid &g, vector<double> &q, vector<double> &qtPrev, vector <double> &sigma, vector <int> &Ix, elementsMesh &eM)
{
   auto &x = g.x, &b = A.b;
   auto &M = eM.M;
   const int xSize = x.size();
   int l = 0;

   for (int p = 0; p < xSize - 1; p++)
   {
      int numArea = returnNumberSubAreaOfMesh(Ix, M, p, l) - 1;
      double hk = x[p + 1] - x[p], coefG = (lambda[numArea](q[p]) + lambda[numArea](q[p + 1])) / (2.0 * hk),
                  coefC = sigma[numArea] * hk / (6.0 * deltaT), coefb1 = hk / 6.0;
      vector <int> globalNum = { p, p + 1 };
      for (int i = 0; i < 2; i++)
      {
         vector <double> f = { F[numArea](x[p], tValue), F[numArea](x[p + 1], tValue) };
         vector <double> tempqt = { qtPrev[p], qtPrev[p + 1] };
         double sumbi = 0;
         for (int j = 0; j < 2; j++)
         {
            double elem = coefG * G[i][j] + coefC * C[i][j];
            addElementInMatrix(A, elem, globalNum[i], globalNum[j]);
            sumbi += coefb1 * C[i][j] * f[j] + coefC * C[i][j] * tempqt[j];
         }
         b[globalNum[i]] += sumbi;
      }
   }
}

void calcGlobalMatrix¿ForMN(double v, double tValue, double deltaT, matrix &A, grid &g, vector<double> &q, vector<double> &qtPrev, vector <double> &sigma, vector <int> &Ix, elementsMesh &eM)
{
   auto &x = g.x, &b = A.b;
   auto &M = eM.M;
   const int xSize = x.size();
   int l = 0;

   for (int p = 0; p < xSize - 1; p++)
   {
      int numArea = returnNumberSubAreaOfMesh(Ix, M, p, l) - 1;
      double hk = x[p + 1] - x[p], coefG = (lambda[numArea](q[p]) + lambda[numArea](q[p + 1])) / (2.0 * hk),
         coefC = sigma[numArea] * hk / (6.0 * deltaT), coefb1 = hk / 6.0,
         diffLambda1 = diffLambda[numArea](q[p]) / (2.0 * hk), 
         diffLambda2 = diffLambda[numArea](q[p + 1]) / (2.0 * hk);
      vector <vector <vector<double>>> matrixDiff = { {{ diffLambda1, diffLambda2 },
                                                       { -diffLambda1, -diffLambda2 }},
                                                      {{ -diffLambda1, -diffLambda2 },
                                                       { diffLambda1, diffLambda2 }} };
      vector <int> globalNum = { p, p + 1 };
      for (int i = 0; i < 2; i++)
      {
         vector <double> f = { F[numArea](x[p], tValue), F[numArea](x[p + 1], tValue) };
         vector <double> tempq = { q[p], q[p + 1] };
         vector <double> tempqt = { qtPrev[p], qtPrev[p + 1] };
         double sumbiMain = 0;
         for (int j = 0; j < 2; j++)
         {
            double elemA = 0;
            double sumbi2 = 0, elemAdiff = 0;
            for (int r = 0; r < 2; r++)
            {
               elemAdiff += matrixDiff[i][r][j] * tempq[r];
               sumbi2 += matrixDiff[i][j][r] * tempq[r];
            }
            elemA += v * elemAdiff;
            sumbiMain += v * sumbi2 * tempq[j];
            addElementInMatrix(A, elemA, globalNum[i], globalNum[j]);
         }
         b[globalNum[i]] += sumbiMain;
      }
   }
}

void boundaryCondition(double tValue, matrix &A, vector <vector <int>> &bC, grid &g)
{
   auto &ggl = A.ggl, &ggu = A.ggu;
   auto &di = A.di, &b = A.b, &x = g.x;
   const int countCond = 2, xLast = x.size() - 1;

   for (int k = 0; k < countCond; k++)
   {
      int typeCond = bC[k][0], numFunc = bC[k][1] - 1;
      if (typeCond == 1)
         if (k == 0)
         {
            ggu[1][0] = 0;
            di[0] = 1;
            b[0] = ug[numFunc](tValue);
         }
         else
         {
            ggl[xLast][0] = 0;
            di[xLast] = 1;
            b[xLast] = ug[numFunc](tValue);
      }
   }
}

double calcDiscrepancy(matrix &A, vector <double> &q)
{
   auto &ggl = A.ggl, &ggu = A.ggu;
   auto &di = A.di, &b = A.b;
   const int size = di.size();
   double discrepancy = 0, normb = 0;
  
   for (int i = 0; i < size; i++)
   {
      normb += b[i] * b[i];
      int j1 = i - 1, j2 = i + 1;
      double sum = di[i] * q[i];
      if (j1 >= 0) sum += ggl[i][0] * q[j1];
      if (j2 < size) sum += ggu[j2][0] * q[j2];
      discrepancy += (b[i] - sum) * (b[i] - sum);
   }
   discrepancy = sqrt(discrepancy / normb);
   return discrepancy;
}

int calcMPI(double tValue, double deltaT, matrix &A, grid &g, vector<double> &q, vector<double> &q0, vector<vector <int>> &bC, vector <double> &sigma, vector <int> &Ix, elementsMesh &eM)
{
   calcGlobalMatrix(tValue, deltaT, A, g, q0, q0, sigma, Ix, eM);
   boundaryCondition(tValue, A, bC, g);
   double discrepancy = calcDiscrepancy(A, q0), eps = 1e-14, w = 1;
   vector <double> qTemp{ };
   qTemp.resize(g.x.size());
   vecEqualVec(qTemp, q0);
   int k = 0;

   if (discrepancy < eps)
   {
      vecEqualVec(q, q0);
      return 0;
   }
   for (k = 1; discrepancy > eps && k < 10000; k++)
   {
      vector <double> y{ };
      calcLU(A);
      calcY(A, y);
      calcQ(A, y, q);
      calcRelaxation(qTemp, q, w);
      vecEqualVec(qTemp, q);
      clearMatrix(A);
      calcGlobalMatrix(tValue, deltaT, A, g, q, q0, sigma, Ix, eM);
      boundaryCondition(tValue, A, bC, g);
      discrepancy = calcDiscrepancy(A, q);
      printf_s("   k: %d discrepancy: %.15lf\n", k, discrepancy);
   }

   return k - 1;
}

int calcMN(double tValue, double deltaT, matrix &A, grid &g, vector<double> &q, vector<double> &q0, vector<vector <int>> &bC, vector <double> &sigma, vector <int> &Ix, elementsMesh &eM)
{
   calcGlobalMatrix(tValue, deltaT, A, g, q0, q0, sigma, Ix, eM);
   boundaryCondition(tValue, A, bC, g);
   double discrepancy = calcDiscrepancy(A, q0), eps = 1e-14, w = 1;
   vector <double> qTemp{ };
   qTemp.resize(g.x.size());
   vecEqualVec(qTemp, q0);
   vecEqualVec(q, q0);
   int k = 0;
   if (discrepancy < eps)
      return 0;

   for (k = 1; discrepancy > eps && k < 10000; k++)
   {
      vector <double> y{ };
      double v = 1;
      calcGlobalMatrix¿ForMN(v, tValue, deltaT, A, g, q, q0, sigma, Ix, eM);
      boundaryCondition(tValue, A, bC, g);
      calcLU(A);
      while (!matrixPositiveDefinite(A))
      {
         v /= 2;
         clearMatrix(A);
         calcGlobalMatrix(tValue, deltaT, A, g, q, q0, sigma, Ix, eM);
         calcGlobalMatrix¿ForMN(v, tValue, deltaT, A, g, q, q0, sigma, Ix, eM);
         boundaryCondition(tValue, A, bC, g);
         calcLU(A);
      }
      calcY(A, y);
      calcQ(A, y, q);
      calcRelaxation(qTemp, q, w);
      vecEqualVec(qTemp, q);
      clearMatrix(A);
      calcGlobalMatrix(tValue, deltaT, A, g, q, q0, sigma, Ix, eM);
      boundaryCondition(tValue, A, bC, g);
      discrepancy = calcDiscrepancy(A, q);
      printf_s("   k: %d discrepancy: %.15lf\n", k, discrepancy);
   }

   return k - 1;
}

bool matrixPositiveDefinite(matrix &A)
{
   auto &di = A.di;
   const int size = di.size();

   for (int k = 0; k < size; k++)
   {
      double det = 1;
      for (int i = 0; i < k + 1; i++)
         det *= di[i];
      if (det < 0)
         return false;
   }
   return true;
}

void calcRelaxation(vector <double> &q0, vector <double> &q1, double w)
{
   const int n = q1.size();
   for (int i = 0; i < n; i++)
      q1[i] = w * q1[i] + (1.0 - w) * q0[i];
}

vector <int> calcSolutionMN(matrix &A, grid &g, vector<vector<double>> &q, vector<vector <int>> &bC, vector <double> &sigma, vector <int> &Ix, elementsMesh &eM)
{
   auto &t = g.t;
   const int tSize = t.size();
   vector <int> countIter{ };
   countIter.resize(tSize);

   for (int ti = 1; ti < tSize; ti++)
   {
      double deltaT = t[ti] - t[ti - 1];
      printf_s("t[%d] = %.5lf\n", ti, t[ti]);
      clearMatrix(A);
      countIter[ti] = calcMN(t[ti], deltaT, A, g, q[ti], q[ti - 1], bC, sigma, Ix, eM);
   }
   return countIter;
}

vector <int> calcSolution(matrix &A, grid &g, vector<vector<double>> &q, vector<vector <int>> &bC, vector <double> &sigma, vector <int> &Ix, elementsMesh &eM)
{
   auto &t = g.t;
   const int tSize = t.size();
   vector <int> countIter{ };
   countIter.resize(tSize);

   for (int ti = 1; ti < tSize; ti++)
   {
      double deltaT = t[ti] - t[ti - 1];
      printf_s("t[%d] = %.5lf\n", ti, t[ti]);
      clearMatrix(A);
      countIter[ti] = calcMPI(t[ti], deltaT, A, g, q[ti], q[ti - 1], bC, sigma, Ix, eM);
   }
   return countIter;
}

void clearMatrix(matrix &A)
{
   auto &ggl = A.ggl, &ggu = A.ggu;
   auto &di = A.di, &b = A.b;
   const int size = di.size();
   for (int i = 0; i < size; i++)
   {
      ggl[i][0] = 0;
      ggu[i][0] = 0;
      di[i] = 0;
      b[i] = 0;
   }
}

void vecEqualVec(vector <double> &q1, vector <double> &q2)
{
   const int VSize = q1.size();
   for (int i = 0; i < VSize; i++)
      q1[i] = q2[i];
}


