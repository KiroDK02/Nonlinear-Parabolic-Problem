#include "functions.h"

void calcLU(matrix &A)
{
   auto &ggl = A.ggl, &ggu = A.ggu;
   auto &di = A.di;
   const int ASize = di.size(), K = 1;

   for (int i = 0; i < ASize; i++)
   {
      double sumDi = 0;
      int j = i - K;
      for (int jl = 0; jl < K; jl++, j++)
      {
         double sum1 = 0, sum2 = 0;
         if (j < 0) continue;
         int mj = i - j;
         for (int mi = 0; mi < jl; mi++, mj++)
         {
            sum1 += ggu[i][mi] * ggl[j][mj];
            sum2 += ggu[j][mj] * ggl[i][mi];
         }
         ggu[i][jl] -= sum1;
         ggl[i][jl] = (ggl[i][jl] - sum2) / di[j];
         sumDi += ggu[i][jl] * ggl[i][jl];
      }
      di[i] -= sumDi;
   }
}

void calcY(matrix &A, vector <double> &y)
{
   auto &ggl = A.ggl;
   auto &b = A.b;
   const int ASize = b.size(), K = 1;
   y.resize(ASize);

   for (int i = 0; i < ASize; i++)
   {
      double sum = 0;
      int j = i - K;
      for (int jl = 0; jl < K; jl++, j++)
      {
         if (j < 0) continue;
         sum += ggl[i][jl] * y[j];
      }
      y[i] = b[i] - sum;
   }
}

void calcQ(matrix &A, vector <double> &y, vector <double> &q)
{
   auto &ggu = A.ggu;
   auto &di = A.di;
   const int ASize = y.size(), K = 1;

   for (int i = ASize - 1; i >= 0; i--)
   {
      q[i] = y[i] / di[i];
      int j = i - K;
      for (int jl = 0; jl < K; jl++, j++)
      {
         if (j < 0) continue;
         y[j] -= ggu[i][jl] * q[i];
      }
   }
   y.clear();
}