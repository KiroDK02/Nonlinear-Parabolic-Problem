#include "functions.h"

void splitMesh(grid &g, elementsMesh &eM, vector <int> &Ix)
{
   auto &x = g.x, &t = g.t, &xM = eM.xM, &tM = eM.tM;
   const int xMSize = xM.size() - 1, tMSize = tM.size() - 1;
   Ix.resize(xMSize + 1);
   std::ifstream grid("grid.txt");

   int nXk = 0;
   x.resize(1, xM[0]);
   for (int i = 0, j = 1; i < xMSize; i++, j++)
   {
      int countInterval = 0;
      double q = 0, step = 0;
      grid >> countInterval >> q;
      nXk += countInterval;
      x.resize(nXk + 1);

      if (q != 1)
      {
         double sumProgression = (pow(q, countInterval) - 1.0) / (q - 1.0);
         step = (xM[i + 1] - xM[i]) / sumProgression;
         int jk = 1;
         for (j; j < nXk; j++, jk++)
            x[j] = xM[i] + step * (pow(q, jk) - 1.0) / (q - 1.0);
      }
      else
      {
         step = (xM[i + 1] - xM[i]) / countInterval;
         int jk = 1;
         for (j; j < nXk; j++, jk++)
            x[j] = xM[i] + step * jk;
      }
      x[j] = xM[i + 1];
      Ix[i + 1] = j;
   }

   int nTk = 0;
   t.resize(1, tM[0]);
   for (int i = 0, j = 1; i < tMSize; i++, j++)
   {
      int countInterval = 0;
      double q = 0, step = 0;
      grid >> countInterval >> q;
      nTk += countInterval;
      t.resize(nTk + 1);

      if (q != 1)
      {
         double sumProgression = (pow(q, countInterval) - 1.0) / (q - 1.0);
         step = (tM[i + 1] - tM[i]) / sumProgression;
         int jk = 1;
         for (j; j < nTk; j++, jk++)
            t[j] = tM[i] + step * (pow(q, jk) - 1.0) / (q - 1.0);
      }
      else
      {
         step = (tM[i + 1] - tM[i]) / countInterval;
         int jk = 1;
         for (j; j < nTk; j++, jk++)
            t[j] = tM[i] + step * jk;
      }
      t[j] = tM[i + 1];
   }
   grid.close();
}

int returnNumberSubAreaOfMesh(vector <int> &Ix, vector <vector<int>> &M, int p, int &l)
{
   const int L = M.size(), L1 = l;

   for (l; l < L; l++)
   {
      int mx0 = Ix[M[l][1]], mx1 = Ix[M[l][2]], m = M[l][0];
      if (mx0 <= p && p <= mx1 && mx0 <= (p + 1) && (p + 1) <= mx1)
         return m;
   }

   for (l = 0; l < L1; l++)
   {
      int mx0 = Ix[M[l][1]], mx1 = Ix[M[l][2]], m = M[l][0];
      if (mx0 <= p && p <= mx1 && mx0 <= (p + 1) && (p + 1) <= mx1)
         return m;
   }
   return -1;
}