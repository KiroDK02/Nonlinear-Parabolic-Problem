#include "functions.h"

   std::function<double(double, double)> uAn = [](double x, double t) { return x; };

int main()
{
   elementsMesh eM{ };
   vector <double> sigma{ };
   vector <vector <int>> bC{ };
   grid g{ };
   vector <int> Ix{ };
   matrix A{ };
   vector <vector <double>> q{ };
   int flag = -1;

   readMesh(eM, sigma, bC, flag);
   splitMesh(g, eM, Ix);
   createPortrait(A, g);
   begApproximation(q, g);
   
   vector <int> countIter{ };
   switch (flag)
   {
      case 1:
      {
         countIter = calcSolution(A, g, q, bC, sigma, Ix, eM);
         break;
      }
      case 2:
      {
         countIter = calcSolutionMN(A, g, q, bC, sigma, Ix, eM);
         break;
      }
      default:
         break;
   }

   vector <vector <double>> qAn{ };
   const int xSize = g.x.size(), tSize = g.t.size();
   qAn.resize(tSize);
   auto &x = g.x, &t = g.t;

   for (int ti = 0; ti < tSize; ti++)
   {
      qAn[ti].resize(xSize);
      for (int p = 0; p < xSize; p++)
         qAn[ti][p] = uAn(x[p], t[ti]);
   }
   
   FILE *fout = NULL;
   fopen_s(&fout, "out.txt", "w");
   for (int ti = 0; ti < tSize; ti += 1)
   {
      double normL2 = sqrt(NIGaussThree(t[ti], g, q[ti]));
      if (ti == 0)
#pragma warning (disable: 6387)
         fprintf_s(fout, "Time layer %d: t = %.5lf, countNodes = %d ||u* - u||L2 = %.4e\n", ti, t[ti], xSize, normL2);
      else
      {
         double ht = t[ti] - t[ti - 1];
         fprintf_s(fout, "Time layer %d: t = %.5lf, countNodes = %d, ht = %.5lf, countIteration = %d, ||u* - u||L2 = %.4e\n", ti, t[ti], xSize, ht, countIter[ti], normL2);
      }
   }
   fclose(fout);
#pragma warning (default: 6387)
   return 0;
}

double calcValue(grid &g, double x, vector <double> &qNum)
{
   auto &gX = g.x;
   int begX = 0, endX = g.x.size() - 1;

   if (x < gX[0] || x > gX[endX])
   {
      printf_s("The point does not belong to the area.");
      return -1;
   }

   while (!(gX[begX] <= x && x <= gX[begX + 1]))
   {
      int indX = (begX + endX) / 2;
      if (gX[indX] < x)
         begX = indX;
      else
         endX = indX;
   }

   vector <int> globalNums = { begX, begX + 1 };
   double x0 = gX[begX], x1 = gX[begX + 1], hx = x1 - x0;
   double valueNumerical = qNum[begX] * (x1 - x) / hx + 
                           qNum[begX + 1] * (x - x0) / hx;
   return valueNumerical;
}

double NIGaussThree(double tValue, grid &g, vector <double> &qNum)
{
   auto &x = g.x;
   const int xSize = x.size();
   double integralValue = 0;

   for (int p = 0; p < xSize - 1; p += 8)
   {
      double hx = x[p + 1] - x[p], mid = (x[p + 1] + x[p]) / 2.0,
             gaussNode = sqrt(0.6);
      double x1 = mid, x2 = mid + gaussNode * hx / 2.0, x3 = mid - gaussNode * hx / 2.0;
      double uNum1 = calcValue(g, x1, qNum), uNum2 = calcValue(g, x2, qNum), uNum3 = calcValue(g, x3, qNum);
      double value1 = uAn(x1, tValue) - uNum1, value2 = uAn(x2, tValue) - uNum2, value3 = uAn(x3, tValue) - uNum3;
      integralValue += hx * (8.0 * value1 * value1 + 5.0 * (value2 * value2 + value3 * value3));
   }
   integralValue /= 18.0;
   return integralValue;
}