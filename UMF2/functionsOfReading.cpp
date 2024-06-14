#include "functions.h"

void readMesh(elementsMesh &eM, vector <double> &sigma, vector <vector <int>> &bC, int &flag)
{
   auto &xM = eM.xM, &tM = eM.tM;
   auto &M = eM.M;
   int xSize{ }, tSize{ }, countSubArea{ };

   std::ifstream mesh("mesh.txt");
   mesh >> xSize;
   xM.resize(xSize);
   for (int i = 0; i < xSize; i++)
      mesh >> xM[i];

   mesh >> tSize;
   tM.resize(tSize);
   for (int i = 0; i < tSize; i++)
      mesh >> tM[i];

   mesh >> countSubArea;
   M.resize(countSubArea);
   for (int element = 0; element < countSubArea; element++)
   {
      M[element].resize(3);
      mesh >> M[element][0];
      for (int i = 1; i < 3; i++)
      {
         int temp{ };
         mesh >> temp;
         M[element][i] = temp - 1;
      }
   }
   mesh >> flag;
   mesh.close();

   std::ifstream parametres("parametres.txt");
   sigma.resize(countSubArea);
   for (int i = 0; i < countSubArea; i++)
      parametres >> sigma[i];
   parametres.close();

   std::ifstream boundaryCondition("boundaryCondition.txt");
   bC.resize(2);
   for (int i = 0; i < 2; i++)
   {
      bC[i].resize(2);
      boundaryCondition >> bC[i][0] >> bC[i][1];
   }
   boundaryCondition.close();
}