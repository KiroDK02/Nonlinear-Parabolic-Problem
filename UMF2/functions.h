#include <iostream>
#include <fstream>
#include <vector>
#include <functional>

using std::vector;

struct elementsMesh
{
   vector <double> xM{ }, tM{ };
   vector <vector <int>> M{ };
};

struct grid
{
   vector <double> x{ }, t{ };
};

struct matrix
{
   vector <vector<double>> ggl{ }, ggu{ };
   vector <double> di{ }, b{ };
};

void readMesh(elementsMesh &eM, vector <double> &sigma, vector <vector <int>> &bC, int &flag);
void splitMesh(grid &g, elementsMesh &eM, vector <int> &Ix);
int returnNumberSubAreaOfMesh(vector <int> &Ix, vector <vector<int>> &M, int p, int &l);
void createPortrait(matrix &A, grid &g);
void addElementInMatrix(matrix &A, double elem, int i, int j);
void begApproximation(vector <vector<double>> &q, grid &g);
vector <int> calcSolution(matrix &A, grid &g, vector<vector<double>> &q, vector<vector <int>> &bC, vector <double> &sigma, vector <int> &Ix, elementsMesh &eM);
int calcMPI(double tValue, double deltaT, matrix &A, grid &g, vector<double> &q, vector<double> &q0, vector<vector <int>> &bC, vector <double> &sigma, vector <int> &Ix, elementsMesh &eM);
void calcGlobalMatrix(double tValue, double deltaT, matrix &A, grid &g, vector<double> &q, vector<double> &qtPrev, vector <double> &sigma, vector <int> &Ix, elementsMesh &eM);
vector <int> calcSolutionMN(matrix &A, grid &g, vector<vector<double>> &q, vector<vector <int>> &bC, vector <double> &sigma, vector <int> &Ix, elementsMesh &eM);
int calcMN(double tValue, double deltaT, matrix &A, grid &g, vector<double> &q, vector<double> &q0, vector<vector <int>> &bC, vector <double> &sigma, vector <int> &Ix, elementsMesh &eM);
void calcGlobalMatrix¿ForMN(double v, double tValue, double deltaT, matrix &A, grid &g, vector<double> &q, vector<double> &qtPrev, vector <double> &sigma, vector <int> &Ix, elementsMesh &eM);
void boundaryCondition(double tValue, matrix &A, vector <vector <int>> &bC, grid &g);
void calcLU(matrix &A);
void calcY(matrix &A, vector <double> &y);
void calcQ(matrix &A, vector <double> &y, vector <double> &q);
double calcDiscrepancy(matrix &A, vector <double> &q);
void clearMatrix(matrix &A);
void vecEqualVec(vector <double> &q1, vector <double> &q2);
double calcValue(grid &g, double x, vector <double> &qNum);
double NIGaussThree(double tValue, grid &g, vector <double> &qNum);
void calcRelaxation(vector <double> &q0, vector <double> &q1, double w);
bool matrixPositiveDefinite(matrix &A);