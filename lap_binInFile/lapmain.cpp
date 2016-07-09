/************************************************************************
*
*  lapmain.cpp
   version 1.0 - 4 September 1996
   author: Roy Jonker @ MagicLogic Optimization Inc.
		   Adrian M. Peter - modified front-end interface to take in
		   input file.
   main program file to run and check Jonker-Volgenant LAP code
*
*************************************************************************/

#include <fstream>
#include "system.h"
#include "gnrl.h"
#include "lap.h"

using namespace std;
  
int main(int argc, char *argv[])
{
  #define COSTRANGE 1000.0
  #define PRINTCOST 0

  int dim=0;
  int startdim=0;
  int enddim=0;
  char dimBuff[4];
  char costBuff[8];
  cost **assigncost, *u, *v, lapcost;
  row i, *colsol;
  col j, *rowsol;
  double runtime;

  //cout << "D: " << sizeof(double) 
//	  << "I: " << sizeof(long int)
//	  << "F: " << sizeof(float) << endl;

  // Open file for input.
  ifstream costMatrixFile(argv[1],ios::in|ios::binary);
  if (!costMatrixFile)
  {
        cerr << "Opening " << argv[1] << " failed" << endl;
		exit(-1);
  }
  

  // Read in cost matrix dimension
  //costMatrixFile >> startdim;
  //costMatrixFile >> enddim;
  costMatrixFile.read(dimBuff,sizeof(int));
  startdim = *(reinterpret_cast<int*>(dimBuff));
  costMatrixFile.read(dimBuff,sizeof(int));
  enddim   = *(reinterpret_cast<int*>(dimBuff));
  

  // Allocate memory for cost matrix.
  assigncost = new cost*[enddim];
  for (i = 0; i < enddim; i++)
    assigncost[i] = new cost[enddim];

  // The row solution treats the x vector in (xy^T) as the reference 
  // and list the indicies that each element of the x vector maps to in the
  // y vector.  For example suppose x=[1 0 0] and y=[0 0 1].  At the end of
  // the algorithm, rowsol=[2 0 1].  This means that the first element of x
  // correspondes to the 3 element of y (need to +1 since rowsol values are in 
  // C++ array indicies), 2nd element maps to the 1st and last element maps to
  // the 2nd element.
  rowsol = new col[enddim];
  colsol = new row[enddim];
  u = new cost[enddim];
  v = new cost[enddim];
  // Notice that this for loop only runs once since in our case 
  // startdim is equal to enddim.
  for (dim = startdim; dim <= enddim; dim++)
  {
    // Read in cost matrix from file.
    for (j = 0; j < dim; j++)
		for (i = 0; i < dim; i++)
		{
			costMatrixFile.read(costBuff,sizeof(double));
			assigncost[i][j]=*(reinterpret_cast<double*>(costBuff));
			//costMatrixFile >> assigncost[i][j];
		}

#if (PRINTCOST) 
	double min=0, max=0;
    for (i = 0; i < dim; i++)
    {
      printf("\n");
      for (j = 0; j < dim; j++)
	  {
        //printf("%8.12f \n", assigncost[i][j]);
		if(assigncost[i][j] < min)
			min = assigncost[i][j];
		if(assigncost[i][j] > max)
			max = assigncost[i][j];
	  }
    }
	cout << "Max value: " << max << endl 
		<< "Min value: " << min << endl;
#endif
    
    printf("\nstart\n");
    runtime = seconds();
    lapcost = lap(dim, assigncost, rowsol, colsol, u, v);
    runtime = seconds() - runtime;
    printf("\n\ndim  %4d - lap cost %5d - runtime %6.3f\n", dim, lapcost, runtime);
  
    checklap(dim, assigncost, rowsol, colsol, u, v);

#if (PRINTCOST) 
	cout << "Row sol:\n";
	for(i = 0; i < enddim; ++i)
		printf("%4d ", rowsol[i]);
	cout << "\nCol sol: \n";
	for(j = 0; j < enddim; ++j)
		printf("%4d ", colsol[j]);
	cout << "\nu array:\n";
	for(i = 0; i < enddim; ++i)
		printf("%4d ", u[i]);
	cout << "\nv array: \n";
	for(j = 0; j < enddim; ++j)
		printf("%4d ", v[j]);
#endif
  }
  // Open file for output.
  ofstream assignSolFile(argv[2]);
  if (!assignSolFile)
  {
        cerr << "Opening " << argv[2] << " failed" << endl;
		exit(-1);
  }
  // Write out solution to linear assigment.
  // Currently only writing out row solution.
  for(i = 0; i < enddim; ++i)
		assignSolFile << rowsol[i] << endl;

  delete[] assigncost;
  delete[] rowsol;
  delete[] colsol;
  delete[] u;
  delete[] v;
  costMatrixFile.close();	
  assignSolFile.close();

  /*printf("\n\npress key\n");
  char c;
  cin >> c;*/
  return 0;
}

