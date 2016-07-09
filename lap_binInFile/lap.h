/************************************************************************
*
*  lap.h
   version 1.0 - 21 june 1996
   author  Roy Jonker, MagicLogic Optimization Inc.
   
   header file for LAP
*
**************************************************************************/

/*************** CONSTANTS  *******************/

  #define BIG 100000

/*************** TYPES      *******************/
/* Changed cost to double  */
  typedef int row;
  typedef int col;
  typedef double cost;

/*************** FUNCTIONS  *******************/
/* Changed assigncost to double  
   Changed u and v to double  */
extern float lap(int dim, double **assigncost,
               int *rowsol, int *colsol, double *u, double *v);
/* Changed assigncost to double  
   Changed u and v to double  */
extern void checklap(int dim, double **assigncost,
                     int *rowsol, int *colsol, double *u, double *v);

