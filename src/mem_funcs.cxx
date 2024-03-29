#include "../include/all.hxx"
#include <cstring>

/*#########################################################################*/
/*#########################################################################*/
/*#########################################################################*/
/***************************************************************************/
//
//   MMM     MMM    MMMMMMMM    MMM     MMM    OOO     RRRRRRR   YYYY   YYYY
//   MM M   M MM    MM          MM M   M MM  000 000   RR   RRR   YYY   YYY
//   MM  M M  MM    MM          MM  M M  MM OO     OO  TT    RR    YYY YYY 
//   MM   M   MM    MM          MM   M   MM OO     OO  TT  RRR       YYY
//   MM       MM    MM          MM       MM OO     OO  TTRRR         YYY
//   MM       MM    MMMMMM      MM       MM OO     OO  TT   R        YYY
//   MM       MM    MM          MM       MM OO     OO  TT    R       YYY
//   MM       MM    MM          MM       MM  OOO OOO   RR     R      YYY
//   MM       MM    MMMMMMMMMM  MM       MM    OOO     RR     RR    YYYYY  
//
//       allocation funcions, 1-D and 2-D integer and double_precision
//--------------------------------------------------------------------------


/*####################################################################*/
/* ALLOCATE MEMORY FOR A 2-D INTEGER ARRAY AND SET ITS VALUES TO ZERO */
/*####################################################################*/
void itwo_alloc(int ***array,int rows, int cols)
{
int  i,frows,fcols;
//int error=0;

if ((rows==0)||(cols==0))
  {
  printf("Error: Attempting to allocate array of size 0\n");
  exit(0);
  }

frows=rows+1;  /* added one for FORTRAN numbering */
fcols=cols+1;  /* added one for FORTRAN numbering */

*array=(int **)malloc(frows*sizeof(int *));
if (*array) 
  {
  memset((*array), 0, frows*sizeof(int*));
  for (i=0; i<frows; i++)
    {
    (*array)[i] =(int *)malloc(fcols*sizeof(int ));
    if ((*array)[i] == NULL)
      {
	//error = 1;
      i = frows;
      }
     else memset((*array)[i], 0, fcols*sizeof(int )); 
     }
   }
return;
}


/*#############################################################################*/
/* ALLOCATE MEMORY FOR A 2-D DOUBLE PRECISION ARRAY AND SET ITS VALUES TO ZERO */
/*#############################################################################*/
void dtwo_alloc(double ***array,int rows, int cols)
{
int  i,frows,fcols;
//int error=0;

if ((rows==0)||(cols==0))
  {
  printf("Error: Attempting to allocate array of size 0\n");
  exit(0);
  }

frows=rows+1;  /* added one for FORTRAN numbering */
fcols=cols+1;  /* added one for FORTRAN numbering */

*array=(double **)malloc(frows*sizeof(double *));
if (*array) 
  {
  memset((*array), 0, frows*sizeof(double *));
  for (i=0; i<frows; i++)
    {
    (*array)[i] =(double *)malloc(fcols*sizeof(double ));
    if ((*array)[i] == NULL)
      {
	//error = 1;
      i = frows;
      }
     else memset((*array)[i], 0, fcols*sizeof(double )); 
     }
   }
return;
}

/*#############################################################################*/
/* ALLOCATE MEMORY FOR A 1-D DOUBLE PRECISION ARRAY AND SET ITS VALUES TO ZERO */
/*#############################################################################*/
void d_alloc(double **var,int size)
{
  size++;  /* just for safety */

   *var = (double *)malloc(size * sizeof(double));
   if (*var == NULL)
      {
      printf("Problem allocating memory for array in d_alloc.\n");
      return;
      }
   else memset(*var,0,size*sizeof(double));
   return;
}

/*####################################################################*/
/* ALLOCATE MEMORY FOR A 1-D INTEGER ARRAY AND SET ITS VALUES TO ZERO */
/*#####################################################################*/
void i_alloc(int **var,int size)
{
   size++;  /* just for safety */

   *var = (int *)malloc(size * sizeof(int));
   if (*var == NULL)
      {
      printf("Problem allocating memory in i_alloc\n");
      return; 
      }
   else memset(*var,0,size*sizeof(int));
   return;
}

/*###########################################################################*/
/* ALLOCATE MEMORY FOR A 1-D FLOATING POINT ARRAY AND SET ITS VALUES TO ZERO */
/*###########################################################################*/
void f_alloc(float **var,int size)
{
  size++;  /* just for safety */

   *var = (float *)malloc(size * sizeof(float));
   if (*var == NULL)
      {
      printf("Problem allocating memory for array in f_alloc.\n");
      return;
      }
   else memset(*var,0,size*sizeof(float));
   return;
}
