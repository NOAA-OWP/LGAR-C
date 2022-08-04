#ifndef _X_STUFF_H

#define _X_STUFF_H

/******************INCLUDED HEADER FILES *******************************/
#include "X11/Xlib.h"
#include "X11/Xatom.h"
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
/******************END OF INCLUDED HEADER FILES ************************/

/******************MY DEFINES ******************************************/
/******************END OF MY DEFINES************************************/

/****** SUBS ****/
extern int rescale();
extern void xshow();
                  
/******************XLIB SPECIFIC NEEDS**********************************/
extern int errno;
char *ProgramName;
/* char *malloc();  Comment out on DEC alpha */
Window XCreateWindow();
Window myWin;
Display *theDisplay;
Pixmap pmap1;
GC gengc[10];
GC image_gc;
/*****************END OF XLIB SPECIFIC NEEDS****************************/


   /* needed for pixmap creation */
   unsigned int pmap_depth; /* Bits per pixel value */
 
   XSegment xtics[50];
   XSegment ytics[20];
   XPoint profile[300];
   XPoint bprofile[300];

  
   static int smalrectw;

   int x_textx[50],x_texty[50],y_textx[20],y_texty[20],textlen;
   int r_textx[17],r_texty[17];
   int yaxx,yaxy,xaxx,xaxy;
   char topstring[60];
   char btext1[19],btext2[19],btext3[19],btext4[19];
   int btextx,btexty;
   char xscale_text[50][12], yscale_text[20][12];
   char xtext[50],ytext1[16],ytext2[16];
   char gename[30][15];
   int bed1x,bed1y,bed2x,bed2y;


/* GLOBAL X-WINDOWS VARIABLES */
   int winw ; /* was 770   controls width of window */
   int winh ; /* was 455   controls height of window */
   int xorg;    /* controls x origin coord of hydrograph box */
   int yorg;   /* controls y origin coord of hydrograph box */
   int xline;  /* controls length of x axis in hydrograph box*/
   int yline;  /* controls length of y axis in hydrograph box */
   float yaxdiff; /* first guess of maximum depth in channel (m) */
   int ticlen;
   int ptextlen;
   int rescaleflag;


#endif  // _X_STUFF_H
