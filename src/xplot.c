#include "X_stuff.h"


extern void   xplot(int shade,int num_fronts_to_plot, int xorg, int yorg, int *topleftx, int *toplefty,
                    int *fwidth, int *fheight, int refresh_flag)

{
   int j,dcnt;
   int black=0;
   int white=1;
   int blue=6;
   int green=5;
   int violet=7; 
   int x,y,w,h;

   dcnt=0;

   while(dcnt<10)  
      {       
      XEvent pe;
      XExposeEvent *ee;
      XConfigureEvent *ce;
      j++;

      if(refresh_flag==1 )
         {
         XFillRectangle(theDisplay,pmap1,gengc[white],0,0,winw,winh);
         }

      XDrawLine(theDisplay,pmap1,gengc[black],xorg,yorg+yline,xorg,yorg);
      XDrawLine(theDisplay,pmap1,gengc[black],xorg,yorg,xorg+xline,yorg);

   
      /*  Fill rectangles for all the wetting fronts */
      for(j=0;j<num_fronts_to_plot;j++)
        {
        x=topleftx[j];
        y=toplefty[j];
        w=fwidth[j];
        h=fheight[j];
        XFillRectangle(theDisplay,pmap1,gengc[blue],x,y,w,h);
        XFillRectangle(theDisplay,pmap1,gengc[blue],x,y,w,h);
        }
      xshow();
      dcnt++;
      } 
   return;
}
