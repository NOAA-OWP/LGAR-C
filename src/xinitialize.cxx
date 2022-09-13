#include "X_stuff.h"

extern int xinitialize(int winw, int winh)
{
   int i, j ;
   char *geom = NULL;
   int winx, winy;
   register int xdir, ydir;
   register int xoff, yoff;
   register int centerX, centerY;
   XGCValues xgcv;
   XSetWindowAttributes xswa;
   Colormap map;
   unsigned int bw = 1;
   char *display = NULL;
   Status status;
   char *fg = NULL;
   char *bg = NULL;
   char *bd = NULL;
   char bgname[12];
   int fg_pix, bg_pix, bd_pix;
   XColor fg_def, fg_exact, bg_def, bg_exact, bd_def, bd_exact;
   int bs = NotUseful;
   Visual visual;

/* Added Variable List */
   int k,igc,eofcheck;
   unsigned line_width= 0;
   int line_style=LineSolid;
   int cap_style=CapButt;
   int join_style=JoinMiter;
   int x_text_off=20;
   int y_text_off=4;
   int rec_wid=2,rec_widx,rec_widy;
   int line_draw=1;
   int numdat=1,numpts;
   float zdisc,zcol;
   int ii,numvals;
   float axscale;
   float maxx,maxy,minx,miny,xscale,yscale,zdiff,maxz,minz;
   int scalef;

   XSegment yax;
   XSegment xax;

   strcpy(gename[0],"black");
   strcpy(gename[1],"white");
   strcpy(gename[2],"red");
   strcpy(gename[3],"orange");
   strcpy(gename[4],"yellow");
   strcpy(gename[5],"green");
   strcpy(gename[6],"blue");
   strcpy(gename[7],"violet");
   strcpy(gename[8],"grey30");
   strcpy(gename[9],"grey70");


//   if(rescaleflag==1) xrescale();  /* catches the first one */
//   rescaleflag=0;  /* turns off rescaleflag after first one is caught */

   sprintf(xtext,"Infiltrated Depth");
   sprintf(ytext1,"Stage");
   sprintf(ytext2," (m)");
   btextx=400;
   btexty=50;

   yaxx=10;
   yaxy=475;
   xaxx=300;
   xaxy=630;


   /* lay out x and y axis tic marks */

   if (!(theDisplay = XOpenDisplay(display)))
      {
      perror("Cannot open display\n");
      exit(-1);
      }

   map = XDefaultColormap(theDisplay,DefaultScreen(theDisplay));
   if (fg) 
      {
      status = XAllocNamedColor(theDisplay, map, fg, &fg_def, &fg_exact);
      fg_pix = status ? fg_def.pixel : WhitePixel(theDisplay, 
					     DefaultScreen(theDisplay));
      } 
   else
      fg_pix = WhitePixel(theDisplay, DefaultScreen(theDisplay));

   strcpy(bgname,"grey40");
   bg=bgname;
   if (bg)  
      {
      status = XAllocNamedColor(theDisplay, map, bg, &bg_def, &bg_exact);
      bg_pix = status ? bg_def.pixel : BlackPixel(theDisplay, 
					     DefaultScreen(theDisplay));
      } 
   else
      bg_pix = BlackPixel(theDisplay, DefaultScreen(theDisplay));

   if (bd) 
      {
      status = XAllocNamedColor(theDisplay, map, bd, &bd_def, &bd_exact);
      bd_pix = status ? bd_def.pixel : WhitePixel(theDisplay, 
					     DefaultScreen(theDisplay));
      }
   else
      bd_pix = WhitePixel(theDisplay, DefaultScreen(theDisplay));



   xswa.backing_store = bs;
   xswa.event_mask = ExposureMask | StructureNotifyMask;
   xswa.background_pixel = bg_pix;
   xswa.border_pixel = bd_pix;
   visual.visualid = CopyFromParent;
   myWin = XCreateWindow(theDisplay,
		RootWindow(theDisplay, DefaultScreen(theDisplay)),
		winx, winy, winw, winh, bw, 
		DefaultDepth(theDisplay, DefaultScreen(theDisplay)), 
		InputOutput,&visual,CWEventMask | CWBackingStore | 
		CWBorderPixel | CWBackPixel, &xswa);

   XChangeProperty(theDisplay, myWin, XA_WM_NAME, XA_STRING, 8,
		PropModeReplace, (unsigned char *)"LGAR Infilt", 11);
   XMapWindow(theDisplay, myWin);                   

   /* create a Pixmap, pmap1 */
   pmap_depth=DefaultDepth(theDisplay,DefaultScreen(theDisplay));
   pmap1=XCreatePixmap(theDisplay,myWin,winw,winh,pmap_depth);
   image_gc=XCreateGC(theDisplay,myWin,0,0);
   strcpy(topstring,"Annual Rainfall (mm)");

   for(igc=0;igc<10;igc++)
      {
      status = XAllocNamedColor(theDisplay,map,gename[igc],&fg_def, &fg_exact);
      fg_pix = status ? fg_def.pixel : WhitePixel(theDisplay, 
			   DefaultScreen(theDisplay));
      xgcv.foreground = fg_pix;
      xgcv.function = GXcopy;    /* was GXinvert */
      xgcv.plane_mask = AllPlanes;
      xgcv.fill_style = FillSolid;
      xgcv.graphics_exposures = False; 
      gengc[igc] = XCreateGC(theDisplay, myWin,
		GCForeground | GCFunction | GCPlaneMask | GCFillStyle |
		GCGraphicsExposures, &xgcv);
      XSetLineAttributes(theDisplay,gengc[igc],line_width,line_style,cap_style,
   	        join_style);
      }

   /* fill the pixmap with the specified background color */
   XFillRectangle(theDisplay,pmap1,gengc[1],0,0,winw,winh);
   return(0);
}

