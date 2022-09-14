#ifndef XSHOW_H_INCLUDED
#define XSHOW_H_INCLUDED

#include "../include/X_stuff.hxx"

extern void xshow()
{
   XCopyArea(theDisplay,pmap1,myWin,image_gc,0,0,winw,winh,0,0);
   XFlush(theDisplay);  
   return;
}

#endif
