#include "X_stuff.h"

extern void xshow()
{
   XCopyArea(theDisplay,pmap1,myWin,image_gc,0,0,winw,winh,0,0);
   XFlush(theDisplay);  
   return;
}
