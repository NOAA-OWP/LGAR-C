
#include "all.h"


/* find the corners of the rectangles to draw */

extern int  xcategorize(int xmin, int xmax, int ymax, int num_soil_layers, int *soil_type,
                        double *cum_layer_thickness_cm, struct soil_properties_ *soil_properties, 
                        int *topleftx, int *toplefty, int *fwidth, int *fheight)
                               //
                               //
{                              //                   
int i=0;                       //             <-----------------winw--------------->     
double z,wc;                   //             ^         xorg                   xmax
int layer,soil;                //             |           |         xline        |
double the,thr;                //             |   yorg   >----------------------->
int x,y,w,h;                   //             |           |
struct wetting_front *current; //             |           |     x,y
                               //             |           |      *--------| -  
                               //             w         y |      |        | ^
current = head;                //             i         l |      |        | |
                               //             n         i |      |        | |      
                               //             h         n |      |        | h
                               //             |         e |      |        | |
// NOTES:                                     |           |      |        | |
// xmin = smallest theta_r in all layers      |           |      |        | | 
// xmax = largest theta_e in all layers       |           |      |        | v
// ymin = z=0                                 |           |      ---------- -
// ymax = total depth of soil                 |   ymax -  V      
// x,y,h,w are in units of pixels (int)       V                  |<--w--->|
// 


do    // loop through the fronts
  {
  layer = current->layer_num;
  soil = soil_type[layer];
  the = soil_properties[soil].theta_e;
  thr = soil_properties[soil].theta_r;

  // note 1: we can draw all rectangles starting from x=theta_r.
  x = (int)((double)xline*(thr-xmin)/(the-xmin));

  // note 2: all fronts are in contact with the top of the layer they live in, so:
  y = (int)((double)yline*(cum_layer_thickness_cm[layer-1]/ymax));  

  // note 3: iff we use theta_r as the left side, then the right side x coord is common ffor all rects. too.
  w = (int)((double)xline*(current->theta-thr)/(the-xmin));

  if(current->to_bottom == TRUE)  // this front touches the bottom of the current layer
    {                             // so its height is the number of pixels that spans that layer

    h = (int)((double)yline*(cum_layer_thickness_cm[layer] - cum_layer_thickness_cm[layer-1])/
                             cum_layer_thickness_cm[num_soil_layers] + cum_layer_thickness_cm[layer-1]);
    }
  else  // front is connected to top of layer it lives in, but not to the bottom.
    {
    h = (int)((double)xline*(current->depth_cm - cum_layer_thickness_cm[layer-1])/
                                                                  cum_layer_thickness_cm[num_soil_layers]);
    }
    
  topleftx[i] = x + xorg;
  toplefty[i] = y + yorg; 
  fwidth  [i] = w;
  fheight [i] = h; 
  i++;
  if(i > MAX_NUM_WETTING_FRONTS)
    {
    fprintf(stdout,"More wetting fronts than array elements in xcategorize().\n");
    fprintf(stdout,"Increase MAX_NUM_WETTING_FRONTS.  Program stopped.\n");
    exit(-1);
    }
  current=current->next;
  } while(current != NULL);
  
return(i);  // this returns the number of rectangles to draw
}
