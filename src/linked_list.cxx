#include "../include/all.hxx"

//#####################################################################################
/* authors : Fred Ogden and Ahmad Jan and Peter La Follette
   year    : 2022
   email   : ahmad.jan@noaa.gov
   - The file constains linked list functionality for the wetting front
   - Originally written by Fred Ogden (most of it), and modified/extended by Ahmad Jan and Peter La Follette*/
//#####################################################################################

/*#########################################################*/
/*#########################################################*/
/*#########################################################*/
//---------------------------------------------------------
//
// linked   LL           II     SSS      TTTTTTTTTTT
//          LL           II   SSS  SSS        TT
//          LL           II  SS      SS       TT
//          LL           II  SS               TT
//          LL           II   SSSS            TT
//          LL           II       SSS         TT
//          LL           II         SSS       TT
//          LL           II  SS      SS       TT
//          LL           II   SSS  SSS        TT
//          LLLLLLLLLL   II     SSSS          TT  functions
//___________________________________________________________


/*###########################################################*/
/* listDelete() - deletes memory allocated to a linked list  */
/* This function must be called on any list to deallocate    */
/* the dynamic memory used in creating and manipultating the */
/* list. (added by NJF)                                                     */
/*###########################################################*/
extern void listDelete(struct wetting_front* head)
{
  while (head != NULL) {
    struct wetting_front *next = head->next;
    free( head );
    head = next;
  }
}

/*#########################################################*/
/* listPrint() - prints a linked list to screen             */
/*#########################################################*/
extern void listPrint(struct wetting_front* head)
{
  struct wetting_front *current = head;
  printf("\n[ ");

  //start from the beginning
  while(current != NULL) {
    if (current->next == NULL)
      printf("(%lf,%6.14f,%d,%d,%d, %e, %lf, %6.14f, %d) ] \n",current->depth_cm, current->theta, current->layer_num,
	     current->front_num, current->to_bottom, current->dzdt_cm_per_h, current->K_cm_per_h, current->psi_cm, current->is_WF_GW);
    else
      printf("(%lf,%6.14f,%d,%d,%d, %e, %lf, %6.14f, %d)\n",current->depth_cm, current->theta, current->layer_num,
	     current->front_num, current->to_bottom, current->dzdt_cm_per_h, current->K_cm_per_h, current->psi_cm, current->is_WF_GW);
    current = current->next;
  }

}

/*###########################################################*/
/* listCopy() - copies a linked list to another linked list  */
/*###########################################################*/
extern struct wetting_front* listCopy(struct wetting_front* current, struct wetting_front* state_previous)
{
  if (current == NULL) {
    return NULL;
  }
  else {

    struct wetting_front* wf = (struct wetting_front*)malloc(sizeof(struct wetting_front));

    wf->depth_cm = current->depth_cm;
    wf->theta = current->theta;
    wf->psi_cm = current->psi_cm;
    wf->K_cm_per_h = current->K_cm_per_h;
    wf->layer_num = current->layer_num;
    wf->front_num = current->front_num;
    wf->to_bottom = current->to_bottom;
    wf->dzdt_cm_per_h = current->dzdt_cm_per_h;
    wf->is_WF_GW = current->is_WF_GW;
    wf->next = listCopy(current->next, NULL);

    if (state_previous == NULL)
      state_previous = wf;
    return wf;
  }
}


/*#######################################################*/
/* listInsertFirst - adds a list entry to start of list  */
/*#######################################################*/
extern void listInsertFirst(double depth, double theta, int front_num, int layer_num, bool bottom_flag, struct wetting_front** head, bool is_WF_GW)
{

  //create a link
  struct wetting_front *link = (struct wetting_front*) malloc(sizeof(struct wetting_front));

  link->depth_cm = depth;
  link->theta = theta;
  link->front_num = front_num;
  link->layer_num = layer_num;
  link->to_bottom = bottom_flag;
  link->dzdt_cm_per_h = (double)0.0;
  link->is_WF_GW = is_WF_GW; 
  
  //point it to old first wetting_front
  link->next = *head;

 
  //point head to new first wetting_front  FORGET THIS CODE WON'T WORK
  *head = link;
  link = (*head)->next;
  
  while (link != NULL) {
    link->front_num++;
    link=link->next;
  }

}


/*#######################################################*/
/* listDeleteFirst -deletes the first entry of a linked list */
/*#######################################################*/
extern struct wetting_front* listDeleteFirst(struct wetting_front* head)
{

  //save reference to first link
  struct wetting_front *tempLink = head;

  //mark next to first link as first
  head = head->next;

  //return the deleted link
  return tempLink;
}



/*#######################################################*/
/* listIsEmpty - checks to see if the list is empty          */
/*#######################################################*/
bool listIsEmpty(struct wetting_front* head)
{
  return head == NULL;
}



/*#######################################################*/
/* listLength - counts how many items are in the list    */
/*#######################################################*/
int listLength(struct wetting_front* head)
{
  int listLength = 0;
  struct wetting_front *current;

  for(current = head; current != NULL; current = current->next) {
    listLength++;
  }

  return listLength;
}

/*#######################################################*/
/* listLength_surface - counts how many surface WFs are in the list    */
/*#######################################################*/
int listLength_surface(struct wetting_front* head)
{
  int listLength = 0;
  struct wetting_front *current;

  for(current = head; current != NULL; current = current->next) {
    if (!(current->is_WF_GW)){
      listLength++;
      }
  }
  return listLength;
}

/*#######################################################*/
/* listLength_TO_WFs_above_surface_WFs - counts how many TO WFs are above surface WFs    */
/*#######################################################*/
int listLength_TO_WFs_above_surface_WFs(struct wetting_front* head)
{
  int listLength = 0;
  struct wetting_front *current;

  if (listLength_surface(head)==0){
    listLength = 0;
  }
  else{
    for(current = head; current != NULL; current = current->next) {
      if (current->is_WF_GW){
        listLength++;
        }
      else{
        break;
      }
    }
  }
  return listLength;
}

int listLength_TO_WFs_in_rz(double rzd, struct wetting_front* head)
{
  int listLength = 0;
  struct wetting_front *current;

  for(current = head; current != NULL; current = current->next) {
    if ( (current->is_WF_GW) && (current->to_bottom==0) && (current->depth_cm<rzd) ){
      listLength++;
      }
  }
  return listLength;
}

int listLength_TO_WFs_in_rz_nonzero_depth(double rzd, struct wetting_front* head)
{
  int listLength = 0;
  struct wetting_front *current;

  for(current = head; current != NULL; current = current->next) {
    if ( (current->is_WF_GW) && (current->to_bottom==0) && (current->depth_cm!=0.0) && (current->depth_cm<rzd) ){
      listLength++;
      }
  }
  return listLength;
}


/*#######################################################################################*/
/* listFindFront -searches linked list for a particular front number and returns that struct */
/*#######################################################################################*/
extern struct wetting_front* listFindFront(int key, struct wetting_front* head, struct wetting_front* head_old)
{

  //start from the first link
  struct wetting_front* current=NULL;
  if (head_old == NULL)
    current = head;
  else
    current = head_old;

  //iff list is empty
  if (head == NULL) {
    return NULL;
  }

  //list not empty, so search through list
  while(current->front_num != key) { // as soon as front_num==key, this loop stops
    //if it is last wetting_front
    if(current->next == NULL) {
      return NULL;
    }
    else {
      //go to next link
      current = current->next;
    }
  }

  //if data found, return the current Link
  return current;
}


/*##############################################################*/
/* listDeleteFront -delete the front with a particular front number */
/*##############################################################*/
extern struct wetting_front* listDeleteFront(int front_num, struct wetting_front** head, int *soil_type, struct soil_properties_ *soil_properties)
{
  //start from the first link
  struct wetting_front* current = *head;
  struct wetting_front* previous = NULL;

  //if list is empty
  // debugging
  bool front_found = false;
  while(current != NULL) {
    if (current->front_num == front_num) {
      front_found = true;
    }
    current = current->next;
  }

  if (!front_found) abort();

  current = *head;
  if(*head == NULL) {
    return NULL; // can't do anything, there is no list.
  }

  //navigate through list
  while(current->front_num != front_num) {
    //if it is last wetting_front
    if(current->next == NULL) {
      return NULL;
    }
    else {
      //store reference to current link
      previous = current;
      //move to next link
      current = current->next;
    }
  }

  //found a match, update the link
  if(current == *head) {
    //change first to point to next link
    *head = (*head)->next;
    previous = *head;
  }
  else {
    //bypass the current link
    previous->next = current->next;
    previous = current->next;

  }
  
  current = previous;

  while(previous != NULL) { // decrement all front numbers
    previous->front_num--;
    previous = previous->next;
  }

  int front_num_return = current->front_num;

  if (front_num!=1){
    for (int wf = listLength(*head)-1; wf != 0; wf--) {
      struct wetting_front *current_temp = listFindFront(wf, *head, NULL);
      struct wetting_front *next_temp = current_temp->next;
      if ( (current_temp->to_bottom==TRUE) ){
        current_temp->is_WF_GW = next_temp->is_WF_GW;
        current_temp->psi_cm = next_temp->psi_cm;

        int soil_num_k1 = soil_type[current_temp->layer_num]; 
        double theta_e_k   = soil_properties[soil_num_k1].theta_e;
        double theta_r_k   = soil_properties[soil_num_k1].theta_r;
        double vg_a_k      = soil_properties[soil_num_k1].vg_alpha_per_cm;
        double vg_m_k      = soil_properties[soil_num_k1].vg_m;
        double vg_n_k      = soil_properties[soil_num_k1].vg_n;
        current_temp->theta = calc_theta_from_h(current_temp->psi_cm, vg_a_k, vg_m_k, vg_n_k,theta_e_k,theta_r_k);
      }
    }
  }

  current = listFindFront(front_num_return, *head, NULL);

  return current;
}


/*####################################################################*/
/* listInsertFront -creates a new front at specified position in list */
/* and increase all front numbers greater than or equal to the        */
/* new front number by 1.                                             */
/*####################################################################*/
extern struct wetting_front* listInsertFront(double depth, double theta, int new_front_num,
                                             int layer_num, bool bottom_flag, struct wetting_front** head, bool is_WF_GW)
{
  //start from the first link
  struct wetting_front* current = NULL;
  struct wetting_front* previous = NULL;

  if(*head == NULL) { // list is empty
    
    if(new_front_num==1) { // create it
      //create a link
      struct wetting_front *link = (struct wetting_front*) malloc(sizeof(struct wetting_front));

      link->depth_cm = depth;
      link->theta = theta;
      link->front_num = new_front_num;
      link->layer_num = layer_num;
      link->to_bottom = bottom_flag;
      link->dzdt_cm_per_h = (double)(0.0);
      link->is_WF_GW = is_WF_GW;
      link->next = NULL;
      *head=link;
      return link;
    }
    else { // list is empty and desired link number is >1, unable to create
      return NULL;
    }
  }
  
  // list is not empty, so search through list
  previous = *head;
  do {
    if (previous->front_num == new_front_num-1) { // this is where we want to insert it
      //create a new link
      struct wetting_front *link = (struct wetting_front*) malloc(sizeof(struct wetting_front));

      link->depth_cm = depth;
      link->theta = theta;
      link->front_num = new_front_num;
      link->layer_num = layer_num;
      link->to_bottom = bottom_flag;
      link->dzdt_cm_per_h = (double)(0.0);
      link->is_WF_GW = is_WF_GW;
      link->next = previous->next ;
      previous->next = link;

  //     while(link->next != NULL) { // increment all front numbers
	// current = link->next;
	// current->front_num++;
	// previous = current;
	// current = current->next;
  //     }
      current = *head;
      for (int wf = 1; wf != (listLength(*head)+1); wf++){
        current->front_num = wf;
        current = current->next;
      }

      return link;
    }

    previous = previous->next;
  } while(TRUE);

}

/*############################################################################*/
/* listInsertFrontAtDepth -creates a new front at specified depth in the soil */
/* and increase all front numbers greater than or equal to the                */
/* new front number by 1. Determines layer number, front number and           */
/* determine if the new front is at the bottom of the layer                   */
/*############################################################################*/
extern struct wetting_front* listInsertFrontAtDepth(int num_layers, double *cum_layer_thickness,
                                                    double depth, double theta, struct wetting_front* head, bool is_WF_GW)//this fxn has not been completely updated to work with LGARTO but is no longer used in LGAR.cxx
{
  int el_layer=0;
  bool extends_to_bottom_flag = FALSE;
  bool layer_found_flag = FALSE;

  struct wetting_front *link = NULL;

  //start from the first link
  struct wetting_front* current = head;
  struct wetting_front* previous = NULL;

  if(head == NULL) { // list is empty
    // Kinda weird.  Shouldn't call this function to start a list.  Create link in the first position

    link = (struct wetting_front*) malloc(sizeof(struct wetting_front));

    link->depth_cm = depth;
    link->theta = theta;
    link->to_bottom = FALSE;

    layer_found_flag = listFindLayer(link,num_layers,cum_layer_thickness, &el_layer, &extends_to_bottom_flag);

    if(layer_found_flag ==FALSE) {
      return NULL; // this should never happen, unless it asks to create a front not in the soil
    }

    link->layer_num = el_layer;
    link->to_bottom = extends_to_bottom_flag;
    link->dzdt_cm_per_h = (double)(-1.0);
    link->is_WF_GW = is_WF_GW;
    link->next = NULL;   // because this is the first link.
    head = link;  // don't forget this.  Must keep head current with the beginning of the list.
    return link;
  }
  else {
    // list is not empty, so search through list fo find where this front fits in.  Might be first though...
    if(depth < current->depth_cm) {
      // yep.  It's first
      //create a link and put it at the beginning of the list

      link = (struct wetting_front*) malloc(sizeof(struct wetting_front));

      link->depth_cm = depth;
      link->theta = theta;
      link->to_bottom = FALSE;

      layer_found_flag = listFindLayer(link,num_layers,cum_layer_thickness, &el_layer, &extends_to_bottom_flag);

      if(layer_found_flag ==FALSE) {
	return NULL; // this should never happen, unless it asks to create a front not in the soil
      }

      link->layer_num = el_layer;
      link->to_bottom = extends_to_bottom_flag;
      link->dzdt_cm_per_h = (double)(-1.0);
      link->is_WF_GW = is_WF_GW;
      link->next = head;   // point to the old first list entry
      head=link;  // don't forget this.  Must keep head point to the beginning of the list
      return link;
    }
    else {
      // the new one is not the first wetting front in the list sorted by depth, search to see where it fits.

      previous=current;
      current=current->next;
      do {
	if((depth > previous->depth_cm ) && (depth <= current->depth_cm) ) {
	  // it is in this interval

	  // create a link
	  link = (struct wetting_front*) malloc(sizeof(struct wetting_front));

	  link->depth_cm = depth;
	  link->theta = theta;
	  link->to_bottom = FALSE;

	  layer_found_flag = listFindLayer(link,num_layers,cum_layer_thickness,
					   &el_layer, &extends_to_bottom_flag);

	  if(layer_found_flag ==FALSE) {
	    return NULL; // this should never happen, unless it asks to create a front not in the soil
	  }

	  link->layer_num = el_layer;
	  link->to_bottom = extends_to_bottom_flag;
	  link->dzdt_cm_per_h = (double)(-1.0);
    link->is_WF_GW = is_WF_GW;
	  link->front_num = current->front_num;
	  link->next = current;   // point to one replaced
	  previous->next = link;  // make the previous one point to this new one
	  current = previous->next;
	  break;
	}

	previous=current;
	current=previous->next;
      } while ( current != NULL );

      if (layer_found_flag == TRUE && current->next != NULL) {
	//re-number the links that follow the current link
	current=current->next;

	while (current != NULL) { // increment all front numbers
	  //printf("front number: %3d, current->next=%p\n",current->front_num,current->next);
	  current->front_num++;
	  previous = current;
	  current = current->next;
        }

	return link;
      }
      else {
	return NULL;
      }
    }
  }

}


/*##############################################################*/
/* listFindLayer -find what layer a newly created link lives in */
/*###############################################################*/
extern bool listFindLayer( struct wetting_front* link, int num_layers, double *cum_layer_thickness_cm,
                           int *lives_in_layer,bool *extends_to_bottom_flag)
{

  int layer;
  double depth;

  (*lives_in_layer)=0;
  depth=link->depth_cm;
  for (layer=1; layer<=num_layers; layer++) {
    if(depth <= cum_layer_thickness_cm[layer] && depth > cum_layer_thickness_cm[layer-1]) {
      (*lives_in_layer)  = layer;
      if ( is_epsilon_less_than (depth-cum_layer_thickness_cm[layer],1.0E-05) )
	(*extends_to_bottom_flag)=TRUE;
      break;
    }
  }
  if ( (*lives_in_layer) == 0)
    return FALSE;
  else
    return TRUE;
}

/*#################################################################################################*/
/* listSortFrontsByDepth -if fronts get out of order, this routine sorts them back into order by depth */
/*#################################################################################################*/
extern void listSortFrontsByDepth(struct wetting_front *head)
{
  int i, j, k, tempKey;
  double tempData;
  struct wetting_front *current;
  struct wetting_front *next;

  int size = listLength(head);
  k = size ;

  for ( i = 0 ; i < size - 1 ; i++, k-- ) {
    current = head;
    next = head->next;

    for ( j = 1 ; j < k ; j++ ) {
      if ( current->depth_cm > next->depth_cm ) {
	tempData = current->depth_cm;
	current->depth_cm = next->depth_cm;
	next->depth_cm = tempData;

	tempData = current->theta;
	current->theta = next->theta;
	next->theta = tempData;

	tempKey = current->layer_num;
	current->layer_num = next->layer_num;
	next->layer_num = tempKey;

  tempKey = current->to_bottom;
	current->to_bottom = next->to_bottom;
	next->to_bottom = tempKey;

	// tempKey = current->front_num;
	// current->front_num = next->front_num;
	// next->front_num = tempKey;

  tempData = current->psi_cm;
	current->psi_cm = next->psi_cm;
	next->psi_cm = tempData;

  tempData = current->K_cm_per_h;
	current->K_cm_per_h = next->K_cm_per_h;
	next->K_cm_per_h = tempData;

  tempData = current->dzdt_cm_per_h;
	current->dzdt_cm_per_h = next->dzdt_cm_per_h;
	next->dzdt_cm_per_h = tempData;

  tempKey = current->is_WF_GW;
	current->is_WF_GW = next->is_WF_GW;
	next->is_WF_GW = tempKey;

      }

      current = current->next;
      next = next->next;
    }
  }

}

/*#################################################################################################*/
/* listSortFrontsByDepth -if fronts get out of order in terms of capillary head, this routine sorts them back into order by depth */
/*#################################################################################################*/
extern void listSortFrontsByPsi(struct wetting_front *head)
{
  int i, j, k, tempKey;
  double tempData;
  struct wetting_front *current;
  struct wetting_front *next;

  int size = listLength(head);
  k = size ;

  for ( i = 0 ; i < size - 1 ; i++, k-- ) {
    current = head;
    next = head->next;

    for ( j = 1 ; j < k ; j++ ) {
      if ( current->depth_cm < next->depth_cm && listLength_surface(head)==1) {
	tempData = current->depth_cm;
	current->depth_cm = next->depth_cm;
	next->depth_cm = tempData;

	tempData = current->theta;
	current->theta = next->theta;
	next->theta = tempData;

	tempKey = current->layer_num;
	current->layer_num = next->layer_num;
	next->layer_num = tempKey;

  tempKey = current->to_bottom;
	current->to_bottom = next->to_bottom;
	next->to_bottom = tempKey;

	// tempKey = current->front_num;
	// current->front_num = next->front_num;
	// next->front_num = tempKey;

  tempData = current->psi_cm;
	current->psi_cm = next->psi_cm;
	next->psi_cm = tempData;

  tempData = current->K_cm_per_h;
	current->K_cm_per_h = next->K_cm_per_h;
	next->K_cm_per_h = tempData;

  tempData = current->dzdt_cm_per_h;
	current->dzdt_cm_per_h = next->dzdt_cm_per_h;
	next->dzdt_cm_per_h = tempData;

  tempKey = current->is_WF_GW;
	current->is_WF_GW = next->is_WF_GW;
	next->is_WF_GW = tempKey;

      }

      current = current->next;
      next = next->next;
    }
  }

}


/*#################################################################################################*/
/* listSendToTop - Sometimes TO wetting fronts will have their depth set to 0. This function moves all of these wetting fronts to the top of the list. */
/*#################################################################################################*/
extern void listSendToTop(struct wetting_front *head)
{
  int i, j, k, tempKey;
  double tempData;
  struct wetting_front *current;
  struct wetting_front *next;

  int size = listLength(head);
  k = size ;

  for ( i = 0 ; i < size - 1 ; i++, k-- ) {
    current = head;
    next = head->next;

    for ( j = 1 ; j < k ; j++ ) {
      if ( (current->depth_cm > next->depth_cm) && (next->depth_cm == 0) ) {
	tempData = current->depth_cm;
	current->depth_cm = next->depth_cm;
	next->depth_cm = tempData;

	tempData = current->theta;
	current->theta = next->theta;
	next->theta = tempData;

	tempKey = current->layer_num;
	current->layer_num = next->layer_num;
	next->layer_num = tempKey;

  tempKey = current->to_bottom;
	current->to_bottom = next->to_bottom;
	next->to_bottom = tempKey;

  tempData = current->psi_cm;
	current->psi_cm = next->psi_cm;
	next->psi_cm = tempData;

  tempData = current->K_cm_per_h;
	current->K_cm_per_h = next->K_cm_per_h;
	next->K_cm_per_h = tempData;

  tempData = current->dzdt_cm_per_h;
	current->dzdt_cm_per_h = next->dzdt_cm_per_h;
	next->dzdt_cm_per_h = tempData;

  tempKey = current->is_WF_GW;
	current->is_WF_GW = next->is_WF_GW;
	next->is_WF_GW = tempKey; 
      }

      current = current->next;
      next = next->next;
    }
  }

}

extern int lgarto_count_fronts_for_excessive_calc(struct wetting_front **head){
  int listLength = 0;
  struct wetting_front *current;

  for(current = *head; current != NULL; current = current->next) {
    if (current->is_WF_GW==0){
      break;
    }
    if (current->depth_cm!=0.0){
      break;
    }
    listLength++;
  }

  return listLength;
}


/*####################################################################*/
/* listReverseOrder  - reverses the order of the current linked list  */
/*####################################################################*/
// probably not needed, left in as an example of how it's done
extern void listReverseOrder(struct wetting_front** head_ref)
{
  printf("In Reverse order.... \n ");
  struct wetting_front* prev   = NULL;
  struct wetting_front* current = *head_ref;  // notice that it is passed down as a pointer to a pointer, that explains
  struct wetting_front* next;                 // why they call this "head_ref", it is a pointer referencing the
  // pointer struct wetting_front *head.
  while (current != NULL) {
    next  = current->next;
    current->next = prev;
    prev = current;
    current = next;
  }

  *head_ref = prev;
}

/*#################################################################################################*/
/* GW_fronts_among_surf_WFs - Sometimes surface or TO WFs move in such a way that a TO WFs is in between surface WFs, which is not technically possible.
This function returns true if this is the case so it can be corrected. */
/*#################################################################################################*/
extern bool GW_fronts_among_surf_WFs(struct wetting_front *head){
  struct wetting_front *current;
  current = listFindFront(listLength_TO_WFs_above_surface_WFs(head) + 1, head, NULL);
  int surf_count = 0;
  bool surface_WFs_deeper_than_top_mobile_TO_WF = false;
  for (int wf = 1; wf != (listLength(head)); wf++){
    if (current->is_WF_GW==FALSE){
      surf_count = surf_count + 1;
    }
    else{
      if (surf_count!=listLength_surface(head)){
        surface_WFs_deeper_than_top_mobile_TO_WF = true;
        break;
      }
    }
    current = current->next;
    if (current==NULL) {break;}
  }

  return(surface_WFs_deeper_than_top_mobile_TO_WF);
}
