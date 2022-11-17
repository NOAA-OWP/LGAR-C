#include "../include/all.hxx"

//#####################################################################################
/* authors : Fred Ogden and Ahmad Jan
   year    : 2022
   email   : ahmad.jan@noaa.gov
   - The file constains linked list functionality for the wetting front
   - Originally written by Fred Ogden (most of it), and modified/extended by Ahmad Jan */
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



/*#########################################################*/
/* listPrint() - prints a linked list to screen             */
/*#########################################################*/
extern void listPrint()
{
  struct wetting_front *current = head;
  printf("\n[ ");

  //start from the beginning
  while(current != NULL)
    {
      if (current->next == NULL)
	printf("(%lf,%6.14f,%d,%d,%d, %e, %lf %6.14f) ] \n",current->depth_cm, current->theta, current->layer_num,current->front_num, current->to_bottom, current->dzdt_cm_per_s, current->K_cm_per_s, current->psi_cm);
      else
	printf("(%lf,%6.14f,%d,%d,%d, %e, %lf %6.14f)\n",current->depth_cm, current->theta, current->layer_num,current->front_num, current->to_bottom, current->dzdt_cm_per_s, current->K_cm_per_s, current->psi_cm);
      current = current->next;
  }
  
}

/*###########################################################*/
/* listCopy() - copies a linked list to another linked list  */
/*###########################################################*/
extern struct wetting_front* listCopy(struct wetting_front* current)
{
  //  struct wetting_front *current = head;
  if (current == NULL) {
    return NULL;
  }
  else {

    struct wetting_front* wf = (struct wetting_front*)malloc(sizeof(struct wetting_front));
    
    wf->depth_cm = current->depth_cm;
    wf->theta = current->theta;
    wf->psi_cm = current->psi_cm;
    wf->K_cm_per_s = current->K_cm_per_s;
    wf->layer_num = current->layer_num;
    wf->front_num = current->front_num;
    wf->to_bottom = current->to_bottom;
    wf->dzdt_cm_per_s = current->dzdt_cm_per_s;
    wf->next = listCopy(current->next);

    if (state_previous == NULL) 
      state_previous = wf;
    return wf;
  }
}

extern void listPrintC(struct wetting_front pcurrent)
{
/*#########################################################*/
/* listPrint() -prints a linked list to screen             */
/*#########################################################*/
struct wetting_front *current = &pcurrent;
printf("\n[ ");

//start from the beginning
while(current != NULL)
  {
    printf("(%lf,%lf,%d,%d,%d, %e, %lf %lf)\n",current->depth_cm, current->theta, current->layer_num,current->front_num, 
	 current->to_bottom, current->dzdt_cm_per_s, current->K_cm_per_s, current->psi_cm);

  current = current->next;
  }
printf(" ]\n");
}


/*#######################################################*/
/* listInsertFirst - adds a list entry to start of list  */
/*#######################################################*/
extern void listInsertFirst(double depth, double theta, int front_num, int layer_num, bool bottom_flag)
{

  //create a link
  struct wetting_front *link = (struct wetting_front*) malloc(sizeof(struct wetting_front));
  
  link->depth_cm = depth;
  link->theta = theta;  
  link->front_num = front_num;
  link->layer_num = layer_num;
  link->to_bottom = bottom_flag;
  link->dzdt_cm_per_s = (double)0.0;
  
  //point it to old first wetting_front
  link->next = head;
  
  //point head to new first wetting_front  FORGET THIS CODE WON'T WORK
  head = link;
  link = head->next;
  
  while (link != NULL)
    {
      link->front_num++;
      link=link->next;
    }
  
}


/*#######################################################*/
/* listDeleteFirst -deletes the first entry of a linked list */
/*#######################################################*/
extern struct wetting_front* listDeleteFirst() 
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
bool listIsEmpty() 
{  
  return head == NULL;
}



/*#######################################################*/
/* listLength - counts how many items are in the list    */
/*#######################################################*/
int listLength()
{
  int listLength = 0;
  struct wetting_front *current;
  
  for(current = head; current != NULL; current = current->next) {
    listLength++;
  }
  
  return listLength;
}


/*#######################################################################################*/
/* listFindFront -searches linked list for a particular front number and returns that struct */
/*#######################################################################################*/
extern struct wetting_front* listFindFront(int key, struct wetting_front* head_old)
{
  
  //start from the first link
  //struct wetting_front* current = head;
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
extern struct wetting_front* listDeleteFront(int front_num)
{
  //start from the first link
  struct wetting_front* current = head;
  struct wetting_front* previous = NULL;
  
  //iff list is empty
  //listPrint();
  // debugging
  bool front_found = false;
  while(current != NULL) {
    if (current->front_num == front_num) {
      //printf("Deleting front found....! %d \n", front_num);
      front_found = true;
    }
    current = current->next;
  }
  
  if (!front_found) abort();
  
  current = head;
  if(head == NULL) {
    return NULL; // can't ddo anything, there is no list.
  }
  
  //navigate through list
  while(current->front_num != front_num) {
    //iff it is last wetting_front
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
  if(current == head) {
    //change first to point to next link
    head = head->next;
    previous = head;
    // printf("deleting (current == head) .... \n");
  }
  else {
    //bypass the current link
    previous->next = current->next;
    previous = current->next;
    // printf("deleting (current != head) .... %lf %lf \n",current->depth_cm, previous->depth_cm);
    
  }
  
  current = previous;
  //listPrint();
  //  previous = head;
  while(previous != NULL) { // decrement all front numbers
    //current = link->next;
    //printf("shifting front num = %lf %d \n", previous->depth_cm, previous->front_num);
    previous->front_num--;
    previous = previous->next;  
  }
  //listPrint();
  
  //if (current != NULL)
  //printf("deleted .................... %lf \n", current->depth_cm);
  
  return current;
}


/*####################################################################*/
/* listInsertFront -creates a new front at specified position in list */
/* and increase all front numbers greater than or equal to the        */
/* new front number by 1.                                             */
/*####################################################################*/
extern struct wetting_front* listInsertFront(double depth, double theta, int new_front_num, 
                                           int layer_num, bool bottom_flag)
{
  //start from the first link
  struct wetting_front* link = NULL;
  struct wetting_front* current = NULL;
  struct wetting_front* previous = NULL;
  
  if(head == NULL) { // list is empty
    if(new_front_num==1) { // create it
      //create a link
      struct wetting_front *link = (struct wetting_front*) malloc(sizeof(struct wetting_front));
      
      link->depth_cm = depth;
      link->theta = theta;      
      link->front_num = new_front_num;
      link->layer_num = layer_num;
      link->to_bottom = bottom_flag;
      link->dzdt_cm_per_s = (double)(0.0);
      link->next = NULL;
      head=link;
      return link;
    }
    else { // list is empty and desired link number is >1, unable to create
      return NULL;
    }
  }
  
  //printf("inserting %d \n", new_front_num);
  // list is not empty, so search through list
  previous = head;
  do {
    //    printf("insertingC: %lf %d %d \n", previous->theta, previous->front_num, new_front_num);
    if (previous->front_num == new_front_num-1) { // this is where we want to insert it
      //create a new link
      struct wetting_front *link = (struct wetting_front*) malloc(sizeof(struct wetting_front));
      //    printf("insertingA \n");
      link->depth_cm = depth;
      link->theta = theta;
      link->front_num = new_front_num;
      link->layer_num = layer_num;
      link->to_bottom = bottom_flag;
      link->dzdt_cm_per_s = (double)(0.0);
      link->next = previous->next ;
      previous->next = link;
      //    printf("insertingB \n");
      while(link->next != NULL) { // increment all front numbers
	//	printf("insertingD %d \n", link->front_num);
	
	current = link->next;
	current->front_num++;
	previous = current;
	current = current->next;
	/*
	  current = link->next;
	  current->front_num++;
	  link=current;
	  printf("insertingE %d \n", link->front_num);
	*/
      }
      //abort();
      return link;
    }
    previous = previous-> next;
  } while(TRUE);
  
}

/*############################################################################*/
/* listInsertFrontAtDepth -creates a new front at specified depth in the soil */
/* and increase all front numbers greater than or equal to the                */
/* new front number by 1. Determines layer number, front number and           */
/* determine if the new front is at the bottom of the layer                   */
/*############################################################################*/
extern struct wetting_front* listInsertFrontAtDepth(int num_layers, double *cum_layer_thickness,
                                                    double depth, double theta)
{
  int layer;
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
      // this should never happen, unless it asks to create a front not in the soil
      return NULL;
    }
    
    link->layer_num = el_layer;
    link->to_bottom = extends_to_bottom_flag;
    link->dzdt_cm_per_s = (double)(-1.0);
    link->next = NULL;   // because this is the first link.
    head=link;  // don't forget this.  Must keep head current with the beginning of the list.
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
	// this should never happen, unless it asks to create a front not in the soil
	return NULL;
      }
      
      link->layer_num = el_layer;
      link->to_bottom = extends_to_bottom_flag;
      link->dzdt_cm_per_s = (double)(-1.0);
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
	    // this should never happen, unless it asks to create a front not in the soil
	    return NULL;
	  }
	  
	  link->layer_num = el_layer;
	  link->to_bottom = extends_to_bottom_flag;
	  link->dzdt_cm_per_s = (double)(-1.0);
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
extern void listSortFrontsByDepth() 
{
  int i, j, k, tempKey;
  double tempData;
  struct wetting_front *current;
  struct wetting_front *next;
  
  int size = listLength();
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
	
	tempKey = current->front_num;
	current->front_num = next->front_num;
	next->front_num = tempKey;
      }
                      
      current = current->next;
      next = next->next;
    }
  }
  
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

