#include "../include/all.hxx"

/*#########################################################################*/
/*#########################################################################*/
/*#########################################################################*/
/***************************************************************************/
//
//         UU      UU   TTTTTTTTTTTT     II     LL
//         UU      UU        TT          II     LL
//         UU      UU        TT          II     LL
//         UU      UU        TT          II     LL
//         UU      UU        TT          II     LL
//         UU      UU        TT          II     LL
//         UU      UU        TT          II     LL
//          UU    UU         TT          II     LL
//          UUU  UUU         TT          II     LL
//             UU            TT          II     LLLLLLLLL   ities
//
//---------------------------------------------------------------------------

extern bool is_epsilon_less_than(double a, double eps)
{
/***********************************************************************/
/* testing to see if two real numbers are equal hampered by precision. */
/** this is the "close enough function                                **/
/** it returns TRUE if a real number is                               **/
/** close enough (epsilon) to zero.                                   **/
/***********************************************************************/
if(fabs(a)<eps) return TRUE;  // note use of the real number absolute value function.
else            return FALSE;
}
