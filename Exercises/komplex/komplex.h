#ifndef HAVE_KOMPLEX_H /* this is necessary for multiple includes */
#define HAVE_KOMPLEX_H

struct komplex {double re; double im;};
typedef struct komplex komplex;

void    komplex_print    (char* s, komplex z);   /* prints string s and then komplex z */
void    komplex_set      (komplex* z, double x, double y);   /* z is set to x+i*y */
komplex komplex_new      (double x, double y);   /* returns x+i*y */
komplex komplex_add      (komplex a, komplex b); /* returns a+b */
komplex komplex_sub      (komplex a, komplex b); /* returns a-b */


#endif
