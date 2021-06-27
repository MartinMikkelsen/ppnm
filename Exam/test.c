#include<math . h>
#include<s t d l i b . h>
#include<s t d i o . h>
#define FOR(k) for(int k=0;k<dim;k++)
#define RND (double)rand()/RANDMAX
double strata(
int dim , double f(int dim, double∗ x),
double∗ a, double∗ b,
double acc, double eps,
int n_reuse , double mean_reuse)
{
int N=16∗dim;
double V=1; FOR(k) V∗=b[k]−a[k];
int n_left[dim], n_right[dim];
double x[dim] , mean_left[dim], mean_right[dim], mean=0;
FOR(k){ mean_left[k]= 0; me an ri gh t [ k ]= 0; n l e f t [ k ]= 0; n r i g h t [ k ]= 0; }
fo r ( int i =0; i<N; i ++){
FOR( k ) x [ k]=a [ k]+RND∗( b [ k]−a [ k ] ) ;
double f x=f ( dim , x ) ;
mean+=f x ;
FOR( k ){
i f ( x [ k] >( a [ k]+b [ k ] ) / 2 ) { n r i g h t [ k]++; me an ri gh t [ k]+= f x ; }
e l s e { n l e f t [ k]++; m e a n l e f t [ k]+= f x ; }
}
}
mean/=N;
FOR( k ){ m e a n l e f t [ k]/= n l e f t [ k ] ; me an ri gh t [ k]/= n r i g h t [ k ] ; }
int kdiv =0; double maxvar=0;
FOR( k ){
double var=f a b s ( me an ri gh t [ k]− m e a n l e f t [ k ] ) ;
i f ( var>maxvar ){ maxvar=var ; kdiv=k ; }
}
double i n t e g =(mean∗N+mean reuse ∗ n r e u s e ) / (N+n r e u s e )∗V;
double e r r o r=f a b s ( mean reuse−mean )∗V;
double t o l e r=acc+f a b s ( i n t e g )∗ ep s ;
i f ( e r r o r <t o l e r ) return i n t e g ;
double a2 [ dim ] , b2 [ dim ] ; FOR( k ) a2 [ k]=a [ k ] ; FOR( k ) b2 [ k]=b [ k ] ;
a2 [ kdiv ]=( a [ kdiv ]+b [ kdiv ] ) / 2 ; b2 [ kdiv ]=( a [ kdiv ]+b [ kdiv ] ) / 2 ;
double i n t e g l e f t=
s t r a t a ( dim , f , a , b2 , acc / s q r t ( 2 ) , eps , n l e f t [ kdiv ] , m e a n l e f t [ kdiv ] ) ;
double i n t e g r i g h t=
s t r a t a ( dim , f , a2 , b , acc / s q r t ( 2 ) , eps , n r i g h t [ kdiv ] , me an ri gh t [ kdiv ] ) ;
retu
