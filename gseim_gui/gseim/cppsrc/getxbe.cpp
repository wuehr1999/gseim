#include "global.h"
#include "xbeusr.h"
#include "xbejac.h"
#include <cstdlib>
#include <iostream>
#include <string>
using namespace std;
void x_abc_to_dq(
 Global &G, XbeUsr &X, XbeJac &J);
void x_abc_to_dq_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_and_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_clock(
 Global &G, XbeUsr &X, XbeJac &J);
void x_cmpr_1_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_cmpr_1_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_cmpr_2_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_cmpr_2_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_cmpr_simple_2_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_cmpr_simple_2_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_cmprh_1_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_cmprh_2_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_const(
 Global &G, XbeUsr &X, XbeJac &J);
void x_cosine(
 Global &G, XbeUsr &X, XbeJac &J);
void x_diff(
 Global &G, XbeUsr &X, XbeJac &J);
void x_dq_to_abc(
 Global &G, XbeUsr &X, XbeJac &J);
void x_dummy_sink(
 Global &G, XbeUsr &X, XbeJac &J);
void x_dummy_source(
 Global &G, XbeUsr &X, XbeJac &J);
void x_filter_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_indmc1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_indmc2a(
 Global &G, XbeUsr &X, XbeJac &J);
void x_indmc2b(
 Global &G, XbeUsr &X, XbeJac &J);
void x_integrator(
 Global &G, XbeUsr &X, XbeJac &J);
void x_integrator_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_lag_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_lag_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_limiter(
 Global &G, XbeUsr &X, XbeJac &J);
void x_modulo(
 Global &G, XbeUsr &X, XbeJac &J);
void x_modulo_twopi(
 Global &G, XbeUsr &X, XbeJac &J);
void x_mult_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_multscl(
 Global &G, XbeUsr &X, XbeJac &J);
void x_not(
 Global &G, XbeUsr &X, XbeJac &J);
void x_or_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_pole_complex_order_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_pole_complex_order_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_pole_real_order_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_pole_real_order_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_pole_real_order_3(
 Global &G, XbeUsr &X, XbeJac &J);
void x_pole_real_order_4(
 Global &G, XbeUsr &X, XbeJac &J);
void x_pole_real_order_5(
 Global &G, XbeUsr &X, XbeJac &J);
void x_pwl10_xy(
 Global &G, XbeUsr &X, XbeJac &J);
void x_pwl20(
 Global &G, XbeUsr &X, XbeJac &J);
void x_sine(
 Global &G, XbeUsr &X, XbeJac &J);
void x_srcac(
 Global &G, XbeUsr &X, XbeJac &J);
void x_sum_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_sum_3(
 Global &G, XbeUsr &X, XbeJac &J);
void x_sum_4(
 Global &G, XbeUsr &X, XbeJac &J);
void x_sum_5(
 Global &G, XbeUsr &X, XbeJac &J);
void x_sum_6(
 Global &G, XbeUsr &X, XbeJac &J);
void x_triangle_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_triangle_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_triangle_3(
 Global &G, XbeUsr &X, XbeJac &J);
void x_user_fn_1_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_user_fn_2_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_user_fn_3_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_user_fn_4_3(
 Global &G, XbeUsr &X, XbeJac &J);
void x_user_fn_5_3(
 Global &G, XbeUsr &X, XbeJac &J);
void x_vsi_3ph_1(
 Global &G, XbeUsr &X, XbeJac &J);
void get_xbe(
 const int i_xbel,
 Global &G,
 XbeUsr &X,
 XbeJac &J) {
   switch (i_xbel) {
     case 0:
       x_abc_to_dq(G,X,J);
       break;
     case 1:
       x_abc_to_dq_2(G,X,J);
       break;
     case 2:
       x_and_2(G,X,J);
       break;
     case 3:
       x_clock(G,X,J);
       break;
     case 4:
       x_cmpr_1_1(G,X,J);
       break;
     case 5:
       x_cmpr_1_2(G,X,J);
       break;
     case 6:
       x_cmpr_2_1(G,X,J);
       break;
     case 7:
       x_cmpr_2_2(G,X,J);
       break;
     case 8:
       x_cmpr_simple_2_1(G,X,J);
       break;
     case 9:
       x_cmpr_simple_2_2(G,X,J);
       break;
     case 10:
       x_cmprh_1_1(G,X,J);
       break;
     case 11:
       x_cmprh_2_1(G,X,J);
       break;
     case 12:
       x_const(G,X,J);
       break;
     case 13:
       x_cosine(G,X,J);
       break;
     case 14:
       x_diff(G,X,J);
       break;
     case 15:
       x_dq_to_abc(G,X,J);
       break;
     case 16:
       x_dummy_sink(G,X,J);
       break;
     case 17:
       x_dummy_source(G,X,J);
       break;
     case 18:
       x_filter_1(G,X,J);
       break;
     case 19:
       x_indmc1(G,X,J);
       break;
     case 20:
       x_indmc2a(G,X,J);
       break;
     case 21:
       x_indmc2b(G,X,J);
       break;
     case 22:
       x_integrator(G,X,J);
       break;
     case 23:
       x_integrator_1(G,X,J);
       break;
     case 24:
       x_lag_1(G,X,J);
       break;
     case 25:
       x_lag_2(G,X,J);
       break;
     case 26:
       x_limiter(G,X,J);
       break;
     case 27:
       x_modulo(G,X,J);
       break;
     case 28:
       x_modulo_twopi(G,X,J);
       break;
     case 29:
       x_mult_2(G,X,J);
       break;
     case 30:
       x_multscl(G,X,J);
       break;
     case 31:
       x_not(G,X,J);
       break;
     case 32:
       x_or_2(G,X,J);
       break;
     case 33:
       x_pole_complex_order_1(G,X,J);
       break;
     case 34:
       x_pole_complex_order_2(G,X,J);
       break;
     case 35:
       x_pole_real_order_1(G,X,J);
       break;
     case 36:
       x_pole_real_order_2(G,X,J);
       break;
     case 37:
       x_pole_real_order_3(G,X,J);
       break;
     case 38:
       x_pole_real_order_4(G,X,J);
       break;
     case 39:
       x_pole_real_order_5(G,X,J);
       break;
     case 40:
       x_pwl10_xy(G,X,J);
       break;
     case 41:
       x_pwl20(G,X,J);
       break;
     case 42:
       x_sine(G,X,J);
       break;
     case 43:
       x_srcac(G,X,J);
       break;
     case 44:
       x_sum_2(G,X,J);
       break;
     case 45:
       x_sum_3(G,X,J);
       break;
     case 46:
       x_sum_4(G,X,J);
       break;
     case 47:
       x_sum_5(G,X,J);
       break;
     case 48:
       x_sum_6(G,X,J);
       break;
     case 49:
       x_triangle_1(G,X,J);
       break;
     case 50:
       x_triangle_2(G,X,J);
       break;
     case 51:
       x_triangle_3(G,X,J);
       break;
     case 52:
       x_user_fn_1_1(G,X,J);
       break;
     case 53:
       x_user_fn_2_1(G,X,J);
       break;
     case 54:
       x_user_fn_3_1(G,X,J);
       break;
     case 55:
       x_user_fn_4_3(G,X,J);
       break;
     case 56:
       x_user_fn_5_3(G,X,J);
       break;
     case 57:
       x_vsi_3ph_1(G,X,J);
       break;
   }
   return;
}
