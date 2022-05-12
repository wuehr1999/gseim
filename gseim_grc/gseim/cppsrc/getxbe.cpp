#include "global.h"
#include "xbeusr.h"
#include "xbejac.h"
#include <cstdlib>
#include <iostream>
#include <string>
using namespace std;
void x_abc_to_alphabeta_3(
 Global &G, XbeUsr &X, XbeJac &J);
void x_abc_to_dq(
 Global &G, XbeUsr &X, XbeJac &J);
void x_abc_to_dq_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_abc_to_dq0_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_abs(
 Global &G, XbeUsr &X, XbeJac &J);
void x_and_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_atan2_rad(
 Global &G, XbeUsr &X, XbeJac &J);
void x_clock(
 Global &G, XbeUsr &X, XbeJac &J);
void x_clock_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_clock_3(
 Global &G, XbeUsr &X, XbeJac &J);
void x_clock_3ph(
 Global &G, XbeUsr &X, XbeJac &J);
void x_clock_thyr(
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
void x_cos(
 Global &G, XbeUsr &X, XbeJac &J);
void x_dead_zone(
 Global &G, XbeUsr &X, XbeJac &J);
void x_delay_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_delay_discrete(
 Global &G, XbeUsr &X, XbeJac &J);
void x_delay_onestep(
 Global &G, XbeUsr &X, XbeJac &J);
void x_diff(
 Global &G, XbeUsr &X, XbeJac &J);
void x_div(
 Global &G, XbeUsr &X, XbeJac &J);
void x_dq_to_abc(
 Global &G, XbeUsr &X, XbeJac &J);
void x_dq0_to_abc_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_dummy_sink(
 Global &G, XbeUsr &X, XbeJac &J);
void x_dummy_source(
 Global &G, XbeUsr &X, XbeJac &J);
void x_edge_delay(
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
void x_linear(
 Global &G, XbeUsr &X, XbeJac &J);
void x_modulo(
 Global &G, XbeUsr &X, XbeJac &J);
void x_modulo_twopi(
 Global &G, XbeUsr &X, XbeJac &J);
void x_monostable_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_multscl(
 Global &G, XbeUsr &X, XbeJac &J);
void x_mult_2(
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
void x_pulse10(
 Global &G, XbeUsr &X, XbeJac &J);
void x_pwl10_xy(
 Global &G, XbeUsr &X, XbeJac &J);
void x_pwl20(
 Global &G, XbeUsr &X, XbeJac &J);
void x_pwm20_1(
 Global &G, XbeUsr &X, XbeJac &J);
void x_sampler(
 Global &G, XbeUsr &X, XbeJac &J);
void x_signum(
 Global &G, XbeUsr &X, XbeJac &J);
void x_sin(
 Global &G, XbeUsr &X, XbeJac &J);
void x_src_ac(
 Global &G, XbeUsr &X, XbeJac &J);
void x_sum_2(
 Global &G, XbeUsr &X, XbeJac &J);
void x_sum_3(
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
void x_xfer_fn(
 Global &G, XbeUsr &X, XbeJac &J);
void get_xbe(
 const int i_xbel,
 Global &G,
 XbeUsr &X,
 XbeJac &J) {
   switch (i_xbel) {
     case 0:
       x_abc_to_alphabeta_3(G,X,J);
       break;
     case 1:
       x_abc_to_dq(G,X,J);
       break;
     case 2:
       x_abc_to_dq_2(G,X,J);
       break;
     case 3:
       x_abc_to_dq0_2(G,X,J);
       break;
     case 4:
       x_abs(G,X,J);
       break;
     case 5:
       x_and_2(G,X,J);
       break;
     case 6:
       x_atan2_rad(G,X,J);
       break;
     case 7:
       x_clock(G,X,J);
       break;
     case 8:
       x_clock_1(G,X,J);
       break;
     case 9:
       x_clock_3(G,X,J);
       break;
     case 10:
       x_clock_3ph(G,X,J);
       break;
     case 11:
       x_clock_thyr(G,X,J);
       break;
     case 12:
       x_cmpr_1_1(G,X,J);
       break;
     case 13:
       x_cmpr_1_2(G,X,J);
       break;
     case 14:
       x_cmpr_2_1(G,X,J);
       break;
     case 15:
       x_cmpr_2_2(G,X,J);
       break;
     case 16:
       x_cmpr_simple_2_1(G,X,J);
       break;
     case 17:
       x_cmpr_simple_2_2(G,X,J);
       break;
     case 18:
       x_cmprh_1_1(G,X,J);
       break;
     case 19:
       x_cmprh_2_1(G,X,J);
       break;
     case 20:
       x_const(G,X,J);
       break;
     case 21:
       x_cos(G,X,J);
       break;
     case 22:
       x_dead_zone(G,X,J);
       break;
     case 23:
       x_delay_1(G,X,J);
       break;
     case 24:
       x_delay_discrete(G,X,J);
       break;
     case 25:
       x_delay_onestep(G,X,J);
       break;
     case 26:
       x_diff(G,X,J);
       break;
     case 27:
       x_div(G,X,J);
       break;
     case 28:
       x_dq_to_abc(G,X,J);
       break;
     case 29:
       x_dq0_to_abc_2(G,X,J);
       break;
     case 30:
       x_dummy_sink(G,X,J);
       break;
     case 31:
       x_dummy_source(G,X,J);
       break;
     case 32:
       x_edge_delay(G,X,J);
       break;
     case 33:
       x_indmc1(G,X,J);
       break;
     case 34:
       x_indmc2a(G,X,J);
       break;
     case 35:
       x_indmc2b(G,X,J);
       break;
     case 36:
       x_integrator(G,X,J);
       break;
     case 37:
       x_integrator_1(G,X,J);
       break;
     case 38:
       x_lag_1(G,X,J);
       break;
     case 39:
       x_lag_2(G,X,J);
       break;
     case 40:
       x_limiter(G,X,J);
       break;
     case 41:
       x_linear(G,X,J);
       break;
     case 42:
       x_modulo(G,X,J);
       break;
     case 43:
       x_modulo_twopi(G,X,J);
       break;
     case 44:
       x_monostable_1(G,X,J);
       break;
     case 45:
       x_multscl(G,X,J);
       break;
     case 46:
       x_mult_2(G,X,J);
       break;
     case 47:
       x_not(G,X,J);
       break;
     case 48:
       x_or_2(G,X,J);
       break;
     case 49:
       x_pole_complex_order_1(G,X,J);
       break;
     case 50:
       x_pole_complex_order_2(G,X,J);
       break;
     case 51:
       x_pole_real_order_1(G,X,J);
       break;
     case 52:
       x_pole_real_order_2(G,X,J);
       break;
     case 53:
       x_pole_real_order_3(G,X,J);
       break;
     case 54:
       x_pole_real_order_4(G,X,J);
       break;
     case 55:
       x_pole_real_order_5(G,X,J);
       break;
     case 56:
       x_pulse10(G,X,J);
       break;
     case 57:
       x_pwl10_xy(G,X,J);
       break;
     case 58:
       x_pwl20(G,X,J);
       break;
     case 59:
       x_pwm20_1(G,X,J);
       break;
     case 60:
       x_sampler(G,X,J);
       break;
     case 61:
       x_signum(G,X,J);
       break;
     case 62:
       x_sin(G,X,J);
       break;
     case 63:
       x_src_ac(G,X,J);
       break;
     case 64:
       x_sum_2(G,X,J);
       break;
     case 65:
       x_sum_3(G,X,J);
       break;
     case 66:
       x_triangle_1(G,X,J);
       break;
     case 67:
       x_triangle_2(G,X,J);
       break;
     case 68:
       x_triangle_3(G,X,J);
       break;
     case 69:
       x_user_fn_1_1(G,X,J);
       break;
     case 70:
       x_user_fn_2_1(G,X,J);
       break;
     case 71:
       x_user_fn_3_1(G,X,J);
       break;
     case 72:
       x_user_fn_4_3(G,X,J);
       break;
     case 73:
       x_user_fn_5_3(G,X,J);
       break;
     case 74:
       x_vsi_3ph_1(G,X,J);
       break;
     case 75:
       x_xfer_fn(G,X,J);
       break;
   }
   return;
}
