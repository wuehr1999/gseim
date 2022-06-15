#include "global.h"
#include "ebeusr.h"
#include "ebejac.h"
#include <cstdlib>
#include <iostream>
#include <string>
using namespace std;
void e_ammeter(
 Global &G, EbeUsr &X, EbeJac &J);
void e_ammeter_fb(
 Global &G, EbeUsr &X, EbeJac &J);
void e_c(
 Global &G, EbeUsr &X, EbeJac &J);
void e_diode_r(
 Global &G, EbeUsr &X, EbeJac &J);
void e_diode_spice(
 Global &G, EbeUsr &X, EbeJac &J);
void e_dummy_e(
 Global &G, EbeUsr &X, EbeJac &J);
void e_ground(
 Global &G, EbeUsr &X, EbeJac &J);
void e_l(
 Global &G, EbeUsr &X, EbeJac &J);
void e_r(
 Global &G, EbeUsr &X, EbeJac &J);
void e_solar_module_1(
 Global &G, EbeUsr &X, EbeJac &J);
void e_solar_module_rs(
 Global &G, EbeUsr &X, EbeJac &J);
void e_switch_1(
 Global &G, EbeUsr &X, EbeJac &J);
void e_thyristor(
 Global &G, EbeUsr &X, EbeJac &J);
void e_voltmeter(
 Global &G, EbeUsr &X, EbeJac &J);
void e_voltmeter_1(
 Global &G, EbeUsr &X, EbeJac &J);
void e_voltmeter_fb(
 Global &G, EbeUsr &X, EbeJac &J);
void e_vsrc_ac(
 Global &G, EbeUsr &X, EbeJac &J);
void e_vsrc_clock(
 Global &G, EbeUsr &X, EbeJac &J);
void e_vsrc_dc(
 Global &G, EbeUsr &X, EbeJac &J);
void e_vsrc_pulse10(
 Global &G, EbeUsr &X, EbeJac &J);
void e_vsrc_x(
 Global &G, EbeUsr &X, EbeJac &J);
void e_xfmr_l1l2(
 Global &G, EbeUsr &X, EbeJac &J);
void e_xfmr_level0_1ph(
 Global &G, EbeUsr &X, EbeJac &J);
void e_xfmr_level2_1ph(
 Global &G, EbeUsr &X, EbeJac &J);
void e_xfmr_level0_1ph_1_2(
 Global &G, EbeUsr &X, EbeJac &J);
void get_ebe(
 const int i_ebel,
 Global &G,
 EbeUsr &X,
 EbeJac &J) {
   switch (i_ebel) {
     case 0:
       e_ammeter(G,X,J);
       break;
     case 1:
       e_ammeter_fb(G,X,J);
       break;
     case 2:
       e_c(G,X,J);
       break;
     case 3:
       e_diode_r(G,X,J);
       break;
     case 4:
       e_diode_spice(G,X,J);
       break;
     case 5:
       e_dummy_e(G,X,J);
       break;
     case 6:
       e_ground(G,X,J);
       break;
     case 7:
       e_l(G,X,J);
       break;
     case 8:
       e_r(G,X,J);
       break;
     case 9:
       e_solar_module_1(G,X,J);
       break;
     case 10:
       e_solar_module_rs(G,X,J);
       break;
     case 11:
       e_switch_1(G,X,J);
       break;
     case 12:
       e_thyristor(G,X,J);
       break;
     case 13:
       e_voltmeter(G,X,J);
       break;
     case 14:
       e_voltmeter_1(G,X,J);
       break;
     case 15:
       e_voltmeter_fb(G,X,J);
       break;
     case 16:
       e_vsrc_ac(G,X,J);
       break;
     case 17:
       e_vsrc_clock(G,X,J);
       break;
     case 18:
       e_vsrc_dc(G,X,J);
       break;
     case 19:
       e_vsrc_pulse10(G,X,J);
       break;
     case 20:
       e_vsrc_x(G,X,J);
       break;
     case 21:
       e_xfmr_l1l2(G,X,J);
       break;
     case 22:
       e_xfmr_level0_1ph(G,X,J);
       break;
     case 23:
       e_xfmr_level2_1ph(G,X,J);
       break;
     case 24:
       e_xfmr_level0_1ph_1_2(G,X,J);
       break;
   }
   return;
}
