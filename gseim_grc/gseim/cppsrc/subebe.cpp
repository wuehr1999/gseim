#include "global.h"
#include "ebeusr.h"
#include "ebejac.h"
#include "utils.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
void e_ammeter(Global &G,EbeUsr &X,EbeJac &J) {
   double cur_p;
   double cur_p_s;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int na_cur_p = 0;
   const int nas_cur_p_s = 0;
   const int no_i = 0;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   if (G.flags[G.i_dc] || G.flags[G.i_trns]) {
     cur_p = X.val_aux[na_cur_p];

     if (G.flags[G.i_function]) {
       X.f[nf_1] =  cur_p;
       X.f[nf_2] = -cur_p;
       X.f[nf_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdaux[nf_1][na_cur_p] =  1.0;
       J.dfdaux[nf_2][na_cur_p] = -1.0;
       J.dfdv  [nf_3][nnd_p   ] =  1.0;
       J.dfdv  [nf_3][nnd_n   ] = -1.0;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     cur_p_s = X.val_auxs[nas_cur_p_s];
     if (G.flags[G.i_function]) {
       X.h[nh_1] =  cur_p_s;
       X.h[nh_2] = -cur_p_s;
       X.h[nh_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
     }
     if (G.flags[G.i_jacobian]) {
       J.dhdauxs[nh_1][nas_cur_p_s] =  1.0;
       J.dhdauxs[nh_2][nas_cur_p_s] = -1.0;
       J.dhdv   [nh_3][nnd_p      ] =  1.0;
       J.dhdv   [nh_3][nnd_n      ] = -1.0;
     }
     return;
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_i] = X.cur_nd[nnd_p];
      return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_nd[nnd_n] = 0.0;
     X.val_nd[nnd_p] = 0.0;
     return;
   }
   return;
}
void e_ammeter_fb(Global &G,EbeUsr &X,EbeJac &J) {
   double c1,c2;
   double i_fb;
   double k_scale;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int nx_i_fb = 0;
   const int nr_k_scale = 0;
   const int no_i = 0;
   const int no_i_fb = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   k_scale = X.rprm[nr_k_scale];

   if (k_scale == 0.0) {
     cout << "ammeter_fb: k_scale = 0 ? Halting..." << endl;
     exit(1);
   }
   if (G.flags[G.i_dc] || G.flags[G.i_trns]) {
     i_fb = X.val_xvr[nx_i_fb];

     if (G.flags[G.i_function]) {
       c1 = i_fb/k_scale;
       X.f[nf_1] =  c1;
       X.f[nf_2] = -c1;
       X.f[nf_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
     }
     if (G.flags[G.i_jacobian]) {
       c2 = 1.0/k_scale;
       J.dfdxvr[nf_1][nx_i_fb] =  c2;
       J.dfdxvr[nf_2][nx_i_fb] = -c2;
       J.dfdv  [nf_3][nnd_p  ] =  1.0;
       J.dfdv  [nf_3][nnd_n  ] = -1.0;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     i_fb = X.val_xvr[nx_i_fb];

     if (G.flags[G.i_function]) {
       c1 = i_fb/k_scale;
       X.h[nh_1] =  c1;
       X.h[nh_2] = -c1;
       X.h[nh_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
     }
     if (G.flags[G.i_jacobian]) {
       c2 = 1.0/k_scale;
       J.dhdxvr[nh_1][nx_i_fb] =  c2;
       J.dhdxvr[nh_2][nx_i_fb] = -c2;
       J.dhdv  [nh_3][nnd_p  ] =  1.0;
       J.dhdv  [nh_3][nnd_n  ] = -1.0;
     }
     return;
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_i] = X.cur_nd[nnd_p];
      X.outprm[no_i_fb] = X.rprm[nr_k_scale]*X.cur_nd[nnd_p];
      return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_xvr[nx_i_fb] = 0.0;
     return;
   }
   return;
}
void e_c(Global &G,EbeUsr &X,EbeJac &J) {
   double c1;
   double qp,qm;
   double cur_p;
   double c,k_scale;
   double v0;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int nstv_qp = 0;
   const int nstv_qm = 1;
   const int nas_cur_p = 0;
   const int nr_c = 0;
   const int nr_k_scale = 1;
   const int nst_v0 = 0;
   const int no_i = 0;
   const int no_v = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   c = X.rprm[nr_c];
   k_scale = X.rprm[nr_k_scale];
   c1 = c*k_scale;

   if (G.flags[G.i_dc] || G.flags[G.i_trns]) {
      if (G.flags[G.i_function]) {
         X.f[nf_1] = 0.0;
         X.f[nf_2] = 0.0;

         X.g[ng_1] = c1*(X.val_nd[nnd_p]-X.val_nd[nnd_n]);
         X.g[ng_2] = -X.g[ng_1];
      }
      if (G.flags[G.i_jacobian]) {
         J.dgdv[ng_1][nnd_p] =  c1;
         J.dgdv[ng_1][nnd_n] = -c1;
         J.dgdv[ng_2][nnd_p] = -c1;
         J.dgdv[ng_2][nnd_n] =  c1;
      }
      X.val_stv[nstv_qp] = c1*(X.val_nd[nnd_p]-X.val_nd[nnd_n]);
      X.val_stv[nstv_qm] = -X.val_stv[nstv_qp];
   }

   if (G.flags[G.i_startup]) {
      v0 = X.stprm[nst_v0];
//    cout << "c.ebe: v0 = " << v0 << endl;
      cur_p = X.val_auxs[nas_cur_p];
      if (G.flags[G.i_function]) {
         X.h[nh_1] =  cur_p;
         X.h[nh_2] = -cur_p;
         X.h[nh_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n]-v0;
      }
      if (G.flags[G.i_jacobian]) {
         J.dhdauxs[nh_1][nas_cur_p] =  1.0;
         J.dhdauxs[nh_2][nas_cur_p] = -1.0;
         J.dhdv   [nh_3][nnd_p    ] =  1.0;
         J.dhdv   [nh_3][nnd_n    ] = -1.0;
      }
      X.val_stv[nstv_qp] = c1*(X.val_nd[nnd_p]-X.val_nd[nnd_n]);
      X.val_stv[nstv_qm] = -X.val_stv[nstv_qp];
      return;
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
      X.outprm[no_i] = X.cur_nd[nnd_p];
      return;
   }
   if (G.flags[G.i_init_guess]) {
      X.val_nd[nnd_p] = 0.0;
      X.val_nd[nnd_n] = 0.0;
      return;
   }
   return;
}
void e_diode_r(Global &G,EbeUsr &X,EbeJac &J) {
   double vp,vn,r,g,v0;
   double r_on,r_off,v_on,v_on_1;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int nr_r_on = 0;
   const int nr_r_off = 1;
   const int nr_v_on = 2;
   const int nr_v_on_1 = 3;
   const int no_i = 0;
   const int no_v = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nh_1 = 0;
   const int nh_2 = 1;
   if (G.flags[G.i_one_time_parms]) {
     r_on  = X.rprm[nr_r_on ];
     r_off = X.rprm[nr_r_off];
     v_on  = X.rprm[nr_v_on ];

     v_on_1 = v_on*r_off/(r_off-r_on);
     X.rprm[nr_v_on_1] = v_on_1;

     return;
   }
   v_on = X.rprm[nr_v_on];

   if (G.flags[G.i_dc] || G.flags[G.i_trns] || G.flags[G.i_startup]) {
     vp = X.val_nd[nnd_p];
     vn = X.val_nd[nnd_n];
     v0 = vp-vn;

     r_on    = X.rprm[nr_r_on  ];
     r_off   = X.rprm[nr_r_off ];
     v_on_1  = X.rprm[nr_v_on_1];

     if (v0 >= v_on_1) {
       r = r_on;
     } else {
       r = r_off;
     }
     if (r < 1.0e-9) {
       cout << "diode_r: r too small!" << endl;
       exit(1);
     }
     g = 1.0/r;
   }
   if (G.flags[G.i_dc] || G.flags[G.i_trns]) {
     if (G.flags[G.i_function]) {
       if (v0 >= v_on_1) {
         X.f[nf_1] = g*(vp-vn)-g*v_on;
       } else {
         X.f[nf_1] = g*(vp-vn);
       }
       X.f[nf_2] = -X.f[nf_1];
     }

     if (G.flags[G.i_jacobian]) {
       J.dfdv[nf_1][nnd_p] =  g;
       J.dfdv[nf_1][nnd_n] = -g;
       J.dfdv[nf_2][nnd_p] = -g;
       J.dfdv[nf_2][nnd_n] =  g;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_function]) {
       if (v0 >= v_on_1) {
         X.h[nh_1] = g*(vp-vn)-g*v_on;
       } else {
         X.h[nh_1] = g*(vp-vn);
       }
       X.h[nh_2] = -X.h[nh_1];
     }

     if (G.flags[G.i_jacobian]) {
       J.dhdv[nh_1][nnd_p] =  g;
       J.dhdv[nh_1][nnd_n] = -g;
       J.dhdv[nh_2][nnd_p] = -g;
       J.dhdv[nh_2][nnd_n] =  g;
     }
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
     X.outprm[no_i] = X.cur_nd[nnd_p];
     return;
   }
   if (G.flags[G.i_init_guess]) {
      X.val_nd[nnd_n] = 0.0;
      X.val_nd[nnd_p] = 0.0;
      return;
   }
   return;
}
void e_dummy_e(Global &G,EbeUsr &X,EbeJac &J) {
   const int nnd_a = 0;
   const int nf_1 = 0;
   return;
}
void e_ground(Global &G,EbeUsr &X,EbeJac &J) {
   const int nnd_g = 0;
   const int nf_1 = 0;
// do nothing
   return;
}
void e_l(Global &G,EbeUsr &X,EbeJac &J) {
   double l1;
   double cur_p;
   double l,k_scale;
   double i0;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int na_cur_p = 0;
   const int nr_l = 0;
   const int nr_k_scale = 1;
   const int nst_i0 = 0;
   const int no_i = 0;
   const int no_v = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nh_1 = 0;
   const int nh_2 = 1;
   if (G.flags[G.i_dc] || G.flags[G.i_trns]) {
     l = X.rprm[nr_l];
     k_scale = X.rprm[nr_k_scale];
     l1 = l*k_scale;

     cur_p = X.val_aux[na_cur_p];
     if (G.flags[G.i_function]) {
       X.f[nf_1] =  cur_p;
       X.f[nf_2] = -cur_p;
       X.f[nf_3] = (X.val_nd[nnd_p]-X.val_nd[nnd_n])/l1;
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdaux[nf_1][na_cur_p] =  1.0;
       J.dfdaux[nf_2][na_cur_p] = -1.0;
       J.dfdv  [nf_3][nnd_p   ] =  1.0/l1;
       J.dfdv  [nf_3][nnd_n   ] = -1.0/l1;
     }
   }

   if (G.flags[G.i_startup]) {
      i0 = X.stprm[nst_i0];
      if (G.flags[G.i_function]) {
         X.h[nh_1] =  i0;
         X.h[nh_2] = -i0;
      }
      X.val_aux[na_cur_p] = i0;
      return;
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
      X.outprm[no_i] = X.cur_nd[nnd_p];
      return;
   }
   if (G.flags[G.i_init_guess]) {
      X.val_nd[nnd_p] = 0.0;
      X.val_nd[nnd_n] = 0.0;
      return;
   }
   return;
}
void e_r(Global &G,EbeUsr &X,EbeJac &J) {
   double vp,vn,r_eff;
   double r,k_scale,g;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int nr_r = 0;
   const int nr_k_scale = 1;
   const int nr_g = 2;
   const int no_i = 0;
   const int no_v = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nh_1 = 0;
   const int nh_2 = 1;
   if (G.flags[G.i_one_time_parms]) {
     r = X.rprm[nr_r];
     k_scale = X.rprm[nr_k_scale];
     if (r*k_scale < 1.0e-9) {
       cout << "r.ebe: r too small!" << endl;
       cout << "r.ebe: r=" << r << endl;
       cout << "r.ebe: k_scale=" << k_scale << endl;
       cout << "r.ebe: Halting..." << endl;
       exit(1);
     }
     r_eff = r*k_scale;

     g = 1.0e0/r_eff;
     X.rprm[nr_g] = g;
     return;
   }

   if (G.flags[G.i_dc] || G.flags[G.i_trns]) {
     vp = X.val_nd[nnd_p];
     vn = X.val_nd[nnd_n];

     g = X.rprm[nr_g];

     if (G.flags[G.i_function]) {
       X.f[nf_1] = g*(vp-vn);
       X.f[nf_2] = -X.f[nf_1];
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdv[nf_1][nnd_p] =  g;
       J.dfdv[nf_1][nnd_n] = -g;
       J.dfdv[nf_2][nnd_p] = -g;
       J.dfdv[nf_2][nnd_n] =  g;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     vp = X.val_nd[nnd_p];
     vn = X.val_nd[nnd_n];

     g = X.rprm[nr_g];

     if (G.flags[G.i_function]) {
       X.h[nh_1] = g*(vp-vn);
       X.h[nh_2] = -X.f[nf_1];
     }
     if (G.flags[G.i_jacobian]) {
       J.dhdv[nh_1][nnd_p] =  g;
       J.dhdv[nh_1][nnd_n] = -g;
       J.dhdv[nh_2][nnd_p] = -g;
       J.dhdv[nh_2][nnd_n] =  g;
     }
     return;
   }
   if (G.flags[G.i_outvar]) {
     g = X.rprm[nr_g];
     X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
     X.outprm[no_i] = g*(X.val_nd[nnd_p]-X.val_nd[nnd_n]);
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_nd[nnd_n] = 0.0;
     X.val_nd[nnd_p] = 0.0;
     return;
   }
   return;
}
void e_switch_1(Global &G,EbeUsr &X,EbeJac &J) {
   double r,g,vp,vn;
   bool l_closed;
   double x;
   double r_on,r_off,v_on,x_high;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int nx_x = 0;
   const int nr_r_on = 0;
   const int nr_r_off = 1;
   const int nr_v_on = 2;
   const int nr_x_high = 3;
   const int no_i = 0;
   const int no_v = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nh_1 = 0;
   const int nh_2 = 1;
   if (G.flags[G.i_one_time_parms]) {
     v_on  = X.rprm[nr_v_on ];
     r_on  = X.rprm[nr_r_on ];
     r_off = X.rprm[nr_r_off];

     if (v_on < 0.0) {
       cout << "switch_g1.ebe: v_on < 0 ?!" << endl;
       cout << "   Halting..." << endl;
       exit(1);
     }
     if (r_on == 0.0) {
       cout << "switch_g1.ebe: r_on  =  0 ?!" << endl;
       cout << "   Halting..." << endl;
       exit(1);
     }
     if (r_off < r_on) {
       cout << "switch_g1.ebe: r_off < r_on ?!" << endl;
       cout << "   Halting..." << endl;
       exit(1);
     }
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
     X.outprm[no_i] = X.cur_nd[nnd_p];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_nd[nnd_p] = 0.0;
     X.val_nd[nnd_n] = 0.0;
     return;
   }

   x_high = X.rprm[nr_x_high];

   v_on  = X.rprm[nr_v_on ];
   r_on  = X.rprm[nr_r_on ];
   r_off = X.rprm[nr_r_off];

   vp = X.val_nd[nnd_p];
   vn = X.val_nd[nnd_n];

   if (G.flags[G.i_trns] || G.flags[G.i_startup] || G.flags[G.i_dc]) {
     l_closed = (X.val_xvr[nx_x] > (0.5*x_high));
     if (l_closed) {
       g = 1.0/r_on;
     } else {
       g = 1.0/r_off;
     }
     if (G.flags[G.i_trns] || G.flags[G.i_dc]) {
       if (G.flags[G.i_function]) {
         X.f[nf_1] = g*(vp-vn-v_on);
         X.f[nf_2] = -X.f[nf_1];
       }
       if (G.flags[G.i_jacobian]) {
         J.dfdv[nf_1][nnd_p] =  g;
         J.dfdv[nf_1][nnd_n] = -g;
         J.dfdv[nf_2][nnd_p] = -g;
         J.dfdv[nf_2][nnd_n] =  g;
       }
     }
     if (G.flags[G.i_startup]) {
       if (G.flags[G.i_function]) {
         X.h[nh_1] = g*(vp-vn-v_on);
         X.h[nh_2] = -X.h[nh_1];
       }
       if (G.flags[G.i_jacobian]) {
         J.dhdv[nh_1][nnd_p] =  g;
         J.dhdv[nh_1][nnd_n] = -g;
         J.dhdv[nh_2][nnd_p] = -g;
         J.dhdv[nh_2][nnd_n] =  g;
       }
     }
     return;
   }
   return;
}
void e_thyristor(Global &G,EbeUsr &X,EbeJac &J) {
   double vp,vn,r,g,v1a;
   bool l_closed;
   double g_in;
   int flag_on,flag1;
   double r_on,r_off,v_on,x_high,xhb2;
   double l_closed_st;
   const int nnd_anode = 0;
   const int nnd_cathode = 1;
   const int nx_g_in = 0;
   const int ni_flag_on = 0;
   const int ni_flag1 = 1;
   const int nr_r_on = 0;
   const int nr_r_off = 1;
   const int nr_v_on = 2;
   const int nr_x_high = 3;
   const int nr_xhb2 = 4;
   const int nst_l_closed_st = 0;
   const int no_i = 0;
   const int no_v = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nh_1 = 0;
   const int nh_2 = 1;
   if (G.flags[G.i_one_time_parms]) {
     x_high  = X.rprm[nr_x_high];
     xhb2 = 0.5*x_high;
     X.rprm[nr_xhb2] = xhb2;
     return;
   }
   if (G.flags[G.i_dc]) {
     cout << "thyristor.ebe: DC not allowed. Halting..." << endl;
     exit(1);
   }
   xhb2 = X.rprm[nr_xhb2];
   r_on = X.rprm[nr_r_on];
   r_off = X.rprm[nr_r_off];
   v_on = X.rprm[nr_v_on];
   vp = X.val_nd[nnd_anode];
   vn = X.val_nd[nnd_cathode];
   v1a = vp-vn;

   if (G.flags[G.i_startup]) {
     l_closed = (X.stprm[nst_l_closed_st] > xhb2);
   }
   if (G.flags[G.i_trns]) {
     l_closed = (X.iprm[ni_flag_on] == 1);
     if (!l_closed) {
       if (X.val_xvr[nx_g_in] > xhb2) {
         if (v1a >= v_on) {
           l_closed = true;
         }
       }
     } else {
       if (v1a <= v_on) {
         l_closed = false;
       }
     }
   }
   if (G.flags[G.i_startup] || G.flags[G.i_trns]) {
     if (l_closed) {
       X.iprm[ni_flag1] = 1;
       r = r_on;
     } else {
       X.iprm[ni_flag1] = 0;
       r = r_off;
     }
     if (r < 1.0e-9) {
       cout << "thyristor.ebe: r: " << r << " is too small!" << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     g = 1.0/r;
   }
   if (G.flags[G.i_trns]) {
     if (G.flags[G.i_function]) {
       if (l_closed) {
         X.f[nf_1] = g*(vp-vn-v_on);
       } else {
         X.f[nf_1] = g*(vp-vn);
       }
       X.f[nf_2] = - X.f[nf_1];
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdv[nf_1][nnd_anode  ] =  g;
       J.dfdv[nf_1][nnd_cathode] = -g;
       J.dfdv[nf_2][nnd_anode  ] = -g;
       J.dfdv[nf_2][nnd_cathode] =  g;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_function]) {
       if (l_closed) {
         X.h[nh_1] = g*(vp-vn-v_on);
       } else {
         X.h[nh_1] = g*(vp-vn);
       }
       X.h[nh_2] = - X.h[nh_1];
     }
     if (G.flags[G.i_jacobian]) {
       J.dhdv[nh_1][nnd_anode  ] =  g;
       J.dhdv[nh_1][nnd_cathode] = -g;
       J.dhdv[nh_2][nnd_anode  ] = -g;
       J.dhdv[nh_2][nnd_cathode] =  g;
     }
     return;
   }
   if (G.flags[G.i_save_history]) {
     X.iprm[ni_flag_on] = X.iprm[ni_flag1];
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_v] = X.val_nd[nnd_anode]-X.val_nd[nnd_cathode];
     X.outprm[no_i] = X.cur_nd[nnd_anode];
     return;
   }
   if (G.flags[G.i_init_guess]) {
//    Assume the switch to be closed
      X.val_nd[nnd_anode] = 0.0;
      X.val_nd[nnd_cathode] = 0.0;
      return;
   }
   return;
}
void e_voltmeter(Global &G,EbeUsr &X,EbeJac &J) {
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int no_v = 0;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nh_1 = 0;
   const int nh_2 = 1;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_dc] || G.flags[G.i_trns]) {
     if (G.flags[G.i_function]) {
       X.f[nf_1] = 0.0;
       X.f[nf_2] = 0.0;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_function]) {
       X.h[nh_1] = 0.0;
       X.h[nh_2] = 0.0;
     }
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
     return;
   }
   return;
}
void e_voltmeter_1(Global &G,EbeUsr &X,EbeJac &J) {
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int no_v = 0;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nh_1 = 0;
   const int nh_2 = 1;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_dc] || G.flags[G.i_trns]) {
     if (G.flags[G.i_function]) {
       X.f[nf_1] = 0.0;
       X.f[nf_2] = 0.0;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_function]) {
       X.h[nh_1] = 0.0;
       X.h[nh_2] = 0.0;
     }
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
     return;
   }
   return;
}
void e_voltmeter_fb(Global &G,EbeUsr &X,EbeJac &J) {
   double v_fb;
   double k_scale;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int nx_v_fb = 0;
   const int nr_k_scale = 0;
   const int no_v = 0;
   const int no_v_fb = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_dc] || G.flags[G.i_trns]) {
     k_scale = X.rprm[nr_k_scale];
     if (G.flags[G.i_function]) {
       X.f[nf_1] = 0.0;
       X.f[nf_2] = 0.0;
       X.f[nf_3] = X.val_xvr[nx_v_fb]
         -k_scale*(X.val_nd[nnd_p]-X.val_nd[nnd_n]);
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdxvr[nf_3][nx_v_fb] =  1.0;
       J.dfdv  [nf_3][nnd_p  ] = -k_scale;
       J.dfdv  [nf_3][nnd_n  ] =  k_scale;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     k_scale = X.rprm[nr_k_scale];
     if (G.flags[G.i_function]) {
       X.h[nh_1] = 0.0;
       X.h[nh_2] = 0.0;
       X.h[nh_3] = X.val_xvr[nx_v_fb]
         -k_scale*(X.val_nd[nnd_p]-X.val_nd[nnd_n]);
     }
     if (G.flags[G.i_jacobian]) {
       J.dhdxvr[nh_3][nx_v_fb] =  1.0;
       J.dhdv  [nh_3][nnd_p  ] = -k_scale;
       J.dhdv  [nh_3][nnd_n  ] =  k_scale;
     }
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
     X.outprm[no_v_fb] = X.val_xvr[nx_v_fb];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_xvr[nx_v_fb] = 0.0;
     return;
   }
   return;
}
void e_vsrc_ac(Global &G,EbeUsr &X,EbeJac &J) {
   double v0;
   double cur_p;
   double cur_p_s;
   double a,f_hz,phi_deg,t0,omega,phi_rad,vdc;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int na_cur_p = 0;
   const int nas_cur_p_s = 0;
   const int nr_a = 0;
   const int nr_f_hz = 1;
   const int nr_phi_deg = 2;
   const int nr_t0 = 3;
   const int nr_omega = 4;
   const int nr_phi_rad = 5;
   const int nr_vdc = 6;
   const int no_i = 0;
   const int no_v = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   if (G.flags[G.i_one_time_parms]) {
     f_hz    = X.rprm[nr_f_hz];
     phi_deg = X.rprm[nr_phi_deg];
     omega   = G.twopi*f_hz;
     phi_rad = G.deg_to_rad*phi_deg;

     X.rprm[nr_omega  ] = omega;
     X.rprm[nr_phi_rad] = phi_rad;
   }
   if (G.flags[G.i_dc]) {
     vdc = X.rprm[nr_vdc];
     v0 = vdc;
   }
   if (G.flags[G.i_startup] || G.flags[G.i_trns] || G.flags[G.i_init_guess]) {
     a       = X.rprm[nr_a];
     t0      = X.rprm[nr_t0];
     vdc     = X.rprm[nr_vdc];
     omega   = X.rprm[nr_omega];
     phi_rad = X.rprm[nr_phi_rad];

     v0 = a*sin(omega*(G.time_given_e-t0)+phi_rad)+vdc;
   }
   if (G.flags[G.i_dc] || G.flags[G.i_trns] ) {
     cur_p = X.val_aux[na_cur_p];
     if (G.flags[G.i_function]) {
       X.f[nf_1] =  cur_p;
       X.f[nf_2] = -cur_p;
       X.f[nf_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n]-v0;
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdaux[nf_1][na_cur_p] =  1.0;
       J.dfdaux[nf_2][na_cur_p] = -1.0;
       J.dfdv  [nf_3][nnd_p   ] =  1.0;
       J.dfdv  [nf_3][nnd_n   ] = -1.0;
     }
   }
   if (G.flags[G.i_startup]) {
     cur_p_s = X.val_auxs[nas_cur_p_s];
     if (G.flags[G.i_function]) {
       X.h[nh_1] =  cur_p_s;
       X.h[nh_2] = -cur_p_s;
       X.h[nh_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n]-v0;
     }
     if (G.flags[G.i_jacobian]) {
       J.dhdauxs[nh_1][nas_cur_p_s] =  1.0;
       J.dhdauxs[nh_2][nas_cur_p_s] = -1.0;
       J.dhdv   [nh_3][nnd_p      ] =  1.0;
       J.dhdv   [nh_3][nnd_n      ] = -1.0;
     }
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
      X.outprm[no_i] = X.cur_nd[nnd_n];
   }
   if (G.flags[G.i_init_guess]) {
     X.val_nd[nnd_n] = 0.0;
     X.val_nd[nnd_p] = 0.0;
   }
   return;
}
void e_vsrc_clock(Global &G,EbeUsr &X,EbeJac &J) {
   double y0,t0_new,delta_min,del1,del2,t_a,t_b,tnext_p;
   int n;
   double cur_p;
   double cur_p_s;
   double T1,T2,L1,L2,t0,delta1,delta2,T,L0,tk1,tk2,tk3,tk4,tk5,slope1,slope2,
     epsl;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int na_cur_p = 0;
   const int nas_cur_p_s = 0;
   const int nr_T1 = 0;
   const int nr_T2 = 1;
   const int nr_L1 = 2;
   const int nr_L2 = 3;
   const int nr_t0 = 4;
   const int nr_delta1 = 5;
   const int nr_delta2 = 6;
   const int nr_T = 7;
   const int nr_L0 = 8;
   const int nr_tk1 = 9;
   const int nr_tk2 = 10;
   const int nr_tk3 = 11;
   const int nr_tk4 = 12;
   const int nr_tk5 = 13;
   const int nr_slope1 = 14;
   const int nr_slope2 = 15;
   const int nr_epsl = 16;
   const int no_i = 0;
   const int no_v = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   if (G.flags[G.i_one_time_parms]) {
     T1 = X.rprm[nr_T1];
     T2 = X.rprm[nr_T2];
     L1 = X.rprm[nr_L1];
     L2 = X.rprm[nr_L2];
     delta1 = X.rprm[nr_delta1];
     delta2 = X.rprm[nr_delta2];

     delta_min = 0.1*min(delta1,delta2);
     del1 = T1-0.5*(delta1+delta2);
     del2 = T2-0.5*(delta1+delta2);

     if (del1 < delta_min) {
       cout << "vsrc_clock.ebe: T1 is too small. Halting..." << endl; exit(1);
     }
     if (del2 < delta_min) {
       cout << "vsrc_clock.ebe: T2 is too small. Halting..." << endl; exit(1);
     }
     T = T1 + T2;
     tk1 = 0.5*delta1;
     tk2 = T1 - 0.5*delta2;
     tk3 = T1 + 0.5*delta2;
     tk4 = T  - 0.5*delta1;
     tk5 = T  + 0.5*delta1;

     slope1 = (L1-L2)/delta1;
     slope2 = (L2-L1)/delta2;
     epsl = min(delta1,delta2)/10.0;
     L0 = 0.5*(L1+L2);

     X.rprm[nr_T] = T;
     X.rprm[nr_tk1] = tk1;
     X.rprm[nr_tk2] = tk2;
     X.rprm[nr_tk3] = tk3;
     X.rprm[nr_tk4] = tk4;
     X.rprm[nr_tk5] = tk5;
     X.rprm[nr_slope1] = slope1;
     X.rprm[nr_slope2] = slope2;
     X.rprm[nr_epsl] = epsl;
     X.rprm[nr_L0] = L0;

     return;
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
      X.outprm[no_i] = X.cur_nd[nnd_n];
      return;
   }
   if (G.flags[G.i_dc]) {
     cout << "vsrc_clock.ebe: DC not allowed. Halting..." << endl; exit(1);
   }

   t0   = X.rprm[nr_t0  ];
   T    = X.rprm[nr_T   ];
   epsl = X.rprm[nr_epsl];
   tk1  = X.rprm[nr_tk1 ];
   tk2  = X.rprm[nr_tk2 ];
   tk3  = X.rprm[nr_tk3 ];
   tk4  = X.rprm[nr_tk4 ];
   tk5  = X.rprm[nr_tk5 ];

   if (G.time_given_e < t0) {
     n = ((t0-G.time_given_e)/T) + 1;
     t0_new = t0-n*T;
   } else {
     t0_new = t0;
   }
   t_a = G.time_given_e-t0_new;
   t_b = fmod(t_a,T);

   if (abs(t_b-T) < epsl) t_b = 0.0;

   if (G.flags[G.i_next_time]) {
     if (t_b < tk1) {
       tnext_p = tk1;
     } else if (t_b < tk2) {
       tnext_p = tk2;
     } else if (t_b < tk3) {
       tnext_p = tk3;
     } else if (t_b < tk4) {
       tnext_p = tk4;
     } else {
       tnext_p = tk5;
     }
     G.time_nextbreak_e = G.time_given_e + (tnext_p-t_b);
     return;
   }

   L0 = X.rprm[nr_L0];
   L1 = X.rprm[nr_L1];
   L2 = X.rprm[nr_L2];
   slope1 = X.rprm[nr_slope1];
   slope2 = X.rprm[nr_slope2];

   if (t_b < tk1) {
     y0 = L0 + slope1*t_b;
   } else if (t_b < tk2) {
     y0 = L1;
   } else if (t_b < tk3) {
     y0 = L1 + slope2*(t_b-tk2);
   } else if (t_b < tk4) {
     y0 = L2;
   } else {
     y0 = L2 + slope1*(t_b-tk4);
   }

   if (G.flags[G.i_trns]) {
     cur_p = X.val_aux[na_cur_p];

     if (G.flags[G.i_function]) {
       X.f[nf_1] =  cur_p;
       X.f[nf_2] = -cur_p;
       X.f[nf_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n]-y0;
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdaux[nf_1][na_cur_p] =  1.0;
       J.dfdaux[nf_2][na_cur_p] = -1.0;
       J.dfdv  [nf_3][nnd_p   ] =  1.0;
       J.dfdv  [nf_3][nnd_n   ] = -1.0;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     cur_p_s = X.val_auxs[nas_cur_p_s];
     if (G.flags[G.i_function]) {
       X.h[nh_1] =  cur_p_s;
       X.h[nh_2] = -cur_p_s;
       X.h[nh_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n]-y0;
     }
     if (G.flags[G.i_jacobian]) {
       J.dhdauxs[nh_1][nas_cur_p_s] =  1.0;
       J.dhdauxs[nh_2][nas_cur_p_s] = -1.0;
       J.dhdv   [nh_3][nnd_p      ] =  1.0;
       J.dhdv   [nh_3][nnd_n      ] = -1.0;
     }
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_nd[nnd_n] = 0.0;
     X.val_nd[nnd_p] = 0.0;
     return;
   }
   return;
}
void e_vsrc_dc(Global &G,EbeUsr &X,EbeJac &J) {
   double v0,t0;
   double cur_p;
   double cur_p_s;
   double vdc,k_scale;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int na_cur_p = 0;
   const int nas_cur_p_s = 0;
   const int nr_vdc = 0;
   const int nr_k_scale = 1;
   const int no_i = 0;
   const int no_v = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_dc] || G.flags[G.i_trns]) {
     vdc      = X.rprm[nr_vdc];
     k_scale  = X.rprm[nr_k_scale];

     v0 = k_scale*vdc;
     cur_p = X.val_aux[na_cur_p];

     if (G.flags[G.i_function]) {
       X.f[nf_1] =  cur_p;
       X.f[nf_2] = -cur_p;
       X.f[nf_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n]-v0;
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdaux[nf_1][na_cur_p] =  1.0;
       J.dfdaux[nf_2][na_cur_p] = -1.0;
       J.dfdv  [nf_3][nnd_p   ] =  1.0;
       J.dfdv  [nf_3][nnd_n   ] = -1.0;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     vdc = X.rprm[nr_vdc];
     k_scale = X.rprm[nr_k_scale];
     cur_p_s = X.val_auxs[nas_cur_p_s];
     if (G.flags[G.i_function]) {
       X.h[nh_1] =  cur_p_s;
       X.h[nh_2] = -cur_p_s;
       X.h[nh_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n]-k_scale*vdc;
     }
     if (G.flags[G.i_jacobian]) {
       J.dhdauxs[nh_1][nas_cur_p_s] =  1.0;
       J.dhdauxs[nh_2][nas_cur_p_s] = -1.0;
       J.dhdv   [nh_3][nnd_p      ] =  1.0;
       J.dhdv   [nh_3][nnd_n      ] = -1.0;
     }
     return;
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
      X.outprm[no_i] = X.cur_nd[nnd_n];
      return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_nd[nnd_n] = 0.0;
     X.val_nd[nnd_p] = 0.0;
     return;
   }
   return;
}
void e_vsrc_pulse10(Global &G,EbeUsr &X,EbeJac &J) {
   int intrvl;
   double time0,g_value,g_1,g_2,delt_trns,delt_1,t_trns;
   bool l_high;
   double cur_p;
   double cur_p_s;
   int i0,n1;
   double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,
     t20,y_low,y_high,t_rise,t_fall,epsl;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int na_cur_p = 0;
   const int nas_cur_p_s = 0;
   const int ni_i0 = 0;
   const int ni_n1 = 1;
   const int nr_t1 = 0;
   const int nr_t2 = 1;
   const int nr_t3 = 2;
   const int nr_t4 = 3;
   const int nr_t5 = 4;
   const int nr_t6 = 5;
   const int nr_t7 = 6;
   const int nr_t8 = 7;
   const int nr_t9 = 8;
   const int nr_t10 = 9;
   const int nr_t11 = 10;
   const int nr_t12 = 11;
   const int nr_t13 = 12;
   const int nr_t14 = 13;
   const int nr_t15 = 14;
   const int nr_t16 = 15;
   const int nr_t17 = 16;
   const int nr_t18 = 17;
   const int nr_t19 = 18;
   const int nr_t20 = 19;
   const int nr_y_low = 20;
   const int nr_y_high = 21;
   const int nr_t_rise = 22;
   const int nr_t_fall = 23;
   const int nr_epsl = 24;
   const int no_i = 0;
   const int no_v = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   if (G.flags[G.i_one_time_parms]) {
     t_rise = X.rprm[nr_t_rise];
     t_fall = X.rprm[nr_t_fall];
     epsl=0.02*min(t_rise,t_fall);
     X.rprm[nr_epsl] = epsl;
     return;
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
      X.outprm[no_i] = X.cur_nd[nnd_n];
      return;
   }
   if (G.flags[G.i_dc]) {
     cout << "vsrc_pulse10.ebe: DC not allowed. Halting..." << endl; exit(1);
   }
   y_low  = X.rprm[nr_y_low ];
   y_high = X.rprm[nr_y_high];
   t_rise = X.rprm[nr_t_rise];
   t_fall = X.rprm[nr_t_fall];
   epsl   = X.rprm[nr_epsl  ];

   i0 = X.iprm[ni_i0];
   n1 = X.iprm[ni_n1];

   time0 = G.time_given_e;

   intrvl = n1+1;
   for (int i=0; i < n1; i++) {
      if (X.rprm[i] >= time0) {
         intrvl = i+1;
         break;
      }
   }
   if (intrvl % 2 == 1) {
      l_high = (i0 != 0);
   } else {
      l_high = (i0 == 0);
   }
   if (l_high) {
      g_1 = y_low;
      g_2 = y_high;
      delt_trns = t_fall;
      if (intrvl > 1) delt_1 = t_rise;
   } else {
      g_1 = y_high;
      g_2 = y_low;
      delt_trns = t_rise;
      if (intrvl > 1) delt_1 = t_fall;
   }

   if (G.flags[G.i_next_time]) {
     G.time_nextbreak_e = G.time_end;
     if (intrvl == 1) {
       if ((X.rprm[0]-time0) < epsl) {
         G.time_nextbreak_e = X.rprm[0] + delt_trns;
       } else {
         G.time_nextbreak_e = X.rprm[0];
       }
       return;
     } else if (intrvl == (n1+1)) {
       if ((X.rprm[intrvl-2] + delt_1-time0) > epsl) {
         G.time_nextbreak_e = X.rprm[intrvl-2] + delt_1;
       } else {
         G.time_nextbreak_e = G.time_end;
       }
     } else {
       t_trns = X.rprm[intrvl];

       if ((X.rprm[intrvl-2] + delt_1-time0) > epsl) {
         G.time_nextbreak_e = X.rprm[intrvl-2] + delt_1;
       } else if ((X.rprm[intrvl-1]-time0) < epsl) {
         G.time_nextbreak_e = X.rprm[intrvl-1] + delt_trns;
       } else {
         G.time_nextbreak_e = X.rprm[intrvl-1];
       }
     }
     return;
   }

   if (intrvl == 1) {
      g_value = g_2;
   } else {
      t_trns = X.rprm[intrvl-2];
      if (time0 < (t_trns+delt_1)) {
         g_value = g_1 + (g_2-g_1)*(time0-t_trns)/delt_1;
      } else {
         g_value = g_2;
      }
   }

   if (G.flags[G.i_trns]) {
     cur_p = X.val_aux[na_cur_p];

     if (G.flags[G.i_function]) {
       X.f[nf_1] =  cur_p;
       X.f[nf_2] = -cur_p;
       X.f[nf_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n]-g_value;
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdaux[nf_1][na_cur_p] =  1.0;
       J.dfdaux[nf_2][na_cur_p] = -1.0;
       J.dfdv  [nf_3][nnd_p   ] =  1.0;
       J.dfdv  [nf_3][nnd_n   ] = -1.0;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     cur_p_s = X.val_auxs[nas_cur_p_s];
     if (G.flags[G.i_function]) {
       X.h[nh_1] =  cur_p_s;
       X.h[nh_2] = -cur_p_s;
       X.h[nh_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n]-g_value;
     }
     if (G.flags[G.i_jacobian]) {
       J.dhdauxs[nh_1][nas_cur_p_s] =  1.0;
       J.dhdauxs[nh_2][nas_cur_p_s] = -1.0;
       J.dhdv   [nh_3][nnd_p      ] =  1.0;
       J.dhdv   [nh_3][nnd_n      ] = -1.0;
     }
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_nd[nnd_n] = 0.0;
     X.val_nd[nnd_p] = 0.0;
     return;
   }
   return;
}
void e_vsrc_x(Global &G,EbeUsr &X,EbeJac &J) {
   double v0,t0;
   double cur_p;
   double cur_p_s;
   double x_in;
   double k_scale;
   const int nnd_p = 0;
   const int nnd_n = 1;
   const int na_cur_p = 0;
   const int nas_cur_p_s = 0;
   const int nx_x_in = 0;
   const int nr_k_scale = 0;
   const int no_i = 0;
   const int no_v = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_dc] || G.flags[G.i_trns]) {
     k_scale = X.rprm[nr_k_scale];
     v0 = k_scale*X.val_xvr[nx_x_in];
     cur_p = X.val_aux[na_cur_p];

     if (G.flags[G.i_function]) {
       X.f[nf_1] =  cur_p;
       X.f[nf_2] = -cur_p;
       X.f[nf_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n]-v0;
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdaux[nf_1][na_cur_p] =  1.0;
       J.dfdaux[nf_2][na_cur_p] = -1.0;
       J.dfdv  [nf_3][nnd_p   ] =  1.0;
       J.dfdv  [nf_3][nnd_n   ] = -1.0;
       J.dfdxvr[nf_3][nx_x_in ] = -k_scale;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     k_scale = X.rprm[nr_k_scale];
     v0 = k_scale*X.val_xvr[nx_x_in];
     cur_p_s = X.val_auxs[nas_cur_p_s];
     if (G.flags[G.i_function]) {
       X.h[nh_1] =  cur_p_s;
       X.h[nh_2] = -cur_p_s;
       X.h[nh_3] = X.val_nd[nnd_p]-X.val_nd[nnd_n]-v0;
     }
     if (G.flags[G.i_jacobian]) {
       J.dhdauxs[nh_1][nas_cur_p_s] =  1.0;
       J.dhdauxs[nh_2][nas_cur_p_s] = -1.0;
       J.dhdv   [nh_3][nnd_p      ] =  1.0;
       J.dhdv   [nh_3][nnd_n      ] = -1.0;
       J.dhdxvr [nh_3][nx_x_in    ] = -k_scale;
     }
     return;
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_v] = X.val_nd[nnd_p]-X.val_nd[nnd_n];
      X.outprm[no_i] = X.cur_nd[nnd_n];
      return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_nd[nnd_n] = 0.0;
     X.val_nd[nnd_p] = 0.0;
     return;
   }
   return;
}
void e_xfmr_l1l2(Global &G,EbeUsr &X,EbeJac &J) {
   double v1a,v2a;
   int nfp_p_p=0;
   int nfp_s_p=1;
   int nfp_p_n=2;
   int nfp_s_n=3;
   double i1,i2,i1d,i2d;
   double l1,l2,k,m;
   double i10,i20;
   const int nnd_p_p = 0;
   const int nnd_s_p = 1;
   const int nnd_p_n = 2;
   const int nnd_s_n = 3;
   const int na_i1 = 0;
   const int na_i2 = 1;
   const int na_i1d = 2;
   const int na_i2d = 3;
   const int nr_l1 = 0;
   const int nr_l2 = 1;
   const int nr_k = 2;
   const int nr_m = 3;
   const int nst_i10 = 0;
   const int nst_i20 = 1;
   const int no_ip = 0;
   const int no_is = 1;
   const int no_vp = 2;
   const int no_vs = 3;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nf_4 = 3;
   const int nf_5 = 4;
   const int nf_6 = 5;
   const int nf_7 = 6;
   const int nf_8 = 7;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   const int nh_4 = 3;
   if (G.flags[G.i_dc]) {
     cout << "xfmr_basic.ebe: dc not implemented." << endl;
     cout << "  Halting..." << endl; exit(1);
   }
   l1 = X.rprm[nr_l1];
   l2 = X.rprm[nr_l2];
   m  = X.rprm[nr_m ];

   if (G.flags[G.i_one_time_parms]) {
     k = X.rprm[nr_k];
     m = k*sqrt(l1*l2);
     X.rprm[nr_m] = m;
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_aux[na_i1 ] = 0.0;
     X.val_aux[na_i2 ] = 0.0;
     X.val_aux[na_i1d] = 0.0;
     X.val_aux[na_i2d] = 0.0;
     return;
   }
   if (G.flags[G.i_trns]) {
     i1  = X.val_aux[na_i1 ];
     i2  = X.val_aux[na_i2 ];
     i1d = X.val_aux[na_i1d];
     i2d = X.val_aux[na_i2d];

     v1a = X.val_nd[nnd_p_p]-X.val_nd[nnd_p_n];
     v2a = X.val_nd[nnd_s_p]-X.val_nd[nnd_s_n];

     if (G.flags[G.i_function]) {
       X.f[nfp_p_p] =  i1;
       X.f[nfp_p_n] = -i1;
       X.f[nfp_s_p] =  i2;
       X.f[nfp_s_n] = -i2;
       X.f[nf_5] = i1d;
       X.f[nf_6] = i2d;
       X.f[nf_7] = v1a - l1*i1d - m*i2d;
       X.f[nf_8] = v2a - l2*i2d - m*i1d;
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdaux[nfp_p_p][na_i1] =  1.0;
       J.dfdaux[nfp_p_n][na_i1] = -1.0;
       J.dfdaux[nfp_s_p][na_i2] =  1.0;
       J.dfdaux[nfp_s_n][na_i2] = -1.0;

       J.dfdaux[nf_5][na_i1d] = 1.0;
       J.dfdaux[nf_6][na_i2d] = 1.0;

       J.dfdv  [nf_7][nnd_p_p] =  1.0;
       J.dfdv  [nf_7][nnd_p_n] = -1.0;
       J.dfdaux[nf_7][na_i1d] = -l1;
       J.dfdaux[nf_7][na_i2d] = -m;

       J.dfdv  [nf_8][nnd_s_p] =  1.0;
       J.dfdv  [nf_8][nnd_s_n] = -1.0;
       J.dfdaux[nf_8][na_i2d] = -l2;
       J.dfdaux[nf_8][na_i1d] = -m;
     }
   }
   if (G.flags[G.i_startup]) {
      i10 = X.stprm[nst_i10];
      i20 = X.stprm[nst_i20];
      cout << "xfmr_l1l2.ebe: startup not implemented. Halting..." << endl;
      exit(1);
      if (G.flags[G.i_function]) {
         X.h[nh_1] =  i10;
         X.h[nh_2] = -i10;
         X.h[nh_3] =  i20;
         X.h[nh_4] = -i20;
      }
      X.val_aux[na_i1] = i10;
      X.val_aux[na_i2] = i20;
      return;
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_vp] = X.val_nd[nnd_p_p]-X.val_nd[nnd_p_n];
      X.outprm[no_vs] = X.val_nd[nnd_s_p]-X.val_nd[nnd_s_n];
      X.outprm[no_ip] = X.cur_nd[nnd_p_p];
      X.outprm[no_is] = X.cur_nd[nnd_s_p];
      return;
   }
   return;
}
void e_xfmr_level0_1ph(Global &G,EbeUsr &X,EbeJac &J) {
   int nfp_p_p=0;
   int nfp_s_p=1;
   int nfp_p_n=2;
   int nfp_s_n=3;
   double cur_p_p,cur_s_p;
   double p_turns,s_turns;
   double ip0,is0;
   const int nnd_p_p = 0;
   const int nnd_s_p = 1;
   const int nnd_p_n = 2;
   const int nnd_s_n = 3;
   const int na_cur_p_p = 0;
   const int na_cur_s_p = 1;
   const int nr_p_turns = 0;
   const int nr_s_turns = 1;
   const int nst_ip0 = 0;
   const int nst_is0 = 1;
   const int no_ip = 0;
   const int no_is = 1;
   const int no_vp = 2;
   const int no_vs = 3;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nf_4 = 3;
   const int nf_5 = 4;
   const int nf_6 = 5;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   const int nh_4 = 3;
   if (G.flags[G.i_dc]) {
     cout << "xfmr_level1_1ph.ebe: dc not implemented." << endl;
     cout << "  Halting..." << endl; exit(1);
   }
   if (G.flags[G.i_init_guess]) {
     X.val_aux[na_cur_p_p] = 0.0;
     X.val_aux[na_cur_s_p] = 0.0;
     return;
   }
   if (G.flags[G.i_trns]) {
     p_turns = X.rprm[nr_p_turns];
     s_turns = X.rprm[nr_s_turns];

     cur_p_p = X.val_aux[na_cur_p_p];
     cur_s_p = X.val_aux[na_cur_s_p];

     if (G.flags[G.i_function]) {
       X.f[nfp_p_p] =  cur_p_p;
       X.f[nfp_p_n] = -cur_p_p;
       X.f[nfp_s_p] =  cur_s_p;
       X.f[nfp_s_n] = -cur_s_p;

       X.f[nf_5] = (X.val_nd[nnd_p_p]-X.val_nd[nnd_p_n])/p_turns -
                 (X.val_nd[nnd_s_p]-X.val_nd[nnd_s_n])/s_turns;
       X.f[nf_6] = p_turns*cur_p_p + s_turns*cur_s_p;
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdaux[nfp_p_p][na_cur_p_p] =  1.0;
       J.dfdaux[nfp_p_n][na_cur_p_p] = -1.0;
       J.dfdaux[nfp_s_p][na_cur_s_p] =  1.0;
       J.dfdaux[nfp_s_n][na_cur_s_p] = -1.0;

       J.dfdv[nf_5][nnd_p_p] =  1.0/p_turns;
       J.dfdv[nf_5][nnd_p_n] = -1.0/p_turns;
       J.dfdv[nf_5][nnd_s_p] = -1.0/s_turns;
       J.dfdv[nf_5][nnd_s_n] =  1.0/s_turns;

       J.dfdaux[nf_6][na_cur_p_p] = p_turns;
       J.dfdaux[nf_6][na_cur_s_p] = s_turns;
     }
   }
   if (G.flags[G.i_startup]) {
      ip0 = X.stprm[nst_ip0];
      is0 = X.stprm[nst_is0];
      if (G.flags[G.i_function]) {
         X.h[nh_1] =  ip0;
         X.h[nh_2] = -ip0;
         X.h[nh_3] =  is0;
         X.h[nh_4] = -is0;
      }
      cout << "xfmr_level0_1ph.ebe: startup not implemented. Halting..." << endl;
      exit(1);
      X.val_aux[na_cur_p_p] = ip0;
      X.val_aux[na_cur_s_p] = is0;
      return;
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_vp] = X.val_nd[nnd_p_p]-X.val_nd[nnd_p_n];
      X.outprm[no_vs] = X.val_nd[nnd_s_p]-X.val_nd[nnd_s_n];
      X.outprm[no_ip] = X.cur_nd[nnd_p_p];
      X.outprm[no_is] = X.cur_nd[nnd_s_p];
      return;
   }
   return;
}
void e_xfmr_level2_1ph(Global &G,EbeUsr &X,EbeJac &J) {
   int nfp_p_p=0;
   int nfp_s_p=1;
   int nfp_p_n=2;
   int nfp_s_n=3;
   double cur_p_p,cur_s_p,im;
   double p_turns,s_turns,lm;
   double ip0,is0;
   const int nnd_p_p = 0;
   const int nnd_s_p = 1;
   const int nnd_p_n = 2;
   const int nnd_s_n = 3;
   const int na_cur_p_p = 0;
   const int na_cur_s_p = 1;
   const int na_im = 2;
   const int nr_p_turns = 0;
   const int nr_s_turns = 1;
   const int nr_lm = 2;
   const int nst_ip0 = 0;
   const int nst_is0 = 1;
   const int no_ip = 0;
   const int no_is = 1;
   const int no_vp = 2;
   const int no_vs = 3;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nf_4 = 3;
   const int nf_5 = 4;
   const int nf_6 = 5;
   const int nf_7 = 6;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   const int nh_4 = 3;
   if (G.flags[G.i_dc]) {
     cout << "xfmr_level2_1ph.ebe: dc not implemented." << endl;
     cout << "  Halting..." << endl; exit(1);
   }
   if (G.flags[G.i_init_guess]) {
     X.val_aux[na_cur_p_p] = 0.0;
     X.val_aux[na_cur_s_p] = 0.0;
     X.val_aux[na_im     ] = 0.0;
     return;
   }
   if (G.flags[G.i_trns]) {
     p_turns = X.rprm[nr_p_turns];
     s_turns = X.rprm[nr_s_turns];
     lm      = X.rprm[nr_lm     ];

     cur_p_p = X.val_aux[na_cur_p_p];
     cur_s_p = X.val_aux[na_cur_s_p];
     im      = X.val_aux[na_im     ];

     if (G.flags[G.i_function]) {
       X.f[nfp_p_p] =  cur_p_p;
       X.f[nfp_p_n] = -cur_p_p;
       X.f[nfp_s_p] =  cur_s_p;
       X.f[nfp_s_n] = -cur_s_p;

       X.f[nf_5] = (X.val_nd[nnd_p_p]-X.val_nd[nnd_p_n])/p_turns -
                   (X.val_nd[nnd_s_p]-X.val_nd[nnd_s_n])/s_turns;
       X.f[nf_6] = p_turns*(cur_p_p-im) + s_turns*cur_s_p;
       X.f[nf_7] = (X.val_nd[nnd_p_p]-X.val_nd[nnd_p_n])/lm;
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdaux[nfp_p_p][na_cur_p_p] =  1.0;
       J.dfdaux[nfp_p_n][na_cur_p_p] = -1.0;
       J.dfdaux[nfp_s_p][na_cur_s_p] =  1.0;
       J.dfdaux[nfp_s_n][na_cur_s_p] = -1.0;

       J.dfdv[nf_5][nnd_p_p] =  1.0/p_turns;
       J.dfdv[nf_5][nnd_p_n] = -1.0/p_turns;
       J.dfdv[nf_5][nnd_s_p] = -1.0/s_turns;
       J.dfdv[nf_5][nnd_s_n] =  1.0/s_turns;

       J.dfdaux[nf_6][na_cur_p_p] =  p_turns;
       J.dfdaux[nf_6][na_im     ] = -p_turns;
       J.dfdaux[nf_6][na_cur_s_p] =  s_turns;

       J.dfdv[nf_7][nnd_p_p] =  1.0/lm;
       J.dfdv[nf_7][nnd_p_n] = -1.0/lm;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
      ip0 = X.stprm[nst_ip0];
      is0 = X.stprm[nst_is0];
      cout << "xfmr_level2_1ph.ebe: startup not implemented. Halting..." << endl;
      exit(1);
      if (G.flags[G.i_function]) {
         X.h[nh_1] =  ip0;
         X.h[nh_2] = -ip0;
         X.h[nh_3] =  is0;
         X.h[nh_4] = -is0;
      }
      X.val_aux[na_cur_p_p] = ip0;
      X.val_aux[na_cur_s_p] = is0;
      return;
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_vp] = X.val_nd[nnd_p_p]-X.val_nd[nnd_p_n];
      X.outprm[no_vs] = X.val_nd[nnd_s_p]-X.val_nd[nnd_s_n];
      X.outprm[no_ip] = X.cur_nd[nnd_p_p];
      X.outprm[no_is] = X.cur_nd[nnd_s_p];
      return;
   }
   return;
}
void e_xfmr_level0_1ph_1_2(Global &G,EbeUsr &X,EbeJac &J) {
   int nfp_s1_n=0;
   int nfp_s2_p=1;
   int nfp_p_p =2;
   int nfp_s1_p=3;
   int nfp_p_n =4;
   int nfp_s2_n=5;
   double cur_p_p,cur_s1_p,cur_s2_p;
   double p_turns,s1_turns,s2_turns;
   double ip0,is10,is20;
   const int nnd_s1_n = 0;
   const int nnd_s2_p = 1;
   const int nnd_p_p = 2;
   const int nnd_s1_p = 3;
   const int nnd_p_n = 4;
   const int nnd_s2_n = 5;
   const int na_cur_p_p = 0;
   const int na_cur_s1_p = 1;
   const int na_cur_s2_p = 2;
   const int nr_p_turns = 0;
   const int nr_s1_turns = 1;
   const int nr_s2_turns = 2;
   const int nst_ip0 = 0;
   const int nst_is10 = 1;
   const int nst_is20 = 2;
   const int no_ip = 0;
   const int no_is1 = 1;
   const int no_is2 = 2;
   const int no_vp = 3;
   const int no_vs1 = 4;
   const int no_vs2 = 5;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nf_4 = 3;
   const int nf_5 = 4;
   const int nf_6 = 5;
   const int nf_7 = 6;
   const int nf_8 = 7;
   const int nf_9 = 8;
   const int nh_1 = 0;
   const int nh_2 = 1;
   const int nh_3 = 2;
   const int nh_4 = 3;
   const int nh_5 = 4;
   const int nh_6 = 5;
   if (G.flags[G.i_dc]) {
     cout << "xfmr_level0_1ph_1_2.ebe: dc not implemented." << endl;
     cout << "  Halting..." << endl; exit(1);
   }
   if (G.flags[G.i_init_guess]) {
     X.val_aux[na_cur_p_p ] = 0.0;
     X.val_aux[na_cur_s1_p] = 0.0;
     X.val_aux[na_cur_s2_p] = 0.0;
     return;
   }
   if (G.flags[G.i_trns]) {
     p_turns  = X.rprm[nr_p_turns ];
     s1_turns = X.rprm[nr_s1_turns];
     s2_turns = X.rprm[nr_s2_turns];

     cur_p_p  = X.val_aux[na_cur_p_p ];
     cur_s1_p = X.val_aux[na_cur_s1_p];
     cur_s2_p = X.val_aux[na_cur_s2_p];

     if (G.flags[G.i_function]) {
       X.f[nfp_p_p ] =  cur_p_p;
       X.f[nfp_p_n ] = -cur_p_p;
       X.f[nfp_s1_p] =  cur_s1_p;
       X.f[nfp_s1_n] = -cur_s1_p;
       X.f[nfp_s2_p] =  cur_s2_p;
       X.f[nfp_s2_n] = -cur_s2_p;

       X.f[nf_7] = (X.val_nd[nnd_p_p]-X.val_nd[nnd_p_n])/p_turns -
                   (X.val_nd[nnd_s1_p]-X.val_nd[nnd_s1_n])/s1_turns;
       X.f[nf_8] = (X.val_nd[nnd_p_p]-X.val_nd[nnd_p_n])/p_turns -
                   (X.val_nd[nnd_s2_p]-X.val_nd[nnd_s2_n])/s2_turns;
       X.f[nf_9] = p_turns*cur_p_p +
                   s1_turns*cur_s1_p +
                   s2_turns*cur_s2_p;
     }
     if (G.flags[G.i_jacobian]) {
       J.dfdaux[nfp_p_p ][na_cur_p_p ] =  1.0;
       J.dfdaux[nfp_p_n ][na_cur_p_p ] = -1.0;
       J.dfdaux[nfp_s1_p][na_cur_s1_p] =  1.0;
       J.dfdaux[nfp_s1_n][na_cur_s1_p] = -1.0;
       J.dfdaux[nfp_s2_p][na_cur_s2_p] =  1.0;
       J.dfdaux[nfp_s2_n][na_cur_s2_p] = -1.0;

       J.dfdv[nf_7][nnd_p_p] =  1.0/p_turns;
       J.dfdv[nf_7][nnd_p_n] = -1.0/p_turns;
       J.dfdv[nf_7][nnd_s1_p] = -1.0/s1_turns;
       J.dfdv[nf_7][nnd_s1_n] =  1.0/s1_turns;

       J.dfdv[nf_8][nnd_p_p] =  1.0/p_turns;
       J.dfdv[nf_8][nnd_p_n] = -1.0/p_turns;
       J.dfdv[nf_8][nnd_s2_p] = -1.0/s2_turns;
       J.dfdv[nf_8][nnd_s2_n] =  1.0/s2_turns;

       J.dfdaux[nf_9][na_cur_p_p] = p_turns;
       J.dfdaux[nf_9][na_cur_s1_p] = s1_turns;
       J.dfdaux[nf_9][na_cur_s2_p] = s2_turns;
     }
   }
   if (G.flags[G.i_startup]) {
      ip0  = X.stprm[nst_ip0 ];
      is10 = X.stprm[nst_is10];
      is20 = X.stprm[nst_is20];

      cout << "xfmr_level0_1ph_1_2.ebe: startup not implemented. Halting..." << endl;
      exit(1);
      if (G.flags[G.i_function]) {
         X.h[nh_1] =  ip0;
         X.h[nh_2] = -ip0;
         X.h[nh_3] =  is10;
         X.h[nh_4] = -is10;
         X.h[nh_5] =  is20;
         X.h[nh_6] = -is20;
      }
      X.val_aux[na_cur_p_p ] = ip0;
      X.val_aux[na_cur_s1_p] = is10;
      X.val_aux[na_cur_s2_p] = is20;
      return;
   }
   if (G.flags[G.i_outvar]) {
      X.outprm[no_vp ] = X.val_nd[nnd_p_p ]-X.val_nd[nnd_p_n ];
      X.outprm[no_vs1] = X.val_nd[nnd_s1_p]-X.val_nd[nnd_s1_n];
      X.outprm[no_vs2] = X.val_nd[nnd_s2_p]-X.val_nd[nnd_s2_n];
      X.outprm[no_ip ] = X.cur_nd[nnd_p_p ];
      X.outprm[no_is1] = X.cur_nd[nnd_s1_p];
      X.outprm[no_is2] = X.cur_nd[nnd_s2_p];
      return;
   }
   return;
}
