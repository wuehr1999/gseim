/*
Copyright (C) 2022 - Mahesh Patil <mbpatil@ee.iitb.ac.in>
This file is part of GSEIM.

GSEIM is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <vector>

using namespace std;

void user_function(
   const int index_fn,
   const double time0,
   double* x,
   double* y,
   vector<int> &iprm,
   vector<double> &rprm);
#include "global.h"
#include "xbeusr.h"
#include "xbejac.h"
#include "utils.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
void x_abc_to_alphabeta_3(Global &G,XbeUsr &X,XbeJac &J) {
   static double k1=0.0, k2=0.0, k3=0.0;
   double a,b,c;
   double alpha,beta;
   const int nvr_a = 0;
   const int nvr_b = 1;
   const int nvr_c = 2;
   const int nvr_alpha = 3;
   const int nvr_beta = 4;
   const int no_a = 0;
   const int no_b = 1;
   const int no_c = 2;
   const int no_alpha = 3;
   const int no_beta = 4;
   const int ng_1 = 0;
   const int ng_2 = 1;
// cout << "abc_to_alphabeta_1.xbe" << endl;

   if (G.flags[G.i_one_time_parms]) {
     k1 = sqrt(2.0/3.0);
     k2 = sqrt(1.0/6.0);
     k3 = 1.0/(sqrt(2.0));
     return;
   }
   if (G.flags[G.i_init_guess]) {
     a = X.val_vr[nvr_a];
     b = X.val_vr[nvr_b];
     c = X.val_vr[nvr_c];

     X.val_vr[nvr_alpha] = k1*a - k2*(b+c);
     X.val_vr[nvr_beta ] = k3*(b-c);

     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_a] = X.val_vr[nvr_a];
     X.outprm[no_b] = X.val_vr[nvr_b];
     X.outprm[no_c] = X.val_vr[nvr_c];
     X.outprm[no_alpha] = X.val_vr[nvr_alpha];
     X.outprm[no_beta ] = X.val_vr[nvr_beta ];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     a = X.val_vr[nvr_a];
     b = X.val_vr[nvr_b];
     c = X.val_vr[nvr_c];

     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_alpha] = k1*a - k2*(b+c);
       X.val_vr[nvr_beta ] = k3*(b-c);
     } else if (G.flags[G.i_implicit]) {
       alpha = X.val_vr[nvr_alpha];
       beta  = X.val_vr[nvr_beta ];
       if (G.flags[G.i_function]) {
         X.g[ng_1] = alpha - k1*a + k2*(b+c);
         X.g[ng_2] = beta - k3*(b-c);
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_alpha] = 1.0;
         J.dgdvr[ng_1][nvr_a    ] = -k1;
         J.dgdvr[ng_1][nvr_b    ] =  k2;
         J.dgdvr[ng_1][nvr_c    ] =  k2;

         J.dgdvr[ng_2][nvr_beta] = 1.0;
         J.dgdvr[ng_2][nvr_b   ] = -k3;
         J.dgdvr[ng_2][nvr_c   ] =  k3;
       }
     }
     return;
   }
   return;
}
void x_abc_to_dq(Global &G,XbeUsr &X,XbeJac &J) {
   static double k1=0,k2=0,k3=0;
   double van,vbn,vcn;
   double vqs,vds;
   const int nvr_van = 0;
   const int nvr_vbn = 1;
   const int nvr_vcn = 2;
   const int nvr_vqs = 3;
   const int nvr_vds = 4;
   const int no_van = 0;
   const int no_vbn = 1;
   const int no_vcn = 2;
   const int no_vqs = 3;
   const int no_vds = 4;
   const int ng_1 = 0;
   const int ng_2 = 1;
// cout << "abc_to_dq.xbe" << endl;

   if (G.flags[G.i_one_time_parms]) {
     k1 = 2.0/3.0;
     k2 = 1.0/3.0;
     k3 = 1.0/(sqrt(3.0));
     return;
   }
   if (G.flags[G.i_init_guess]) {
     van = X.val_vr[nvr_van];
     vbn = X.val_vr[nvr_vbn];
     vcn = X.val_vr[nvr_vcn];

     X.val_vr[nvr_vqs] = k1*van - k2*(vbn + vcn);
     X.val_vr[nvr_vds] = - k3*(vbn - vcn);
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_van] = X.val_vr[nvr_van];
     X.outprm[no_vbn] = X.val_vr[nvr_vbn];
     X.outprm[no_vcn] = X.val_vr[nvr_vcn];
     X.outprm[no_vqs] = X.val_vr[nvr_vqs];
     X.outprm[no_vds] = X.val_vr[nvr_vds];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     van = X.val_vr[nvr_van];
     vbn = X.val_vr[nvr_vbn];
     vcn = X.val_vr[nvr_vcn];
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_vqs] = k1*van - k2*(vbn + vcn);
       X.val_vr[nvr_vds] = - k3*(vbn - vcn);
     } else if (G.flags[G.i_implicit]) {
       vqs = X.val_vr[nvr_vqs];
       vds = X.val_vr[nvr_vds];
       if (G.flags[G.i_function]) {
         X.g[ng_1] = vqs -(k1*van - k2*(vbn + vcn));
         X.g[ng_2] = vds -(- k3*(vbn - vcn));
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_vqs] = 1.0;
         J.dgdvr[ng_1][nvr_van] = -k1;
         J.dgdvr[ng_1][nvr_vbn] =  k2;
         J.dgdvr[ng_1][nvr_vcn] =  k2;

         J.dgdvr[ng_2][nvr_vds] = 1.0;
         J.dgdvr[ng_2][nvr_vbn] =  k3;
         J.dgdvr[ng_2][nvr_vcn] = -k3;
       }
     }
     return;
   }
   return;
}
void x_abc_to_dq0_2(Global &G,XbeUsr &X,XbeJac &J) {
   static double k1=0,k2=0,k3=0;
   double c1,s1;
   double xa,xb,xc,theta;
   double xd,xq,x0;
   double c,s;
   const int nvr_xa = 0;
   const int nvr_xb = 1;
   const int nvr_xc = 2;
   const int nvr_theta = 3;
   const int nvr_xd = 4;
   const int nvr_xq = 5;
   const int nvr_x0 = 6;
   const int no_xa = 0;
   const int no_xb = 1;
   const int no_xc = 2;
   const int no_xd = 3;
   const int no_xq = 4;
   const int no_x0 = 5;
   const int na_c = 0;
   const int na_s = 1;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   const int ng_4 = 3;
   const int ng_5 = 4;
   if (G.flags[G.i_one_time_parms]) {
     k1 = 1.0/3.0;
     k2 = 1.0/(sqrt(3.0));
     k3 = 2.0/3.0;
     return;
   }
   if (G.flags[G.i_init_guess]) {
     xa    = X.val_vr[nvr_xa   ];
     xb    = X.val_vr[nvr_xb   ];
     xc    = X.val_vr[nvr_xc   ];
     theta = X.val_vr[nvr_theta];

     c = cos(theta);
     s = sin(theta);

     X.val_aux[na_c] = c;
     X.val_aux[na_s] = s;

     X.val_vr[nvr_xd] = - xa*(-k3*c) - xb*(k1*c-k2*s) - xc*( k1*c+k2*s);
     X.val_vr[nvr_xq] = - xa*(-k3*s) - xb*(k2*c+k1*s) - xc*(-k2*c+k1*s);
     X.val_vr[nvr_x0] = k1*(xa + xb + xc);

     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_xa] = X.val_vr[nvr_xa];
     X.outprm[no_xb] = X.val_vr[nvr_xb];
     X.outprm[no_xc] = X.val_vr[nvr_xc];
     X.outprm[no_xd] = X.val_vr[nvr_xd];
     X.outprm[no_xq] = X.val_vr[nvr_xq];
     X.outprm[no_x0] = X.val_vr[nvr_x0];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     xa    = X.val_vr[nvr_xa   ];
     xb    = X.val_vr[nvr_xb   ];
     xc    = X.val_vr[nvr_xc   ];
     theta = X.val_vr[nvr_theta];

     c1 = cos(theta);
     s1 = sin(theta);

     if (G.flags[G.i_explicit]) {
       X.val_aux[na_c] = c1;
       X.val_aux[na_s] = s1;

       X.val_vr[nvr_xd] = - xa*(-k3*c) - xb*(k1*c-k2*s) - xc*( k1*c+k2*s);
       X.val_vr[nvr_xq] = - xa*(-k3*s) - xb*(k2*c+k1*s) - xc*(-k2*c+k1*s);
       X.val_vr[nvr_x0] = k1*(xa + xb + xc);
     } else if (G.flags[G.i_implicit]) {
       xd = X.val_vr[nvr_xd];
       xq = X.val_vr[nvr_xq];
       x0 = X.val_vr[nvr_x0];

       c = X.val_aux[na_c];
       s = X.val_aux[na_s];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = xd + xa*(-k3*c) + xb*(k1*c-k2*s) +xc*( k1*c+k2*s);
         X.g[ng_2] = xq + xa*(-k3*s) + xb*(k2*c+k1*s) +xc*(-k2*c+k1*s);
         X.g[ng_3] = x0 + xa*(-k1) + xb*(-k1) +xc*(-k1);
         X.g[ng_4] = c - c1;
         X.g[ng_5] = s - s1;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr [ng_1][nvr_xd] = 1.0;
         J.dgdvr [ng_1][nvr_xa] = -k3*c;
         J.dgdvr [ng_1][nvr_xb] = k1*c-k2*s;
         J.dgdvr [ng_1][nvr_xc] = k1*c+k2*s;
         J.dgdaux[ng_1][na_c  ] = -xa*k3 + (xb+xc)*k1;
         J.dgdaux[ng_1][na_s  ] = -xb*k2 + xc*k2;

         J.dgdvr [ng_2][nvr_xq] = 1.0;
         J.dgdvr [ng_2][nvr_xa] = -k3*s;
         J.dgdvr [ng_2][nvr_xb] = k2*c+k1*s;
         J.dgdvr [ng_2][nvr_xc] = -k2*c+k1*s;
         J.dgdaux[ng_2][na_c  ] = (xb-xc)*k2;
         J.dgdaux[ng_2][na_s  ] = -xa*k3 + (xb+xc)*k1;

         J.dgdvr [ng_3][nvr_x0] = 1.0;
         J.dgdvr [ng_3][nvr_xa] = -k1;
         J.dgdvr [ng_3][nvr_xb] = -k1;
         J.dgdvr [ng_3][nvr_xc] = -k1;

         J.dgdaux[ng_4][na_c     ] = 1.0;
         J.dgdvr [ng_4][nvr_theta] =  s1;

         J.dgdaux[ng_5][na_s     ] = 1.0;
         J.dgdvr [ng_5][nvr_theta] = -c1;
       }
     }
     return;
   }
   return;
}
void x_abc_to_dq_2(Global &G,XbeUsr &X,XbeJac &J) {
   static double k4=0;
   double a,b,c,cost,sint;
   double d,q;
   double alpha,beta;
   const int nvr_a = 0;
   const int nvr_b = 1;
   const int nvr_c = 2;
   const int nvr_cost = 3;
   const int nvr_sint = 4;
   const int nvr_d = 5;
   const int nvr_q = 6;
   const int no_a = 0;
   const int no_b = 1;
   const int no_c = 2;
   const int no_d = 3;
   const int no_q = 4;
   const int na_alpha = 0;
   const int na_beta = 1;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   const int ng_4 = 3;
   if (G.flags[G.i_one_time_parms]) {
     k4 = 0.5*(sqrt(3.0));
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_d] = X.val_vr[nvr_d];
     X.outprm[no_q] = X.val_vr[nvr_q];
     X.outprm[no_a] = X.val_vr[nvr_a];
     X.outprm[no_b] = X.val_vr[nvr_b];
     X.outprm[no_c] = X.val_vr[nvr_c];
     return;
   }
   a = X.val_vr[nvr_a];
   b = X.val_vr[nvr_b];
   c = X.val_vr[nvr_c];
   cost = X.val_vr[nvr_cost];
   sint = X.val_vr[nvr_sint];

   if (G.flags[G.i_init_guess]) {
     alpha = 1.5*a;
     beta  = k4*(b-c);
 
     X.val_aux[na_alpha] = alpha;
     X.val_aux[na_beta ] = beta;

     X.val_vr[nvr_q] = -sint*alpha + cost*beta;
     X.val_vr[nvr_d] =  cost*alpha + sint*beta;

     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       alpha = 1.5*a;
       beta  = k4*(b-c);
       X.val_vr[nvr_q] = -sint*alpha + cost*beta;
       X.val_vr[nvr_d] =  cost*alpha + sint*beta;
     } else if (G.flags[G.i_implicit]) {
       alpha = X.val_aux[na_alpha];
       beta  = X.val_aux[na_beta ];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = alpha - 1.5*a;
         X.g[ng_2] = beta - k4*(b-c);

         X.g[ng_3] = q + sint*alpha - cost*beta;
         X.g[ng_4] = d - cost*alpha - sint*beta;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdaux[ng_1][na_alpha] =  1.0;
         J.dgdvr [ng_1][nvr_a   ] = -1.5;

         J.dgdaux[ng_2][na_beta] =  1.0;
         J.dgdvr [ng_2][nvr_b  ] = -k4;
         J.dgdvr [ng_2][nvr_c  ] =  k4; 

         J.dgdvr [ng_3][nvr_q   ] =  1.0;
         J.dgdvr [ng_3][nvr_sint] =  alpha;
         J.dgdvr [ng_3][nvr_cost] = -beta;
         J.dgdaux[ng_3][na_alpha] =  sint;
         J.dgdaux[ng_3][na_beta ] = -cost;

         J.dgdvr [ng_4][nvr_d   ] =  1.0;
         J.dgdvr [ng_4][nvr_sint] = -beta;
         J.dgdvr [ng_4][nvr_cost] = -alpha;
         J.dgdaux[ng_4][na_alpha] = -cost;
         J.dgdaux[ng_4][na_beta ] = -sint;
       }
     }
     return;
   }
   return;
}
void x_abs(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double y;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = abs(X.val_vr[nvr_x]);
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = abs(X.val_vr[nvr_x]);
     } else if (G.flags[G.i_implicit]) {
       x = X.val_vr[nvr_x];
       y = X.val_vr[nvr_y];
       if (x >= 0) {
         if (G.flags[G.i_function]) {
           X.g[ng_1] = y - x;
         }
         if (G.flags[G.i_jacobian]) {
           J.dgdvr[ng_1][nvr_y] =  1.0;
           J.dgdvr[ng_1][nvr_x] = -1.0;
         }
       } else {
         if (G.flags[G.i_function]) {
           X.g[ng_1] = y + x;
         }
         if (G.flags[G.i_jacobian]) {
           J.dgdvr[ng_1][nvr_y] = 1.0;
           J.dgdvr[ng_1][nvr_x] = 1.0;
         }
       }
     }
     return;
   }
   return;
}
void x_and_2(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   bool X1,X2,Y;
   double x1,x2;
   double y;
   double y_high,hb2;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int nr_y_high = 0;
   const int nr_hb2 = 1;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     y_high = X.rprm[nr_y_high];
     hb2 = 0.5*y_high;
     X.rprm[nr_hb2] = hb2;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   y_high = X.rprm[nr_y_high];
   hb2 = X.rprm[nr_hb2];

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];
   X1 = x1 > hb2;
   X2 = x2 > hb2;
   Y = X1 && X2;

   if (Y) {
     y0 = y_high;
   } else {
     y0 = 0.0;
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_and_3(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   bool X1,X2,X3,Y;
   double x1,x2,x3;
   double y;
   double y_high,hb2;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_x3 = 2;
   const int nvr_y = 3;
   const int nr_y_high = 0;
   const int nr_hb2 = 1;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_x3 = 2;
   const int no_y = 3;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     y_high = X.rprm[nr_y_high];
     hb2 = 0.5*y_high;
     X.rprm[nr_hb2] = hb2;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_x3] = X.val_vr[nvr_x3];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   y_high = X.rprm[nr_y_high];
   hb2 = X.rprm[nr_hb2];

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];
   x3 = X.val_vr[nvr_x3];
   X1 = x1 > hb2;
   X2 = x2 > hb2;
   X3 = x3 > hb2;
   Y = (X1 && X2) && X3;

   if (Y) {
     y0 = y_high;
   } else {
     y0 = 0.0;
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_atan2_rad(Global &G,XbeUsr &X,XbeJac &J) {
   double x,y;
   double theta;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nvr_theta = 2;
   const int no_x = 0;
   const int no_y = 1;
   const int no_theta = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_init_guess]) {
     x = X.val_vr[nvr_x];
     y = X.val_vr[nvr_y];

     X.val_vr[nvr_theta] = atan2(y,x);
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x    ] = X.val_vr[nvr_x    ];
     X.outprm[no_y    ] = X.val_vr[nvr_y    ];
     X.outprm[no_theta] = X.val_vr[nvr_theta];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     x = X.val_vr[nvr_x];
     y = X.val_vr[nvr_y];

     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_theta] = atan2(y,x);
     } else if (G.flags[G.i_implicit]) {
       theta = X.val_vr[nvr_theta];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = theta - atan2(y,x);
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_theta] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_average_mv_1(Global &G,XbeUsr &X,XbeJac &J) {
   int flag_active_edge;
   double time0,delt,delx_sum;
   double clk,x;
   double y;
   int active_pos_edge,active_neg_edge;
   double clk_high,dt,clk_cross,clk_prev,t_lapsed,x_sum,t_prev,x_prev,y0;
   const int nvr_clk = 0;
   const int nvr_x = 1;
   const int nvr_y = 2;
   const int ni_active_pos_edge = 0;
   const int ni_active_neg_edge = 1;
   const int nr_clk_high = 0;
   const int nr_dt = 1;
   const int nr_clk_cross = 2;
   const int nr_clk_prev = 3;
   const int nr_t_lapsed = 4;
   const int nr_x_sum = 5;
   const int nr_t_prev = 6;
   const int nr_x_prev = 7;
   const int nr_y0 = 8;
   const int no_clk = 0;
   const int no_x = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     active_pos_edge = X.iprm[ni_active_pos_edge];
     active_neg_edge = X.iprm[ni_active_neg_edge];

     if (active_pos_edge == 0) {
       if (active_neg_edge == 0) {
         cout << "average_mv_1.xbe: one of pos/neg edge must be 1. Halting..." << endl;
         exit(1);
       }
     } else {
       if (active_neg_edge != 0) {
         cout << "average_mv_1.xbe: only one of pos/neg edge must be 1. Halting..." << endl;
         exit(1);
       }
     }
     clk_high = X.rprm[nr_clk_high];

     if (clk_high < 0.0) {
       cout << "average_mv_1.xbe: check clk_high. Halting..." << endl;
       exit(1);
     }

     clk_cross = 0.5*clk_high;
     X.rprm[nr_clk_cross] = clk_cross;

     X.rprm[nr_clk_prev] = 0.0;
     X.rprm[nr_t_lapsed] = 0.0;
     X.rprm[nr_x_sum   ] = 0.0;
     X.rprm[nr_t_prev  ] = 0.0;
     X.rprm[nr_x_prev  ] = 0.0;
     X.rprm[nr_y0      ] = 0.0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x  ] = X.val_vr[nvr_x  ];
     X.outprm[no_y  ] = X.val_vr[nvr_y  ];
     X.outprm[no_clk] = X.val_vr[nvr_clk];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_next_time]) {
     time0 = G.time_given_x;

     active_pos_edge = X.iprm[ni_active_pos_edge];

     clk_prev  = X.rprm[nr_clk_prev ];
     clk_cross = X.rprm[nr_clk_cross];

     clk = X.val_vr[nvr_clk];

     if (active_pos_edge == 1) {
       if ((clk_prev <= clk_cross) && (clk >= clk_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     } else {
       if ((clk_prev >= clk_cross) && (clk <= clk_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     }
     if (flag_active_edge == 1) {
       G.time_nextbreak_x = time0 + X.rprm[nr_dt];
     } else {
       G.time_nextbreak_x = G.time_end;
     }
     return;
   }
   if (G.flags[G.i_save_history]) {
     time0 = G.time_given_x;
     active_pos_edge = X.iprm[ni_active_pos_edge];

     clk_prev  = X.rprm[nr_clk_prev ];
     clk_cross = X.rprm[nr_clk_cross];

     clk = X.val_vr[nvr_clk];

     if (active_pos_edge == 1) {
       if ((clk_prev <= clk_cross) && (clk >= clk_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     } else {
       if ((clk_prev >= clk_cross) && (clk <= clk_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     }
     t_prev = X.rprm[nr_t_prev];
     x_prev = X.rprm[nr_x_prev];
     delt = time0-t_prev;
     delx_sum = 0.5*delt*(X.val_vr[nvr_x]+x_prev);

     if (flag_active_edge == 1) {
       t_lapsed = X.rprm[nr_t_lapsed];
       x_sum = X.rprm[nr_x_sum] + delx_sum;
       if (t_lapsed > 1.0e-10) {
         X.rprm[nr_y0] = x_sum/t_lapsed;
       }
       X.rprm[nr_t_lapsed] = 0.0;
       X.rprm[nr_x_sum   ] = 0.0;
     } else {
       X.rprm[nr_t_lapsed] += delt;
       X.rprm[nr_x_sum   ] += delx_sum;
     }
     X.rprm[nr_clk_prev] = clk;
     X.rprm[nr_t_prev  ] = time0;
     X.rprm[nr_x_prev  ] = X.val_vr[nvr_x];

     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = 0.0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     y0 = X.rprm[nr_y0];
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
   }
   return;
}
void x_average_mv_2(Global &G,XbeUsr &X,XbeJac &J) {
   double time0,t0,t1,tp;
   int n_t,n_x;
   int n_popped;
   double x;
   double y;
   double T,y0,epsl;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_T = 0;
   const int nr_y0 = 1;
   const int nr_epsl = 2;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     X.que1 = queue<double>();
     X.que2 = queue<double>();
     X.rprm[nr_y0] = 0.0;
     X.rprm[nr_epsl] = X.rprm[nr_T]/1.0e3;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x  ] = X.val_vr[nvr_x  ];
     X.outprm[no_y  ] = X.val_vr[nvr_y  ];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_save_history]) {
     time0 = G.time_given_x;
     x = X.val_vr[nvr_x];
     T = X.rprm[nr_T];
     epsl = X.rprm[nr_epsl];

     n_t = X.que1.size();
     n_x = X.que2.size();

     if (n_t != n_x) {
       cout << "average_mv_2.xbe: n_t: " << n_t << " n_x: " << n_x << endl;
       cout << "  are not equal. Halting..." << endl; exit(1);
     }

     X.que1.push(time0);
     X.que2.push(x);

     t1 = time0 - 1.1*T;

     if (n_t > 2) {
       n_popped = 0;

       while(true) {
         if (X.que1.front() <= t1) {
           X.que1.pop();
           X.que2.pop();
           n_popped++;
         } else {
           break;
         }
       }
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = 0.0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }

   time0 = G.time_given_x;
   T = X.rprm[nr_T];
   epsl = X.rprm[nr_epsl];

   n_t = X.que1.size();
   n_x = X.que2.size();

   if (n_t != n_x) {
     cout << "average_mv_2.xbe: n_t: " << n_t << " n_x: " << n_x << endl;
     cout << "  are not equal. Halting..." << endl; exit(1);
   }
   if (n_t > 2) {
     x = X.val_vr[nvr_x];
     t0 = X.que1.front();
     tp = t0 + T;
     if (time0 >= tp) {
       X.rprm[nr_y0] = calc_avg_1(X.que1,X.que2,time0,x,T);
     }
   }
   if (G.flags[G.i_trns]) {
     y0 = X.rprm[nr_y0];
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
   }
   return;
}
void x_clock(Global &G,XbeUsr &X,XbeJac &J) {
   double y0,t0_new,delta_min,del1,del2,t_a,t_b,tnext_p;
   int n;
   double y;
   double T1,T2,L1,L2,t0,delta1,delta2,T,L0,tk1,tk2,tk3,tk4,tk5,slope1,slope2,
     epsl;
   const int nvr_y = 0;
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
   const int no_y = 0;
   const int ng_1 = 0;
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
       cout << "clock.xbe: T1 is too small. Check delta1, delta2. Halting..." << endl; exit(1);
     }
     if (del2 < delta_min) {
       cout << "clock.xbe: T2 is too small. Check delta1, delta2. Halting..." << endl; exit(1);
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
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   t0   = X.rprm[nr_t0  ];
   T    = X.rprm[nr_T   ];
   epsl = X.rprm[nr_epsl];
   tk1  = X.rprm[nr_tk1 ];
   tk2  = X.rprm[nr_tk2 ];
   tk3  = X.rprm[nr_tk3 ];
   tk4  = X.rprm[nr_tk4 ];
   tk5  = X.rprm[nr_tk5 ];

   if (G.time_given_x < t0) {
     n = ((t0-G.time_given_x)/T) + 1;
     t0_new = t0-n*T;
   } else {
     t0_new = t0;
   }
   t_a = G.time_given_x-t0_new;
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
     G.time_nextbreak_x = G.time_given_x + (tnext_p-t_b);
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

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_clock_1(Global &G,XbeUsr &X,XbeJac &J) {
   double y0,t0_new,delta_min,del1,del2,t_a,t_b,tnext_p;
   double L1,L2;
   int n;
   double y;
   double f_hz,D,y_high,delta1,delta2,t0,T1,T2,T,L0,tk1,tk2,tk3,tk4,tk5,slope1,
     slope2,epsl;
   const int nvr_y = 0;
   const int nr_f_hz = 0;
   const int nr_D = 1;
   const int nr_y_high = 2;
   const int nr_delta1 = 3;
   const int nr_delta2 = 4;
   const int nr_t0 = 5;
   const int nr_T1 = 6;
   const int nr_T2 = 7;
   const int nr_T = 8;
   const int nr_L0 = 9;
   const int nr_tk1 = 10;
   const int nr_tk2 = 11;
   const int nr_tk3 = 12;
   const int nr_tk4 = 13;
   const int nr_tk5 = 14;
   const int nr_slope1 = 15;
   const int nr_slope2 = 16;
   const int nr_epsl = 17;
   const int no_y = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     f_hz = X.rprm[nr_f_hz];
     D    = X.rprm[nr_D   ];

     if (D >= 1.0) {
       cout << "clock_1.xbe: D.ge.1.0 ? Halting..." << endl;
       exit(1);
     }
     if (D <= 0.0) {
       cout << "clock_1.xbe: D.le.0.0 ? Halting..." << endl;
       exit(1);
     }
     T = 1.0/f_hz;
     T1 = D*T;
     T2 = T - T1;

     X.rprm[nr_T1] = T1;
     X.rprm[nr_T2] = T2;
     X.rprm[nr_T ] = T;

     X.rprm[nr_T1] = T1;
     X.rprm[nr_T2] = T2;
     X.rprm[nr_T] = T;

     L1 = X.rprm[nr_y_high];
     L2 = 0.0;

     delta1 = X.rprm[nr_delta1];
     delta2 = X.rprm[nr_delta2];

     delta_min = 0.1*min(delta1,delta2);
     del1 = T1-0.5*(delta1+delta2);
     del2 = T2-0.5*(delta1+delta2);

     if (del1 < delta_min) {
       cout << "clock_1.xbe: T1 is too small. Check delta1, delta2. Halting..." << endl; exit(1);
     }
     if (del2 < delta_min) {
       cout << "clock_1.xbe: T2 is too small. Check delta1, delta2. Halting..." << endl; exit(1);
     }
     tk1 = 0.5*delta1;
     tk2 = T1 - 0.5*delta2;
     tk3 = T1 + 0.5*delta2;
     tk4 = T  - 0.5*delta1;
     tk5 = T  + 0.5*delta1;

     slope1 = (L1-L2)/delta1;
     slope2 = (L2-L1)/delta2;
     epsl = min(delta1,delta2)/10.0;
     L0 = 0.5*(L1+L2);

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
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   t0   = X.rprm[nr_t0  ];
   T    = X.rprm[nr_T   ];
   epsl = X.rprm[nr_epsl];
   tk1  = X.rprm[nr_tk1 ];
   tk2  = X.rprm[nr_tk2 ];
   tk3  = X.rprm[nr_tk3 ];
   tk4  = X.rprm[nr_tk4 ];
   tk5  = X.rprm[nr_tk5 ];

   if (G.time_given_x < t0) {
     n = ((t0-G.time_given_x)/T) + 1;
     t0_new = t0-n*T;
   } else {
     t0_new = t0;
   }
   t_a = G.time_given_x-t0_new;
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
     G.time_nextbreak_x = G.time_given_x + (tnext_p-t_b);
     return;
   }

   L0 = X.rprm[nr_L0];
   L1 = X.rprm[nr_y_high];
   L2 = 0.0;
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

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_clock_1a(Global &G,XbeUsr &X,XbeJac &J) {
   double y0,t0_new,delta_min,del1,del2,t_a,t_b,tnext_p;
   double L1,L2;
   int n;
   double y;
   double f_hz,D,y_low,y_high,delta1,delta2,t0,T1,T2,T,L0,tk1,tk2,tk3,tk4,tk5,
     slope1,slope2,epsl;
   const int nvr_y = 0;
   const int nr_f_hz = 0;
   const int nr_D = 1;
   const int nr_y_low = 2;
   const int nr_y_high = 3;
   const int nr_delta1 = 4;
   const int nr_delta2 = 5;
   const int nr_t0 = 6;
   const int nr_T1 = 7;
   const int nr_T2 = 8;
   const int nr_T = 9;
   const int nr_L0 = 10;
   const int nr_tk1 = 11;
   const int nr_tk2 = 12;
   const int nr_tk3 = 13;
   const int nr_tk4 = 14;
   const int nr_tk5 = 15;
   const int nr_slope1 = 16;
   const int nr_slope2 = 17;
   const int nr_epsl = 18;
   const int no_y = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     f_hz = X.rprm[nr_f_hz];
     D    = X.rprm[nr_D   ];

     if (D >= 1.0) {
       cout << "clock_1a.xbe: D.ge.1.0 ? Halting..." << endl;
       exit(1);
     }
     if (D <= 0.0) {
       cout << "clock_1a.xbe: D.le.0.0 ? Halting..." << endl;
       exit(1);
     }
     T = 1.0/f_hz;
     T1 = D*T;
     T2 = T - T1;

     X.rprm[nr_T1] = T1;
     X.rprm[nr_T2] = T2;
     X.rprm[nr_T ] = T;

     X.rprm[nr_T1] = T1;
     X.rprm[nr_T2] = T2;
     X.rprm[nr_T] = T;

     L1 = X.rprm[nr_y_high];
     L2 = X.rprm[nr_y_low ];

     delta1 = X.rprm[nr_delta1];
     delta2 = X.rprm[nr_delta2];

     delta_min = 0.1*min(delta1,delta2);
     del1 = T1-0.5*(delta1+delta2);
     del2 = T2-0.5*(delta1+delta2);

     if (del1 < delta_min) {
       cout << "clock_1a.xbe: T1 is too small. Check delta1, delta2. Halting..." << endl; exit(1);
     }
     if (del2 < delta_min) {
       cout << "clock_1a.xbe: T2 is too small. Check delta1, delta2. Halting..." << endl; exit(1);
     }
     tk1 = 0.5*delta1;
     tk2 = T1 - 0.5*delta2;
     tk3 = T1 + 0.5*delta2;
     tk4 = T  - 0.5*delta1;
     tk5 = T  + 0.5*delta1;

     slope1 = (L1-L2)/delta1;
     slope2 = (L2-L1)/delta2;
     epsl = min(delta1,delta2)/10.0;
     L0 = 0.5*(L1+L2);

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
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   t0   = X.rprm[nr_t0  ];
   T    = X.rprm[nr_T   ];
   epsl = X.rprm[nr_epsl];
   tk1  = X.rprm[nr_tk1 ];
   tk2  = X.rprm[nr_tk2 ];
   tk3  = X.rprm[nr_tk3 ];
   tk4  = X.rprm[nr_tk4 ];
   tk5  = X.rprm[nr_tk5 ];

   if (G.time_given_x < t0) {
     n = ((t0-G.time_given_x)/T) + 1;
     t0_new = t0-n*T;
   } else {
     t0_new = t0;
   }
   t_a = G.time_given_x-t0_new;
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
     G.time_nextbreak_x = G.time_given_x + (tnext_p-t_b);
     return;
   }

   L0 = X.rprm[nr_L0];
   L1 = X.rprm[nr_y_high];
   L2 = X.rprm[nr_y_low];
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

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_clock_3(Global &G,XbeUsr &X,XbeJac &J) {
   double y0,t0_new,delta_min,del1,del2,t_a,t_b,tnext_p;
   double t10,t20,t30,t40;
   int n;
   double y;
   int indx;
   double alpha,beta,T,y_high,delta1,delta2,T1,T2,t0,L1,L2,L0,tk1,tk2,tk3,tk4,
     tk5,slope1,slope2,epsl;
   const int nvr_y = 0;
   const int ni_indx = 0;
   const int nr_alpha = 0;
   const int nr_beta = 1;
   const int nr_T = 2;
   const int nr_y_high = 3;
   const int nr_delta1 = 4;
   const int nr_delta2 = 5;
   const int nr_T1 = 6;
   const int nr_T2 = 7;
   const int nr_t0 = 8;
   const int nr_L1 = 9;
   const int nr_L2 = 10;
   const int nr_L0 = 11;
   const int nr_tk1 = 12;
   const int nr_tk2 = 13;
   const int nr_tk3 = 14;
   const int nr_tk4 = 15;
   const int nr_tk5 = 16;
   const int nr_slope1 = 17;
   const int nr_slope2 = 18;
   const int nr_epsl = 19;
   const int no_y = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     y_high = X.rprm[nr_y_high];
     if (y_high < 1.0) {
       cout << "clock_3: y_high < 1.0? Halting..." << endl; exit(1);
     }

     T     = X.rprm[nr_T    ];
     alpha = X.rprm[nr_alpha];
     beta  = X.rprm[nr_beta ];

     delta1 = X.rprm[nr_delta1];
     delta2 = X.rprm[nr_delta2];

     t10 = ((180.0-2*alpha)/360.0)*T;
     t20 = T-t10;
     t30 = (alpha/360.0)*T;
     t40 = (beta/360.0)*T;
     epsl = min(delta1,delta2)/100.0;

     indx = X.iprm[ni_indx];

     if (indx == 1) {
       T1 = t10;
       T2 = t20;
       t0 = t30+t40;
       L1 = y_high; L2 = 0.0;
     } else if (indx == 2) {
       T1 = t10;
       T2 = t20;
       t0 = t30+t40;
       L1 = 0.0; L2 = y_high;
     } else if (indx == 3) {
       T1 = t10;
       T2 = t20;
       t0 = t30+t40+(0.5*T);
       L1 = y_high; L2 = 0.0;
     } else if (indx == 4) {
       T1 = t10;
       T2 = t20;
       t0 = t30+t40+(0.5*T);
       L1 = 0.0; L2 = y_high;
     } else {
       cout << "clock_3: wrong indx value: " << indx << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     delta_min = 0.1*min(delta1,delta2);
     del1 = T1-0.5*(delta1+delta2);
     del2 = T2-0.5*(delta1+delta2);

     if (del1 < delta_min) {
       cout << "clock.xbe: T1 is too small. Check delta1, delta2. Halting..." << endl; exit(1);
     }
     if (del2 < delta_min) {
       cout << "clock.xbe: T2 is too small. Check delta1, delta2. Halting..." << endl; exit(1);
     }
     tk1 = 0.5*delta1;
     tk2 = T1 - 0.5*delta2;
     tk3 = T1 + 0.5*delta2;
     tk4 = T  - 0.5*delta1;
     tk5 = T  + 0.5*delta1;

     slope1 = (L1-L2)/delta1;
     slope2 = (L2-L1)/delta2;
     epsl = min(delta1,delta2)/10.0;
     L0 = 0.5*(L1+L2);

     X.rprm[nr_T1] = T1;
     X.rprm[nr_T2] = T2;
     X.rprm[nr_L1] = L1;
     X.rprm[nr_L2] = L2;
     X.rprm[nr_t0] = t0;
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
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   t0   = X.rprm[nr_t0  ];
   T    = X.rprm[nr_T   ];
   epsl = X.rprm[nr_epsl];
   tk1  = X.rprm[nr_tk1 ];
   tk2  = X.rprm[nr_tk2 ];
   tk3  = X.rprm[nr_tk3 ];
   tk4  = X.rprm[nr_tk4 ];
   tk5  = X.rprm[nr_tk5 ];

   if (G.time_given_x < t0) {
     n = ((t0-G.time_given_x)/T) + 1;
     t0_new = t0-n*T;
   } else {
     t0_new = t0;
   }
   t_a = G.time_given_x-t0_new;
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
     G.time_nextbreak_x = G.time_given_x + (tnext_p-t_b);
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

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_clock_3ph(Global &G,XbeUsr &X,XbeJac &J) {
   double y0,t0_new,dt_min,del1,del2,t_a,t_b,tnext_p;
   int n;
   double y;
   int index1,flag_frequency,flag_period;
   double x_low,x_high,frequency,T,D,alpha,dt,T1,T2,t0,dt1,dt2,L0,L1,L2,tk1,tk2,
     tk3,tk4,tk5,slope1,slope2,epsl;
   const int nvr_y = 0;
   const int ni_index1 = 0;
   const int ni_flag_frequency = 1;
   const int ni_flag_period = 2;
   const int nr_x_low = 0;
   const int nr_x_high = 1;
   const int nr_frequency = 2;
   const int nr_T = 3;
   const int nr_D = 4;
   const int nr_alpha = 5;
   const int nr_dt = 6;
   const int nr_T1 = 7;
   const int nr_T2 = 8;
   const int nr_t0 = 9;
   const int nr_dt1 = 10;
   const int nr_dt2 = 11;
   const int nr_L0 = 12;
   const int nr_L1 = 13;
   const int nr_L2 = 14;
   const int nr_tk1 = 15;
   const int nr_tk2 = 16;
   const int nr_tk3 = 17;
   const int nr_tk4 = 18;
   const int nr_tk5 = 19;
   const int nr_slope1 = 20;
   const int nr_slope2 = 21;
   const int nr_epsl = 22;
   const int no_y = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     flag_frequency = X.iprm[ni_flag_frequency];
     flag_period = X.iprm[ni_flag_period];

     if ((flag_frequency == 0) && (flag_period == 0)) {
       cout << "clock_3ph.xbe: check flag_frequency and flag_period" << endl;
       cout << "  Both cannot be zero." << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     if ((flag_frequency != 0) && (flag_period != 0)) {
       cout << "clock_3ph.xbe: check flag_frequency and flag_period" << endl;
       cout << "  Both cannot be non-zero." << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     index1 = X.iprm[ni_index1];
     if ((index1 < 1) || (index1 > 6)) {
       cout << "clock_3ph.xbe: index1 must be in the range 1 to 6" << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     if (flag_frequency != 0) {
       frequency = X.rprm[nr_frequency];
       T = 1.0/frequency;
       X.rprm[nr_T] = T;
     }
     if (flag_period != 0) {
       T = X.rprm[nr_T];
       frequency = 1.0/T;
       X.rprm[nr_frequency] = frequency;
     }
     x_low  = X.rprm[nr_x_low ];
     x_high = X.rprm[nr_x_high];
     T      = X.rprm[nr_T     ];
     dt     = X.rprm[nr_dt    ];
     D      = X.rprm[nr_D     ];
     alpha  = X.rprm[nr_alpha ];

     T1 = D*T;
     T2 = T - T1;
     t0 = (alpha/360.0)*T + ((double)(index1-1))*(T/6.0);
     dt1 = dt;
     dt2 = dt;

     X.rprm[nr_T1 ] = T1;
     X.rprm[nr_T2 ] = T2;
     X.rprm[nr_t0 ] = t0;
     X.rprm[nr_dt1] = dt;
     X.rprm[nr_dt2] = dt;

     L1 = x_high;
     L2 = x_low;

     slope1 = (L1-L2)/dt1;
     slope2 = (L2-L1)/dt2;
     epsl = min(dt1,dt2)/10.0;
     L0 = 0.5*(L1+L2);

     X.rprm[nr_slope1] = slope1;
     X.rprm[nr_slope2] = slope2;
     X.rprm[nr_L1] = L1;
     X.rprm[nr_L2] = L2;
     X.rprm[nr_L0] = L0;
     X.rprm[nr_epsl] = epsl;

     dt_min = 0.1*min(dt1,dt2);
     del1 = T1-0.5*(dt1+dt2);
     del2 = T2-0.5*(dt1+dt2);

     if (del1 < dt_min) {
       cout << "clock_3ph.xbe: T1 is too small. Check dt1, dt2. Halting..." << endl; exit(1);
     }
     if (del2 < dt_min) {
       cout << "clock_3ph.xbe: T2 is too small. Check dt1, dt2. Halting..." << endl; exit(1);
     }
     tk1 = 0.5*dt1;
     tk2 = T1 - 0.5*dt2;
     tk3 = T1 + 0.5*dt2;
     tk4 = T  - 0.5*dt1;
     tk5 = T  + 0.5*dt1;

     X.rprm[nr_tk1] = tk1;
     X.rprm[nr_tk2] = tk2;
     X.rprm[nr_tk3] = tk3;
     X.rprm[nr_tk4] = tk4;
     X.rprm[nr_tk5] = tk5;

     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   t0   = X.rprm[nr_t0  ];
   T    = X.rprm[nr_T   ];
   epsl = X.rprm[nr_epsl];
   tk1  = X.rprm[nr_tk1 ];
   tk2  = X.rprm[nr_tk2 ];
   tk3  = X.rprm[nr_tk3 ];
   tk4  = X.rprm[nr_tk4 ];
   tk5  = X.rprm[nr_tk5 ];

   if (G.time_given_x < t0) {
     n = ((t0-G.time_given_x)/T) + 1;
     t0_new = t0-n*T;
   } else {
     t0_new = t0;
   }
   t_a = G.time_given_x-t0_new;
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
     G.time_nextbreak_x = G.time_given_x + (tnext_p-t_b);
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

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_clock_thyr(Global &G,XbeUsr &X,XbeJac &J) {
   double y0,t0_new,dt_min,del1,del2,t_a,t_b,tnext_p;
   int n;
   double y;
   int flag_frequency,flag_period,flag_tw_degrees;
   double x_low,x_high,frequency,T,tw_deg,tw,alpha,beta,dt,T1,T2,t0,dt1,dt2,L0,
     L1,L2,tk1,tk2,tk3,tk4,tk5,slope1,slope2,epsl;
   const int nvr_y = 0;
   const int ni_flag_frequency = 0;
   const int ni_flag_period = 1;
   const int ni_flag_tw_degrees = 2;
   const int nr_x_low = 0;
   const int nr_x_high = 1;
   const int nr_frequency = 2;
   const int nr_T = 3;
   const int nr_tw_deg = 4;
   const int nr_tw = 5;
   const int nr_alpha = 6;
   const int nr_beta = 7;
   const int nr_dt = 8;
   const int nr_T1 = 9;
   const int nr_T2 = 10;
   const int nr_t0 = 11;
   const int nr_dt1 = 12;
   const int nr_dt2 = 13;
   const int nr_L0 = 14;
   const int nr_L1 = 15;
   const int nr_L2 = 16;
   const int nr_tk1 = 17;
   const int nr_tk2 = 18;
   const int nr_tk3 = 19;
   const int nr_tk4 = 20;
   const int nr_tk5 = 21;
   const int nr_slope1 = 22;
   const int nr_slope2 = 23;
   const int nr_epsl = 24;
   const int no_y = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     flag_frequency = X.iprm[ni_flag_frequency];
     flag_period = X.iprm[ni_flag_period];

     if ((flag_frequency == 0) && (flag_period == 0)) {
       cout << "clock_thyr.xbe: check flag_frequency and flag_period" << endl;
       cout << "  Both cannot be zero." << endl;
       cout << "  Halting..." << endl;
       exit(1);
     }
     if ((flag_frequency != 0) && (flag_period != 0)) {
       cout << "clock_thyr.xbe: check flag_frequency and flag_period" << endl;
       cout << "  Both cannot be non-zero." << endl;
       cout << "  Halting..." << endl;
       exit(1);
     }
     if (flag_frequency != 0) {
       frequency = X.rprm[nr_frequency];
       T = 1.0/frequency;
       X.rprm[nr_T] = T;
     }
     if (flag_period != 0) {
       T = X.rprm[nr_T];
       frequency = 1.0/T;
       X.rprm[nr_frequency] = frequency;
     }
     flag_tw_degrees = X.iprm[ni_flag_tw_degrees];
     if (flag_tw_degrees != 0) {
       tw_deg = X.rprm[nr_tw_deg];
       tw = T*(tw_deg/360.0);
       X.rprm[nr_tw] = tw;
     }

     x_low  = X.rprm[nr_x_low ];
     x_high = X.rprm[nr_x_high];
     T      = X.rprm[nr_T     ];
     alpha  = X.rprm[nr_alpha ];
     beta   = X.rprm[nr_beta  ];
     tw     = X.rprm[nr_tw    ];
     dt     = X.rprm[nr_dt    ];

     T1 = tw;
     T2 = T - T1;
     t0 = ((alpha+beta)/360.0)*T;
     dt1 = dt;
     dt2 = dt;

     X.rprm[nr_T1 ] = T1;
     X.rprm[nr_T2 ] = T2;
     X.rprm[nr_t0 ] = t0;
     X.rprm[nr_dt1] = dt;
     X.rprm[nr_dt2] = dt;

     L1 = x_high;
     L2 = x_low;

     slope1 = (L1-L2)/dt1;
     slope2 = (L2-L1)/dt2;
     epsl = min(dt1,dt2)/10.0;
     L0 = 0.5*(L1+L2);

     X.rprm[nr_slope1] = slope1;
     X.rprm[nr_slope2] = slope2;
     X.rprm[nr_L1] = L1;
     X.rprm[nr_L2] = L2;
     X.rprm[nr_L0] = L0;
     X.rprm[nr_epsl] = epsl;

     dt_min = 0.1*min(dt1,dt2);
     del1 = T1-0.5*(dt1+dt2);
     del2 = T2-0.5*(dt1+dt2);

     if (del1 < dt_min) {
       cout << "clock_thyr.xbe: T1 is too small. Check dt1, dt2." << endl; exit(1);
       cout << "  T1: " << T1 << ". Halting..." << endl;
     }
     if (del2 < dt_min) {
       cout << "clock_thyr.xbe: T2 is too small. Check dt1, dt2." << endl; exit(1);
       cout << "  T2: " << T2 << ". Halting..." << endl;
     }
     tk1 = 0.5*dt1;
     tk2 = T1 - 0.5*dt2;
     tk3 = T1 + 0.5*dt2;
     tk4 = T  - 0.5*dt1;
     tk5 = T  + 0.5*dt1;

     X.rprm[nr_tk1] = tk1;
     X.rprm[nr_tk2] = tk2;
     X.rprm[nr_tk3] = tk3;
     X.rprm[nr_tk4] = tk4;
     X.rprm[nr_tk5] = tk5;

     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   t0   = X.rprm[nr_t0  ];
   T    = X.rprm[nr_T   ];
   epsl = X.rprm[nr_epsl];
   tk1  = X.rprm[nr_tk1 ];
   tk2  = X.rprm[nr_tk2 ];
   tk3  = X.rprm[nr_tk3 ];
   tk4  = X.rprm[nr_tk4 ];
   tk5  = X.rprm[nr_tk5 ];

   if (G.time_given_x < t0) {
     n = ((t0-G.time_given_x)/T) + 1;
     t0_new = t0-n*T;
   } else {
     t0_new = t0;
   }
   t_a = G.time_given_x-t0_new;
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
     G.time_nextbreak_x = G.time_given_x + (tnext_p-t_b);
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

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_cmpr_1_1(Global &G,XbeUsr &X,XbeJac &J) {
   int iter;
   double y0,time0;
   double x;
   double y;
   int flag_invert,flag_quad;
   double x0,y_low,y_high,x_1,x_2,t_1,t_2,epsl,delt_min,delt_nrml;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int ni_flag_invert = 0;
   const int ni_flag_quad = 1;
   const int nr_x0 = 0;
   const int nr_y_low = 1;
   const int nr_y_high = 2;
   const int nr_x_1 = 3;
   const int nr_x_2 = 4;
   const int nr_t_1 = 5;
   const int nr_t_2 = 6;
   const int nr_epsl = 7;
   const int nr_delt_min = 8;
   const int nr_delt_nrml = 9;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   x0 = X.rprm[nr_x0];

   flag_invert = X.iprm[ni_flag_invert];
   flag_quad = X.iprm[ni_flag_quad];

   x = X.val_vr[nvr_x];

   if (G.flags[G.i_save_history]) {

     x_1 = X.rprm[nr_x_1];
     t_1 = X.rprm[nr_t_1];

     x_2 = x_1;
     t_2 = t_1;

     x_1 = x;
     t_1 = G.time_given_x;

     X.rprm[nr_x_1] = x_1;
     X.rprm[nr_x_2] = x_2;
     X.rprm[nr_t_1] = t_1;
     X.rprm[nr_t_2] = t_2;
     return;
   }
   if (G.flags[G.i_next_time]) {
     iter = G.iter_trns_x;
     time0 = G.time_given_x;

     x_1 = X.rprm[nr_x_1];
     x_2 = X.rprm[nr_x_2];
     t_1 = X.rprm[nr_t_1];
     t_2 = X.rprm[nr_t_2];
     epsl      = X.rprm[nr_epsl     ];
     delt_min  = X.rprm[nr_delt_min ];
     delt_nrml = X.rprm[nr_delt_nrml];

     G.time_nextbreak_x = get_tnext(iter,flag_quad,
       t_1,t_2,(x_1-x0),(x_2-x0),
       time0,(x-x0),epsl,delt_min,delt_nrml);

     return;
   }

   if (flag_invert == 0) {
     if (x >= x0) {
       y0 = X.rprm[nr_y_high];
     } else {
       y0 = X.rprm[nr_y_low];
     }   
   } else {
     if (x >= x0) {
       y0 = X.rprm[nr_y_low];
     } else {
       y0 = X.rprm[nr_y_high];
     }   
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_cmpr_1_2(Global &G,XbeUsr &X,XbeJac &J) {
   int iter;
   double y10,y20,time0;
   double x;
   double y1,y2;
   int flag_invert,flag_quad;
   double x0,y_low,y_high,x_1,x_2,t_1,t_2,epsl,delt_min,delt_nrml;
   const int nvr_x = 0;
   const int nvr_y1 = 1;
   const int nvr_y2 = 2;
   const int ni_flag_invert = 0;
   const int ni_flag_quad = 1;
   const int nr_x0 = 0;
   const int nr_y_low = 1;
   const int nr_y_high = 2;
   const int nr_x_1 = 3;
   const int nr_x_2 = 4;
   const int nr_t_1 = 5;
   const int nr_t_2 = 6;
   const int nr_epsl = 7;
   const int nr_delt_min = 8;
   const int nr_delt_nrml = 9;
   const int no_x = 0;
   const int no_y1 = 1;
   const int no_y2 = 2;
   const int ng_1 = 0;
   const int ng_2 = 1;
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y1] = X.val_vr[nvr_y1];
     X.outprm[no_y2] = X.val_vr[nvr_y2];
     return;
   }

   x0 = X.rprm[nr_x0];

   flag_invert = X.iprm[ni_flag_invert];
   flag_quad = X.iprm[ni_flag_quad];

   x = X.val_vr[nvr_x];

   if (G.flags[G.i_save_history]) {

     x_1 = X.rprm[nr_x_1];
     t_1 = X.rprm[nr_t_1];

     x_2 = x_1;
     t_2 = t_1;

     x_1 = x;
     t_1 = G.time_given_x;

     X.rprm[nr_x_1] = x_1;
     X.rprm[nr_x_2] = x_2;
     X.rprm[nr_t_1] = t_1;
     X.rprm[nr_t_2] = t_2;
     return;
   }
   if (G.flags[G.i_next_time]) {
     iter = G.iter_trns_x;
     time0 = G.time_given_x;

     x_1 = X.rprm[nr_x_1];
     x_2 = X.rprm[nr_x_2];
     t_1 = X.rprm[nr_t_1];
     t_2 = X.rprm[nr_t_2];
     epsl      = X.rprm[nr_epsl     ];
     delt_min  = X.rprm[nr_delt_min ];
     delt_nrml = X.rprm[nr_delt_nrml];

     G.time_nextbreak_x = get_tnext(iter,flag_quad,
       t_1,t_2,(x_1-x0),(x_2-x0),
       time0,(x-x0),epsl,delt_min,delt_nrml);

     return;
   }

   if (flag_invert == 0) {
     if (x >= x0) {
       y10 = X.rprm[nr_y_high];
       y20 = X.rprm[nr_y_low];
     } else {
       y10 = X.rprm[nr_y_low];
       y20 = X.rprm[nr_y_high];
     }   
   } else {
     if (x >= x0) {
       y10 = X.rprm[nr_y_low];
       y20 = X.rprm[nr_y_high];
     } else {
       y10 = X.rprm[nr_y_high];
       y20 = X.rprm[nr_y_low];
     }   
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y1] = y10;
     X.val_vr[nvr_y2] = y20;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y1] = y10;
       X.val_vr[nvr_y2] = y20;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y1] - y10;
         X.g[ng_2] = X.val_vr[nvr_y2] - y20;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y1] = 1.0;
         J.dgdvr[ng_2][nvr_y2] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_cmpr_2_1(Global &G,XbeUsr &X,XbeJac &J) {
   int iter;
   double y0,time0;
   double x1,x2;
   double y;
   int flag_invert,flag_quad;
   double y_low,y_high,x1_1,x1_2,x2_1,x2_2,t_1,t_2,epsl,delt_min,delt_nrml;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int ni_flag_invert = 0;
   const int ni_flag_quad = 1;
   const int nr_y_low = 0;
   const int nr_y_high = 1;
   const int nr_x1_1 = 2;
   const int nr_x1_2 = 3;
   const int nr_x2_1 = 4;
   const int nr_x2_2 = 5;
   const int nr_t_1 = 6;
   const int nr_t_2 = 7;
   const int nr_epsl = 8;
   const int nr_delt_min = 9;
   const int nr_delt_nrml = 10;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   flag_invert = X.iprm[ni_flag_invert];
   flag_quad = X.iprm[ni_flag_quad];

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];

   if (G.flags[G.i_save_history]) {

     x1_1 = X.rprm[nr_x1_1];
     x2_1 = X.rprm[nr_x2_1];
     t_1 = X.rprm[nr_t_1];

     x1_2 = x1_1;
     x2_2 = x2_1;
     t_2 = t_1;

     x1_1 = x1;
     x2_1 = x2;
     t_1 = G.time_given_x;

     X.rprm[nr_x1_1] = x1_1;
     X.rprm[nr_x1_2] = x1_2;
     X.rprm[nr_x2_1] = x2_1;
     X.rprm[nr_x2_2] = x2_2;
     X.rprm[nr_t_1] = t_1;
     X.rprm[nr_t_2] = t_2;
     return;
   }
   if (G.flags[G.i_next_time]) {
     iter = G.iter_trns_x;
     time0 = G.time_given_x;

     x1_1 = X.rprm[nr_x1_1];
     x1_2 = X.rprm[nr_x1_2];
     x2_1 = X.rprm[nr_x2_1];
     x2_2 = X.rprm[nr_x2_2];
     t_1  = X.rprm[nr_t_1 ];
     t_2  = X.rprm[nr_t_2 ];
     epsl      = X.rprm[nr_epsl     ];
     delt_min  = X.rprm[nr_delt_min ];
     delt_nrml = X.rprm[nr_delt_nrml];

     G.time_nextbreak_x = get_tnext(iter,flag_quad,
       t_1,t_2,(x1_1-x2_1),(x1_2-x2_2),
       time0,(x1-x2),epsl,delt_min,delt_nrml);

     return;
   }

   if (flag_invert == 0) {
     if (x1 >= x2) {
       y0 = X.rprm[nr_y_high];
     } else {
       y0 = X.rprm[nr_y_low];
     }   
   } else {
     if (x1 >= x2) {
       y0 = X.rprm[nr_y_low];
     } else {
       y0 = X.rprm[nr_y_high];
     }   
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_cmpr_2_2(Global &G,XbeUsr &X,XbeJac &J) {
   int iter;
   double y10,y20,time0;
   double x1,x2;
   double y1,y2;
   int flag_invert,flag_quad;
   double y_low,y_high,x1_1,x1_2,x2_1,x2_2,t_1,t_2,epsl,delt_min,delt_nrml;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y1 = 2;
   const int nvr_y2 = 3;
   const int ni_flag_invert = 0;
   const int ni_flag_quad = 1;
   const int nr_y_low = 0;
   const int nr_y_high = 1;
   const int nr_x1_1 = 2;
   const int nr_x1_2 = 3;
   const int nr_x2_1 = 4;
   const int nr_x2_2 = 5;
   const int nr_t_1 = 6;
   const int nr_t_2 = 7;
   const int nr_epsl = 8;
   const int nr_delt_min = 9;
   const int nr_delt_nrml = 10;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y1 = 2;
   const int no_y2 = 3;
   const int ng_1 = 0;
   const int ng_2 = 1;
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y1] = X.val_vr[nvr_y1];
     X.outprm[no_y2] = X.val_vr[nvr_y2];
     return;
   }

   flag_invert = X.iprm[ni_flag_invert];
   flag_quad = X.iprm[ni_flag_quad];

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];

   if (G.flags[G.i_save_history]) {

     x1_1 = X.rprm[nr_x1_1];
     x2_1 = X.rprm[nr_x2_1];
     t_1 = X.rprm[nr_t_1];

     x1_2 = x1_1;
     x2_2 = x2_1;
     t_2 = t_1;

     x1_1 = x1;
     x2_1 = x2;
     t_1 = G.time_given_x;

     X.rprm[nr_x1_1] = x1_1;
     X.rprm[nr_x1_2] = x1_2;
     X.rprm[nr_x2_1] = x2_1;
     X.rprm[nr_x2_2] = x2_2;
     X.rprm[nr_t_1] = t_1;
     X.rprm[nr_t_2] = t_2;
     return;
   }
   if (G.flags[G.i_next_time]) {
     iter = G.iter_trns_x;
     time0 = G.time_given_x;

     x1_1 = X.rprm[nr_x1_1];
     x1_2 = X.rprm[nr_x1_2];
     x2_1 = X.rprm[nr_x2_1];
     x2_2 = X.rprm[nr_x2_2];
     t_1  = X.rprm[nr_t_1 ];
     t_2  = X.rprm[nr_t_2 ];
     epsl      = X.rprm[nr_epsl     ];
     delt_min  = X.rprm[nr_delt_min ];
     delt_nrml = X.rprm[nr_delt_nrml];

     G.time_nextbreak_x = get_tnext(iter,flag_quad,
       t_1,t_2,(x1_1-x2_1),(x1_2-x2_2),
       time0,(x1-x2),epsl,delt_min,delt_nrml);

     return;
   }

   if (flag_invert == 0) {
     if (x1 >= x2) {
       y10 = X.rprm[nr_y_high];
       y20 = X.rprm[nr_y_low];
     } else {
       y10 = X.rprm[nr_y_low];
       y20 = X.rprm[nr_y_high];
     }   
   } else {
     if (x1 >= x2) {
       y10 = X.rprm[nr_y_low];
       y20 = X.rprm[nr_y_high];
     } else {
       y10 = X.rprm[nr_y_high];
       y20 = X.rprm[nr_y_low];
     }   
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y1] = y10;
     X.val_vr[nvr_y2] = y20;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y1] = y10;
       X.val_vr[nvr_y2] = y20;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y1] - y10;
         X.g[ng_2] = X.val_vr[nvr_y2] - y20;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y1] = 1.0;
         J.dgdvr[ng_2][nvr_y2] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_cmpr_simple_2_1(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   double x1,x2;
   double y;
   int flag_invert;
   double y_low,y_high;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int ni_flag_invert = 0;
   const int nr_y_low = 0;
   const int nr_y_high = 1;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   flag_invert = X.iprm[ni_flag_invert];

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];

   if (flag_invert == 0) {
     if (x1 >= x2) {
       y0 = X.rprm[nr_y_high];
     } else {
       y0 = X.rprm[nr_y_low];
     }   
   } else {
     if (x1 >= x2) {
       y0 = X.rprm[nr_y_low];
     } else {
       y0 = X.rprm[nr_y_high];
     }   
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }

   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_cmpr_simple_2_2(Global &G,XbeUsr &X,XbeJac &J) {
   double y10,y20;
   double x1,x2;
   double y1,y2;
   int flag_invert;
   double y_low,y_high;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y1 = 2;
   const int nvr_y2 = 3;
   const int ni_flag_invert = 0;
   const int nr_y_low = 0;
   const int nr_y_high = 1;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y1 = 2;
   const int no_y2 = 3;
   const int ng_1 = 0;
   const int ng_2 = 1;
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y1] = X.val_vr[nvr_y1];
     X.outprm[no_y2] = X.val_vr[nvr_y2];
     return;
   }

   flag_invert = X.iprm[ni_flag_invert];

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];

   if (flag_invert == 0) {
     if (x1 >= x2) {
       y10 = X.rprm[nr_y_high];
       y20 = X.rprm[nr_y_low];
     } else {
       y10 = X.rprm[nr_y_low];
       y20 = X.rprm[nr_y_high];
     }   
   } else {
     if (x1 >= x2) {
       y10 = X.rprm[nr_y_low];
       y20 = X.rprm[nr_y_high];
     } else {
       y10 = X.rprm[nr_y_high];
       y20 = X.rprm[nr_y_low];
     }   
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y1] = y10;
     X.val_vr[nvr_y2] = y20;
     return;
   }

   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y1] = y10;
       X.val_vr[nvr_y2] = y20;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y1] - y10;
         X.g[ng_2] = X.val_vr[nvr_y2] - y20;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y1] = 1.0;
         J.dgdvr[ng_2][nvr_y2] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_cmprh_1_1(Global &G,XbeUsr &X,XbeJac &J) {
   int iter;
   double y0,delx,time0,offset1;
   double x;
   double y;
   int flag_invert,flag_quad;
   double x0,y_low,y_high,h,x_1,x_2,t_1,t_2,epsl,delt_min,delt_nrml,hby2,y_old,
     y_half;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int ni_flag_invert = 0;
   const int ni_flag_quad = 1;
   const int nr_x0 = 0;
   const int nr_y_low = 1;
   const int nr_y_high = 2;
   const int nr_h = 3;
   const int nr_x_1 = 4;
   const int nr_x_2 = 5;
   const int nr_t_1 = 6;
   const int nr_t_2 = 7;
   const int nr_epsl = 8;
   const int nr_delt_min = 9;
   const int nr_delt_nrml = 10;
   const int nr_hby2 = 11;
   const int nr_y_old = 12;
   const int nr_y_half = 13;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     h = X.rprm[nr_h];
     hby2 = 0.5*h;
     X.rprm[nr_hby2] = hby2;

     y_low  = X.rprm[nr_y_low ];
     y_high = X.rprm[nr_y_high];

     if (y_low >= y_high) {
       cout << "cmprh_1_1.xbe: y_low <= y_high?" << endl;
       cout << "  y_low=" << y_low << endl;
       cout << "  y_high=" << y_high << endl;
       cout << "  Halting..." << endl;
       exit (1);
     }
     y_half = 0.5*(y_high+y_low);
     X.rprm[nr_y_half] = y_half;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   flag_invert = X.iprm[ni_flag_invert];
   flag_quad = X.iprm[ni_flag_quad];

   x = X.val_vr[nvr_x];

   if (G.flags[G.i_save_history]) {

     x_1 = X.rprm[nr_x_1];
     t_1 = X.rprm[nr_t_1];

     x_2 = x_1;
     t_2 = t_1;

     x_1 = x;
     t_1 = G.time_given_x;

     X.rprm[nr_x_1] = x_1;
     X.rprm[nr_x_2] = x_2;
     X.rprm[nr_t_1] = t_1;
     X.rprm[nr_t_2] = t_2;

     X.rprm[nr_y_old] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_next_time]) {
     iter = G.iter_trns_x;
     time0 = G.time_given_x;

     x0       = X.rprm[nr_x0    ];
     hby2     = X.rprm[nr_hby2  ];
     y_old    = X.rprm[nr_y_old ];
     y_half   = X.rprm[nr_y_half];

     if (flag_invert == 0) {
       if (y_old > y_half) {
         offset1 = x0 - hby2;
       } else {
         offset1 = x0 + hby2;
       }
     } else {
       if (y_old < y_half) {
         offset1 = x0 - hby2;
       } else {
         offset1 = x0 + hby2;
       }
     }
     x_1 = X.rprm[nr_x_1];
     x_2 = X.rprm[nr_x_2];
     t_1 = X.rprm[nr_t_1];
     t_2 = X.rprm[nr_t_2];

     epsl      = X.rprm[nr_epsl     ];
     delt_min  = X.rprm[nr_delt_min ];
     delt_nrml = X.rprm[nr_delt_nrml];

     G.time_nextbreak_x = get_tnext(iter,flag_quad,
       t_1,t_2,(x_1-offset1),(x_2-offset1),
       time0,(x-offset1),epsl,delt_min,delt_nrml);

     return;
   }
   hby2   = X.rprm[nr_hby2  ];
   y_low  = X.rprm[nr_y_low ];
   y_high = X.rprm[nr_y_high];
   y_half = X.rprm[nr_y_half];
   y_old  = X.rprm[nr_y_old ];
   x0     = X.rprm[nr_x0    ];

   if (G.flags[G.i_init_guess]) {
     delx = x-x0;

     if (flag_invert == 0) {
       if (delx > hby2) {
         y0 = y_high;
       } else if (delx < (-hby2)) {
         y0 = y_low;
       } else {
//       set arbitrarily to low
         y0 = y_low;
       }
     } else {
       if (delx > hby2) {
         y0 = y_low;
       } else if (delx < (-hby2)) {
         y0 = y_high;
       } else {
//       set arbitrarily to low
         y0 = y_low;
       }
     }
     X.val_vr[nvr_y] = y0;
     return;
   }

   delx = x-x0;

   if (flag_invert == 0) {
     if (y_old > y_half) {
       if (delx < (-hby2)) {
         y0 = y_low;
       } else {
         y0 = y_high;
       }
     } else {
       if (delx > hby2) {
         y0 = y_high;
       } else {
         y0 = y_low;
       }
     }
   } else {
     if (y_old < y_half) {
       if (delx < (-hby2)) {
         y0 = y_high;
       } else {
         y0 = y_low;
       }
     } else {
       if (delx > hby2) {
         y0 = y_low;
       } else {
         y0 = y_high;
       }
     }
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_cmprh_2_1(Global &G,XbeUsr &X,XbeJac &J) {
   int iter;
   double y0,delx,time0,offset1;
   double x1,x2;
   double y;
   int flag_invert,flag_quad;
   double x0,y_low,y_high,h,x1_1,x1_2,x2_1,x2_2,t_1,t_2,epsl,delt_min,delt_nrml,
     hby2,y_old,y_half;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int ni_flag_invert = 0;
   const int ni_flag_quad = 1;
   const int nr_x0 = 0;
   const int nr_y_low = 1;
   const int nr_y_high = 2;
   const int nr_h = 3;
   const int nr_x1_1 = 4;
   const int nr_x1_2 = 5;
   const int nr_x2_1 = 6;
   const int nr_x2_2 = 7;
   const int nr_t_1 = 8;
   const int nr_t_2 = 9;
   const int nr_epsl = 10;
   const int nr_delt_min = 11;
   const int nr_delt_nrml = 12;
   const int nr_hby2 = 13;
   const int nr_y_old = 14;
   const int nr_y_half = 15;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     h = X.rprm[nr_h];
     hby2 = 0.5*h;
     X.rprm[nr_hby2] = hby2;

     y_low  = X.rprm[nr_y_low ];
     y_high = X.rprm[nr_y_high];

     if (y_low >= y_high) {
       cout << "cmprh_2_1.xbe: y_low <= y_high?" << endl;
       cout << "  y_low=" << y_low << endl;
       cout << "  y_high=" << y_high << endl;
       cout << "  Halting..." << endl;
       exit (1);
     }
     y_half = 0.5*(y_high+y_low);
     X.rprm[nr_y_half] = y_half;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   flag_invert = X.iprm[ni_flag_invert];
   flag_quad = X.iprm[ni_flag_quad];

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];

   if (G.flags[G.i_save_history]) {

     x1_1 = X.rprm[nr_x1_1];
     x2_1 = X.rprm[nr_x2_1];
     t_1 = X.rprm[nr_t_1];

     x1_2 = x1_1;
     x2_2 = x2_1;
     t_2 = t_1;

     x1_1 = x1;
     x2_1 = x2;
     t_1 = G.time_given_x;

     X.rprm[nr_x1_1] = x1_1;
     X.rprm[nr_x1_2] = x1_2;
     X.rprm[nr_x2_1] = x2_1;
     X.rprm[nr_x2_2] = x2_2;
     X.rprm[nr_t_1] = t_1;
     X.rprm[nr_t_2] = t_2;

     X.rprm[nr_y_old] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_next_time]) {
     iter = G.iter_trns_x;
     time0 = G.time_given_x;

     x0       = X.rprm[nr_x0    ];
     hby2     = X.rprm[nr_hby2  ];
     y_old    = X.rprm[nr_y_old ];
     y_half   = X.rprm[nr_y_half];

     if (flag_invert == 0) {
       if (y_old > y_half) {
         offset1 = x0 - hby2;
       } else {
         offset1 = x0 + hby2;
       }
     } else {
       if (y_old < y_half) {
         offset1 = x0 - hby2;
       } else {
         offset1 = x0 + hby2;
       }
     }
     x1_1 = X.rprm[nr_x1_1];
     x1_2 = X.rprm[nr_x1_2];
     x2_1 = X.rprm[nr_x2_1];
     x2_2 = X.rprm[nr_x2_2];
     t_1 = X.rprm[nr_t_1];
     t_2 = X.rprm[nr_t_2];

     epsl      = X.rprm[nr_epsl     ];
     delt_min  = X.rprm[nr_delt_min ];
     delt_nrml = X.rprm[nr_delt_nrml];

     G.time_nextbreak_x = get_tnext(iter,flag_quad,
       t_1,t_2,(x1_1-x2_1-offset1),(x1_2-x2_2-offset1),
       time0,(x1-x2-offset1),epsl,delt_min,delt_nrml);

     return;
   }
   hby2   = X.rprm[nr_hby2  ];
   y_low  = X.rprm[nr_y_low ];
   y_high = X.rprm[nr_y_high];
   y_half = X.rprm[nr_y_half];
   y_old  = X.rprm[nr_y_old ];
   x0     = X.rprm[nr_x0    ];

   if (G.flags[G.i_init_guess]) {
     delx = x1-x2-x0;

     if (flag_invert == 0) {
       if (delx > hby2) {
         y0 = y_high;
       } else if (delx < (-hby2)) {
         y0 = y_low;
       } else {
//       set arbitrarily to low
         y0 = y_low;
       }
     } else {
       if (delx > hby2) {
         y0 = y_low;
       } else if (delx < (-hby2)) {
         y0 = y_high;
       } else {
//       set arbitrarily to low
         y0 = y_low;
       }
     }
     X.val_vr[nvr_y] = y0;
     return;
   }

   delx = x1-x2-x0;

   if (flag_invert == 0) {
     if (y_old > y_half) {
       if (delx < (-hby2)) {
         y0 = y_low;
       } else {
         y0 = y_high;
       }
     } else {
       if (delx > hby2) {
         y0 = y_high;
       } else {
         y0 = y_low;
       }
     }
   } else {
     if (y_old < y_half) {
       if (delx < (-hby2)) {
         y0 = y_high;
       } else {
         y0 = y_low;
       }
     } else {
       if (delx > hby2) {
         y0 = y_low;
       } else {
         y0 = y_high;
       }
     }
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_const(Global &G,XbeUsr &X,XbeJac &J) {
   double k0;
   double y;
   double k,k_scale;
   const int nvr_y = 0;
   const int nr_k = 0;
   const int nr_k_scale = 1;
   const int no_y = 0;
   const int ng_1 = 0;
// cout << "const.xbe" << endl;
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   k0 = X.rprm[nr_k_scale]*X.rprm[nr_k];
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = k0;
//   cout << "const.xbe: k0 = " << k0 << endl;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = k0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - k0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_cos(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   double x;
   double y;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = cos(X.val_vr[nvr_x]);
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = cos(X.val_vr[nvr_x]);
     } else if (G.flags[G.i_implicit]) {
       x = X.val_vr[nvr_x];
       y = X.val_vr[nvr_y];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = y - cos(x);
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] =  1.0;
         J.dgdvr[ng_1][nvr_x] = sin(x);
       }
     }
     return;
   }
   return;
}
void x_dead_zone(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   double x;
   double y;
   double xmin,xmax;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_xmin = 0;
   const int nr_xmax = 1;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   xmin = X.rprm[nr_xmin];
   xmax = X.rprm[nr_xmax];

   x = X.val_vr[nvr_x];
   if (x >= xmax) {
     y0 = x - xmax;
   } else if (x <= xmin) {
     y0 = x - xmin;
   } else {
     y0 = 0.0;
   }   

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       y = X.val_vr[nvr_y];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = y - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] =  1.0;
         if (x >= xmax) {
           J.dgdvr[ng_1][nvr_x] = 0.0;
         } else if (x <= xmin) {
           J.dgdvr[ng_1][nvr_x] = 0.0;
         } else {
           J.dgdvr[ng_1][nvr_x] = -1.0;
         }   
       }
     }
     return;
   }
   return;
}
void x_decoder_2_4(Global &G,XbeUsr &X,XbeJac &J) {
   const int NX=2;
   const int NY=4;
   int k1,XIN[NX];
   double Y0[NY];
   double x0,x1;
   double y0,y1,y2,y3;
   double y00,y01,y02,y03,y10,y11,y12,y13,y20,y21,y22,y23,y30,y31,y32,y33,
     x_high,hb2;
   const int nvr_x0 = 0;
   const int nvr_x1 = 1;
   const int nvr_y0 = 2;
   const int nvr_y1 = 3;
   const int nvr_y2 = 4;
   const int nvr_y3 = 5;
   const int nr_y00 = 0;
   const int nr_y01 = 1;
   const int nr_y02 = 2;
   const int nr_y03 = 3;
   const int nr_y10 = 4;
   const int nr_y11 = 5;
   const int nr_y12 = 6;
   const int nr_y13 = 7;
   const int nr_y20 = 8;
   const int nr_y21 = 9;
   const int nr_y22 = 10;
   const int nr_y23 = 11;
   const int nr_y30 = 12;
   const int nr_y31 = 13;
   const int nr_y32 = 14;
   const int nr_y33 = 15;
   const int nr_x_high = 16;
   const int nr_hb2 = 17;
   const int no_x0 = 0;
   const int no_x1 = 1;
   const int no_y0 = 2;
   const int no_y1 = 3;
   const int no_y2 = 4;
   const int no_y3 = 5;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   const int ng_4 = 3;
   if (G.flags[G.i_one_time_parms]) {
     x_high = X.rprm[nr_x_high];
     hb2 = 0.5*x_high;
     X.rprm[nr_hb2] = hb2;

     X.vec3d_1.resize(NX);
     for (int i = 0; i < NX; ++i) {
       X.vec3d_1[i].resize(NX);
       for (int j = 0; j < NX; ++j) {
         X.vec3d_1[i][j].resize(NY);
       }
     }
     k1 = 0;
     for (int i = 0; i < NX; ++i) {
       for (int j = 0; j < NX; ++j) {
         for (int k = 0; k < NY; ++k) {
           X.vec3d_1[i][j][k] = X.rprm[k1];
           k1++;
         }
       }
     }
     for (int i = 0; i < NX; ++i) {
       for (int j = 0; j < NX; ++j) {
         for (int k = 0; k < NY; ++k) {
           cout
             << " i: " << i
             << " j: " << j
             << " k: " << k
             << " vec3d_1: " << X.vec3d_1[i][j][k] << endl;
         }
       }
     }
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x0] = X.val_vr[nvr_x0];
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_y0] = X.val_vr[nvr_y0];
     X.outprm[no_y1] = X.val_vr[nvr_y1];
     X.outprm[no_y2] = X.val_vr[nvr_y2];
     X.outprm[no_y3] = X.val_vr[nvr_y3];
     return;
   }
   hb2 = X.rprm[nr_hb2];

   for (int i = 0; i < NX; ++i) {
     if (X.val_vr[i] > hb2) {
       XIN[i] = 1;
     } else {
       XIN[i] = 0;
     }
   }
   for (int i = 0; i < NY; ++i) {
     Y0[i] = X.vec3d_1[XIN[1]][XIN[0]][i];
   }

   if (G.flags[G.i_init_guess]) {
     for (int i = 0; i < NY; ++i) {
       X.val_vr[i + NX] = Y0[i];
     }
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       for (int i = 0; i < NY; ++i) {
         X.val_vr[i + NX] = Y0[i];
       }
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         for (int i = 0; i < NY; ++i) {
           X.g[i] = X.val_vr[i + NX] - Y0[i];
         }
       }
       if (G.flags[G.i_jacobian]) {
         for (int i = 0; i < NY; ++i) {
           J.dgdvr[i][i + NX] = 1.0;
         }
       }
     }
     return;
   }
   return;
}
void x_delay_discrete(Global &G,XbeUsr &X,XbeJac &J) {
   double time0,rsd1,rsd2;
   double t0_new,t_a,t_b,t_c,t_d,y0;
   int n;
   int l_cross;
   double x;
   double y;
   int n_delay;
   double T,t0,dt,y_current,y_old_1,y_old_2,y_old_3,y_old_4,epsl1,epsl2;
   double y_st;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int ni_n_delay = 0;
   const int nr_T = 0;
   const int nr_t0 = 1;
   const int nr_dt = 2;
   const int nr_y_current = 3;
   const int nr_y_old_1 = 4;
   const int nr_y_old_2 = 5;
   const int nr_y_old_3 = 6;
   const int nr_y_old_4 = 7;
   const int nr_epsl1 = 8;
   const int nr_epsl2 = 9;
   const int nst_y_st = 0;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     n_delay = X.iprm[ni_n_delay];
     if ((n_delay < 1) || (n_delay > 4)) {
       cout << "delay_discrete: n_delay must be 1, 2, or 3. Halting.." << endl;
       exit (1);
     }
     dt = X.rprm[nr_dt];
     epsl1 = dt/10.0;
     epsl2 = dt/100.0;
     X.rprm[nr_epsl1] = epsl1;
     X.rprm[nr_epsl2] = epsl2;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
       }   
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }   
     }   
     return;
   }   
   if (G.flags[G.i_trns] || G.flags[G.i_save_history]) {

     time0 = G.time_given_x;
     T = X.rprm[nr_T];
     y_current = X.rprm[nr_y_current];

     t0 = X.rprm[nr_t0];
     dt = X.rprm[nr_dt];

     epsl1 = X.rprm[nr_epsl1];
     epsl2 = X.rprm[nr_epsl2];

     n_delay = X.iprm[ni_n_delay];

     if (time0 < t0) {
       n = ((int)((t0-time0)/T)) + 1;
       t0_new = t0 - n*T;
     } else {
       t0_new = t0;
     }
               
     t_a = (time0-t0_new);
     t_b = fmod(t_a,T);
     if (abs(t_b-T) < epsl2) t_b = 0.0;
     t_c = T - dt;

     t_d = abs(t_b-T);
     l_cross = 0;
     if ((t_d < epsl1) || (abs(t_b) < epsl1)) {
       l_cross = 1;
     }
     if (l_cross == 1) {
       if (n_delay == 1) {
         y0 = X.rprm[nr_y_current];
       } else if (n_delay == 2) {
         y0 = X.rprm[nr_y_old_1];
       } else if (n_delay == 3) {
         y0 = X.rprm[nr_y_old_2];
       } else if (n_delay == 4) {
         y0 = X.rprm[nr_y_old_3];
       }
     } else {
       if (n_delay == 1) {
         y0 = X.rprm[nr_y_old_1];
       } else if (n_delay == 2) {
         y0 = X.rprm[nr_y_old_2];
       } else if (n_delay == 3) {
         y0 = X.rprm[nr_y_old_3];
       } else if (n_delay == 4) {
         y0 = X.rprm[nr_y_old_4];
       }
     }
   }
   if (G.flags[G.i_save_history]) {
     if (l_cross == 1) {
       X.rprm[nr_y_old_4] = X.rprm[nr_y_old_3];
       X.rprm[nr_y_old_3] = X.rprm[nr_y_old_2];
       X.rprm[nr_y_old_2] = X.rprm[nr_y_old_1];
       X.rprm[nr_y_old_1] = X.rprm[nr_y_current];

       X.rprm[nr_y_current] = X.val_vr[nvr_x];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_delay_discrete_1(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   int n;
   double x;
   double y;
   int n_delay,sampler_index;
   double y_current,y_old_1,y_old_2,y_old_3,y_old_4;
   double y_st;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int ni_n_delay = 0;
   const int ni_sampler_index = 1;
   const int nr_y_current = 0;
   const int nr_y_old_1 = 1;
   const int nr_y_old_2 = 2;
   const int nr_y_old_3 = 3;
   const int nr_y_old_4 = 4;
   const int nst_y_st = 0;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;

   if (G.flags[G.i_one_time_parms]) {
     n_delay = X.iprm[ni_n_delay];
     if ((n_delay < 1) || (n_delay > 4)) {
       cout << "delay_discrete_1: n_delay must be 1, 2, or 3. Halting.." << endl;
       exit (1);
     }
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
       }   
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }   
     }   
     return;
   }   
   if (G.flags[G.i_save_history]) {
     sampler_index = X.iprm[ni_sampler_index];

     if (G.sampler_flag[sampler_index] == 1) {
       X.rprm[nr_y_old_4] = X.rprm[nr_y_old_3];
       X.rprm[nr_y_old_3] = X.rprm[nr_y_old_2];
       X.rprm[nr_y_old_2] = X.rprm[nr_y_old_1];
       X.rprm[nr_y_old_1] = X.rprm[nr_y_current];

       X.rprm[nr_y_current] = X.val_vr[nvr_x];
     }
     return;
   }

   if (G.flags[G.i_trns]) {
     sampler_index = X.iprm[ni_sampler_index];

     n_delay = X.iprm[ni_n_delay];

     if (G.sampler_flag[sampler_index] == 1) {
       if (n_delay == 1) {
         y0 = X.rprm[nr_y_current];
       } else if (n_delay == 2) {
         y0 = X.rprm[nr_y_old_1];
       } else if (n_delay == 3) {
         y0 = X.rprm[nr_y_old_2];
       } else if (n_delay == 4) {
         y0 = X.rprm[nr_y_old_3];
       }
     } else {
       if (n_delay == 1) {
         y0 = X.rprm[nr_y_old_1];
       } else if (n_delay == 2) {
         y0 = X.rprm[nr_y_old_2];
       } else if (n_delay == 3) {
         y0 = X.rprm[nr_y_old_3];
       } else if (n_delay == 4) {
         y0 = X.rprm[nr_y_old_4];
       }
     }

     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_delay_onestep(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double y;
   double x_last;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_x_last = 0;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     X.rprm[nr_x_last] = 0.0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.val_vr[nvr_x];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - X.val_vr[nvr_x];
       }   
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }   
     }   
     return;
   }   
   if (G.flags[G.i_save_history]) {
     X.rprm[nr_x_last] = X.val_vr[nvr_x];
     return;
   }
   if (G.flags[G.i_trns]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.rprm[nr_x_last];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - X.rprm[nr_x_last];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_delay_onestep_1(Global &G,XbeUsr &X,XbeJac &J) {
   double time0;
   bool flag_edge;
   double x;
   double y;
   double x_high,dt,x_last,x_cross;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_x_high = 0;
   const int nr_dt = 1;
   const int nr_x_last = 2;
   const int nr_x_cross = 3;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     x_high = X.rprm[nr_x_high];

     if (x_high < 0.0) {
       cout << "delay_onestep_1.xbe: check x_high. Halting..." << endl;
       exit(1);
     }
     X.rprm[nr_x_cross] = 0.5*x_high;
     X.rprm[nr_x_last] = 0.0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_next_time]) {
     time0 = G.time_given_x;

     x_last  = X.rprm[nr_x_last ];
     x_cross = X.rprm[nr_x_cross];

     x = X.val_vr[nvr_x];

     flag_edge = false;

     if ((x_last <= x_cross) && (x > x_cross)) {
       flag_edge = true;
     } else if ((x_last >= x_cross) && (x < x_cross)) {
       flag_edge = true;
     }
     if (flag_edge) {
       G.time_nextbreak_x = time0 + X.rprm[nr_dt];
     } else {
       G.time_nextbreak_x = G.time_end;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.val_vr[nvr_x];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - X.val_vr[nvr_x];
       }   
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }   
     }   
     return;
   }   
   if (G.flags[G.i_save_history]) {
     X.rprm[nr_x_last] = X.val_vr[nvr_x];
     return;
   }
   if (G.flags[G.i_trns]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.rprm[nr_x_last];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - X.rprm[nr_x_last];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_diff(Global &G,XbeUsr &X,XbeJac &J) {
   double x1,x2;
   double y;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = X.val_vr[nvr_x1] - X.val_vr[nvr_x2];
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.val_vr[nvr_x1] - X.val_vr[nvr_x2];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - X.val_vr[nvr_x1] + X.val_vr[nvr_x2];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
         J.dgdvr[ng_1][nvr_x1] = -1.0;
         J.dgdvr[ng_1][nvr_x2] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_div(Global &G,XbeUsr &X,XbeJac &J) {
   double x20,x2inv;
   double x1,x2;
   double y;
   double k,delta;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int nr_k = 0;
   const int nr_delta = 1;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   k = X.rprm[nr_k];
   delta = X.rprm[nr_delta];

   if (G.flags[G.i_init_guess]) {
     x1 = X.val_vr[nvr_x1];
     x2 = X.val_vr[nvr_x2];
     if (abs(x2) <= delta) {
       if (x2 < 0) {
         x20 = -delta;
       } else {
         x20 = delta;
       }
     } else {
       x20 = x2;
     }
     X.val_vr[nvr_y] = k*x1/x20;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];
   if (abs(x2) <= delta) {
     if (x2 < 0) {
       x20 = -delta;
     } else {
       x20 = delta;
     }
   } else {
     x20 = x2;
   }
   x2inv = 1.0/x20;

   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = k*x1*x2inv;
     } else if (G.flags[G.i_implicit]) {
       y = X.val_vr[nvr_y];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = y - k*x1*x2inv;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
         J.dgdvr[ng_1][nvr_x1] = -k*x2inv;
         J.dgdvr[ng_1][nvr_x2] =  k*x1*x2inv*x2inv;
       }
     }
     return;
   }
   return;
}
void x_dq0_to_abc_2(Global &G,XbeUsr &X,XbeJac &J) {
   static double k1=0.5,k2=0;
   double c1,s1;
   double xd,xq,x0,theta;
   double xa,xb,xc;
   double c,s;
   const int nvr_xd = 0;
   const int nvr_xq = 1;
   const int nvr_x0 = 2;
   const int nvr_theta = 3;
   const int nvr_xa = 4;
   const int nvr_xb = 5;
   const int nvr_xc = 6;
   const int no_xd = 0;
   const int no_xq = 1;
   const int no_x0 = 2;
   const int no_xa = 3;
   const int no_xb = 4;
   const int no_xc = 5;
   const int na_c = 0;
   const int na_s = 1;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   const int ng_4 = 3;
   const int ng_5 = 4;
   if (G.flags[G.i_one_time_parms]) {
     k2 = 0.5*sqrt(3.0);
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_xd] = X.val_vr[nvr_xd];
     X.outprm[no_xq] = X.val_vr[nvr_xq];
     X.outprm[no_x0] = X.val_vr[nvr_x0];
     X.outprm[no_xa] = X.val_vr[nvr_xa];
     X.outprm[no_xb] = X.val_vr[nvr_xb];
     X.outprm[no_xc] = X.val_vr[nvr_xc];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     xd = X.val_vr[nvr_xd];
     xq = X.val_vr[nvr_xq];
     x0 = X.val_vr[nvr_x0];

     theta = X.val_vr[nvr_theta];

     c = cos(theta);
     s = sin(theta);

     X.val_vr [nvr_xa] = xd*c + xq*s + x0;
     X.val_vr [nvr_xb] = - xd*(k1*c-k2*s) - xq*(k2*c+k1*s) + x0;
     X.val_vr [nvr_xc] = - xd*(k1*c+k2*s) - xq*(-k2*c+k1*s) + x0;
     X.val_aux[na_c  ] = c;
     X.val_aux[na_s  ] = s;

     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     xd = X.val_vr[nvr_xd];
     xq = X.val_vr[nvr_xq];
     x0 = X.val_vr[nvr_x0];

     theta = X.val_vr[nvr_theta];

     c1 = cos(theta);
     s1 = sin(theta);

     if (G.flags[G.i_explicit]) {
       X.val_vr [nvr_xa] = xd*c1 + xq*s1 + x0;
       X.val_vr [nvr_xb] = - xd*(k1*c1-k2*s1) - xq*(k2*c1+k1*s1) + x0;
       X.val_vr [nvr_xc] = - xd*(k1*c1+k2*s1) - xq*(-k2*c1+k1*s1) + x0;
       X.val_aux[na_c  ] = c1;
       X.val_aux[na_s  ] = s1;
     } else if (G.flags[G.i_implicit]) {
       xa = X.val_vr [nvr_xa];
       xb = X.val_vr [nvr_xb];
       xc = X.val_vr [nvr_xc];
       c  = X.val_aux[na_c  ];
       s  = X.val_aux[na_s  ];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = xa + xd*(-c) + xq*(-s) - x0;
         X.g[ng_2] = xb + xd*(k1*c-k2*s) + xq*(k2*c+k1*s) - x0;
         X.g[ng_3] = xc + xd*(k1*c+k2*s) + xq*(-k2*c+k1*s) - x0;
         X.g[ng_4] = c - c1;
         X.g[ng_5] = s - s1;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr [ng_1][nvr_xa] = 1.0;
         J.dgdvr [ng_1][nvr_xd] = -c;
         J.dgdvr [ng_1][nvr_xq] = -s;
         J.dgdvr [ng_1][nvr_x0] = -1.0;
         J.dgdaux[ng_1][na_c  ] = -xd;
         J.dgdaux[ng_1][na_s  ] = -xq;

         J.dgdvr [ng_2][nvr_xb] = 1.0;
         J.dgdvr [ng_2][nvr_xd] = k1*c-k2*s;
         J.dgdvr [ng_2][nvr_xq] = k2*c+k1*s;
         J.dgdvr [ng_2][nvr_x0] = -1.0;
         J.dgdaux[ng_2][na_c  ] = xd*k1+xq*k2;
         J.dgdaux[ng_2][na_s  ] = -xd*k2+xq*k1;

         J.dgdvr [ng_3][nvr_xc] = 1.0;
         J.dgdvr [ng_3][nvr_xd] = k1*c+k2*s;
         J.dgdvr [ng_3][nvr_xq] = -k2*c+k1*s;
         J.dgdvr [ng_3][nvr_x0] = -1.0;
         J.dgdaux[ng_3][na_c  ] = xd*k1-xq*k2;
         J.dgdaux[ng_3][na_s  ] = xd*k2+xq*k1;

         J.dgdaux[ng_4][na_c     ] = 1.0;
         J.dgdvr [ng_4][nvr_theta] =  s1;

         J.dgdaux[ng_5][na_s     ] = 1.0;
         J.dgdvr [ng_5][nvr_theta] = -c1;
       }
     }
     return;
   }
   return;
}
void x_dq_to_abc(Global &G,XbeUsr &X,XbeJac &J) {
   static double k1=0,k2=0,k3=0;
   double d,q,cost,sint;
   double a,b,c;
   double alpha,beta;
   const int nvr_d = 0;
   const int nvr_q = 1;
   const int nvr_cost = 2;
   const int nvr_sint = 3;
   const int nvr_a = 4;
   const int nvr_b = 5;
   const int nvr_c = 6;
   const int no_d = 0;
   const int no_q = 1;
   const int no_a = 2;
   const int no_b = 3;
   const int no_c = 4;
   const int na_alpha = 0;
   const int na_beta = 1;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   const int ng_4 = 3;
   const int ng_5 = 4;
   if (G.flags[G.i_one_time_parms]) {
     k1 = 2.0/3.0;
     k2 = 1.0/3.0;
     k3 = 1.0/(sqrt(3.0));
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_d] = X.val_vr[nvr_d];
     X.outprm[no_q] = X.val_vr[nvr_q];
     X.outprm[no_a] = X.val_vr[nvr_a];
     X.outprm[no_b] = X.val_vr[nvr_b];
     X.outprm[no_c] = X.val_vr[nvr_c];
     return;
   }
   cost = X.val_vr[nvr_cost];
   sint = X.val_vr[nvr_sint];
   d = X.val_vr[nvr_d];
   q = X.val_vr[nvr_q];

   if (G.flags[G.i_init_guess]) {
     alpha = cost*d - sint*q;
     beta  = sint*d + cost*q;

     X.val_aux[na_alpha] = alpha;
     X.val_aux[na_beta ] = beta;
     X.val_vr[nvr_a] =  k1*alpha;
     X.val_vr[nvr_b] = -k2*alpha + k3*beta;
     X.val_vr[nvr_c] = -k2*alpha - k3*beta;

     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       alpha = cost*d - sint*q;
       beta  = sint*d + cost*q;

       X.val_aux[na_alpha] = alpha;
       X.val_aux[na_beta ] = beta;

       X.val_vr[nvr_a] =  k1*alpha;
       X.val_vr[nvr_b] = -k2*alpha + k3*beta;
       X.val_vr[nvr_c] = -k2*alpha - k3*beta;
     } else if (G.flags[G.i_implicit]) {
       alpha = X.val_aux[na_alpha];
       beta  = X.val_aux[na_beta ];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = alpha - cost*d + sint*q;
         X.g[ng_2] = beta - sint*d - cost*q;

         X.g[ng_3] = a - k1*alpha;
         X.g[ng_4] = b + k2*alpha - k3*beta;
         X.g[ng_5] = c + k2*alpha + k3*beta;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdaux[ng_1][na_alpha] =  1.0;
         J.dgdvr [ng_1][nvr_cost] = -d;
         J.dgdvr [ng_1][nvr_sint] =  q;
         J.dgdvr [ng_1][nvr_d   ] = -cost;
         J.dgdvr [ng_1][nvr_q   ] =  sint;

         J.dgdaux[ng_2][na_beta ] =  1.0;
         J.dgdvr [ng_2][nvr_cost] = -q;
         J.dgdvr [ng_2][nvr_sint] = -d;
         J.dgdvr [ng_2][nvr_d   ] = -sint;
         J.dgdvr [ng_2][nvr_q   ] = -cost;

         J.dgdvr [ng_3][nvr_a   ] =  1.0;
         J.dgdaux[ng_3][na_alpha] = -k1;

         J.dgdvr [ng_4][nvr_b   ] =  1.0;
         J.dgdaux[ng_4][na_alpha] =  k2;
         J.dgdaux[ng_4][na_beta ] = -k3;

         J.dgdvr [ng_5][nvr_c   ] =  1.0;
         J.dgdaux[ng_5][na_alpha] =  k2;
         J.dgdaux[ng_5][na_beta ] =  k3;
       }
     }
     return;
   }
   return;
}
void x_dummy_sink(Global &G,XbeUsr &X,XbeJac &J) {
   double y;
   const int nvr_y = 0;
   const int no_y = 0;
   return;
}
void x_dummy_source(Global &G,XbeUsr &X,XbeJac &J) {
   double y;
   const int nvr_y = 0;
   const int no_y = 0;
   return;
}
void x_edge_delay(Global &G,XbeUsr &X,XbeJac &J) {
   const int flag_empty=0;
   const int flag_ltoh =1;
   const int flag_htol =2;

   double t0,t1,time0,delta,tmin;
   double y0,t_next_1,t_next_2;
   int flag1,flag2,nloc,nextloc,top,l;
   double x;
   double y;
   int n_delay,flag_frequency,flag_period,flag_zero_delay;
   double x_low,x_high,frequency,T,theta_delay,theta_delay_1,t_delay,x_last,
     t_low_to_high,t_high_to_low,x_cross,epsl2,epsl3;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int ni_n_delay = 0;
   const int ni_flag_frequency = 1;
   const int ni_flag_period = 2;
   const int ni_flag_zero_delay = 3;
   const int nr_x_low = 0;
   const int nr_x_high = 1;
   const int nr_frequency = 2;
   const int nr_T = 3;
   const int nr_theta_delay = 4;
   const int nr_theta_delay_1 = 5;
   const int nr_t_delay = 6;
   const int nr_x_last = 7;
   const int nr_t_low_to_high = 8;
   const int nr_t_high_to_low = 9;
   const int nr_x_cross = 10;
   const int nr_epsl2 = 11;
   const int nr_epsl3 = 12;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     X.edge_delay_nmax = 50;
     X.edge_delay_nloc = X.edge_delay_nmax;
     X.edge_delay_tchange.resize(X.edge_delay_nmax);
     X.edge_delay_flag.resize(X.edge_delay_nmax);

     flag_frequency = X.iprm[ni_flag_frequency];
     flag_period    = X.iprm[ni_flag_period   ];

     if ((flag_frequency == 0) && (flag_period == 0)) {
       cout << "edge_delay.xbe: check flag_frequency and flag_period" << endl;
       cout << "  Both cannot be zero." << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     if ((flag_frequency != 0) && (flag_period != 0)) {
       cout << "edge_delay.xbe: check flag_frequency and flag_period" << endl;
       cout << "  Both cannot be non-zero." << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     if (flag_frequency != 0) {
       frequency = X.rprm[nr_frequency];
       T = 1.0/frequency;
       X.rprm[nr_T] = T;
     }
     if (flag_period != 0) {
       T = X.rprm[nr_T];
       frequency = 1.0/T;
       X.rprm[nr_frequency] = frequency;
     }

     x_high = X.rprm[nr_x_high];
     x_low  = X.rprm[nr_x_low ];
     if (x_high < x_low) {
       cout << "edge_delay.xbe: x_high and x_low are not defined correctly." << endl;
       cout << "  x_high=" << x_high << " x_low=" << x_low << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     X.rprm[nr_x_cross] = 0.5*(x_high+x_low);

     n_delay = X.iprm[ni_n_delay];
     theta_delay = X.rprm[nr_theta_delay];
     theta_delay_1 = X.rprm[nr_theta_delay_1];
     t_delay = (T/360.0)*(theta_delay*((double)(n_delay))+theta_delay_1);
     X.rprm[nr_t_delay] = t_delay;
     if (t_delay < 1.0e-8) {
       X.iprm[ni_flag_zero_delay] = 1;
     } else {
       X.iprm[ni_flag_zero_delay] = 0;
     }

     epsl2 = 0.001*T/360.0;
     X.rprm[nr_epsl2] = epsl2;
     X.rprm[nr_epsl3] = 1.1*epsl2;

     X.rprm[nr_t_low_to_high] = G.time_end;
     X.rprm[nr_t_high_to_low] = G.time_end;

     X.rprm[nr_x_last] = x_low;

//   initialise

     X.edge_delay_nextloc = 0;
     X.edge_delay_top     = 0;

     nloc = X.edge_delay_nloc;
     for (int i=0; i < nloc; i++) {
       X.edge_delay_flag[i] = flag_empty;
     }
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = 0.0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }

   flag_zero_delay = X.iprm[ni_flag_zero_delay];

   t_delay = X.rprm[nr_t_delay];
   nloc    = X.edge_delay_nloc;
   nextloc = X.edge_delay_nextloc;
   top     = X.edge_delay_top;

   time0 = G.time_given_x;

   if (G.flags[G.i_save_history]) {
     if (flag_zero_delay == 0) {
       x_last   = X.rprm[nr_x_last ];
       x = X.val_vr[nvr_x];
       x_cross = X.rprm[nr_x_cross];
       if ((x >= x_cross) && (x_last <= x_cross)) {
         t0 = time0 + t_delay;

         X.edge_delay_tchange[nextloc] = t0;
         X.edge_delay_flag   [nextloc] = flag_ltoh;

         for (int i=0; i < nloc; i++) {
           if (X.edge_delay_flag[i] == flag_empty) {
             X.edge_delay_nextloc = i;
             if (i > top) X.edge_delay_top = i;
             goto jump0;
           }
         }
         cout << "edge_delay: empty location not found (1). Halting..." << endl;
         exit(1);
         jump0:;
       }
       if ((x <= x_cross) && (x_last >= x_cross)) {
         t0 = time0 + t_delay;
         X.edge_delay_tchange[nextloc] = t0;
         X.edge_delay_flag   [nextloc] = flag_htol;

         for (int i=0; i < nloc; i++) {
           if (X.edge_delay_flag[i] == flag_empty) {
             X.edge_delay_nextloc = i;
             if (i > top) X.edge_delay_top = i;
             goto jump1;
           }
         }
         cout << "edge_delay: empty location not found (2). Halting..." << endl;
         exit(1);
         jump1:;
       }
       X.rprm[nr_x_last] = x;

       top = X.edge_delay_top;
       for (int i=0; i <= top; i++) {
         if (time0 > X.edge_delay_tchange[i]) {
           X.edge_delay_flag[i] = flag_empty;
         }
       }
       for (int i=0; i < nloc; i++) {
         if (X.edge_delay_flag[i] == flag_empty) {
           X.edge_delay_nextloc = i;
           if (i > top) X.edge_delay_top = i;
           goto jump2;
         }
       }
       cout << "edge_delay: empty location not found (3). Halting..." << endl;
       exit(1);
       jump2:;
     }
     return;
   }

   if (G.flags[G.i_next_time]) {
     if (flag_zero_delay == 0) {
       epsl2 = X.rprm[nr_epsl2];
       epsl3 = X.rprm[nr_epsl3];

       tmin = G.time_end;
       l = -1;

       for (int i=0; i <= top; i++) {
         flag1 = X.edge_delay_flag[i];
         if (flag1 != flag_empty) {
           t1 = X.edge_delay_tchange[i];
           if (t1 < tmin) {
             tmin = t1;
             l = i;
           }
         }
       }

       if (l == -1) {
         t_next_1 = G.time_end;
       } else {
         if (time0 < tmin) {
           delta = tmin - time0;
           if (delta < epsl3) {
             t_next_1 = tmin + epsl2;
           } else {
             t_next_1 = tmin - epsl2;
           }
         }
       }
       G.time_nextbreak_x = t_next_1;
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     if (flag_zero_delay == 1) {
       y0 = X.val_vr[nvr_x];
     } else {
       tmin = G.time_end;
       l = -1;

       for (int i=0; i <= top; i++) {
         flag1 = X.edge_delay_flag[i];
         if (flag1 != flag_empty) {
           t1 = X.edge_delay_tchange[i];
           if (t1 < tmin) {
             tmin = t1;
             l = i;
           }
         }
       }
       if (l == -1) {
         y0 = X.val_vr[nvr_y];
       } else {
         if (time0 >= tmin) {
           flag2 = X.edge_delay_flag[l];
           if (flag2 == flag_ltoh) {
             y0 = X.rprm[nr_x_high];
           } else if (flag2 == flag_htol) {
             y0 = X.rprm[nr_x_low];
           } else {
             cout << "edge_delay: flag2=" << flag2 << " is not expected." << endl;
             cout << "   Halting..." << endl; exit(1);
           }
         } else {
           y0 = X.val_vr[nvr_y];
         }
       }
     }
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_edge_delay_1(Global &G,XbeUsr &X,XbeJac &J) {
   int flag_active_edge;
   double time0,y0;
   double x;
   double y;
   double x_high,t_delay,t_rise,x_prev,x_cross,t2,epsl1;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_x_high = 0;
   const int nr_t_delay = 1;
   const int nr_t_rise = 2;
   const int nr_x_prev = 3;
   const int nr_x_cross = 4;
   const int nr_t2 = 5;
   const int nr_epsl1 = 6;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     x_high  = X.rprm[nr_x_high ];
     t_delay = X.rprm[nr_t_delay];

     if (x_high < 0.0) {
       cout << "edge_delay_1.xbe: x_high < 0? Halting..." << endl;
       exit(1);
     }
     x_cross = 0.5*x_high;
     X.rprm[nr_x_cross] = x_cross;

     t_rise = X.rprm[nr_t_rise];
     epsl1 = 1.1*t_rise;

     X.rprm[nr_epsl1] = epsl1;

     X.rprm[nr_x_prev] = 0.0;
     X.rprm[nr_t2    ] = G.time_begin;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_save_history]) {
     time0 = G.time_given_x;

     x_prev  = X.rprm[nr_x_prev ];
     x_cross = X.rprm[nr_x_cross];

     x = X.val_vr[nvr_x];

     if ((x_prev <= x_cross) && (x >= x_cross)) {
       flag_active_edge = 1;
     } else {
       flag_active_edge = 0;
     }
     if (flag_active_edge == 1) {
       X.rprm[nr_t2] = time0 + X.rprm[nr_t_delay];
     }
     X.rprm[nr_x_prev] = x;

     return;
   }
   if (G.flags[G.i_next_time]) {
     time0 = G.time_given_x;
     t2 = X.rprm[nr_t2];

     x_prev  = X.rprm[nr_x_prev ];
     x_cross = X.rprm[nr_x_cross];
     x = X.val_vr[nvr_x];

     if ((x_prev <= x_cross) && (x >= x_cross)) {
       flag_active_edge = 1;
     } else {
       flag_active_edge = 0;
     }
     if (flag_active_edge == 1) {
       G.time_nextbreak_x = time0 + X.rprm[nr_t_delay];
       return;
     }
     if (time0 < t2) {
       t_rise = X.rprm[nr_t_rise];
       epsl1 = X.rprm[nr_epsl1];

       if ((t2-time0) <= epsl1) {
         G.time_nextbreak_x = t2 + t_rise;
       } else {
         G.time_nextbreak_x = t2 - t_rise;
       }
     } else {
       G.time_nextbreak_x = G.time_end;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = 0.0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     time0 = G.time_given_x;
     x = X.val_vr[nvr_x];

     x_cross = X.rprm[nr_x_cross];
     x_prev  = X.rprm[nr_x_prev ];
     t2      = X.rprm[nr_t2     ];

     if (time0 < t2) {
       y0 = 0.0;
     } else {
       if ((x < x_cross) || (x_prev < x_cross)) {
         y0 = 0.0;
         X.rprm[nr_t2] = G.time_begin;
       } else {
         y0 = X.rprm[nr_x_high];
       }
     }

     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
   }
   return;
}
void x_indmc1(Global &G,XbeUsr &X,XbeJac &J) {
   double ids,iqs,idr,iqr,wr,tem0;
   double p;
   static double k4=0;
   double ids_psids,ids_psidr,iqs_psiqs,iqs_psiqr;
   double idr_psids,idr_psidr,iqr_psiqs,iqr_psiqr;
   double tem0_iqs,tem0_idr,tem0_ids,tem0_iqr;
   double tem0_psids,tem0_psidr,tem0_psiqs,tem0_psiqr;
   double wr_wrm;
   double vqs,vds,tl;
   double wrm;
   double psids,psidr,psiqs,psiqr;
   int poles;
   double rs,lls,lm,llr,rr,j,ls,lr,le,l1,l2,l3,x1,x2;
   double psids0,psiqs0,psidr0,psiqr0,wrm0;
   const int nvr_vqs = 0;
   const int nvr_vds = 1;
   const int nvr_tl = 2;
   const int nvr_wrm = 3;
   const int ni_poles = 0;
   const int nr_rs = 0;
   const int nr_lls = 1;
   const int nr_lm = 2;
   const int nr_llr = 3;
   const int nr_rr = 4;
   const int nr_j = 5;
   const int nr_ls = 6;
   const int nr_lr = 7;
   const int nr_le = 8;
   const int nr_l1 = 9;
   const int nr_l2 = 10;
   const int nr_l3 = 11;
   const int nr_x1 = 12;
   const int nr_x2 = 13;
   const int nst_psids0 = 0;
   const int nst_psiqs0 = 1;
   const int nst_psidr0 = 2;
   const int nst_psiqr0 = 3;
   const int nst_wrm0 = 4;
   const int no_wrm = 0;
   const int no_tem = 1;
   const int no_vds = 2;
   const int no_vqs = 3;
   const int no_ia = 4;
   const int no_ib = 5;
   const int no_ic = 6;
   const int na_psids = 0;
   const int na_psidr = 1;
   const int na_psiqs = 2;
   const int na_psiqr = 3;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nf_4 = 3;
   const int nf_5 = 4;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   const int ng_4 = 3;
   const int ng_5 = 4;
// cout << "indmc1.xbe" << endl;

   if (G.flags[G.i_one_time_parms]) {
     k4 = 0.5*(sqrt(3.0));

     lls = X.rprm[nr_lls];
     lm  = X.rprm[nr_lm ];
     llr = X.rprm[nr_llr];

     ls = lls + lm;
     lr = llr + lm;
     le = (ls*lr/lm) - lm;
     l1 = lr/(lm*le);
     l2 = 1.0 + (lls/lm);
     l3 = lls/lm;

     X.rprm[nr_ls] = ls;
     X.rprm[nr_lr] = lr;
     X.rprm[nr_le] = le;
     X.rprm[nr_l1] = l1;
     X.rprm[nr_l2] = l2;
     X.rprm[nr_l3] = l3;

     poles = X.iprm[ni_poles];
     p = (double)(poles);
     x1 = 0.75*p*lm;
     x2 = 0.5*p;

     X.rprm[nr_x1] = x1;
     X.rprm[nr_x2] = x2;

     return;
   }
   if (G.flags[G.i_init_guess]) {
//   for init guess, use start-up parameters (for this element)
//   can add separate igparms for this purpose later.
     X.val_aux[na_psids] = X.stprm[nst_psids0];
     X.val_aux[na_psiqs] = X.stprm[nst_psiqs0];
     X.val_aux[na_psidr] = X.stprm[nst_psidr0];
     X.val_aux[na_psiqr] = X.stprm[nst_psiqr0];
     X.val_vr [nvr_wrm ] = X.stprm[nst_wrm0  ];
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_wrm] = X.val_vr[nvr_wrm];
     X.outprm[no_vds] = X.val_vr[nvr_vds];
     X.outprm[no_vqs] = X.val_vr[nvr_vqs];

     psids = X.val_aux[na_psids];
     psidr = X.val_aux[na_psidr];
     psiqs = X.val_aux[na_psiqs];
     psiqr = X.val_aux[na_psiqr];

     le  = X.rprm[nr_le];
     lm  = X.rprm[nr_lm];
     l1  = X.rprm[nr_l1];
     l2  = X.rprm[nr_l2];

     ids = (l1*psids) - (psidr/le);
     iqs = (l1*psiqs) - (psiqr/le);
     idr = (psids/lm) - (l2*ids);
     iqr = (psiqs/lm) - (l2*iqs);

     X.outprm[no_ia] = iqs;
     X.outprm[no_ib] = -0.5*iqs-k4*ids;
     X.outprm[no_ic] = -0.5*iqs+k4*ids;

     x1 = X.rprm[nr_x1];
     tem0 = x1*(iqs*idr-ids*iqr);
     X.outprm[no_tem] = tem0;

     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_aux[na_psids] = X.stprm[nst_psids0];
       X.val_aux[na_psiqs] = X.stprm[nst_psiqs0];
       X.val_aux[na_psidr] = X.stprm[nst_psidr0];
       X.val_aux[na_psiqr] = X.stprm[nst_psiqr0];
       X.val_vr [nvr_wrm ] = X.stprm[nst_wrm0  ];
     } else {
       X.h[nf_1] = X.val_aux[na_psids] - X.stprm[nst_psids0];
       X.h[nf_2] = X.val_aux[na_psiqs] - X.stprm[nst_psiqs0];
       X.h[nf_3] = X.val_aux[na_psidr] - X.stprm[nst_psidr0];
       X.h[nf_4] = X.val_aux[na_psiqr] - X.stprm[nst_psiqr0];
       X.h[nf_5] = X.val_vr [nvr_wrm ] - X.stprm[nst_wrm0  ];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = X.val_aux[na_psids] - X.val_aux_u[na_psids];
         X.h[nf_2] = X.val_aux[na_psiqs] - X.val_aux_u[na_psiqs];
         X.h[nf_3] = X.val_aux[na_psidr] - X.val_aux_u[na_psidr];
         X.h[nf_4] = X.val_aux[na_psiqr] - X.val_aux_u[na_psiqr];
         X.h[nf_5] = X.val_vr [nvr_wrm ] - X.val_vr_u [nvr_wrm ];
       } else {
         rs  = X.rprm[nr_rs ];
         lls = X.rprm[nr_lls];
         lm  = X.rprm[nr_lm ];
         rr  = X.rprm[nr_rr ];
         j   = X.rprm[nr_j  ];
         le  = X.rprm[nr_le ];
         l1  = X.rprm[nr_l1 ];
         l2  = X.rprm[nr_l2 ];
         l3  = X.rprm[nr_l3 ];
         x1  = X.rprm[nr_x1 ];
         x2  = X.rprm[nr_x2 ];

         vqs = X.val_vr[nvr_vqs];
         vds = X.val_vr[nvr_vds];
         wrm = X.val_vr[nvr_wrm];
         tl  = X.val_vr[nvr_tl ];

         psids = X.val_aux[na_psids];
         psidr = X.val_aux[na_psidr];
         psiqs = X.val_aux[na_psiqs];
         psiqr = X.val_aux[na_psiqr];

         ids = (l1*psids) - (psidr/le);
         iqs = (l1*psiqs) - (psiqr/le);

         idr = (psids/lm) - (l2*ids);
         iqr = (psiqs/lm) - (l2*iqs);

         tem0 = x1*(iqs*idr-ids*iqr);
         wr   = x2*wrm;

         X.f[nf_1] = vds-rs*ids;
         X.f[nf_2] = vqs-rs*iqs;
         X.f[nf_3] = (-wr)*psiqr-rr*idr;
         X.f[nf_4] = ( wr)*psidr-rr*iqr;
         X.f[nf_5] = (tem0-tl)/j;
       }
     } else {
//     There is some repetition here; never mind.
       rs  = X.rprm[nr_rs ];
       lls = X.rprm[nr_lls];
       lm  = X.rprm[nr_lm ];
       rr  = X.rprm[nr_rr ];
       j   = X.rprm[nr_j  ];
       le  = X.rprm[nr_le ];
       l1  = X.rprm[nr_l1 ];
       l2  = X.rprm[nr_l2 ];
       l3  = X.rprm[nr_l3 ];
       x1  = X.rprm[nr_x1 ];
       x2  = X.rprm[nr_x2 ];

       vqs = X.val_vr[nvr_vqs];
       vds = X.val_vr[nvr_vds];
       wrm = X.val_vr[nvr_wrm];
       tl  = X.val_vr[nvr_tl ];

       psids = X.val_aux[na_psids];
       psidr = X.val_aux[na_psidr];
       psiqs = X.val_aux[na_psiqs];
       psiqr = X.val_aux[na_psiqr];

       if (G.flags[G.i_function] || G.flags[G.i_jacobian]) {
         ids = (l1*psids) - (psidr/le);
         iqs = (l1*psiqs) - (psiqr/le);
         idr = (psids/lm) - (l2*ids);
         iqr = (psiqs/lm) - (l2*iqs);
         tem0 = x1*(iqs*idr-ids*iqr);
         wr = x2*wrm;

         if (G.flags[G.i_function]) {
           X.g[ng_1] = vds-rs*ids;
           X.g[ng_2] = vqs-rs*iqs;
           X.g[ng_3] = (-wr)*psiqr-rr*idr;
           X.g[ng_4] = ( wr)*psidr-rr*iqr;
           X.g[ng_5] = (tem0-tl)/j;
         }
       }
       if (G.flags[G.i_jacobian]) {
//       ids = (l1*psids) - (psidr/le);
         ids_psids = l1;
         ids_psidr = -1.0/le;

//       iqs = (l1*psiqs) - (psiqr/le);
         iqs_psiqs = l1;
         iqs_psiqr = -1.0/le;

//       idr = (psids/lm) - (l2*ids);
         idr_psids = (1.0/lm) - (l2*ids_psids);
         idr_psidr =          - (l2*ids_psidr);

//       iqr = (psiqs/lm) - (l2*iqs);
         iqr_psiqs = (1.0/lm) - (l2*iqs_psiqs);
         iqr_psiqr =          - (l2*iqs_psiqr);

//       tem0 = x1*(iqs*idr-ids*iqr);
         tem0_iqs =  x1*idr;
         tem0_idr =  x1*iqs;
         tem0_ids = -x1*iqr;
         tem0_iqr = -x1*ids;

         tem0_psids =
           tem0_idr*idr_psids +
           tem0_ids*ids_psids;

         tem0_psidr =
           tem0_idr*idr_psidr +
           tem0_ids*ids_psidr;

         tem0_psiqs =
           tem0_iqs*iqs_psiqs +
           tem0_iqr*iqr_psiqs;

         tem0_psiqr =
           tem0_iqs*iqs_psiqr +
           tem0_iqr*iqr_psiqr;

//       wr = x2*wrm;
         wr_wrm = x2;

//       X.g[ng_1] = vds-rs*ids;
         J.dgdvr[ng_1][nvr_vds] = 1.0;
         J.dgdaux[ng_1][na_psids] = -rs*ids_psids;
         J.dgdaux[ng_1][na_psidr] = -rs*ids_psidr;

//       X.g[ng_2] = vqs-rs*iqs;
         J.dgdvr[ng_2][nvr_vqs] = 1.0;
         J.dgdaux[ng_2][na_psiqs] = -rs*iqs_psiqs;
         J.dgdaux[ng_2][na_psiqr] = -rs*iqs_psiqr;

//       X.g[ng_3] = (-wr)*psiqr-rr*idr;
         J.dgdvr[ng_3][nvr_wrm] = (-wr_wrm)*psiqr;
         J.dgdaux[ng_3][na_psiqr] = (-wr);
         J.dgdaux[ng_3][na_psids] = -rr*idr_psids;
         J.dgdaux[ng_3][na_psidr] = -rr*idr_psidr;

//       X.g[ng_4] = ( wr)*psidr-rr*iqr;
         J.dgdvr[ng_4][nvr_wrm] = (wr_wrm)*psidr;
         J.dgdaux[ng_4][na_psidr] = wr;
         J.dgdaux[ng_4][na_psiqs] = -rr*iqr_psiqs;
         J.dgdaux[ng_4][na_psiqr] = -rr*iqr_psiqr;

//       X.g[ng_5] = (tem0-tl)/j;
         J.dgdvr[ng_5][nvr_tl] = -1.0/j;
         J.dgdaux[ng_5][na_psids] = tem0_psids/j;
         J.dgdaux[ng_5][na_psidr] = tem0_psidr/j;
         J.dgdaux[ng_5][na_psiqs] = tem0_psiqs/j;
         J.dgdaux[ng_5][na_psiqr] = tem0_psiqr/j;

//       cout << "J.dgdvr[ng_1][nvr_vds] = " << J.dgdvr[ng_1][nvr_vds] << endl;
//       cout << "J.dgdaux[ng_1][na_psids] = " << J.dgdaux[ng_1][na_psids] << endl;
//       cout << "J.dgdaux[ng_1][na_psidr] = " << J.dgdaux[ng_1][na_psidr] << endl;
//       cout << "J.dgdvr[ng_2][nvr_vqs] = " << J.dgdvr[ng_2][nvr_vqs] << endl;
//       cout << "J.dgdaux[ng_2][na_psiqs] = " << J.dgdaux[ng_2][na_psiqs] << endl;
//       cout << "J.dgdaux[ng_2][na_psiqr] = " << J.dgdaux[ng_2][na_psiqr] << endl;
//       cout << "J.dgdvr[ng_3][nvr_wrm] = " << J.dgdvr[ng_3][nvr_wrm] << endl;
//       cout << "J.dgdaux[ng_3][na_psiqr] = " << J.dgdaux[ng_3][na_psiqr] << endl;
//       cout << "J.dgdaux[ng_3][na_psids] = " << J.dgdaux[ng_3][na_psids] << endl;
//       cout << "J.dgdaux[ng_3][na_psidr] = " << J.dgdaux[ng_3][na_psidr] << endl;
//       cout << "J.dgdvr[ng_4][nvr_wrm] = " << J.dgdvr[ng_4][nvr_wrm] << endl;
//       cout << "J.dgdaux[ng_4][na_psidr] = " << J.dgdaux[ng_4][na_psidr] << endl;
//       cout << "J.dgdaux[ng_4][na_psiqs] = " << J.dgdaux[ng_4][na_psiqs] << endl;
//       cout << "J.dgdaux[ng_4][na_psiqr] = " << J.dgdaux[ng_4][na_psiqr] << endl;
//       cout << "J.dgdvr[ng_5][nvr_tl] = " << J.dgdvr[ng_5][nvr_tl] << endl;
//       cout << "J.dgdaux[ng_5][na_psids] = " << J.dgdaux[ng_5][na_psids] << endl;
//       cout << "J.dgdaux[ng_5][na_psidr] = " << J.dgdaux[ng_5][na_psidr] << endl;
//       cout << "J.dgdaux[ng_5][na_psiqs] = " << J.dgdaux[ng_5][na_psiqs] << endl;
//       cout << "J.dgdaux[ng_5][na_psiqr] = " << J.dgdaux[ng_5][na_psiqr] << endl;
       }
     }
     return;
   }

   return;
}
void x_indmc2a(Global &G,XbeUsr &X,XbeJac &J) {
   double ids,iqs,idr,iqr,wr,tem0;
   double p;
   static double k4=0;
   double ids_psids,ids_psidr,iqs_psiqs,iqs_psiqr;
   double idr_psids,idr_psidr,iqr_psiqs,iqr_psiqr;
   double tem0_iqs,tem0_idr,tem0_ids,tem0_iqr;
   double tem0_psids,tem0_psidr,tem0_psiqs,tem0_psiqr;
   double wr_wrm;
   double vqs,vds,tl;
   double wrm,psids,psidr,psiqs,psiqr;
   int poles;
   double rs,lls,lm,llr,rr,j,ls,lr,le,l1,l2,l3,x1,x2;
   double psids0,psiqs0,psidr0,psiqr0,wrm0;
   const int nvr_vqs = 0;
   const int nvr_vds = 1;
   const int nvr_tl = 2;
   const int nvr_wrm = 3;
   const int nvr_psids = 4;
   const int nvr_psidr = 5;
   const int nvr_psiqs = 6;
   const int nvr_psiqr = 7;
   const int ni_poles = 0;
   const int nr_rs = 0;
   const int nr_lls = 1;
   const int nr_lm = 2;
   const int nr_llr = 3;
   const int nr_rr = 4;
   const int nr_j = 5;
   const int nr_ls = 6;
   const int nr_lr = 7;
   const int nr_le = 8;
   const int nr_l1 = 9;
   const int nr_l2 = 10;
   const int nr_l3 = 11;
   const int nr_x1 = 12;
   const int nr_x2 = 13;
   const int nst_psids0 = 0;
   const int nst_psiqs0 = 1;
   const int nst_psidr0 = 2;
   const int nst_psiqr0 = 3;
   const int nst_wrm0 = 4;
   const int no_wrm = 0;
   const int no_tem = 1;
   const int no_vds = 2;
   const int no_vqs = 3;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nf_4 = 3;
   const int nf_5 = 4;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   const int ng_4 = 3;
   const int ng_5 = 4;
// cout << "indmc1.xbe" << endl;

   if (G.flags[G.i_one_time_parms]) {
     k4 = 0.5*(sqrt(3.0));

     lls = X.rprm[nr_lls];
     lm  = X.rprm[nr_lm ];
     llr = X.rprm[nr_llr];

     ls = lls + lm;
     lr = llr + lm;
     le = (ls*lr/lm) - lm;
     l1 = lr/(lm*le);
     l2 = 1.0 + (lls/lm);
     l3 = lls/lm;

     X.rprm[nr_ls] = ls;
     X.rprm[nr_lr] = lr;
     X.rprm[nr_le] = le;
     X.rprm[nr_l1] = l1;
     X.rprm[nr_l2] = l2;
     X.rprm[nr_l3] = l3;

     poles = X.iprm[ni_poles];
     p = (double)(poles);
     x1 = 0.75*p*lm;
     x2 = 0.5*p;

     X.rprm[nr_x1] = x1;
     X.rprm[nr_x2] = x2;

     return;
   }
   if (G.flags[G.i_init_guess]) {
//   for init guess, use start-up parameters (for this element)
//   can add separate igparms for this purpose later.
     X.val_vr[nvr_psids] = X.stprm[nst_psids0];
     X.val_vr[nvr_psiqs] = X.stprm[nst_psiqs0];
     X.val_vr[nvr_psidr] = X.stprm[nst_psidr0];
     X.val_vr[nvr_psiqr] = X.stprm[nst_psiqr0];
     X.val_vr[nvr_wrm  ] = X.stprm[nst_wrm0  ];
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_wrm] = X.val_vr[nvr_wrm];
     X.outprm[no_vds] = X.val_vr[nvr_vds];
     X.outprm[no_vqs] = X.val_vr[nvr_vqs];

     psids = X.val_vr[nvr_psids];
     psidr = X.val_vr[nvr_psidr];
     psiqs = X.val_vr[nvr_psiqs];
     psiqr = X.val_vr[nvr_psiqr];

     le  = X.rprm[nr_le];
     lm  = X.rprm[nr_lm];
     l1  = X.rprm[nr_l1];
     l2  = X.rprm[nr_l2];

     ids = (l1*psids) - (psidr/le);
     iqs = (l1*psiqs) - (psiqr/le);
     idr = (psids/lm) - (l2*ids);
     iqr = (psiqs/lm) - (l2*iqs);

     x1 = X.rprm[nr_x1];
     tem0 = x1*(iqs*idr-ids*iqr);
     X.outprm[no_tem] = tem0;

     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_psids] = X.stprm[nst_psids0];
       X.val_vr[nvr_psiqs] = X.stprm[nst_psiqs0];
       X.val_vr[nvr_psidr] = X.stprm[nst_psidr0];
       X.val_vr[nvr_psiqr] = X.stprm[nst_psiqr0];
       X.val_vr[nvr_wrm  ] = X.stprm[nst_wrm0  ];
     } else {
       X.h[nf_1] = X.val_vr[nvr_psids] - X.stprm[nst_psids0];
       X.h[nf_2] = X.val_vr[nvr_psiqs] - X.stprm[nst_psiqs0];
       X.h[nf_3] = X.val_vr[nvr_psidr] - X.stprm[nst_psidr0];
       X.h[nf_4] = X.val_vr[nvr_psiqr] - X.stprm[nst_psiqr0];
       X.h[nf_5] = X.val_vr[nvr_wrm  ] - X.stprm[nst_wrm0  ];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = X.val_vr[nvr_psids] - X.val_vr_u[nvr_psids];
         X.h[nf_2] = X.val_vr[nvr_psiqs] - X.val_vr_u[nvr_psiqs];
         X.h[nf_3] = X.val_vr[nvr_psidr] - X.val_vr_u[nvr_psidr];
         X.h[nf_4] = X.val_vr[nvr_psiqr] - X.val_vr_u[nvr_psiqr];
         X.h[nf_5] = X.val_vr[nvr_wrm  ] - X.val_vr_u[nvr_wrm  ];
       } else {
         rs  = X.rprm[nr_rs ];
         lls = X.rprm[nr_lls];
         lm  = X.rprm[nr_lm ];
         rr  = X.rprm[nr_rr ];
         j   = X.rprm[nr_j  ];
         le  = X.rprm[nr_le ];
         l1  = X.rprm[nr_l1 ];
         l2  = X.rprm[nr_l2 ];
         l3  = X.rprm[nr_l3 ];
         x1  = X.rprm[nr_x1 ];
         x2  = X.rprm[nr_x2 ];

         vqs = X.val_vr[nvr_vqs];
         vds = X.val_vr[nvr_vds];
         wrm = X.val_vr[nvr_wrm];
         tl  = X.val_vr[nvr_tl ];

         psids = X.val_vr[nvr_psids];
         psidr = X.val_vr[nvr_psidr];
         psiqs = X.val_vr[nvr_psiqs];
         psiqr = X.val_vr[nvr_psiqr];

         ids = (l1*psids) - (psidr/le);
         iqs = (l1*psiqs) - (psiqr/le);

         idr = (psids/lm) - (l2*ids);
         iqr = (psiqs/lm) - (l2*iqs);

         tem0 = x1*(iqs*idr-ids*iqr);
         wr   = x2*wrm;

         X.f[nf_1] = vds-rs*ids;
         X.f[nf_2] = vqs-rs*iqs;
         X.f[nf_3] = (-wr)*psiqr-rr*idr;
         X.f[nf_4] = ( wr)*psidr-rr*iqr;
         X.f[nf_5] = (tem0-tl)/j;
       }
     } else {
//     There is some repetition here; never mind.
       rs  = X.rprm[nr_rs ];
       lls = X.rprm[nr_lls];
       lm  = X.rprm[nr_lm ];
       rr  = X.rprm[nr_rr ];
       j   = X.rprm[nr_j  ];
       le  = X.rprm[nr_le ];
       l1  = X.rprm[nr_l1 ];
       l2  = X.rprm[nr_l2 ];
       l3  = X.rprm[nr_l3 ];
       x1  = X.rprm[nr_x1 ];
       x2  = X.rprm[nr_x2 ];

       vqs = X.val_vr[nvr_vqs];
       vds = X.val_vr[nvr_vds];
       wrm = X.val_vr[nvr_wrm];
       tl  = X.val_vr[nvr_tl ];

       psids = X.val_vr[nvr_psids];
       psidr = X.val_vr[nvr_psidr];
       psiqs = X.val_vr[nvr_psiqs];
       psiqr = X.val_vr[nvr_psiqr];

       if (G.flags[G.i_function] || G.flags[G.i_jacobian]) {
         ids = (l1*psids) - (psidr/le);
         iqs = (l1*psiqs) - (psiqr/le);
         idr = (psids/lm) - (l2*ids);
         iqr = (psiqs/lm) - (l2*iqs);
         tem0 = x1*(iqs*idr-ids*iqr);
         wr = x2*wrm;

         if (G.flags[G.i_function]) {
           X.g[ng_1] = vds-rs*ids;
           X.g[ng_2] = vqs-rs*iqs;
           X.g[ng_3] = (-wr)*psiqr-rr*idr;
           X.g[ng_4] = ( wr)*psidr-rr*iqr;
           X.g[ng_5] = (tem0-tl)/j;
         }
       }
       if (G.flags[G.i_jacobian]) {
//       ids = (l1*psids) - (psidr/le);
         ids_psids = l1;
         ids_psidr = -1.0/le;

//       iqs = (l1*psiqs) - (psiqr/le);
         iqs_psiqs = l1;
         iqs_psiqr = -1.0/le;

//       idr = (psids/lm) - (l2*ids);
         idr_psids = (1.0/lm) - (l2*ids_psids);
         idr_psidr =          - (l2*ids_psidr);

//       iqr = (psiqs/lm) - (l2*iqs);
         iqr_psiqs = (1.0/lm) - (l2*iqs_psiqs);
         iqr_psiqr =          - (l2*iqs_psiqr);

//       tem0 = x1*(iqs*idr-ids*iqr);
         tem0_iqs =  x1*idr;
         tem0_idr =  x1*iqs;
         tem0_ids = -x1*iqr;
         tem0_iqr = -x1*ids;

         tem0_psids =
           tem0_idr*idr_psids +
           tem0_ids*ids_psids;

         tem0_psidr =
           tem0_idr*idr_psidr +
           tem0_ids*ids_psidr;

         tem0_psiqs =
           tem0_iqs*iqs_psiqs +
           tem0_iqr*iqr_psiqs;

         tem0_psiqr =
           tem0_iqs*iqs_psiqr +
           tem0_iqr*iqr_psiqr;

//       wr = x2*wrm;
         wr_wrm = x2;

//       X.g[ng_1] = vds-rs*ids;
         J.dgdvr[ng_1][nvr_vds] = 1.0;
         J.dgdvr[ng_1][nvr_psids] = -rs*ids_psids;
         J.dgdvr[ng_1][nvr_psidr] = -rs*ids_psidr;

//       X.g[ng_2] = vqs-rs*iqs;
         J.dgdvr[ng_2][nvr_vqs] = 1.0;
         J.dgdvr[ng_2][nvr_psiqs] = -rs*iqs_psiqs;
         J.dgdvr[ng_2][nvr_psiqr] = -rs*iqs_psiqr;

//       X.g[ng_3] = (-wr)*psiqr-rr*idr;
         J.dgdvr[ng_3][nvr_wrm] = (-wr_wrm)*psiqr;
         J.dgdvr[ng_3][nvr_psiqr] = (-wr);
         J.dgdvr[ng_3][nvr_psids] = -rr*idr_psids;
         J.dgdvr[ng_3][nvr_psidr] = -rr*idr_psidr;

//       X.g[ng_4] = ( wr)*psidr-rr*iqr;
         J.dgdvr[ng_4][nvr_wrm] = (wr_wrm)*psidr;
         J.dgdvr[ng_4][nvr_psidr] = wr;
         J.dgdvr[ng_4][nvr_psiqs] = -rr*iqr_psiqs;
         J.dgdvr[ng_4][nvr_psiqr] = -rr*iqr_psiqr;

//       X.g[ng_5] = (tem0-tl)/j;
         J.dgdvr[ng_5][nvr_tl] = -1.0/j;
         J.dgdvr[ng_5][nvr_psids] = tem0_psids/j;
         J.dgdvr[ng_5][nvr_psidr] = tem0_psidr/j;
         J.dgdvr[ng_5][nvr_psiqs] = tem0_psiqs/j;
         J.dgdvr[ng_5][nvr_psiqr] = tem0_psiqr/j;
       }
     }
     return;
   }

   return;
}
void x_indmc2b(Global &G,XbeUsr &X,XbeJac &J) {
   double ids,iqs,idr,iqr,wr;
   double ids_psids,ids_psidr,iqs_psiqs,iqs_psiqr;
   static double k4=0;
   double psids,psidr,psiqs,psiqr;
   double ia,ib,ic;
   double lls,lm,llr,ls,lr,le,l1;
   const int nvr_psids = 0;
   const int nvr_psidr = 1;
   const int nvr_psiqs = 2;
   const int nvr_psiqr = 3;
   const int nvr_ia = 4;
   const int nvr_ib = 5;
   const int nvr_ic = 6;
   const int nr_lls = 0;
   const int nr_lm = 1;
   const int nr_llr = 2;
   const int nr_ls = 3;
   const int nr_lr = 4;
   const int nr_le = 5;
   const int nr_l1 = 6;
   const int no_ia = 0;
   const int no_ib = 1;
   const int no_ic = 2;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   if (G.flags[G.i_one_time_parms]) {
     k4 = 0.5*(sqrt(3.0));

     lls = X.rprm[nr_lls];
     lm  = X.rprm[nr_lm ];
     llr = X.rprm[nr_llr];

     ls = lls + lm;
     lr = llr + lm;
     le = (ls*lr/lm) - lm;
     l1 = lr/(lm*le);

     X.rprm[nr_ls] = ls;
     X.rprm[nr_lr] = lr;
     X.rprm[nr_le] = le;
     X.rprm[nr_l1] = l1;

     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_ia] = 0.0;
     X.val_vr[nvr_ib] = 0.0;
     X.val_vr[nvr_ic] = 0.0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_ia] = X.val_vr[nvr_ia];
     X.outprm[no_ib] = X.val_vr[nvr_ib];
     X.outprm[no_ic] = X.val_vr[nvr_ic];
     return;
   }

   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       le  = X.rprm[nr_le];
       l1  = X.rprm[nr_l1];

       psids = X.val_vr[nvr_psids];
       psidr = X.val_vr[nvr_psidr];
       psiqs = X.val_vr[nvr_psiqs];
       psiqr = X.val_vr[nvr_psiqr];

       ids = (l1*psids) - (psidr/le);
       iqs = (l1*psiqs) - (psiqr/le);

       X.val_vr[nvr_ia] = iqs;
       X.val_vr[nvr_ib] = -0.5*iqs-k4*ids;
       X.val_vr[nvr_ic] = -0.5*iqs+k4*ids;
     } else if (G.flags[G.i_implicit]) {
       le  = X.rprm[nr_le ];
       l1  = X.rprm[nr_l1 ];

       ia = X.val_vr[nvr_ia];
       ib = X.val_vr[nvr_ib];
       ic = X.val_vr[nvr_ic];

       psids = X.val_vr[nvr_psids];
       psidr = X.val_vr[nvr_psidr];
       psiqs = X.val_vr[nvr_psiqs];
       psiqr = X.val_vr[nvr_psiqr];

       ids = (l1*psids) - (psidr/le);
       ids_psids = l1;
       ids_psidr = -1.0/le;

       iqs = (l1*psiqs) - (psiqr/le);
       iqs_psiqs = l1;
       iqs_psiqr = -1.0/le;

       if (G.flags[G.i_function]) {
         X.g[ng_1] = ia - iqs;
         X.g[ng_2] = ib -(-0.5*iqs-k4*ids);
         X.g[ng_3] = ic -(-0.5*iqs+k4*ids);
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_ia] = 1.0;
         J.dgdvr[ng_1][nvr_psiqs] = -iqs_psiqs;
         J.dgdvr[ng_1][nvr_psiqr] = -iqs_psiqr;

         J.dgdvr[ng_2][nvr_ib] = 1.0;
         J.dgdvr[ng_2][nvr_psids] = k4*ids_psids;
         J.dgdvr[ng_2][nvr_psidr] = k4*ids_psidr;
         J.dgdvr[ng_2][nvr_psiqs] = 0.5*iqs_psiqs;
         J.dgdvr[ng_2][nvr_psiqr] = 0.5*iqs_psiqr;

         J.dgdvr[ng_3][nvr_ic] = 1.0;
         J.dgdvr[ng_3][nvr_psids] = -k4*ids_psids;
         J.dgdvr[ng_3][nvr_psidr] = -k4*ids_psidr;
         J.dgdvr[ng_3][nvr_psiqs] = 0.5*iqs_psiqs;
         J.dgdvr[ng_3][nvr_psiqr] = 0.5*iqs_psiqr;
       }
     }
     return;
   }
   return;
}
void x_integrator(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double y;
   double k;
   double y_st;
   double y_ig;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_k = 0;
   const int nst_y_st = 0;
   const int nig_y_ig = 0;
   const int no_x = 0;
   const int no_y = 1;
   const int nf_1 = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr [nvr_y] = X.igprm[nig_y_ig];
     return;
   }
   if (G.flags[G.i_outvar]) {
//   cout << "integrator: outvar" << endl;
//   cout << "integrator: X.val_vr[nvr_x] = " << X.val_vr[nvr_x] << endl;
//   cout << "integrator: X.val_vr[nvr_y] = " << X.val_vr[nvr_y] << endl;
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
//     cout << "integrator: X.val_vr[nvr_y] = " << X.val_vr[nvr_y] << endl;
     } else if (G.flags[G.i_implicit]) {
//     cout << "integrator: startup, implicit" << endl;
//     cout << "integrator: X.val_vr[nvr_x] = " << X.val_vr[nvr_x] << endl;
//     cout << "integrator: X.val_vr[nvr_y] = " << X.val_vr[nvr_y] << endl;
       X.h[nf_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
//       use h[..] for both startup and explicit/alg loop:
         X.h[nf_1] = X.val_vr[nvr_y] - X.val_vr_u[nvr_y];
       } else {
         k = X.rprm[nr_k];
         x = X.val_vr[nvr_x];
         X.f[nf_1] = k*x;
//       cout << "integrator: x = " << x << " k = " << k << endl;
       }
     } else if (G.flags[G.i_implicit]) {
       k = X.rprm[nr_k];
       x = X.val_vr[nvr_x];
       if (G.flags[G.i_function]) {
         X.g[ng_1] = k*x;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_x] = k;
       }
     }
     return;
   }
   return;
}
void x_integrator_1(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double y;
   double k,tau,ki;
   double y_st;
   double y_ig;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_k = 0;
   const int nr_tau = 1;
   const int nr_ki = 2;
   const int nst_y_st = 0;
   const int nig_y_ig = 0;
   const int no_x = 0;
   const int no_y = 1;
   const int nf_1 = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
       k = X.rprm[nr_k];
       tau = X.rprm[nr_tau];
       ki = k/tau;
       X.rprm[nr_ki] = ki;
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr [nvr_y] = X.igprm[nig_y_ig];
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
     } else if (G.flags[G.i_implicit]) {
       X.h[nf_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = X.val_vr[nvr_y] - X.val_vr_u[nvr_y];
       } else {
         ki = X.rprm[nr_ki];
         x = X.val_vr[nvr_x];
         X.f[nf_1] = ki*x;
       }
     } else if (G.flags[G.i_implicit]) {
       ki = X.rprm[nr_ki];
       x = X.val_vr[nvr_x];
       if (G.flags[G.i_function]) {
         X.g[ng_1] = ki*x;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_x] = ki;
       }
     }
     return;
   }
   return;
}
void x_integrator_reset_1(Global &G,XbeUsr &X,XbeJac &J) {
   int flag_active_edge;
   double time0;
   double x,r;
   double y;
   int active_pos_edge,active_neg_edge;
   double k,r_high,delt_min,y_reset,r_prev,r_cross;
   double y_st;
   double y_ig;
   const int nvr_x = 0;
   const int nvr_r = 1;
   const int nvr_y = 2;
   const int ni_active_pos_edge = 0;
   const int ni_active_neg_edge = 1;
   const int nr_k = 0;
   const int nr_r_high = 1;
   const int nr_delt_min = 2;
   const int nr_y_reset = 3;
   const int nr_r_prev = 4;
   const int nr_r_cross = 5;
   const int nst_y_st = 0;
   const int nig_y_ig = 0;
   const int no_x = 0;
   const int no_y = 1;
   const int no_r = 2;
   const int nf_1 = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     active_pos_edge = X.iprm[ni_active_pos_edge];
     active_neg_edge = X.iprm[ni_active_neg_edge];

     if (active_pos_edge == 0) {
       if (active_neg_edge == 0) {
         cout << "integrator_reset_1.xbe: one of pos/neg edge must be 1. Halting..." << endl;
         exit(1);
       }
     } else {
       if (active_neg_edge != 0) {
         cout << "integrator_reset_1.xbe: only one of pos/neg edge must be 1. Halting..." << endl;
         exit(1);
       }
     }
     r_high = X.rprm[nr_r_high];

     if (r_high < 0.0) {
       cout << "integrator_reset_1.xbe: check r_high. Halting..." << endl;
       exit(1);
     }
     r_cross = 0.5*r_high;
     X.rprm[nr_r_cross] = r_cross;
     X.rprm[nr_r_prev] = 0.0;
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr [nvr_y] = X.igprm[nig_y_ig];
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     X.outprm[no_r] = X.val_vr[nvr_r];
     return;
   }
   if (G.flags[G.i_reset_x]) {
     time0 = G.time_given_x;
     active_pos_edge = X.iprm[ni_active_pos_edge];

     r_prev  = X.rprm[nr_r_prev ];
     r_cross = X.rprm[nr_r_cross];

     r = X.val_vr[nvr_r];

     if (active_pos_edge == 1) {
       if ((r_prev <= r_cross) && (r >= r_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     } else {
       if ((r_prev >= r_cross) && (r <= r_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     }
     if (flag_active_edge == 1) {
       X.val_vr[nvr_y] = X.rprm[nr_y_reset];
     }
     return;
   }
   if (G.flags[G.i_save_history]) {
     X.rprm[nr_r_prev] = X.val_vr[nvr_r];
     return;
   }
   if (G.flags[G.i_next_time]) {
     time0 = G.time_given_x;

     active_pos_edge = X.iprm[ni_active_pos_edge];
     r_prev  = X.rprm[nr_r_prev ];
     r_cross = X.rprm[nr_r_cross];
     r = X.val_vr[nvr_r];

     if (active_pos_edge == 1) {
       if ((r_prev <= r_cross) && (r >= r_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     } else {
       if ((r_prev >= r_cross) && (r <= r_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     }
     if (flag_active_edge == 1) {
       G.time_nextbreak_x = time0 + X.rprm[nr_delt_min];
     } else {
       G.time_nextbreak_x = G.time_end;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
     } else if (G.flags[G.i_implicit]) {
       X.h[nf_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = X.val_vr[nvr_y] - X.val_vr_u[nvr_y];
       } else {
         k = X.rprm[nr_k];
         x = X.val_vr[nvr_x];
         X.f[nf_1] = k*x;
       }
     } else if (G.flags[G.i_implicit]) {
       k = X.rprm[nr_k];
       x = X.val_vr[nvr_x];
       if (G.flags[G.i_function]) {
         X.g[ng_1] = k*x;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_x] = k;
       }
     }
     return;
   }
   return;
}
void x_jkff(Global &G,XbeUsr &X,XbeJac &J) {
   int flag_active_edge;
   double time0,q0,q0bar;
   double clk,j,k;
   double q,qbar;
   int active_pos_edge,active_neg_edge;
   double x_high,dt,clk_prev,j_prev,k_prev,q_next,x_cross;
   double q_st;
   const int nvr_clk = 0;
   const int nvr_j = 1;
   const int nvr_k = 2;
   const int nvr_q = 3;
   const int nvr_qbar = 4;
   const int ni_active_pos_edge = 0;
   const int ni_active_neg_edge = 1;
   const int nr_x_high = 0;
   const int nr_dt = 1;
   const int nr_clk_prev = 2;
   const int nr_j_prev = 3;
   const int nr_k_prev = 4;
   const int nr_q_next = 5;
   const int nr_x_cross = 6;
   const int nst_q_st = 0;
   const int no_clk = 0;
   const int no_j = 1;
   const int no_k = 2;
   const int no_q = 3;
   const int ng_1 = 0;
   const int ng_2 = 1;
   if (G.flags[G.i_one_time_parms]) {
     active_pos_edge = X.iprm[ni_active_pos_edge];
     active_neg_edge = X.iprm[ni_active_neg_edge];

     if (active_pos_edge == 0) {
       if (active_neg_edge == 0) {
         cout << "jkff.xbe: one of pos/neg edge must be 1. Halting..." << endl;
         exit(1);
       }
     } else {
       if (active_neg_edge != 0) {
         cout << "jkff.xbe: only one of pos/neg edge must be 1. Halting..." << endl;
         exit(1);
       }
     }
     x_high = X.rprm[nr_x_high];

     if (x_high < 0.0) {
       cout << "jkff.xbe: check x_high. Halting..." << endl;
       exit(1);
     }
     x_cross = 0.5*x_high;
     X.rprm[nr_x_cross] = x_cross;

     X.rprm[nr_clk_prev] = 0.0;
     X.rprm[nr_j_prev] = 0.0;
     X.rprm[nr_k_prev] = 0.0;
     X.rprm[nr_q_next] = 0.0;
     return;
   }

   if (G.flags[G.i_outvar]) {
     X.outprm[no_clk] = X.val_vr[nvr_clk];
     X.outprm[no_j  ] = X.val_vr[nvr_j  ];
     X.outprm[no_k  ] = X.val_vr[nvr_k  ];
     X.outprm[no_q  ] = X.val_vr[nvr_q  ];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_q   ] = 0.0;
     X.val_vr[nvr_qbar] = X.rprm[nr_x_high];
     return;
   }
   if (G.flags[G.i_save_history]) {
     time0 = G.time_given_x;
     active_pos_edge = X.iprm[ni_active_pos_edge];

     clk_prev = X.rprm[nr_clk_prev];
     x_cross  = X.rprm[nr_x_cross ];
     x_high   = X.rprm[nr_x_high  ];

     clk = X.val_vr[nvr_clk];
     j   = X.val_vr[nvr_j  ];
     k   = X.val_vr[nvr_k  ];
     q   = X.val_vr[nvr_q  ];

     if (active_pos_edge == 1) {
       if ((clk_prev <= x_cross) && (clk >= x_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     } else {
       if ((clk_prev >= x_cross) && (clk <= x_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     }
     if (flag_active_edge == 1) {
       j_prev = X.rprm[nr_j_prev];
       k_prev = X.rprm[nr_k_prev];

       if (j_prev > x_cross) {
         if (k_prev > x_cross) {
           if (q > x_cross) {
             q_next = 0.0;
           } else {
             q_next = x_high;
           }
         } else {
           q_next = x_high;
         }
       } else {
         if (k_prev > x_cross) {
           q_next = 0.0;
         } else {
           q_next = q;
         }
       }
     } else {
       q_next = q;
     }

     X.rprm[nr_clk_prev] = clk;
     X.rprm[nr_j_prev  ] = j;
     X.rprm[nr_k_prev  ] = k;
     X.rprm[nr_q_next  ] = q_next;

     return;
   }
   if (G.flags[G.i_next_time]) {
     time0 = G.time_given_x;

     active_pos_edge = X.iprm[ni_active_pos_edge];

     clk_prev = X.rprm[nr_clk_prev];
     x_cross  = X.rprm[nr_x_cross ];

     clk = X.val_vr[nvr_clk];

     if (active_pos_edge == 1) {
       if ((clk_prev <= x_cross) && (clk >= x_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     } else {
       if ((clk_prev >= x_cross) && (clk <= x_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     }
     if (flag_active_edge == 1) {
       G.time_nextbreak_x = time0 + X.rprm[nr_dt];
     } else {
       G.time_nextbreak_x = G.time_end;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       q_st = X.stprm[nst_q_st];

       x_cross = X.rprm[nr_x_cross];
       x_high  = X.rprm[nr_x_high ];

       if (q_st > x_cross) {
         X.val_vr[nvr_q   ] = x_high;
         X.val_vr[nvr_qbar] = 0.0;
       } else {
         X.val_vr[nvr_q   ] = 0.0;
         X.val_vr[nvr_qbar] = x_high;
       }
     } else if (G.flags[G.i_implicit]) {
       q_st = X.stprm[nst_q_st];

       x_cross = X.rprm[nr_x_cross];
       x_high  = X.rprm[nr_x_high ];

       if (q_st > x_cross) {
         q0    = x_high;
         q0bar = 0.0;
       } else {
         q0    = 0.0;
         q0bar = x_high;
       }
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_q   ] - q0;
         X.g[ng_2] = X.val_vr[nvr_qbar] - q0bar;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_q   ] = 1.0;
         J.dgdvr[ng_2][nvr_qbar] = 1.0;
       }
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     time0 = G.time_given_x;

     if (G.iter_trns_x == 0) {
       if (G.flags[G.i_slv_init]) {
         q0 = X.rprm[nr_q_next];
       } else {
         q0 = X.val_vr[nvr_q];
       }
     } else {
       q0 = X.rprm[nr_q_next];
     }

     x_cross = X.rprm[nr_x_cross];
     x_high  = X.rprm[nr_x_high ];

     if (q0 > x_cross) {
       q0bar = 0.0;
     } else {
       q0bar = x_high;
     }
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_q   ] = q0;
       X.val_vr[nvr_qbar] = q0bar;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_q   ] - q0;
         X.g[ng_2] = X.val_vr[nvr_qbar] - q0bar;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_q   ] = 1.0;
         J.dgdvr[ng_2][nvr_qbar] = 1.0;
       }
     }
   }
   return;
}
void x_lag_1(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   double x;
   double y;
   double tr,k;
   double y_st;
   double y_ig;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_tr = 0;
   const int nr_k = 1;
   const int nst_y_st = 0;
   const int nig_y_ig = 0;
   const int no_x = 0;
   const int no_y = 1;
   const int nf_1 = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     tr = X.rprm[nr_tr];
     k = 1.0/tr;
     X.rprm[nr_k] = k;
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr [nvr_y] = X.igprm[nig_y_ig];
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
     } else if (G.flags[G.i_implicit]) {
       X.h[nf_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = X.val_vr[nvr_y] - X.val_vr_u[nvr_y];
       } else {
         k = X.rprm[nr_k];
         x = X.val_vr[nvr_x];
         y = X.val_vr[nvr_y];
         X.f[nf_1] = k*(-y+x);
       }
     } else if (G.flags[G.i_implicit]) {
       k = X.rprm[nr_k];
       x = X.val_vr[nvr_x];
       y = X.val_vr[nvr_y];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = k*(-y+x);
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_x] =  k;
         J.dgdvr[ng_1][nvr_y] = -k;
       }
     }
     return;
   }
   return;
}
void x_lag_2(Global &G,XbeUsr &X,XbeJac &J) {
   double tr,y0;
   double x;
   double y;
   double t_delay,k;
   double y_st;
   double y_ig;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_t_delay = 0;
   const int nr_k = 1;
   const int nst_y_st = 0;
   const int nig_y_ig = 0;
   const int no_x = 0;
   const int no_y = 1;
   const int nf_1 = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     tr = X.rprm[nr_t_delay]/log(2.0);
     k = 1.0/tr;
     X.rprm[nr_k] = k;
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr [nvr_y] = X.igprm[nig_y_ig];
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
     } else if (G.flags[G.i_implicit]) {
       X.h[nf_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = X.val_vr[nvr_y] - X.val_vr_u[nvr_y];
       } else {
         k = X.rprm[nr_k];
         x = X.val_vr[nvr_x];
         y = X.val_vr[nvr_y];
         X.f[nf_1] = k*(-y+x);
       }
     } else if (G.flags[G.i_implicit]) {
       k = X.rprm[nr_k];
       x = X.val_vr[nvr_x];
       y = X.val_vr[nvr_y];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = k*(-y+x);
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_x] =  k;
         J.dgdvr[ng_1][nvr_y] = -k;
       }
     }
     return;
   }
   return;
}
void x_limiter(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   double x;
   double y;
   double xmin,xmax;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_xmin = 0;
   const int nr_xmax = 1;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   xmin = X.rprm[nr_xmin];
   xmax = X.rprm[nr_xmax];

   x = X.val_vr[nvr_x];
   if (x >= xmax) {
     y0 = xmax;
   } else if (x <= xmin) {
     y0 = xmin;
   } else {
     y0 = x;
   }   

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       y = X.val_vr[nvr_y];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = y - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] =  1.0;
         if (x >= xmax) {
           J.dgdvr[ng_1][nvr_x] = 0.0;
         } else if (x <= xmin) {
           J.dgdvr[ng_1][nvr_x] = 0.0;
         } else {
           J.dgdvr[ng_1][nvr_x] = -1.0;
         }   
       }
     }
     return;
   }
   return;
}
void x_linear(Global &G,XbeUsr &X,XbeJac &J) {
   double eps;
   double x;
   double y;
   double a1,a2,b1,b2,c1,c2,a,b,c;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_a1 = 0;
   const int nr_a2 = 1;
   const int nr_b1 = 2;
   const int nr_b2 = 3;
   const int nr_c1 = 4;
   const int nr_c2 = 5;
   const int nr_a = 6;
   const int nr_b = 7;
   const int nr_c = 8;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     a1 = X.rprm[nr_a1];
     a2 = X.rprm[nr_a2];
     b1 = X.rprm[nr_b1];
     b2 = X.rprm[nr_b2];
     c1 = X.rprm[nr_c1];
     c2 = X.rprm[nr_c2];

     eps = 1.0e-10;

     if (abs(a2) < eps) {
       cout << "linear.xbe: a2 is too small. a2=" << a2 << endl;
       cout << "  halting..." << endl; exit(1);
     }
     if (abs(b2) < eps) {
       cout << "linear.xbe: b2 is too small. b2=" << b2 << endl;
       cout << "  halting..." << endl; exit(1);
     }
     if (abs(c2) < eps) {
       cout << "linear.xbe: c2 is too small. c2=" << c2 << endl;
       cout << "  halting..." << endl; exit(1);
     }
     a = a1/a2;
     b = b1/b2;
     c = c1/c2;

     if (abs(a) < eps) {
       cout << "linear.xbe: a is too small. a=" << a << endl;
       cout << "  halting..." << endl; exit(1);
     }
     X.rprm[nr_a] = a;
     X.rprm[nr_b] = b;
     X.rprm[nr_c] = c;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   a = X.rprm[nr_a];
   b = X.rprm[nr_b];
   c = X.rprm[nr_c];

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = (b/a)*X.val_vr[nvr_x] + (c/a);
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = (b/a)*X.val_vr[nvr_x] + (c/a);
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = a*X.val_vr[nvr_y] - b*X.val_vr[nvr_x] - c;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] =  a;
         J.dgdvr[ng_1][nvr_x] = -b;
       }
     }
     return;
   }
   return;
}
void x_max(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   double x1,x2;
   double y;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];
   y0 = max(x1,x2);

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_min(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   double x1,x2;
   double y;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];
   y0 = min(x1,x2);

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_modulo(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double x1,x2;
   const int nvr_x = 0;
   const int nr_x1 = 0;
   const int nr_x2 = 1;
   const int no_x = 0;
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     return;
   }
   x1 = X.rprm[nr_x1];
   x2 = X.rprm[nr_x2];

   if (G.flags[G.i_init_guess]) {
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     return;
   }
   return;
}
void x_modulo_twopi(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double x1,x2;
   const int nvr_x = 0;
   const int nr_x1 = 0;
   const int nr_x2 = 1;
   const int no_x = 0;
   if (G.flags[G.i_one_time_parms]) {
     X.rprm[nr_x1] = -G.pi;
     X.rprm[nr_x2] =  G.pi;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     return;
   }
   x1 = X.rprm[nr_x1];
   x2 = X.rprm[nr_x2];

   if (G.flags[G.i_init_guess]) {
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     return;
   }
   return;
}
void x_monostable_1(Global &G,XbeUsr &X,XbeJac &J) {
   int flag_active_edge;
   double time0,y0;
   double x;
   double y;
   int active_pos_edge,active_neg_edge;
   double x_low,x_high,y_low,y_high,T,x_prev,t2,x_cross,epsl,epsl1,y_half;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int ni_active_pos_edge = 0;
   const int ni_active_neg_edge = 1;
   const int nr_x_low = 0;
   const int nr_x_high = 1;
   const int nr_y_low = 2;
   const int nr_y_high = 3;
   const int nr_T = 4;
   const int nr_x_prev = 5;
   const int nr_t2 = 6;
   const int nr_x_cross = 7;
   const int nr_epsl = 8;
   const int nr_epsl1 = 9;
   const int nr_y_half = 10;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     active_pos_edge = X.iprm[ni_active_pos_edge];
     active_neg_edge = X.iprm[ni_active_neg_edge];

     if (active_pos_edge == 0) {
       if (active_neg_edge == 0) {
         cout << "monostable_1.xbe: one of pos/neg edge must be 1. Halting..." << endl;
         exit(1);
       }
     } else {
       if (active_neg_edge != 0) {
         cout << "monostable_1.xbe: only one of pos/neg edge must be 1. Halting..." << endl;
         exit(1);
       }
     }
     x_low  = X.rprm[nr_x_low ];
     x_high = X.rprm[nr_x_high];
     y_low  = X.rprm[nr_y_low ];
     y_high = X.rprm[nr_y_high];
     T      = X.rprm[nr_T     ];

     if (x_low > x_high) {
       cout << "monostable_1.xbe: check x_low/x_high. Halting..." << endl;
       exit(1);
     }
     if (y_low > y_high) {
       cout << "monostable_1.xbe: check y_low/y_high. Halting..." << endl;
       exit(1);
     }

     x_cross = x_low + 0.5*(x_high-x_low);
     X.rprm[nr_x_cross] = x_cross;

     y_half = y_low + 0.5*(y_high-y_low);

     epsl = T/100.0;
     epsl1 = 1.1*epsl;

     X.rprm[nr_epsl  ] = epsl;
     X.rprm[nr_epsl1 ] = epsl1;
     X.rprm[nr_y_half] = y_half;

     X.rprm[nr_x_prev] = 0.0;
     X.rprm[nr_t2    ] = G.time_begin;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_save_history]) {
     time0 = G.time_given_x;
     active_pos_edge = X.iprm[ni_active_pos_edge];

     x_prev  = X.rprm[nr_x_prev ];
     x_cross = X.rprm[nr_x_cross];

     x = X.val_vr[nvr_x];

     if (active_pos_edge == 1) {
       if ((x_prev <= x_cross) && (x >= x_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     } else {
       if ((x_prev >= x_cross) && (x <= x_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     }

     if (flag_active_edge == 1) {
       if (X.val_vr[nvr_y] < X.rprm[nr_y_half]) {
         X.rprm[nr_t2] = time0 + X.rprm[nr_T];
       }
     }

     X.rprm[nr_x_prev] = x;

     return;
   }
   if (G.flags[G.i_next_time]) {
     time0 = G.time_given_x;
     t2 = X.rprm[nr_t2];

     active_pos_edge = X.iprm[ni_active_pos_edge];
     x_prev  = X.rprm[nr_x_prev ];
     x_cross = X.rprm[nr_x_cross];
     x = X.val_vr[nvr_x];

     if (active_pos_edge == 1) {
       if ((x_prev <= x_cross) && (x >= x_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     } else {
       if ((x_prev >= x_cross) && (x <= x_cross)) {
         flag_active_edge = 1;
       } else {
         flag_active_edge = 0;
       }
     }
     if (flag_active_edge == 1) {
       G.time_nextbreak_x = time0 + X.rprm[nr_epsl];
       return;
     }

     if (time0 < t2) {
       epsl  = X.rprm[nr_epsl ];
       epsl1 = X.rprm[nr_epsl1];

       if ((t2-time0) <= epsl1) {
         G.time_nextbreak_x = t2 + epsl;
       } else {
         G.time_nextbreak_x = t2 - epsl;
       }
     } else {
       G.time_nextbreak_x = G.time_end;
     }
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = 0.0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     t2 = X.rprm[nr_t2];
     time0 = G.time_given_x;

     if (time0 < t2) {
       y0 = X.rprm[nr_y_high];
     } else {
       y0 = X.rprm[nr_y_low];
     }
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
   }
   return;
}
void x_mult_2(Global &G,XbeUsr &X,XbeJac &J) {
   double x1,x2;
   double y;
   double k;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int nr_k = 0;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   k = X.rprm[nr_k];
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = k*X.val_vr[nvr_x1]*X.val_vr[nvr_x2];
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = k*X.val_vr[nvr_x1]*X.val_vr[nvr_x2];
     } else if (G.flags[G.i_implicit]) {
       y = X.val_vr[nvr_y];
       x1 = X.val_vr[nvr_x1];
       x2 = X.val_vr[nvr_x2];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = y - k*x1*x2;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
         J.dgdvr[ng_1][nvr_x1] = -k*x2;
         J.dgdvr[ng_1][nvr_x2] = -k*x1;
       }
     }
     return;
   }
   return;
}
void x_multscl(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double y;
   double k;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_k = 0;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   k = X.rprm[nr_k];
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = k*X.val_vr[nvr_x];
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = k*X.val_vr[nvr_x];
//     cout << "multscl.xbe: x = " << X.val_vr[nvr_x]
//      << " y = " << X.val_vr[nvr_y] << endl;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - k*X.val_vr[nvr_x];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
         J.dgdvr[ng_1][nvr_x] = -k;
       }
     }
     return;
   }
   return;
}
void x_nand_2(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   bool X1,X2,Y;
   double x1,x2;
   double y;
   double y_high,hb2;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int nr_y_high = 0;
   const int nr_hb2 = 1;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     y_high = X.rprm[nr_y_high];
     hb2 = 0.5*y_high;
     X.rprm[nr_hb2] = hb2;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   y_high = X.rprm[nr_y_high];
   hb2 = X.rprm[nr_hb2];

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];
   X1 = x1 > hb2;
   X2 = x2 > hb2;
   Y = !(X1 && X2);

   if (Y) {
     y0 = y_high;
   } else {
     y0 = 0.0;
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_nor_2(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   bool X1,X2,Y;
   double x1,x2;
   double y;
   double y_high,hb2;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int nr_y_high = 0;
   const int nr_hb2 = 1;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     y_high = X.rprm[nr_y_high];
     hb2 = 0.5*y_high;
     X.rprm[nr_hb2] = hb2;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   y_high = X.rprm[nr_y_high];
   hb2 = X.rprm[nr_hb2];

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];
   X1 = x1 > hb2;
   X2 = x2 > hb2;
   Y = !(X1 || X2);

   if (Y) {
     y0 = y_high;
   } else {
     y0 = 0.0;
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_not(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   bool X1,Y;
   double x;
   double y;
   double y_high,hb2;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_y_high = 0;
   const int nr_hb2 = 1;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     y_high = X.rprm[nr_y_high];
     hb2 = 0.5*y_high;
     X.rprm[nr_hb2] = hb2;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   y_high = X.rprm[nr_y_high];
   hb2 = X.rprm[nr_hb2];

   x = X.val_vr[nvr_x];
   X1 = x > hb2;
   Y = !X1;

   if (Y) {
     y0 = y_high;
   } else {
     y0 = 0.0;
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_or_2(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   bool X1,X2,Y;
   double x1,x2;
   double y;
   double y_high,hb2;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int nr_y_high = 0;
   const int nr_hb2 = 1;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     y_high = X.rprm[nr_y_high];
     hb2 = 0.5*y_high;
     X.rprm[nr_hb2] = hb2;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   y_high = X.rprm[nr_y_high];
   hb2 = X.rprm[nr_hb2];

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];
   X1 = x1 > hb2;
   X2 = x2 > hb2;
   Y = X1 || X2;

   if (Y) {
     y0 = y_high;
   } else {
     y0 = 0.0;
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_or_3(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   bool X1,X2,X3,Y;
   double x1,x2,x3;
   double y;
   double y_high,hb2;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_x3 = 2;
   const int nvr_y = 3;
   const int nr_y_high = 0;
   const int nr_hb2 = 1;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_x3 = 2;
   const int no_y = 3;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     y_high = X.rprm[nr_y_high];
     hb2 = 0.5*y_high;
     X.rprm[nr_hb2] = hb2;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_x3] = X.val_vr[nvr_x3];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   y_high = X.rprm[nr_y_high];
   hb2 = X.rprm[nr_hb2];

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];
   x3 = X.val_vr[nvr_x3];
   X1 = x1 > hb2;
   X2 = x2 > hb2;
   X3 = x3 > hb2;
   Y = (X1 || X2) || X3;

   if (Y) {
     y0 = y_high;
   } else {
     y0 = 0.0;
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_pole_complex_order_1(Global &G,XbeUsr &X,XbeJac &J) {
   double a2,b2;
   double x;
   double y;
   double z1;
   double a,b,alpha,beta;
   double y_st;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_a = 0;
   const int nr_b = 1;
   const int nr_alpha = 2;
   const int nr_beta = 3;
   const int nst_y_st = 0;
   const int na_z1 = 0;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int ng_1 = 0;
   const int ng_2 = 1;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr [nvr_y] = 0.0;
     X.val_aux[na_z1] = 0.0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
       X.val_aux[na_z1] = 0.0;
     } else if (G.flags[G.i_implicit]) {
       X.h[nf_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
       X.h[nf_2] = X.val_aux[na_z1];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     a = X.rprm[nr_a];
     b = X.rprm[nr_b];
     alpha = X.rprm[nr_alpha];
     beta = X.rprm[nr_beta];
     a2 = a + a;
     b2 = b + b;

     x = X.val_vr[nvr_x];
     y = X.val_vr[nvr_y];
     z1 = X.val_aux[na_z1];

     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = y - X.val_vr_u[nvr_y];
         X.h[nf_2] = z1 - X.val_aux_u[na_z1];
       } else {
         X.f[nf_1] = alpha*y + a2*x - beta*z1;
         X.f[nf_2] = alpha*z1 + b2*x + beta*y;
       }
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = alpha*y + a2*x - beta*z1;
         X.g[ng_2] = alpha*z1 + b2*x + beta*y;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr [ng_1][nvr_x] = a2;
         J.dgdvr [ng_1][nvr_y] = alpha;
         J.dgdaux[ng_1][na_z1] = -beta;

         J.dgdvr [ng_2][nvr_x] = b2;
         J.dgdvr [ng_2][nvr_y] = beta;
         J.dgdaux[ng_2][na_z1] = alpha;
       }
     }
     return;
   }
   return;
}
void x_pole_complex_order_2(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double y;
   double z1,z2,z3;
   double a,b,alpha,beta;
   double y_st;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_a = 0;
   const int nr_b = 1;
   const int nr_alpha = 2;
   const int nr_beta = 3;
   const int nst_y_st = 0;
   const int na_z1 = 0;
   const int na_z2 = 1;
   const int na_z3 = 2;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nf_4 = 3;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   const int ng_4 = 3;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr [nvr_y] = 0.0;
     X.val_aux[na_z1] = 0.0;
     X.val_aux[na_z2] = 0.0;
     X.val_aux[na_z3] = 0.0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
       X.val_aux[na_z1] = 0.0;
       X.val_aux[na_z2] = 0.0;
       X.val_aux[na_z3] = 0.0;
     } else if (G.flags[G.i_implicit]) {
       X.h[nf_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
       X.h[nf_2] = X.val_aux[na_z1];
       X.h[nf_3] = X.val_aux[na_z2];
       X.h[nf_4] = X.val_aux[na_z3];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     a = X.rprm[nr_a];
     b = X.rprm[nr_b];
     alpha = X.rprm[nr_alpha];
     beta = X.rprm[nr_beta];

     x = X.val_vr[nvr_x];
     y = X.val_vr[nvr_y];
     z1 = X.val_aux[na_z1];
     z2 = X.val_aux[na_z2];
     z3 = X.val_aux[na_z3];

     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = y - X.val_vr_u[nvr_y];
         X.h[nf_2] = z1 - X.val_aux_u[na_z1];
         X.h[nf_3] = z2 - X.val_aux_u[na_z2];
         X.h[nf_4] = z3 - X.val_aux_u[na_z3];
       } else {
         X.f[nf_1] = alpha*y - beta*z1 + 2.0*z2;
         X.f[nf_2] = alpha*z1 + beta*y + z3;
         X.f[nf_3] = alpha*z2 - 0.5*beta*z3 + a*x;
         X.f[nf_4] = alpha*z3 + 2.0*beta*z2 + 2.0*b*x;
       }
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = alpha*y - beta*z1 + 2.0*z2;
         X.g[ng_2] = alpha*z1 + beta*y + z3;
         X.g[ng_3] = alpha*z2 - 0.5*beta*z3 + a*x;
         X.g[ng_4] = alpha*z3 + 2.0*beta*z2 + 2.0*b*x;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr [ng_1][nvr_y] = alpha;
         J.dgdaux[ng_1][na_z1] = -beta;
         J.dgdaux[ng_1][na_z2] = 2.0;

         J.dgdvr [ng_2][nvr_y] = beta;
         J.dgdaux[ng_2][na_z1] = alpha;
         J.dgdaux[ng_2][na_z3] = 1.0;

         J.dgdvr [ng_3][nvr_x] = a;
         J.dgdaux[ng_3][na_z2] = alpha;
         J.dgdaux[ng_3][na_z3] = -0.5*beta;

         J.dgdvr [ng_4][nvr_x] = 2.0*b;
         J.dgdaux[ng_4][na_z2] = 2.0*beta;
         J.dgdaux[ng_4][na_z3] = alpha;
       }
     }
     return;
   }
   return;
}
void x_pole_real_order_1(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double y;
   double a,alpha;
   double y_st;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_a = 0;
   const int nr_alpha = 1;
   const int nst_y_st = 0;
   const int nf_1 = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr [nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
     } else if (G.flags[G.i_implicit]) {
       X.h[nf_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     a = X.rprm[nr_a];
     alpha = X.rprm[nr_alpha];
     x = X.val_vr[nvr_x];
     y = X.val_vr[nvr_y];

     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = y - X.val_vr_u[nvr_y];
       } else {
         X.f[nf_1] = alpha*y + a*x;
       }
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = alpha*y + a*x;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_x] = a;
         J.dgdvr[ng_1][nvr_y] = alpha;
       }
     }
     return;
   }
   return;
}
void x_pole_real_order_2(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double y;
   double z1;
   double a,alpha;
   double y_st;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_a = 0;
   const int nr_alpha = 1;
   const int nst_y_st = 0;
   const int na_z1 = 0;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int ng_1 = 0;
   const int ng_2 = 1;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr [nvr_y] = 0.0;
     X.val_aux[na_z1] = 0.0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
       X.val_aux[na_z1] = 0.0;
     } else if (G.flags[G.i_implicit]) {
       X.h[nf_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
       X.h[nf_2] = X.val_aux[na_z1];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     a = X.rprm[nr_a];
     alpha = X.rprm[nr_alpha];
     x = X.val_vr[nvr_x];
     y = X.val_vr[nvr_y];
     z1 = X.val_aux[na_z1];

     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = y - X.val_vr_u[nvr_y];
         X.h[nf_2] = z1 - X.val_aux_u[na_z1];
       } else {
         X.f[nf_1] = alpha*y + a*x;
         X.f[nf_2] = alpha*z1 + x;
       }
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = alpha*y + a*x;
         X.g[ng_2] = alpha*z1 + x;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_x] = a;
         J.dgdvr[ng_1][nvr_y] = alpha;

         J.dgdvr [ng_2][nvr_x] = 1.0;
         J.dgdaux[ng_2][na_z1] = alpha;
       }
     }
     return;
   }
   return;
}
void x_pole_real_order_3(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double y;
   double z1,z2;
   double a,alpha;
   double y_st;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_a = 0;
   const int nr_alpha = 1;
   const int nst_y_st = 0;
   const int na_z1 = 0;
   const int na_z2 = 1;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr [nvr_y] = 0.0;
     X.val_aux[na_z1] = 0.0;
     X.val_aux[na_z2] = 0.0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
       X.val_aux[na_z1] = 0.0;
       X.val_aux[na_z2] = 0.0;
     } else if (G.flags[G.i_implicit]) {
       X.h[nf_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
       X.h[nf_2] = X.val_aux[na_z1];
       X.h[nf_3] = X.val_aux[na_z2];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     a = X.rprm[nr_a];
     alpha = X.rprm[nr_alpha];
     x = X.val_vr[nvr_x];
     y = X.val_vr[nvr_y];
     z1 = X.val_aux[na_z1];
     z2 = X.val_aux[na_z2];

     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = y - X.val_vr_u[nvr_y];
         X.h[nf_2] = z1 - X.val_aux_u[na_z1];
         X.h[nf_3] = z2 - X.val_aux_u[na_z2];
       } else {
         X.f[nf_1] = alpha*y + a*z1;
         X.f[nf_2] = alpha*z1 + z2;
         X.f[nf_3] = alpha*z2 + x;
       }
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = alpha*y + a*z1;
         X.g[ng_2] = alpha*z1 + z2;
         X.g[ng_3] = alpha*z2 + x;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr [ng_1][nvr_y] = alpha;
         J.dgdaux[ng_1][na_z1] = a;

         J.dgdaux[ng_2][na_z1] = alpha;
         J.dgdaux[ng_2][na_z2] = 1.0;

         J.dgdaux[ng_3][na_z2] = alpha;
         J.dgdvr [ng_3][nvr_x] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_pole_real_order_4(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double y;
   double z1,z2,z3;
   double a,alpha;
   double y_st;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_a = 0;
   const int nr_alpha = 1;
   const int nst_y_st = 0;
   const int na_z1 = 0;
   const int na_z2 = 1;
   const int na_z3 = 2;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nf_4 = 3;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   const int ng_4 = 3;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr [nvr_y] = 0.0;
     X.val_aux[na_z1] = 0.0;
     X.val_aux[na_z2] = 0.0;
     X.val_aux[na_z3] = 0.0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
       X.val_aux[na_z1] = 0.0;
       X.val_aux[na_z2] = 0.0;
       X.val_aux[na_z3] = 0.0;
     } else if (G.flags[G.i_implicit]) {
       X.h[nf_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
       X.h[nf_2] = X.val_aux[na_z1];
       X.h[nf_3] = X.val_aux[na_z2];
       X.h[nf_4] = X.val_aux[na_z3];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     a = X.rprm[nr_a];
     alpha = X.rprm[nr_alpha];
     x = X.val_vr[nvr_x];
     y = X.val_vr[nvr_y];
     z1 = X.val_aux[na_z1];
     z2 = X.val_aux[na_z2];
     z3 = X.val_aux[na_z3];

     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = y - X.val_vr_u[nvr_y];
         X.h[nf_2] = z1 - X.val_aux_u[na_z1];
         X.h[nf_3] = z2 - X.val_aux_u[na_z2];
         X.h[nf_4] = z3 - X.val_aux_u[na_z3];
       } else {
         X.f[nf_1] = alpha*y + a*z1;
         X.f[nf_2] = alpha*z1 + z2;
         X.f[nf_3] = alpha*z2 + z3;
         X.f[nf_4] = alpha*z3 + x;
       }
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = alpha*y + a*z1;
         X.g[ng_2] = alpha*z1 + z2;
         X.g[ng_3] = alpha*z2 + z3;
         X.g[ng_4] = alpha*z3 + x;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr [ng_1][nvr_y] = alpha;
         J.dgdaux[ng_1][na_z1] = a;

         J.dgdaux[ng_2][na_z1] = alpha;
         J.dgdaux[ng_2][na_z2] = 1.0;

         J.dgdaux[ng_3][na_z2] = alpha;
         J.dgdaux[ng_3][na_z3] = 1.0;

         J.dgdaux[ng_4][na_z3] = alpha;
         J.dgdvr [ng_4][nvr_x] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_pole_real_order_5(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double y;
   double z1,z2,z3,z4;
   double a,alpha;
   double y_st;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int nr_a = 0;
   const int nr_alpha = 1;
   const int nst_y_st = 0;
   const int na_z1 = 0;
   const int na_z2 = 1;
   const int na_z3 = 2;
   const int na_z4 = 3;
   const int nf_1 = 0;
   const int nf_2 = 1;
   const int nf_3 = 2;
   const int nf_4 = 3;
   const int nf_5 = 4;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   const int ng_4 = 3;
   const int ng_5 = 4;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr [nvr_y] = 0.0;
     X.val_aux[na_z1] = 0.0;
     X.val_aux[na_z2] = 0.0;
     X.val_aux[na_z3] = 0.0;
     X.val_aux[na_z4] = 0.0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
       X.val_aux[na_z1] = 0.0;
       X.val_aux[na_z2] = 0.0;
       X.val_aux[na_z3] = 0.0;
       X.val_aux[na_z4] = 0.0;
     } else if (G.flags[G.i_implicit]) {
       X.h[nf_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
       X.h[nf_2] = X.val_aux[na_z1];
       X.h[nf_3] = X.val_aux[na_z2];
       X.h[nf_4] = X.val_aux[na_z3];
       X.h[nf_5] = X.val_aux[na_z4];
     }
     return;
   }
   if (G.flags[G.i_trns]) {
     a = X.rprm[nr_a];
     alpha = X.rprm[nr_alpha];
     x = X.val_vr[nvr_x];
     y = X.val_vr[nvr_y];
     z1 = X.val_aux[na_z1];
     z2 = X.val_aux[na_z2];
     z3 = X.val_aux[na_z3];
     z4 = X.val_aux[na_z4];

     if (G.flags[G.i_explicit]) {
       if (G.flags[G.i_alg_loop]) {
         X.h[nf_1] = y - X.val_vr_u[nvr_y];
         X.h[nf_2] = z1 - X.val_aux_u[na_z1];
         X.h[nf_3] = z2 - X.val_aux_u[na_z2];
         X.h[nf_4] = z3 - X.val_aux_u[na_z3];
         X.h[nf_5] = z4 - X.val_aux_u[na_z4];
       } else {
         X.f[nf_1] = alpha*y + a*z1;
         X.f[nf_2] = alpha*z1 + z2;
         X.f[nf_3] = alpha*z2 + z3;
         X.f[nf_4] = alpha*z3 + z4;
         X.f[nf_5] = alpha*z4 + x;
       }
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = alpha*y + a*z1;
         X.g[ng_2] = alpha*z1 + z2;
         X.g[ng_3] = alpha*z2 + z3;
         X.g[ng_4] = alpha*z3 + z4;
         X.g[ng_5] = alpha*z4 + x;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr [ng_1][nvr_y] = alpha;
         J.dgdaux[ng_1][na_z1] = a;

         J.dgdaux[ng_2][na_z1] = alpha;
         J.dgdaux[ng_2][na_z2] = 1.0;

         J.dgdaux[ng_3][na_z2] = alpha;
         J.dgdaux[ng_3][na_z3] = 1.0;

         J.dgdaux[ng_4][na_z3] = alpha;
         J.dgdaux[ng_4][na_z4] = 1.0;

         J.dgdaux[ng_5][na_z4] = alpha;
         J.dgdvr [ng_5][nvr_x] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_pulse10(Global &G,XbeUsr &X,XbeJac &J) {
   int intrvl;
   double time0,g_value,g_1,g_2,delt_trns,delt_1,t_trns;
   bool l_high;
   double y;
   int i0,n1;
   double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,
     t20,y_low,y_high,t_rise,t_fall,epsl;
   const int nvr_y = 0;
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
   const int no_y = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     t_rise = X.rprm[nr_t_rise];
     t_fall = X.rprm[nr_t_fall];
     epsl=0.02*min(t_rise,t_fall);
     X.rprm[nr_epsl] = epsl;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   y_low  = X.rprm[nr_y_low ];
   y_high = X.rprm[nr_y_high];
   t_rise = X.rprm[nr_t_rise];
   t_fall = X.rprm[nr_t_fall];
   epsl   = X.rprm[nr_epsl  ];

   i0 = X.iprm[ni_i0];
   n1 = X.iprm[ni_n1];

   time0 = G.time_given_x;

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
     G.time_nextbreak_x = G.time_end;
     if (intrvl == 1) {
       if ((X.rprm[0]-time0) < epsl) {
         G.time_nextbreak_x = X.rprm[0] + delt_trns;
       } else {
         G.time_nextbreak_x = X.rprm[0];
       }
       return;
     } else if (intrvl == (n1+1)) {
       if ((X.rprm[intrvl-2] + delt_1-time0) > epsl) {
         G.time_nextbreak_x = X.rprm[intrvl-2] + delt_1;
       } else {
         G.time_nextbreak_x = G.time_end;
       }
     } else {
       t_trns = X.rprm[intrvl];

       if ((X.rprm[intrvl-2] + delt_1-time0) > epsl) {
         G.time_nextbreak_x = X.rprm[intrvl-2] + delt_1;
       } else if ((X.rprm[intrvl-1]-time0) < epsl) {
         G.time_nextbreak_x = X.rprm[intrvl-1] + delt_trns;
       } else {
         G.time_nextbreak_x = X.rprm[intrvl-1];
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

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = g_value;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = g_value;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - g_value;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_pwl10_xy(Global &G,XbeUsr &X,XbeJac &J) {
   int indx,i;
// i_offset is the total number of points (=10 for pwl10)
   const int i_offset=10;
   double x1a,x2a,y1a,y2a,slope,y0;
   double x;
   double y;
   int n;
   double x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int ni_n = 0;
   const int nr_x1 = 0;
   const int nr_x2 = 1;
   const int nr_x3 = 2;
   const int nr_x4 = 3;
   const int nr_x5 = 4;
   const int nr_x6 = 5;
   const int nr_x7 = 6;
   const int nr_x8 = 7;
   const int nr_x9 = 8;
   const int nr_x10 = 9;
   const int nr_y1 = 10;
   const int nr_y2 = 11;
   const int nr_y3 = 12;
   const int nr_y4 = 13;
   const int nr_y5 = 14;
   const int nr_y6 = 15;
   const int nr_y7 = 16;
   const int nr_y8 = 17;
   const int nr_y9 = 18;
   const int nr_y10 = 19;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
// cout << "pwl10_xy.xbe" << endl;

   if (G.flags[G.i_one_time_parms]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   n = X.iprm[ni_n];
   x = X.val_vr[nvr_x];
   y = X.val_vr[nvr_y];

   indx = n-1;
   for (i=0; i < n; i++) {
     if (x <= X.rprm[i]) {
       indx = i-1;
       break;
     }
   }
   if (indx == (-1)) {
     y0 = X.rprm[i_offset];
   } else if (indx == (n-1)) {
     y0 = X.rprm[i_offset+n-1];
   } else {
     x1a = X.rprm[indx];
     x2a = X.rprm[indx+1];
     y1a = X.rprm[i_offset+indx];
     y2a = X.rprm[i_offset+indx+1];

     if (abs(x2a-x1a) < 1.0e-12) {
       cout << "pwl10_xy.xbe: x1 and x2 are too close!" << endl;
       cout << "   Halting.." << endl;
       exit (1);
     }
     slope = (y2a-y1a)/(x2a-x1a);
     y0 = y1a + slope*(x-x1a);
   }

   if (G.flags[G.i_init_guess]) {
//   cout << "pwl10_xy: init_guess: y0 =" << y0 << endl;
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_pwl20(Global &G,XbeUsr &X,XbeJac &J) {
   const int nmax=20;
   int i,intrvl;
   double time0,y0,slp1;
   const double epsl=1.0e-15;
   double y;
   int n;
   double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,
     t20,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v20;
   const int nvr_y = 0;
   const int ni_n = 0;
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
   const int nr_v1 = 20;
   const int nr_v2 = 21;
   const int nr_v3 = 22;
   const int nr_v4 = 23;
   const int nr_v5 = 24;
   const int nr_v6 = 25;
   const int nr_v7 = 26;
   const int nr_v8 = 27;
   const int nr_v9 = 28;
   const int nr_v10 = 29;
   const int nr_v11 = 30;
   const int nr_v12 = 31;
   const int nr_v13 = 32;
   const int nr_v14 = 33;
   const int nr_v15 = 34;
   const int nr_v16 = 35;
   const int nr_v17 = 36;
   const int nr_v18 = 37;
   const int nr_v19 = 38;
   const int nr_v20 = 39;
   const int no_y = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   n = X.iprm[ni_n];
   time0 = G.time_given_x;

   if (G.flags[G.i_next_time]) {
     if (time0 >= X.rprm[n-1]) {
       G.time_nextbreak_x = G.time_end;
     } else if (time0 < X.rprm[0]) {
       G.time_nextbreak_x = X.rprm[0];
     } else {
       for (int i=1; i < n; i++) {
         if (time0 < X.rprm[i]) {
           if (abs(X.rprm[i]-time0) < epsl) {
             intrvl = i;
           } else {
             intrvl = i-1;
           }
           break;
         }
       }
       G.time_nextbreak_x = X.rprm[intrvl+1];
     }
     return;
   }
   if (time0 >= X.rprm[n-1]) {
     y0 = X.rprm[nmax+n-1];
   } else if (time0 < X.rprm[0]) {
     y0 = X.rprm[nmax];
   } else {
     for (int i=1; i < n; i++) {
       if (time0 < X.rprm[i]) {
         if (abs(X.rprm[i]-time0) < epsl) {
           intrvl = i;
         } else {
           intrvl = i-1;
         }
         break;
       }
     }
     slp1 = (X.rprm[nmax+intrvl+1]-X.rprm[nmax+intrvl])/
            (X.rprm[intrvl+1]-X.rprm[intrvl]);
     y0 = X.rprm[nmax+intrvl]+
             slp1*(time0-X.rprm[intrvl]);
   }

   if (G.flags[G.i_init_guess]) {
     cout << "pwl20: init_guess: y0 =" << y0 << endl;
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_pwm20_1(Global &G,XbeUsr &X,XbeJac &J) {
   int k,indx0,n_periods;
   double t_add,tp,tp0,tp1,tp_last,time0,time0p,time0pp;
   double t_diff,t_next_1,y0;
   double y;
   int ndata,index_last,level_0minus;
   double t_1,t_2,t_3,t_4,t_5,t_6,t_7,t_8,t_9,t_10,t_11,t_12,t_13,t_14,t_15,
     t_16,t_17,t_18,t_19,t_20,theta_1,theta_2,theta_3,theta_4,theta_5,theta_6,
     theta_7,theta_8,theta_9,theta_10,theta_11,theta_12,theta_13,theta_14,
     theta_15,theta_16,theta_17,theta_18,theta_19,theta_20,frequency,y_low,
     y_high,theta_delay,t_delay,t_period,epsl1,epsl2,epsl3;
   const int nvr_y = 0;
   const int ni_ndata = 0;
   const int ni_index_last = 1;
   const int ni_level_0minus = 2;
   const int nr_t_1 = 0;
   const int nr_t_2 = 1;
   const int nr_t_3 = 2;
   const int nr_t_4 = 3;
   const int nr_t_5 = 4;
   const int nr_t_6 = 5;
   const int nr_t_7 = 6;
   const int nr_t_8 = 7;
   const int nr_t_9 = 8;
   const int nr_t_10 = 9;
   const int nr_t_11 = 10;
   const int nr_t_12 = 11;
   const int nr_t_13 = 12;
   const int nr_t_14 = 13;
   const int nr_t_15 = 14;
   const int nr_t_16 = 15;
   const int nr_t_17 = 16;
   const int nr_t_18 = 17;
   const int nr_t_19 = 18;
   const int nr_t_20 = 19;
   const int nr_theta_1 = 20;
   const int nr_theta_2 = 21;
   const int nr_theta_3 = 22;
   const int nr_theta_4 = 23;
   const int nr_theta_5 = 24;
   const int nr_theta_6 = 25;
   const int nr_theta_7 = 26;
   const int nr_theta_8 = 27;
   const int nr_theta_9 = 28;
   const int nr_theta_10 = 29;
   const int nr_theta_11 = 30;
   const int nr_theta_12 = 31;
   const int nr_theta_13 = 32;
   const int nr_theta_14 = 33;
   const int nr_theta_15 = 34;
   const int nr_theta_16 = 35;
   const int nr_theta_17 = 36;
   const int nr_theta_18 = 37;
   const int nr_theta_19 = 38;
   const int nr_theta_20 = 39;
   const int nr_frequency = 40;
   const int nr_y_low = 41;
   const int nr_y_high = 42;
   const int nr_theta_delay = 43;
   const int nr_t_delay = 44;
   const int nr_t_period = 45;
   const int nr_epsl1 = 46;
   const int nr_epsl2 = 47;
   const int nr_epsl3 = 48;
   const int no_y = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     frequency = X.rprm[nr_frequency];
     t_period = 1.0/frequency;
     X.rprm[nr_t_period] = t_period;

     ndata = X.iprm[ni_ndata];

     if ((ndata % 2) != 0) {
       cout << "pwm20_1.xbe: ndata=" << ndata << " is not allowed" << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     y_high = X.rprm[nr_y_high];
     y_low  = X.rprm[nr_y_low ];
     if (y_high < y_low) {
       cout << "pwm20_1.xbe: y_high and y_low are not defined correctly." << endl;
       cout << "  y_high=" << y_high << " y_low=" << y_low << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     for (int i = 0; i < ndata; i++) {
       X.rprm[i] = (t_period/360.0)*X.rprm[i+20];
     }
     if ((abs(X.rprm[ndata-1+20]-360.0) < 0.1) && (abs(X.rprm[0+20]) < 0.1)) {
       cout << "pwm20_1.xbe: time points at both 0 and 360 are" << endl;
       cout << "  not allowed. Halting..." << endl; exit(1);
     }
     if (abs(X.rprm[ndata-1+20]-360.0) < 0.1) {
       for (int i = (ndata-2); i >= 0; i--) {
         X.rprm[i+1] = X.rprm[i];
       }
       X.rprm[0] = 0.0;
     }

     X.iprm[ni_index_last] = -1;

     X.rprm[nr_t_delay] = (t_period/360.0)*X.rprm[nr_theta_delay];
     X.rprm[nr_epsl1] = (t_period/1000.0);
     epsl2 = 0.01*t_period/360.0;
     X.rprm[nr_epsl2] = epsl2;
     X.rprm[nr_epsl3] = 1.1*epsl2;

     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = 0.0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   index_last = X.iprm[ni_index_last];

   t_delay  = X.rprm[nr_t_delay ];
   t_period = X.rprm[nr_t_period];
   epsl1    = X.rprm[nr_epsl1   ];
   epsl2    = X.rprm[nr_epsl2   ];
   epsl3    = X.rprm[nr_epsl3   ];

   ndata = X.iprm[ni_ndata];
   time0 = G.time_given_x;

   if (time0 < 0.0) {
     cout << "pwm20_1.gbe: time0=" << time0 << " is not allowed." << endl;
     cout << "  Halting..." << endl; exit(1);
   }

   time0p = time0 - t_delay;
   if (time0p < 0) time0p = time0p + t_period;
   time0pp = fmod(time0p,t_period);
   if (time0pp < 0) time0pp = 0.0;

   if (G.flags[G.i_save_history]) {
     tp_last = X.rprm[ndata-1];
     indx0 = X.iprm[ni_index_last];

     if (time0pp > tp_last) {
       index_last = -1;
     } else {
       for (int i=(indx0+1);i < ndata;i++) {
         tp = X.rprm[i];
         if (tp > time0pp) {
           index_last = i-2;
           if (index_last < -1) index_last = -1;
           break;
         }
       }
     }
     X.iprm[ni_index_last] = index_last;
   }

   if (G.flags[G.i_next_time]) {
     t_diff = (time0-t_delay)/t_period;
     n_periods = (int)(t_diff);
     t_add = n_periods*t_period + t_delay;
     tp_last = X.rprm[ndata-1];

     if (time0pp > tp_last) {
       if (abs(X.rprm[0]) < epsl1) {
         if ((t_period-time0pp) < epsl3) {
           t_next_1 = t_period + epsl2;
         } else {
           t_next_1 = t_period - epsl2;
         }
       } else {
         t_next_1 = t_period + X.rprm[0] - epsl2;
       }
     } else {
       for (int i = (index_last+1); i < ndata; i++) {
         tp1 = X.rprm[i];
         if (time0pp < tp1) {
           if ((tp1-time0pp) < epsl3) {
             t_next_1 = tp1 + epsl2;
           } else {
             t_next_1 = tp1 - epsl2;
           }
           break;
         }
       }
     }
     G.time_nextbreak_x = t_next_1 + t_add;
     return;
   }

   if (G.flags[G.i_trns]) {
     level_0minus = X.iprm[ni_level_0minus];

     for (int i = (index_last+1); i < ndata; i++) {
       tp1 = X.rprm[i];
       if (time0pp < tp1) {
         k = i % 2;
         if (level_0minus == 0) {
           if (k == 0) {
             y0 = X.rprm[nr_y_low];
           } else {
             y0 = X.rprm[nr_y_high];
           }
         } else {
           if (k == 0) {
             y0 = X.rprm[nr_y_high];
           } else {
             y0 = X.rprm[nr_y_low];
           }
         }
         goto jump1;
       }
     }
     if (level_0minus == 0) {
       y0 = X.rprm[nr_y_low];
     } else {
       y0 = X.rprm[nr_y_high];
     }
     jump1:;
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_sampler(Global &G,XbeUsr &X,XbeJac &J) {
   double time0;
   int n;
   double t0_new,t_a,t_b,t_c,t_d,y0;
   double x;
   double y;
   int index;
   double T,t0,v_previous,dt,epsl1,epsl2;
   double y_st;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int ni_index = 0;
   const int nr_T = 0;
   const int nr_t0 = 1;
   const int nr_v_previous = 2;
   const int nr_dt = 3;
   const int nr_epsl1 = 4;
   const int nr_epsl2 = 5;
   const int nst_y_st = 0;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     dt = X.rprm[nr_dt];
     X.rprm[nr_epsl1] = dt/10.0;
     X.rprm[nr_epsl2] = dt/100.0;

     index = X.iprm[ni_index];
     G.sampler_flag[index] = 0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
       }   
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }   
     }   
     return;
   }   
   if (G.flags[G.i_save_history]) {
     X.rprm[nr_v_previous] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_next_time]) {
     T          = X.rprm[nr_T         ];
     t0         = X.rprm[nr_t0        ];
     v_previous = X.rprm[nr_v_previous];
     dt         = X.rprm[nr_dt        ];

     time0 = G.time_given_x;

     if (time0 < t0) {
       n = ((t0-time0)/T)+1;
       t0_new = t0 - n*T;
     } else {
       t0_new = t0;
     }

     epsl1 = X.rprm[nr_epsl1];
     epsl2 = X.rprm[nr_epsl2];
     t_a = time0 - t0_new;
     t_b = fmod(t_a,T);
     if (abs(t_b-T) < epsl2) t_b = 0.0;
     t_c = T - dt;

     if (G.flags[G.i_next_time]) {
       if (t_b == 0.0) {
         G.time_nextbreak_x = time0 + t_c;
       } else if (abs(t_b-t_c) <= epsl2) {
         G.time_nextbreak_x = time0 + dt;
       } else if (t_b < t_c) {
         G.time_nextbreak_x = time0 - t_b + t_c;
       } else {
         G.time_nextbreak_x = time0 - t_b + T;
       }
       return;
     } else {
       t_d = abs(t_b-T);
       index = X.iprm[ni_index];
       if ((t_d < epsl1) || (abs(t_b) < epsl1)) {
         y0 = X.val_vr[nvr_x];
         G.sampler_flag[index] = 1;
       } else {
         y0 = v_previous;
         G.sampler_flag[index] = 0;
       }
       cout << "sampler: time0: " << time0 << " sampler_flag: " << G.sampler_flag[index] << endl;
       if (G.flags[G.i_explicit]) {
         X.val_vr[nvr_y] = y0;
       } else if (G.flags[G.i_implicit]) {
         if (G.flags[G.i_function]) {
           X.g[ng_1] = X.val_vr[nvr_y] - y0;
         }
         if (G.flags[G.i_jacobian]) {
           J.dgdvr[ng_1][nvr_y] = 1.0;
         }
       }
     }
     return;
   }
   return;
}
void x_sampler_1(Global &G,XbeUsr &X,XbeJac &J) {
   double time0;
   int n;
   double t0_new,t_a,t_b,t_c,t_d,y0;
   double x;
   double y;
   int index;
   double T,t0,v_previous,dt,epsl1,epsl2;
   double y_st;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int ni_index = 0;
   const int nr_T = 0;
   const int nr_t0 = 1;
   const int nr_v_previous = 2;
   const int nr_dt = 3;
   const int nr_epsl1 = 4;
   const int nr_epsl2 = 5;
   const int nst_y_st = 0;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     dt = X.rprm[nr_dt];
     X.rprm[nr_epsl1] = dt/10.0;
     X.rprm[nr_epsl2] = dt/100.0;

     index = X.iprm[ni_index];
     G.sampler_flag[index] = 0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = 0.0;
     return;
   }
   if (G.flags[G.i_time_parms] || G.flags[G.i_next_time]) {
     T  = X.rprm[nr_T ];
     t0 = X.rprm[nr_t0];
     dt = X.rprm[nr_dt];

     time0 = G.time_given_x;

     if (time0 < t0) {
       n = ((t0-time0)/T)+1;
       t0_new = t0 - n*T;
     } else {
       t0_new = t0;
     }

     epsl1 = X.rprm[nr_epsl1];
     epsl2 = X.rprm[nr_epsl2];
     t_a = time0 - t0_new;
     t_b = fmod(t_a,T);
     if (abs(t_b-T) < epsl2) t_b = 0.0;
     t_c = T - dt;

     if (G.flags[G.i_next_time]) {
       if (t_b == 0.0) {
         G.time_nextbreak_x = time0 + t_c;
       } else if (abs(t_b-t_c) <= epsl2) {
         G.time_nextbreak_x = time0 + dt;
       } else if (t_b < t_c) {
         G.time_nextbreak_x = time0 - t_b + t_c;
       } else {
         G.time_nextbreak_x = time0 - t_b + T;
       }
       return;
     } else if (G.flags[G.i_time_parms]) {
       t_d = abs(t_b-T);
       index = X.iprm[ni_index];
       if ((t_d < epsl1) || (abs(t_b) < epsl1)) {
         G.sampler_flag[index] = 1;
       } else {
         G.sampler_flag[index] = 0;
       }
       return;
     }
   }
   if (G.flags[G.i_save_history]) {
     X.rprm[nr_v_previous] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = X.stprm[nst_y_st];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - X.stprm[nst_y_st];
       }   
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }   
     }   
     return;
   }   
   if (G.flags[G.i_trns]) {
     index = X.iprm[ni_index];
     if (G.sampler_flag[index] == 1) {
       y0 = X.val_vr[nvr_x];
     } else {
       y0 = X.rprm[nr_v_previous];
     }
//   cout << "sampler_1: time0: " << time0 << " sampler_flag: " << G.sampler_flag[index] << endl;
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_signal_switch(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   double s,x1,x2;
   double y;
   double s_high,s_cross;
   const int nvr_s = 0;
   const int nvr_x1 = 1;
   const int nvr_x2 = 2;
   const int nvr_y = 3;
   const int nr_s_high = 0;
   const int nr_s_cross = 1;
   const int no_s = 0;
   const int no_x1 = 1;
   const int no_x2 = 2;
   const int no_y = 3;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     s_high = X.rprm[nr_s_high];
     s_cross = 0.5*s_high;
     X.rprm[nr_s_cross] = s_cross;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_s ] = X.val_vr[nvr_s ];
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   s_high = X.rprm[nr_s_high];
   s_cross = X.rprm[nr_s_cross];

   s = X.val_vr[nvr_s];
   x2 = X.val_vr[nvr_x2];

   if (s > s_cross) {
     y0 = X.val_vr[nvr_x1];
   } else {
     y0 = X.val_vr[nvr_x2];
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_signum(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   double x;
   double y;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   x = X.val_vr[nvr_x];
   if (x == 0.0) {
     y0 = 0.0;
   } else if (x > 0.0) {
     y0 = 1.0;
   } else {
     y0 = -1.0;
   }   

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       y = X.val_vr[nvr_y];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = y - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_sin(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   double x;
   double y;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int no_x = 0;
   const int no_y = 1;
   const int ng_1 = 0;
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = sin(X.val_vr[nvr_x]);
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = sin(X.val_vr[nvr_x]);
     } else if (G.flags[G.i_implicit]) {
       x = X.val_vr[nvr_x];
       y = X.val_vr[nvr_y];

       if (G.flags[G.i_function]) {
         X.g[ng_1] = y - sin(x);
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] =  1.0;
         J.dgdvr[ng_1][nvr_x] = -cos(x);
       }
     }
     return;
   }
   return;
}
void x_src_ac(Global &G,XbeUsr &X,XbeJac &J) {
   double v0;
   double y;
   double a,f_hz,phi_deg,t0,dc,omega,phi_rad;
   const int nvr_y = 0;
   const int nr_a = 0;
   const int nr_f_hz = 1;
   const int nr_phi_deg = 2;
   const int nr_t0 = 3;
   const int nr_dc = 4;
   const int nr_omega = 5;
   const int nr_phi_rad = 6;
   const int no_y = 0;
   const int ng_1 = 0;

   if (G.flags[G.i_one_time_parms]) {
     f_hz    = X.rprm[nr_f_hz   ];
     phi_deg = X.rprm[nr_phi_deg];

     omega   = G.twopi*f_hz;
     phi_rad = G.deg_to_rad*phi_deg;

     X.rprm[nr_omega  ] = omega;
     X.rprm[nr_phi_rad] = phi_rad;

     return;
   }
   if (G.flags[G.i_init_guess]) {
     a       = X.rprm[nr_a      ];
     t0      = X.rprm[nr_t0     ];
     dc      = X.rprm[nr_dc     ];
     omega   = X.rprm[nr_omega  ];
     phi_rad = X.rprm[nr_phi_rad];

     v0 = a*sin(omega*(G.time_given_x-t0)+phi_rad) + dc;

     X.val_vr[nvr_y] = v0;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     a       = X.rprm[nr_a      ];
     t0      = X.rprm[nr_t0     ];
     dc      = X.rprm[nr_dc     ];
     omega   = X.rprm[nr_omega  ];
     phi_rad = X.rprm[nr_phi_rad];

     v0 = a*sin(omega*(G.time_given_x-t0)+phi_rad) + dc;
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = v0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - v0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_srff_nand(Global &G,XbeUsr &X,XbeJac &J) {
   double time0,q0,q0bar;
   double s,r;
   double q,qbar;
   double x_high,x_cross,q_prev;
   double q_st;
   const int nvr_s = 0;
   const int nvr_r = 1;
   const int nvr_q = 2;
   const int nvr_qbar = 3;
   const int nr_x_high = 0;
   const int nr_x_cross = 1;
   const int nr_q_prev = 2;
   const int nst_q_st = 0;
   const int no_s = 0;
   const int no_r = 1;
   const int no_q = 2;
   const int ng_1 = 0;
   const int ng_2 = 1;
   if (G.flags[G.i_one_time_parms]) {
     x_high = X.rprm[nr_x_high];

     if (x_high < 0.0) {
       cout << "srff_nand.xbe: check x_high. Halting..." << endl;
       exit(1);
     }
     x_cross = 0.5*x_high;
     X.rprm[nr_x_cross] = x_cross;
     X.rprm[nr_q_prev] = 0.0;
     return;
   }

   if (G.flags[G.i_outvar]) {
     X.outprm[no_s] = X.val_vr[nvr_s];
     X.outprm[no_r] = X.val_vr[nvr_r];
     X.outprm[no_q] = X.val_vr[nvr_q];
     return;
   }
   if (G.flags[G.i_save_history]) {
     X.rprm[nr_q_prev] = X.val_vr[nvr_q];
     return;
   }
   if (G.flags[G.i_startup]) {
     q_st = X.stprm[nst_q_st];

     x_cross = X.rprm[nr_x_cross];
     x_high  = X.rprm[nr_x_high ];

     if (G.flags[G.i_explicit]) {
       if (q_st > x_cross) {
         X.val_vr[nvr_q   ] = x_high;
         X.val_vr[nvr_qbar] = 0.0;
       } else {
         X.val_vr[nvr_q   ] = 0.0;
         X.val_vr[nvr_qbar] = x_high;
       }
     } else if (G.flags[G.i_implicit]) {
       if (q_st > x_cross) {
         q0    = x_high;
         q0bar = 0.0;
       } else {
         q0    = 0.0;
         q0bar = x_high;
       }
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_q   ] - q0;
         X.g[ng_2] = X.val_vr[nvr_qbar] - q0bar;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_q   ] = 1.0;
         J.dgdvr[ng_2][nvr_qbar] = 1.0;
       }
     }
     return;
   }
   x_cross = X.rprm[nr_x_cross];
   x_high  = X.rprm[nr_x_high ];
   q_prev  = X.rprm[nr_q_prev ];

   s = X.val_vr[nvr_s];
   r = X.val_vr[nvr_r];

   if (s < x_cross) {
     if (r < x_cross) {
//      Both outputs zero, not expected generally.
        q0 = x_high;
     } else {
        q0 = 0.0;
     }
   } else {
     if (r > x_cross) {
       q0= 0.0;
     } else {
       q0 = q_prev;
     }
   }
   if (G.flags[G.i_init_guess]) {
     if (q0 > x_cross) {
       q0bar = 0.0;
     } else {
       q0bar = x_high;
     }
     X.val_vr[nvr_q   ] = q0;
     X.val_vr[nvr_qbar] = q0bar;
     return;
   }
   if (G.flags[G.i_trns]) {

     if (G.iter_trns_x == 0) {
       if ((G.flags[G.i_slv_previous]) || (G.flags[G.i_slv_readfile])) {
         q0 = X.val_vr[nvr_q];
       }
     }

     x_cross = X.rprm[nr_x_cross];
     x_high  = X.rprm[nr_x_high ];

     if (q0 > x_cross) {
       q0bar = 0.0;
     } else {
       q0bar = x_high;
     }
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_q   ] = q0;
       X.val_vr[nvr_qbar] = q0bar;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_q   ] - q0;
         X.g[ng_2] = X.val_vr[nvr_qbar] - q0bar;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_q   ] = 1.0;
         J.dgdvr[ng_2][nvr_qbar] = 1.0;
       }
     }
   }
   return;
}
void x_srff_nor(Global &G,XbeUsr &X,XbeJac &J) {
   double time0,q0,q0bar;
   double s,r;
   double q,qbar;
   double x_high,x_cross,q_prev;
   double q_st;
   const int nvr_s = 0;
   const int nvr_r = 1;
   const int nvr_q = 2;
   const int nvr_qbar = 3;
   const int nr_x_high = 0;
   const int nr_x_cross = 1;
   const int nr_q_prev = 2;
   const int nst_q_st = 0;
   const int no_s = 0;
   const int no_r = 1;
   const int no_q = 2;
   const int ng_1 = 0;
   const int ng_2 = 1;
   if (G.flags[G.i_one_time_parms]) {
     x_high = X.rprm[nr_x_high];

     if (x_high < 0.0) {
       cout << "srff_nor.xbe: check x_high. Halting..." << endl;
       exit(1);
     }
     x_cross = 0.5*x_high;
     X.rprm[nr_x_cross] = x_cross;
     X.rprm[nr_q_prev] = 0.0;
     return;
   }

   if (G.flags[G.i_outvar]) {
     X.outprm[no_s] = X.val_vr[nvr_s];
     X.outprm[no_r] = X.val_vr[nvr_r];
     X.outprm[no_q] = X.val_vr[nvr_q];
     return;
   }
   if (G.flags[G.i_save_history]) {
     X.rprm[nr_q_prev] = X.val_vr[nvr_q];
     return;
   }
   if (G.flags[G.i_startup]) {
     q_st = X.stprm[nst_q_st];

     x_cross = X.rprm[nr_x_cross];
     x_high  = X.rprm[nr_x_high ];

     if (G.flags[G.i_explicit]) {
       if (q_st > x_cross) {
         X.val_vr[nvr_q   ] = x_high;
         X.val_vr[nvr_qbar] = 0.0;
       } else {
         X.val_vr[nvr_q   ] = 0.0;
         X.val_vr[nvr_qbar] = x_high;
       }
     } else if (G.flags[G.i_implicit]) {
       if (q_st > x_cross) {
         q0    = x_high;
         q0bar = 0.0;
       } else {
         q0    = 0.0;
         q0bar = x_high;
       }
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_q   ] - q0;
         X.g[ng_2] = X.val_vr[nvr_qbar] - q0bar;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_q   ] = 1.0;
         J.dgdvr[ng_2][nvr_qbar] = 1.0;
       }
     }
     return;
   }
   x_cross = X.rprm[nr_x_cross];
   x_high  = X.rprm[nr_x_high ];
   q_prev  = X.rprm[nr_q_prev ];

   s = X.val_vr[nvr_s];
   r = X.val_vr[nvr_r];

   if (s > x_cross) {
     if (r > x_cross) {
//      Both outputs zero, not expected generally.
        q0 = 0.0;
     } else {
        q0 = x_high;
     }
   } else {
     if (r > x_cross) {
       q0= 0.0;
     } else {
       q0 = q_prev;
     }
   }
   if (G.flags[G.i_init_guess]) {
     if (q0 > x_cross) {
       q0bar = 0.0;
     } else {
       q0bar = x_high;
     }
     X.val_vr[nvr_q   ] = q0;
     X.val_vr[nvr_qbar] = q0bar;
     return;
   }
   if (G.flags[G.i_trns]) {

     if (G.iter_trns_x == 0) {
       if ((G.flags[G.i_slv_previous]) || (G.flags[G.i_slv_readfile])) {
         q0 = X.val_vr[nvr_q];
       }
     }

     x_cross = X.rprm[nr_x_cross];
     x_high  = X.rprm[nr_x_high ];

     if (q0 > x_cross) {
       q0bar = 0.0;
     } else {
       q0bar = x_high;
     }
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_q   ] = q0;
       X.val_vr[nvr_qbar] = q0bar;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_q   ] - q0;
         X.g[ng_2] = X.val_vr[nvr_qbar] - q0bar;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_q   ] = 1.0;
         J.dgdvr[ng_2][nvr_qbar] = 1.0;
       }
     }
   }
   return;
}
void x_sum_2(Global &G,XbeUsr &X,XbeJac &J) {
   double x1,x2;
   double y;
   double k1,k2;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int nr_k1 = 0;
   const int nr_k2 = 1;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   k1 = X.rprm[nr_k1];
   k2 = X.rprm[nr_k2];
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = k1*X.val_vr[nvr_x1] + k2*X.val_vr[nvr_x2];
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = k1*X.val_vr[nvr_x1] + k2*X.val_vr[nvr_x2];
//     cout << "sum_2: k1 = " << k1 << " k2 = " << k2
//       << " X.val_vr[nvr_x1] = " << X.val_vr[nvr_x1]
//       << " X.val_vr[nvr_x2] = " << X.val_vr[nvr_x2]
//       << " X.val_vr[nvr_y] = " << X.val_vr[nvr_y]
//       << endl;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - k1*X.val_vr[nvr_x1] - k2*X.val_vr[nvr_x2];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
         J.dgdvr[ng_1][nvr_x1] = -k1;
         J.dgdvr[ng_1][nvr_x2] = -k2;
       }
     }
     return;
   }
   return;
}
void x_sum_3(Global &G,XbeUsr &X,XbeJac &J) {
   double x1,x2,x3;
   double y;
   double k1,k2,k3;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_x3 = 2;
   const int nvr_y = 3;
   const int nr_k1 = 0;
   const int nr_k2 = 1;
   const int nr_k3 = 2;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_x3 = 2;
   const int no_y = 3;
   const int ng_1 = 0;
   k1 = X.rprm[nr_k1];
   k2 = X.rprm[nr_k2];
   k3 = X.rprm[nr_k3];
   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = k1*X.val_vr[nvr_x1] + k2*X.val_vr[nvr_x2] + k3*X.val_vr[nvr_x3];
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_x3] = X.val_vr[nvr_x3];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = k1*X.val_vr[nvr_x1] + k2*X.val_vr[nvr_x2] + k3*X.val_vr[nvr_x3];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y]
            - k1*X.val_vr[nvr_x1]
            - k2*X.val_vr[nvr_x2]
            - k3*X.val_vr[nvr_x3];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
         J.dgdvr[ng_1][nvr_x1] = -k1;
         J.dgdvr[ng_1][nvr_x2] = -k2;
         J.dgdvr[ng_1][nvr_x3] = -k3;
       }
     }
     return;
   }
   return;
}
void x_triangle_1(Global &G,XbeUsr &X,XbeJac &J) {
   double y0,t0_new,delta_min,t_a,t_b,tnext_p;
   int n;
   double y;
   double T1,T2,L1,L2,t0,T,slope1,slope2,epsl;
   const int nvr_y = 0;
   const int nr_T1 = 0;
   const int nr_T2 = 1;
   const int nr_L1 = 2;
   const int nr_L2 = 3;
   const int nr_t0 = 4;
   const int nr_T = 5;
   const int nr_slope1 = 6;
   const int nr_slope2 = 7;
   const int nr_epsl = 8;
   const int no_y = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     T1 = X.rprm[nr_T1];
     T2 = X.rprm[nr_T2];
     L1 = X.rprm[nr_L1];
     L2 = X.rprm[nr_L2];

     T = T1 + T2;
     slope1 = (L2-L1)/T1;
     slope2 = (L1-L2)/T2;
     epsl = T/1000.0;

     X.rprm[nr_T] = T;
     X.rprm[nr_slope1] = slope1;
     X.rprm[nr_slope2] = slope2;
     X.rprm[nr_epsl] = epsl;

     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   t0   = X.rprm[nr_t0  ];
   T    = X.rprm[nr_T   ];
   T1   = X.rprm[nr_T1  ];
   epsl = X.rprm[nr_epsl];

   if (G.time_given_x < t0) {
     n = ((t0-G.time_given_x)/T) + 1;
     t0_new = t0-n*T;
   } else {
     t0_new = t0;
   }
   t_a = G.time_given_x-t0_new;
   t_b = fmod(t_a,T);

   if (abs(t_b-T) < epsl) t_b = 0.0;

   if (G.flags[G.i_next_time]) {
     if (t_b < T1) {
       tnext_p = T1;
     } else {
       tnext_p = T;
     }
     G.time_nextbreak_x = G.time_given_x + (tnext_p-t_b);
     return;
   }

   L1 = X.rprm[nr_L1];
   L2 = X.rprm[nr_L2];
   slope1 = X.rprm[nr_slope1];
   slope2 = X.rprm[nr_slope2];

   if (t_b < T1) {
     y0 = L1 + slope1*t_b;
   } else {
     y0 = L2 + slope2*(t_b-T1);
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_triangle_2(Global &G,XbeUsr &X,XbeJac &J) {
   double y0,t0_new,delta_min,t_a,t_b,tnext_p;
   int n;
   double y;
   int flag_frequency,flag_period;
   double T,frequency,L1,L2,t0,slope1,slope2,epsl,T1,T2;
   const int nvr_y = 0;
   const int ni_flag_frequency = 0;
   const int ni_flag_period = 1;
   const int nr_T = 0;
   const int nr_frequency = 1;
   const int nr_L1 = 2;
   const int nr_L2 = 3;
   const int nr_t0 = 4;
   const int nr_slope1 = 5;
   const int nr_slope2 = 6;
   const int nr_epsl = 7;
   const int nr_T1 = 8;
   const int nr_T2 = 9;
   const int no_y = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     flag_frequency = X.iprm[ni_flag_frequency];
     flag_period    = X.iprm[ni_flag_period   ];

     if ((flag_frequency == 0) && (flag_period == 0)) {
       cout << "triangle_2.xbe: check flag_frequency and flag_period" << endl;
       cout << "  Both cannot be zero." << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     if ((flag_frequency != 0) && (flag_period != 0)) {
       cout << "triangle_2.xbe: check flag_frequency and flag_period" << endl;
       cout << "  Both cannot be non-zero." << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     if (flag_frequency != 0) {
       frequency = X.rprm[nr_frequency];
       T = 1.0/frequency;
       X.rprm[nr_T] = T;
     }
     if (flag_period != 0) {
       T = X.rprm[nr_T];
       frequency = 1.0/T;
       X.rprm[nr_frequency] = frequency;
     }

     T = X.rprm[nr_T];
     T1 = 0.5*T;
     T2 = T1;

     L1 = X.rprm[nr_L1];
     L2 = X.rprm[nr_L2];

     slope1 = (L2-L1)/T1;
     slope2 = (L1-L2)/T2;
     epsl = T/1000.0;

     X.rprm[nr_slope1] = slope1;
     X.rprm[nr_slope2] = slope2;
     X.rprm[nr_epsl] = epsl;

     X.rprm[nr_T1] = T1;
     X.rprm[nr_T2] = T2;

     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   t0   = X.rprm[nr_t0  ];
   T    = X.rprm[nr_T   ];
   T1   = X.rprm[nr_T1  ];
   epsl = X.rprm[nr_epsl];

   if (G.time_given_x < t0) {
     n = ((t0-G.time_given_x)/T) + 1;
     t0_new = t0-n*T;
   } else {
     t0_new = t0;
   }
   t_a = G.time_given_x-t0_new;
   t_b = fmod(t_a,T);

   if (abs(t_b-T) < epsl) t_b = 0.0;

   if (G.flags[G.i_next_time]) {
     if (t_b < T1) {
       tnext_p = T1;
     } else {
       tnext_p = T;
     }
     G.time_nextbreak_x = G.time_given_x + (tnext_p-t_b);
     return;
   }

   L1 = X.rprm[nr_L1];
   L2 = X.rprm[nr_L2];
   slope1 = X.rprm[nr_slope1];
   slope2 = X.rprm[nr_slope2];

   if (t_b < T1) {
     y0 = L1 + slope1*t_b;
   } else {
     y0 = L2 + slope2*(t_b-T1);
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_triangle_3(Global &G,XbeUsr &X,XbeJac &J) {
   double y0,t0_new,delta_min,t_a,t_b,tnext_p;
   int n;
   double y;
   double frequency,L1,L2,t0,slope1,slope2,epsl,T1,T2,T;
   const int nvr_y = 0;
   const int nr_frequency = 0;
   const int nr_L1 = 1;
   const int nr_L2 = 2;
   const int nr_t0 = 3;
   const int nr_slope1 = 4;
   const int nr_slope2 = 5;
   const int nr_epsl = 6;
   const int nr_T1 = 7;
   const int nr_T2 = 8;
   const int nr_T = 9;
   const int no_y = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     frequency = X.rprm[nr_frequency];
     T = 1.0/frequency;
     T1 = 0.5*T;
     T2 = T1;

     L1 = X.rprm[nr_L1];
     L2 = X.rprm[nr_L2];

     slope1 = (L2-L1)/T1;
     slope2 = (L1-L2)/T2;
     epsl = T/1000.0;

     X.rprm[nr_slope1] = slope1;
     X.rprm[nr_slope2] = slope2;
     X.rprm[nr_epsl] = epsl;

     X.rprm[nr_T ] = T;
     X.rprm[nr_T1] = T1;
     X.rprm[nr_T2] = T2;

     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   t0   = X.rprm[nr_t0  ];
   T    = X.rprm[nr_T   ];
   T1   = X.rprm[nr_T1  ];
   epsl = X.rprm[nr_epsl];

   if (G.time_given_x < t0) {
     n = ((t0-G.time_given_x)/T) + 1;
     t0_new = t0-n*T;
   } else {
     t0_new = t0;
   }
   t_a = G.time_given_x-t0_new;
   t_b = fmod(t_a,T);

   if (abs(t_b-T) < epsl) t_b = 0.0;

   if (G.flags[G.i_next_time]) {
     if (t_b < T1) {
       tnext_p = T1;
     } else {
       tnext_p = T;
     }
     G.time_nextbreak_x = G.time_given_x + (tnext_p-t_b);
     return;
   }

   L1 = X.rprm[nr_L1];
   L2 = X.rprm[nr_L2];
   slope1 = X.rprm[nr_slope1];
   slope2 = X.rprm[nr_slope2];

   if (t_b < T1) {
     y0 = L1 + slope1*t_b;
   } else {
     y0 = L2 + slope2*(t_b-T1);
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_triangle_4(Global &G,XbeUsr &X,XbeJac &J) {
   double y0,t0_new,delta_min,t_a,t_b,tnext_p;
   int n;
   double y;
   int flag_frequency,flag_period;
   double T,frequency,L1,L2,D,t0,slope1,slope2,epsl,T1,T2;
   const int nvr_y = 0;
   const int ni_flag_frequency = 0;
   const int ni_flag_period = 1;
   const int nr_T = 0;
   const int nr_frequency = 1;
   const int nr_L1 = 2;
   const int nr_L2 = 3;
   const int nr_D = 4;
   const int nr_t0 = 5;
   const int nr_slope1 = 6;
   const int nr_slope2 = 7;
   const int nr_epsl = 8;
   const int nr_T1 = 9;
   const int nr_T2 = 10;
   const int no_y = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     flag_frequency = X.iprm[ni_flag_frequency];
     flag_period    = X.iprm[ni_flag_period   ];

     if ((flag_frequency == 0) && (flag_period == 0)) {
       cout << "triangle_2.xbe: check flag_frequency and flag_period" << endl;
       cout << "  Both cannot be zero." << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     if ((flag_frequency != 0) && (flag_period != 0)) {
       cout << "triangle_2.xbe: check flag_frequency and flag_period" << endl;
       cout << "  Both cannot be non-zero." << endl;
       cout << "  Halting..." << endl; exit(1);
     }
     if (flag_frequency != 0) {
       frequency = X.rprm[nr_frequency];
       T = 1.0/frequency;
       X.rprm[nr_T] = T;
     }
     if (flag_period != 0) {
       T = X.rprm[nr_T];
       frequency = 1.0/T;
       X.rprm[nr_frequency] = frequency;
     }

     T = X.rprm[nr_T];
     D = X.rprm[nr_D];
     T1 = D*T;
     T2 = T - T1;

     L1 = X.rprm[nr_L1];
     L2 = X.rprm[nr_L2];

     slope1 = (L2-L1)/T1;
     slope2 = (L1-L2)/T2;
     epsl = T/1000.0;

     X.rprm[nr_slope1] = slope1;
     X.rprm[nr_slope2] = slope2;
     X.rprm[nr_epsl] = epsl;

     X.rprm[nr_T1] = T1;
     X.rprm[nr_T2] = T2;

     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y] = X.val_vr[nvr_y];
     return;
   }

   t0   = X.rprm[nr_t0  ];
   T    = X.rprm[nr_T   ];
   T1   = X.rprm[nr_T1  ];
   epsl = X.rprm[nr_epsl];

   if (G.time_given_x < t0) {
     n = ((t0-G.time_given_x)/T) + 1;
     t0_new = t0-n*T;
   } else {
     t0_new = t0;
   }
   t_a = G.time_given_x-t0_new;
   t_b = fmod(t_a,T);

   if (abs(t_b-T) < epsl) t_b = 0.0;

   if (G.flags[G.i_next_time]) {
     if (t_b < T1) {
       tnext_p = T1;
     } else {
       tnext_p = T;
     }
     G.time_nextbreak_x = G.time_given_x + (tnext_p-t_b);
     return;
   }

   L1 = X.rprm[nr_L1];
   L2 = X.rprm[nr_L2];
   slope1 = X.rprm[nr_slope1];
   slope2 = X.rprm[nr_slope2];

   if (t_b < T1) {
     y0 = L1 + slope1*t_b;
   } else {
     y0 = L2 + slope2*(t_b-T1);
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_user_fn_1_1(Global &G,XbeUsr &X,XbeJac &J) {
// declare large enough size to serve other user_fn_x_x
// elements as well

   double time0;
   double x_uf[20];
   double y_uf[20];
   double x1;
   double y1;
   int iprm1,iprm2,index_fn;
   double rprm1,rprm2,rprm3,rprm4,rprm5,rprm6,rprm7,rprm8,rprm9,rprm10;
   const int nvr_x1 = 0;
   const int nvr_y1 = 1;
   const int ni_iprm1 = 0;
   const int ni_iprm2 = 1;
   const int ni_index_fn = 2;
   const int nr_rprm1 = 0;
   const int nr_rprm2 = 1;
   const int nr_rprm3 = 2;
   const int nr_rprm4 = 3;
   const int nr_rprm5 = 4;
   const int nr_rprm6 = 5;
   const int nr_rprm7 = 6;
   const int nr_rprm8 = 7;
   const int nr_rprm9 = 8;
   const int nr_rprm10 = 9;
   const int no_y1 = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y1] = X.val_vr[nvr_y1];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     return;
   }
   time0 = G.time_given_x;
   index_fn = X.iprm[ni_index_fn];
   x_uf[0] = X.val_vr[nvr_x1];
   user_function(index_fn,time0,x_uf,y_uf,X.iprm,X.rprm);

   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y1] = y_uf[0];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y1] - y_uf[0];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y1] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_user_fn_2_1(Global &G,XbeUsr &X,XbeJac &J) {
// declare large enough size to serve other user_fn_x_x
// elements as well

   double time0;
   double x_uf[20];
   double y_uf[20];
   double x1,x2;
   double y1;
   int iprm1,iprm2,index_fn;
   double rprm1,rprm2,rprm3,rprm4,rprm5,rprm6,rprm7,rprm8,rprm9,rprm10;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y1 = 2;
   const int ni_iprm1 = 0;
   const int ni_iprm2 = 1;
   const int ni_index_fn = 2;
   const int nr_rprm1 = 0;
   const int nr_rprm2 = 1;
   const int nr_rprm3 = 2;
   const int nr_rprm4 = 3;
   const int nr_rprm5 = 4;
   const int nr_rprm6 = 5;
   const int nr_rprm7 = 6;
   const int nr_rprm8 = 7;
   const int nr_rprm9 = 8;
   const int nr_rprm10 = 9;
   const int no_y1 = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y1] = X.val_vr[nvr_y1];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     return;
   }
   time0 = G.time_given_x;
   index_fn = X.iprm[ni_index_fn];
   x_uf[0] = X.val_vr[nvr_x1];
   x_uf[1] = X.val_vr[nvr_x2];
   user_function(index_fn,time0,x_uf,y_uf,X.iprm,X.rprm);

   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y1] = y_uf[0];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y1] - y_uf[0];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y1] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_user_fn_3_1(Global &G,XbeUsr &X,XbeJac &J) {
// declare large enough size to serve other user_fn_x_x
// elements as well

   double time0;
   double x_uf[20];
   double y_uf[20];
   double x1,x2,x3;
   double y1;
   int iprm1,iprm2,index_fn;
   double rprm1,rprm2,rprm3,rprm4,rprm5,rprm6,rprm7,rprm8,rprm9,rprm10;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_x3 = 2;
   const int nvr_y1 = 3;
   const int ni_iprm1 = 0;
   const int ni_iprm2 = 1;
   const int ni_index_fn = 2;
   const int nr_rprm1 = 0;
   const int nr_rprm2 = 1;
   const int nr_rprm3 = 2;
   const int nr_rprm4 = 3;
   const int nr_rprm5 = 4;
   const int nr_rprm6 = 5;
   const int nr_rprm7 = 6;
   const int nr_rprm8 = 7;
   const int nr_rprm9 = 8;
   const int nr_rprm10 = 9;
   const int no_y1 = 0;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y1] = X.val_vr[nvr_y1];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     return;
   }
   time0 = G.time_given_x;
   index_fn = X.iprm[ni_index_fn];
   x_uf[0] = X.val_vr[nvr_x1];
   x_uf[1] = X.val_vr[nvr_x2];
   x_uf[2] = X.val_vr[nvr_x3];
   user_function(index_fn,time0,x_uf,y_uf,X.iprm,X.rprm);

   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y1] = y_uf[0];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y1] - y_uf[0];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y1] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_user_fn_4_3(Global &G,XbeUsr &X,XbeJac &J) {
// declare large enough size to serve other user_fn_x_x
// elements as well

   double time0;
   double x_uf[20];
   double y_uf[20];
   double x1,x2,x3,x4;
   double y1,y2,y3;
   int iprm1,iprm2,index_fn;
   double rprm1,rprm2,rprm3,rprm4,rprm5,rprm6,rprm7,rprm8,rprm9,rprm10;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_x3 = 2;
   const int nvr_x4 = 3;
   const int nvr_y1 = 4;
   const int nvr_y2 = 5;
   const int nvr_y3 = 6;
   const int ni_iprm1 = 0;
   const int ni_iprm2 = 1;
   const int ni_index_fn = 2;
   const int nr_rprm1 = 0;
   const int nr_rprm2 = 1;
   const int nr_rprm3 = 2;
   const int nr_rprm4 = 3;
   const int nr_rprm5 = 4;
   const int nr_rprm6 = 5;
   const int nr_rprm7 = 6;
   const int nr_rprm8 = 7;
   const int nr_rprm9 = 8;
   const int nr_rprm10 = 9;
   const int no_y1 = 0;
   const int no_y2 = 1;
   const int no_y3 = 2;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y1] = X.val_vr[nvr_y1];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     return;
   }
   time0 = G.time_given_x;
   index_fn = X.iprm[ni_index_fn];
   x_uf[0] = X.val_vr[nvr_x1];
   x_uf[1] = X.val_vr[nvr_x2];
   x_uf[2] = X.val_vr[nvr_x3];
   x_uf[3] = X.val_vr[nvr_x4];
   user_function(index_fn,time0,x_uf,y_uf,X.iprm,X.rprm);

   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y1] = y_uf[0];
       X.val_vr[nvr_y2] = y_uf[1];
       X.val_vr[nvr_y3] = y_uf[2];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y1] - y_uf[0];
         X.g[ng_2] = X.val_vr[nvr_y2] - y_uf[1];
         X.g[ng_3] = X.val_vr[nvr_y3] - y_uf[2];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y1] = 1.0;
         J.dgdvr[ng_2][nvr_y2] = 1.0;
         J.dgdvr[ng_3][nvr_y3] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_user_fn_5_3(Global &G,XbeUsr &X,XbeJac &J) {
// declare large enough size to serve other user_fn_x_x
// elements as well

   double time0;
   double x_uf[20];
   double y_uf[20];
   double x1,x2,x3,x4,x5;
   double y1,y2,y3;
   int iprm1,iprm2,index_fn;
   double rprm1,rprm2,rprm3,rprm4,rprm5,rprm6,rprm7,rprm8,rprm9,rprm10;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_x3 = 2;
   const int nvr_x4 = 3;
   const int nvr_x5 = 4;
   const int nvr_y1 = 5;
   const int nvr_y2 = 6;
   const int nvr_y3 = 7;
   const int ni_iprm1 = 0;
   const int ni_iprm2 = 1;
   const int ni_index_fn = 2;
   const int nr_rprm1 = 0;
   const int nr_rprm2 = 1;
   const int nr_rprm3 = 2;
   const int nr_rprm4 = 3;
   const int nr_rprm5 = 4;
   const int nr_rprm6 = 5;
   const int nr_rprm7 = 6;
   const int nr_rprm8 = 7;
   const int nr_rprm9 = 8;
   const int nr_rprm10 = 9;
   const int no_y1 = 0;
   const int no_y2 = 1;
   const int no_y3 = 2;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_y1] = X.val_vr[nvr_y1];
     return;
   }
   if (G.flags[G.i_init_guess]) {
     return;
   }
   time0 = G.time_given_x;
   index_fn = X.iprm[ni_index_fn];
   x_uf[0] = X.val_vr[nvr_x1];
   x_uf[1] = X.val_vr[nvr_x2];
   x_uf[2] = X.val_vr[nvr_x3];
   x_uf[3] = X.val_vr[nvr_x4];
   x_uf[4] = X.val_vr[nvr_x5];
   user_function(index_fn,time0,x_uf,y_uf,X.iprm,X.rprm);

   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y1] = y_uf[0];
       X.val_vr[nvr_y2] = y_uf[1];
       X.val_vr[nvr_y3] = y_uf[2];
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y1] - y_uf[0];
         X.g[ng_2] = X.val_vr[nvr_y2] - y_uf[1];
         X.g[ng_3] = X.val_vr[nvr_y3] - y_uf[2];
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y1] = 1.0;
         J.dgdvr[ng_2][nvr_y2] = 1.0;
         J.dgdvr[ng_3][nvr_y3] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_vsi_3ph_1(Global &G,XbeUsr &X,XbeJac &J) {
   double val_a,val_b,val_c;
   double g1,g2,g3,g4,g5,g6;
   double va,vb,vc;
   double vdc,L,Lby2;
   const int nvr_g1 = 0;
   const int nvr_g2 = 1;
   const int nvr_g3 = 2;
   const int nvr_g4 = 3;
   const int nvr_g5 = 4;
   const int nvr_g6 = 5;
   const int nvr_va = 6;
   const int nvr_vb = 7;
   const int nvr_vc = 8;
   const int nr_vdc = 0;
   const int nr_L = 1;
   const int nr_Lby2 = 2;
   const int no_va = 0;
   const int no_vb = 1;
   const int no_vc = 2;
   const int no_g1 = 3;
   const int no_g2 = 4;
   const int no_g3 = 5;
   const int no_g4 = 6;
   const int no_g5 = 7;
   const int no_g6 = 8;
   const int ng_1 = 0;
   const int ng_2 = 1;
   const int ng_3 = 2;
   if (G.flags[G.i_one_time_parms]) {
     L = X.rprm[nr_L];
     Lby2 = 0.5*L;
     X.rprm[nr_Lby2] = Lby2;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_va] = X.val_vr[nvr_va];
     X.outprm[no_vb] = X.val_vr[nvr_vb];
     X.outprm[no_vc] = X.val_vr[nvr_vc];
     X.outprm[no_g1] = X.val_vr[nvr_g1];
     X.outprm[no_g2] = X.val_vr[nvr_g2];
     X.outprm[no_g3] = X.val_vr[nvr_g3];
     X.outprm[no_g4] = X.val_vr[nvr_g4];
     X.outprm[no_g5] = X.val_vr[nvr_g5];
     X.outprm[no_g6] = X.val_vr[nvr_g6];
     return;
   }

   g1 = X.val_vr[nvr_g1];
   g2 = X.val_vr[nvr_g2];
   g3 = X.val_vr[nvr_g3];
   g4 = X.val_vr[nvr_g4];
   g5 = X.val_vr[nvr_g5];
   g6 = X.val_vr[nvr_g6];

   vdc = X.rprm[nr_vdc];
   Lby2 = X.rprm[nr_Lby2];

   if (g1 > Lby2) {
     if (g4 > Lby2) {
//     val_a = 0.5*vdc;
       val_a = 0.0;
     } else {
//     val_a = vdc;
       val_a = 0.5*vdc;
     }
   } else {
     if (g4 > Lby2) {
//     val_a = 0.0;
       val_a = -0.5*vdc;
     } else {
//     val_a = 0.5*vdc;
       val_a = 0.0;
     }
   }
   if (g3 > Lby2) {
     if (g6 > Lby2) {
//     val_b = 0.5*vdc;
       val_b = 0.0;
     } else {
//     val_b = vdc;
       val_b = 0.5*vdc;
     }
   } else {
     if (g6 > Lby2) {
//     val_b = 0.0;
       val_b = -0.5*vdc;
     } else {
//     val_b = 0.5*vdc;
       val_b = 0.0;
     }
   }
   if (g5 > Lby2) {
     if (g2 > Lby2) {
//     val_c = 0.5*vdc;
       val_c = 0.0;
     } else {
//     val_c = vdc;
       val_c = 0.5*vdc;
     }
   } else {
     if (g2 > Lby2) {
//     val_c = 0.0;
       val_c = -0.5*vdc;
     } else {
//     val_c = 0.5*vdc;
       val_c = 0.0;
     }
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_va] = val_a;
     X.val_vr[nvr_vb] = val_b;
     X.val_vr[nvr_vc] = val_c;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_va] = val_a;
       X.val_vr[nvr_vb] = val_b;
       X.val_vr[nvr_vc] = val_c;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_va] - val_a;
         X.g[ng_2] = X.val_vr[nvr_vb] - val_b;
         X.g[ng_3] = X.val_vr[nvr_vc] - val_c;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_va] = 1.0;
         J.dgdvr[ng_2][nvr_vb] = 1.0;
         J.dgdvr[ng_3][nvr_vc] = 1.0;
       }
     }
     return;
   }
   return;
}
void x_xfer_fn(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   double y;
   int scale_coef;
   double a0,a1,a2,a3,a4,a5,b0,b1,b2,b3,b4,b5,f0;
   double y_sv;
   const int nvr_x = 0;
   const int nvr_y = 1;
   const int ni_scale_coef = 0;
   const int nr_a0 = 0;
   const int nr_a1 = 1;
   const int nr_a2 = 2;
   const int nr_a3 = 3;
   const int nr_a4 = 4;
   const int nr_a5 = 5;
   const int nr_b0 = 6;
   const int nr_b1 = 7;
   const int nr_b2 = 8;
   const int nr_b3 = 9;
   const int nr_b4 = 10;
   const int nr_b5 = 11;
   const int nr_f0 = 12;
   const int nst_y_sv = 0;
   return;
}
void x_xor_2(Global &G,XbeUsr &X,XbeJac &J) {
   double y0;
   bool X1,X2,Y;
   double x1,x2;
   double y;
   double y_high,hb2;
   const int nvr_x1 = 0;
   const int nvr_x2 = 1;
   const int nvr_y = 2;
   const int nr_y_high = 0;
   const int nr_hb2 = 1;
   const int no_x1 = 0;
   const int no_x2 = 1;
   const int no_y = 2;
   const int ng_1 = 0;
   if (G.flags[G.i_one_time_parms]) {
     y_high = X.rprm[nr_y_high];
     hb2 = 0.5*y_high;
     X.rprm[nr_hb2] = hb2;
     return;
   }
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x1] = X.val_vr[nvr_x1];
     X.outprm[no_x2] = X.val_vr[nvr_x2];
     X.outprm[no_y ] = X.val_vr[nvr_y ];
     return;
   }
   y_high = X.rprm[nr_y_high];
   hb2 = X.rprm[nr_hb2];

   x1 = X.val_vr[nvr_x1];
   x2 = X.val_vr[nvr_x2];
   X1 = x1 > hb2;
   X2 = x2 > hb2;

   Y = X1 != X2;

   if (Y) {
     y0 = y_high;
   } else {
     y0 = 0.0;
   }

   if (G.flags[G.i_init_guess]) {
     X.val_vr[nvr_y] = y0;
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     if (G.flags[G.i_explicit]) {
       X.val_vr[nvr_y] = y0;
     } else if (G.flags[G.i_implicit]) {
       if (G.flags[G.i_function]) {
         X.g[ng_1] = X.val_vr[nvr_y] - y0;
       }
       if (G.flags[G.i_jacobian]) {
         J.dgdvr[ng_1][nvr_y ] =  1.0;
       }
     }
     return;
   }
   return;
}
void x_x_dummy(Global &G,XbeUsr &X,XbeJac &J) {
   double x;
   const int nvr_x = 0;
   const int no_x = 0;
   if (G.flags[G.i_outvar]) {
     X.outprm[no_x] = X.val_vr[nvr_x];
     return;
   }
   if (G.flags[G.i_one_time_parms]) {
     return;
   }
   if (G.flags[G.i_init_guess]) {
     return;
   }
   if (G.flags[G.i_trns] || G.flags[G.i_startup]) {
     return;
   }
   if (G.flags[G.i_outvar]) {
     return;
   }
   return;
}
