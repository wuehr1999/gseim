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

using namespace std;

void user_f_1(
   double* x,
   double* y,
   int* iprm,
   double* rprm) {

   double a;

   if (iprm[0] == 0) {
     a = 0.5;
   } else {
     a =2.0;
   }

   y[0] = a*rprm[0]*x[0] + rprm[1];
   return;
}

void user_f_eff_d(
   double* x,
   double* y,
   int* iprm,
   double* rprm) {

   double sigma,ls,vdc;
   double k;
   double wmr,isq,eff_d;

   sigma = rprm[0];
   ls    = rprm[1];  
   vdc   = rprm[2];  

   k = sigma*ls/(0.5*vdc);

// cout << "eff_d: k=" << k
//   << " sigma=" << sigma
//   << " ls=" << ls
//   << " vdc=" << vdc
//   << endl;

   wmr = x[0];
   isq = x[1];

   eff_d = k*wmr*isq;
   y[0] = eff_d;

// cout << "eff_d: y[0]=" << y[0] << endl;
   return;
}

void user_f_eff_q(
   double* x,
   double* y,
   int* iprm,
   double* rprm) {

   double sigma,ls,vdc;
   double k1,k2;
   double eff_q;
   double wmr,isd,imr;

   sigma = rprm[0];
   ls    = rprm[1];  
   vdc   = rprm[2];  

   k1 = sigma*ls/(0.5*vdc);
   k2 = (1.0-sigma)*ls/(0.5*vdc);

// cout << "eff_q: k1=" << k1
//   << " k2=" << k2
//   << " sigma=" << sigma
//   << " ls=" << ls
//   << " vdc=" << vdc
//   << endl;

   wmr = x[0];
   isd = x[1];
   imr = x[2];

   eff_q = (k1*wmr*isd + k2*wmr*imr);
   y[0] = eff_q;
// cout << "eff_q: y[0]=" << y[0] << endl;
   return;
}

void user_f_slip(
   double* x,
   double* y,
   int* iprm,
   double* rprm) {

   int poles;
   double tr;
   double w2,isq,imr,wrm,wmr;

   poles = iprm[0];
   tr    = rprm[0];

   isq = x[0];
   imr = x[1];
   wrm = x[2];

   if (fabs(imr) < 1.0e-2) {
     w2 = 1.0;
   } else {
     w2 = isq/(tr*imr);
   }

   wmr = 0.5*((double)(poles))*wrm + w2;

   y[0] = wmr;
   return;
}

void user_f_mux_1(
   const double time0,
   double* x,
   double* y,
   int* iprm,
   double* rprm) {

   double x1,x2,t0;

   x1 = x[0];
   x2 = x[1];
   t0 = rprm[0];

   if (time0 <= t0) {
     y[0] = x1;
     y[1] = 0.0;
   } else {
     y[0] = 0.0;
     y[1] = x2;
   }

   return;
}

void user_function(
   const int index_fn,
   const double time0,
   double* x,
   double* y,
   int* iprm,
   double* rprm) {

// allow up to 2 iprms and 10 rprms

   double a = 2.0;

   switch(index_fn)
   {
   case 1 :
//    simple example:
//    y[0] = a*rprm[0]*x[0] + rprm[1];

      user_f_1(x,y,iprm,rprm);
      break;
   case 2 :
      user_f_eff_d(x,y,iprm,rprm);
      break;
   case 3 :
      user_f_eff_q(x,y,iprm,rprm);
      break;
   case 4 :
      user_f_slip(x,y,iprm,rprm);
      break;
   case 5 :
      user_f_mux_1(time0,x,y,iprm,rprm);
      break;
   case 6 :
      cout << "user_function: index_fn=6 is not implemented. Halting.." << endl; exit (1);
      break;
   case 7 :
      cout << "user_function: index_fn=7 is not implemented. Halting.." << endl; exit (1);
      break;
   case 8 :
      cout << "user_function: index_fn=8 is not implemented. Halting.." << endl; exit (1);
      break;
   case 9 :
      cout << "user_function: index_fn=9 is not implemented. Halting.." << endl; exit (1);
      break;
   case 10 :
      cout << "user_function: index_fn=10 is not implemented. Halting.." << endl; exit (1);
      break;
   default :
      cout << "user_function: Check index_fn. Halting..." << endl;
      exit (1);
   }
   return;
}
