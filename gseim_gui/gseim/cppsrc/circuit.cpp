/*
Copyright (C) 2021 - Mahesh Patil <mbpatil@ee.iitb.ac.in>
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

#include "circuit.h"

Circuit::Circuit(){}

// -----------------------------------------------------------------------------
Circuit::Circuit(
   const std::string &filename,
   const vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Global &global,
   CctFile &cct_file) {

// filename: name of the circuit file

   std::fstream inf;
   vector<std::string> v1; 
   std::string s1,s2,s3;
   int i_xbel,pos;
   int i_ebel,i_ebeu;
   bool flag_1;
   int xbeu_aux0,ebeu_aux0,ebeu_auxs0,ebeu_stv0;

   cout << "Circuit: filename = " << filename << endl;

   set_flags_default();
   assign_const_1<bool>(cct_file.tick_cf,false);

   n_xbeu = 0;
   n_xbeu_vr = 0;
   xbeu_vr_name.clear();

   inf.open(filename,ios::in|ios::binary);
   for (int i_line=0; i_line < cct_file.line_end_circuit; i_line++) {
     next_line(inf,v1);
//   Note: xelement denotes Xbe
     if (i_line > cct_file.line_begin_circuit) {
       if (v1[0] == "xelement") {
         s1 = assign_string_1(v1,1,2,"type");
         i_xbel = find_name<XbeLib>(xbe_lib,s1);

         for (int i=0; i < xbe_lib[i_xbel].n_vr; i++) {
           s2 = xbe_lib[i_xbel].vr_name[i];
           s3 = assign_string_1(v1,1,2,s2);
           find_word_1(xbeu_vr_name,s3,pos,flag_1);
           if (!flag_1) {
             xbeu_vr_name.push_back(s3);
             n_xbeu_vr++;
           }
         }
         n_xbeu++;
       }
     }
   }
   inf.close();

   flag_x = (n_xbeu > 0);

   cout << "Circuit: n_xbeu_vr = " << n_xbeu_vr << endl;
   map_xbeuvr_to_svec.resize(n_xbeu_vr);

   xbe_usr.resize(n_xbeu);
   xbeu_aux_start.resize(n_xbeu);
   n_xbeu = 0;
   xbeu_aux0 = 0;

   inf.open(filename,ios::in|ios::binary);
   for (int i_line=0; i_line < cct_file.line_end_circuit; i_line++) {
     next_line(inf,v1);
     if (i_line > cct_file.line_begin_circuit) {
       if (v1[0] == "xelement") {
         xbeu_aux_start[n_xbeu] = xbeu_aux0;
         xbe_usr[n_xbeu].set_values_1(n_xbeu,v1,xbeu_vr_name,
           xbeu_aux0,xbe_lib,global,cct_file);
         i_xbel = xbe_usr[n_xbeu].index_xbel;
         n_xbeu++;
         xbeu_aux0 += xbe_lib[i_xbel].n_aux;
       }
     }
   }
   inf.close();

   process_xbeu_1(xbe_lib,xbe_usr);

   cct_file.n_ov_xbe = 0;

   for (int i=0; i < cct_file.n_ov; i++) {
     if (cct_file.ov_flag[i] == global.I_OV_XBE) cct_file.n_ov_xbe++;
   }

   cct_file.n_ov_xvr   = 0;

   for (int i=0; i < cct_file.n_ov; i++) {
     if (cct_file.ovl_name[i] == "xvar") {
       cct_file.ov_flag[i] = global.I_OV_XVR;
       find_word_1(xbeu_vr_name,cct_file.ovr_name[i],pos,flag_1);
       if (!flag_1) {
         cout << "Circuit: xvar <" << cct_file.ovr_name[i] << "> not found" << endl;
         cout << "  in outvar statement. Halting..." << endl; exit(1);
       }
       cct_file.ov1[i] = pos;
       cct_file.ov2[i] = 0;
       cct_file.n_ov_xvr++;
     }
   }

// check that all outvars have been assigned:
   check_vec_const_2<int>(cct_file.ov1,-1,flag_1,pos);
   if (flag_1) {
     cout << "Circuit: ov1[" << pos << "] = -1?" << endl;
     cout << "  Check this outvar in the circuit file: "
       << cct_file.ov_name[pos] << endl;
     cout << "  Halting..." << endl; exit(1);
   }
   check_vec_const_2<int>(cct_file.ov2,-1,flag_1,pos);
   if (flag_1) {
     cout << "Circuit: ov2[" << pos << "] = -1?" << endl;
     cout << "  Check this outvar in the circuit file: "
       << cct_file.ov_name[pos] << endl;
     cout << "  Halting..." << endl; exit(1);
   }

   return;
} //end of Circuit::Circuit
// -----------------------------------------------------------------------------
void Circuit::check_save_history(
   const vector<XbeLib> &xbe_lib,
   const vector<XbeUsr> &xbe_usr) {

   int i_xbeu,i_xbel;

   flag_save_history_x = false;

   for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_savehist) {
       flag_save_history_x = true; break;
     }
   }
   flag_save_history = flag_save_history_x;

   return;
} //end of Circuit::check_save_history
// -----------------------------------------------------------------------------
void Circuit::check_reset_x(
   const vector<XbeLib> &xbe_lib,
   const vector<XbeUsr> &xbe_usr) {

   int i_xbeu,i_xbel;
   flag_reset_x = false;

   for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_reset) {
       flag_reset_x = true;
       break;
     }
   }

   return;
} //end of Circuit::check_reset_x
// -----------------------------------------------------------------------------
void Circuit::check_modulo_x(
   const vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr) {

   int i_xbeu,i_xbel;

   flag_modulo_x = false;
   nttl_xbeu_modulo = 0;
   xbeu_modulo_map.clear();

   for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_modulo) {
       flag_modulo_x = true;
       xbeu_modulo_map.push_back(i_xbeu);
       nttl_xbeu_modulo++;
     }
   }

   return;
} //end of Circuit::check_modulo_x
// -----------------------------------------------------------------------------
void Circuit::check_limit_tstep(
   const vector<XbeLib> &xbe_lib,
   const vector<XbeUsr> &xbe_usr) {

   int i_xbeu,i_xbel;

   flag_limit_tstep_x = false;

   for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_lmttstep) {
       flag_limit_tstep_x = true; break;
     }
   }
   flag_limit_tstep = flag_limit_tstep_x;

   return;
} //end of Circuit::check_limit_tstep
// -----------------------------------------------------------------------------
void Circuit::xbe_map_vr(
   const vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr) {

// Assign map_vr_1,map_vr_2
// used for handling modulo elements

   int i_xbeu,i_xbeu_vr,i_vr,i_xbel;
   vector<int> map_vr_1_1d;
   vector<int> map_vr_2_1d;

   n_map_vr.resize(n_xbeu_vr);
   map_vr_1.clear();
   map_vr_2.clear();

   for (i_xbeu_vr=0; i_xbeu_vr < n_xbeu_vr; i_xbeu_vr++) {
     n_map_vr[i_xbeu_vr] = 0;
     map_vr_1_1d.clear();
     map_vr_2_1d.clear();

     for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       for (i_vr=0; i_vr < xbe_lib[i_xbel].n_vr; i_vr++) {
         if (xbe_usr[i_xbeu].vr[i_vr] == i_xbeu_vr) {
           n_map_vr[i_xbeu_vr]++;
           map_vr_1_1d.push_back(i_xbeu);
           map_vr_2_1d.push_back(i_vr);
         }
       }
     }
     map_vr_1.push_back(map_vr_1_1d);
     map_vr_2.push_back(map_vr_2_1d);
   }
   return;
} // end of Circuit::xbe_map_vr
// -----------------------------------------------------------------------------
void Circuit::process_xbeu_1(
   const vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr) {

   int i_xbel,i_xbeu,i_xbeu_vr;
   bool flag_all_visited;
   int n_begin,n_end;
   int n_pass_max = 20;
   bool flag_1;
   int n_beu_1;

   x_vr_visited.resize(n_xbeu_vr);
   x_beu_visited.resize(n_xbeu);
   x_pass.resize(n_xbeu);

   assign_const_1<bool>(x_vr_visited,false);
   assign_const_1<bool>(x_beu_visited,false);
   assign_const_1<int>(x_pass,-1);

   flag_alg_loop = false;

// Step 1: mark output nodes of all xbeu's of type integrate as visited.

   for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     if (xbe_lib[i_xbel].flag_integrate) {
       n_begin = xbe_lib[i_xbel].n_ipvr;
       n_end = n_begin + xbe_lib[i_xbel].n_opvr;

       for (int i_vr=n_begin; i_vr < n_end; i_vr++) {
         i_xbeu_vr = xbe_usr[i_xbeu].vr[i_vr];
         x_vr_visited[i_xbeu_vr] = true;
       }
     }
   }

// Step 2: mark output nodes of all xbeu's of type delay as visited.

   for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     if (xbe_lib[i_xbel].flag_delay) {
       n_begin = xbe_lib[i_xbel].n_ipvr;
       n_end = n_begin + xbe_lib[i_xbel].n_opvr;

       for (int i_vr=n_begin; i_vr < n_end; i_vr++) {
         i_xbeu_vr = xbe_usr[i_xbeu].vr[i_vr];
         x_vr_visited[i_xbeu_vr] = true;
       }
     }
   }

// Step 3: mark output nodes of all xbeu's of type evaluate-SRC as visited.

   for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     if (xbe_lib[i_xbel].flag_evaluate) {
       if (xbe_lib[i_xbel].flag_source) {
         n_begin = xbe_lib[i_xbel].n_ipvr;
         n_end = n_begin + xbe_lib[i_xbel].n_opvr;

         for (int i_vr=n_begin; i_vr < n_end; i_vr++) {
           i_xbeu_vr = xbe_usr[i_xbeu].vr[i_vr];
           x_vr_visited[i_xbeu_vr] = true;
         }
       }
     }
   }

// Step 4: make passes of xbeu's of type evaluate-NONSRC.
// We will handle source type xbeu's separately. They simply need
// to be processed (in any order) before anything else.

   x_n_pass = 0;

   for (int i_pass=0; i_pass < n_pass_max; i_pass++) {
     for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;

       if (xbe_lib[i_xbel].flag_evaluate) {
         if (!xbe_lib[i_xbel].flag_source) {
           if (!x_beu_visited[i_xbeu]) {
             flag_1 = true;
             n_end = xbe_lib[i_xbel].n_ipvr;
             for (int i_vr=0; i_vr < n_end; i_vr++) {
               i_xbeu_vr = xbe_usr[i_xbeu].vr[i_vr];
               if (!x_vr_visited[i_xbeu_vr]) {
                 flag_1 = false; break;
               }
             }
//           check if this xbeu can be ticked.
             if (flag_1) {
               x_n_pass = i_pass + 1;
               x_beu_visited[i_xbeu] = true;
               x_pass[i_xbeu] = i_pass;
             }
           }
         }
       }
     }

//   The current pass is over.
//   Mark output nodes of elements visited as visited nodes.
//   (This will also mark nodes already marked; never mind)
     for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;

       if (xbe_lib[i_xbel].flag_evaluate) {
         if (!xbe_lib[i_xbel].flag_source) {
           if (x_beu_visited[i_xbeu]) {
             n_begin = xbe_lib[i_xbel].n_ipvr;
             n_end = n_begin + xbe_lib[i_xbel].n_opvr;

             for (int i_vr=n_begin; i_vr < n_end; i_vr++) {
               i_xbeu_vr = xbe_usr[i_xbeu].vr[i_vr];
               x_vr_visited[i_xbeu_vr] = true;
             }
           }
         }
       }
     }

//   Check if all eval/non-source xbeu's have been visited.

     flag_all_visited = true;
     for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;

       if (xbe_lib[i_xbel].flag_evaluate) {
         if (!xbe_lib[i_xbel].flag_source) {
           if (!x_beu_visited[i_xbeu]) {
             flag_all_visited = false; break;
           }
         }
       }
     }
     if (flag_all_visited) break;
   }

   if (!flag_all_visited) {
     cout << "Circuit::process_xbeu_1: ordering failed even after" << endl;
     cout << "  completing " << n_pass_max << " passes." << endl;
     cout << "  -> Algebraic loop." << endl;

     flag_alg_loop = true;

//   In this case, it is not requred to check whether the xbe nodes
//   have been visited or not (that information will not be used anyway),
//   so skip that part:

     goto jump1;
   }

// Check if the input vars of integrator-type xbeu's have been visited.

   for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     if (xbe_lib[i_xbel].flag_integrate) {
       n_end = xbe_lib[i_xbel].n_ipvr;
       for (int i_vr=0; i_vr < n_end; i_vr++) {
         i_xbeu_vr = xbe_usr[i_xbeu].vr[i_vr];

         if (!x_vr_visited[i_xbeu_vr]) {
           cout << "Circuit::process_xbeu_1: an input var of an integrator type" << endl;
           cout << "  xbeu has not been visited." << endl;
           cout << "  i_xbeu_vr=" << i_xbeu_vr << endl;
           flag_alg_loop = true;
           goto jump1;
         }
       }
     }
   }

// Check if the input vars of delay-type xbeu's have been visited.

   for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     if (xbe_lib[i_xbel].flag_delay) {
       n_end = xbe_lib[i_xbel].n_ipvr;
       for (int i_vr=0; i_vr < n_end; i_vr++) {
         i_xbeu_vr = xbe_usr[i_xbeu].vr[i_vr];

         if (!x_vr_visited[i_xbeu_vr]) {
           cout << "Circuit::process_xbeu_1: an input var of a delay type" << endl;
           cout << "  xbeu has not been visited." << endl;
           cout << "  i_xbeu_vr=" << i_xbeu_vr << endl;
           flag_alg_loop = true;
           goto jump1;
         }
       }
     }
   }

// gather elements (eval and non-src) belonging to each pass:

   x_pass_n_beu.resize(x_n_pass);
   x_pass_beu.resize(x_n_pass);

   assign_const_1<int>(x_pass_n_beu,0);

// assign pass_beu[][]:

   for (int i_pass=0; i_pass < x_n_pass; i_pass++) {
     n_beu_1 = 0;
     for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;

       if (xbe_lib[i_xbel].flag_evaluate) {
         if (!xbe_lib[i_xbel].flag_source) {
           if (x_pass[i_xbeu] == i_pass) {
             x_pass_beu[i_pass].push_back(i_xbeu);
             n_beu_1++;
             x_pass_n_beu[i_pass] = n_beu_1;
           }
         }
       }
     }
   }
   jump1:;

   cout << "Circuit::process_xbeu_1: " << endl;
   cout << "  flag_alg_loop = " << flag_alg_loop << endl;

// compute n_intgrtr,n_delay,n_eval_src,n_eval_nonsrc.

   x_n_intgrtr     = 0;
   x_n_delay       = 0;
   x_n_eval_src    = 0;
   x_n_eval_nonsrc = 0;

   x_n_ttlaux = 0;
   x_n_ttlg = 0;

   for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     if (xbe_lib[i_xbel].flag_integrate) {
       x_n_intgrtr++;
     }
     if (xbe_lib[i_xbel].flag_delay) {
       x_n_delay++;
     }
     if (xbe_lib[i_xbel].flag_evaluate) {
       if (xbe_lib[i_xbel].flag_source) {
         x_n_eval_src++;
       } else {
         x_n_eval_nonsrc++;
       }
     }
     x_n_ttlaux += xbe_lib[i_xbel].n_aux;
     x_n_ttlg   += xbe_lib[i_xbel].n_g;
   }

   val_xvr.resize(n_xbeu_vr);
   val_xvr_new.resize(n_xbeu_vr);

   return;
} //end of Circuit::process_xbeu_1
// -----------------------------------------------------------------------------
void Circuit::assign_flag_linear_x(
   const vector<XbeLib> &xbe_lib,
   const vector<XbeUsr> &xbe_usr) {

   int i_xbeu,i_xbel;

   flag_linear_x = true;

   for (i_xbeu=0; i_xbeu < n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (!xbe_lib[i_xbel].flag_jac_const) {
       flag_linear_x = false;
       break;
     }
   }
   return;
} //end of Circuit::assign_flag_linear_x
// -----------------------------------------------------------------------------
void Circuit::set_flags_default() {

   flag_save_history_x = false;
   flag_save_history   = false;
   flag_reset_x        = false;
   flag_modulo_x       = false;
   flag_limit_tstep_x  = false;
   flag_limit_tstep    = false;
   flag_limit_newton_x = false;
   flag_x              = false;
   flag_x_only         = false;
   flag_x_explicit     = false;
   flag_x_matrix       = false;
   flag_linear_x       = false;

   return;
} //end of Circuit::set_flags_default
