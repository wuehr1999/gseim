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

#include "routines1.h"

// -----------------------------------------------------------------------------
void clear_solvec_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SysMat &smat) {

   int i_xbeu,i_xbel,n_vr1,n_aux1;

   for (int i=0; i < cct.n_xbeu_vr; i++) {
    cct.val_xvr[i] = 0.0;
   }

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     n_vr1 = xbe_lib[i_xbel].n_vr;
     for (int i=0; i < n_vr1; i++) {
       xbe_usr[i_xbeu].val_vr[i] = 0.0;
       xbe_usr[i_xbeu].val_vr_new[i] = 0.0;
     }

     n_aux1 = xbe_lib[i_xbel].n_aux;
     for (int i=0; i < n_aux1; i++) {
       xbe_usr[i_xbeu].val_aux[i] = 0.0;
       xbe_usr[i_xbeu].val_aux_new[i] = 0.0;
     }
   }
   for (int i=0; i < smat.n_solvec_x; i++) {
     smat.svec_x[i] = 0.0;
   }

   return;
} // end of clear_solvec_x
// -----------------------------------------------------------------------------
void clear_solvec_x_1(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct) {

   int i_xbeu,i_xbel,n_vr1,n_aux1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     n_vr1 = xbe_lib[i_xbel].n_vr;
     for (int i=0; i < n_vr1; i++) {
       xbe_usr[i_xbeu].val_vr_0[i] = 0.0;
     }
     n_aux1 = xbe_lib[i_xbel].n_aux;
     for (int i=0; i < n_aux1; i++) {
       xbe_usr[i_xbeu].val_aux_0[i] = 0.0;
     }
   }

   return;
} // end of clear_solvec_x_1
// -----------------------------------------------------------------------------
void cct_to_xbe_1(
   const int i_xbeu,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct) {

// assign xbe_usr.val_vr from cct.val_xvr

   int i_xbel,i_xbeu_vr,n_vr1;

   i_xbel = xbe_usr[i_xbeu].index_xbel;

   n_vr1 = xbe_lib[i_xbel].n_vr;
   for (int i=0; i < n_vr1; i++) {
     i_xbeu_vr = xbe_usr[i_xbeu].vr[i];
     xbe_usr[i_xbeu].val_vr[i] = cct.val_xvr[i_xbeu_vr];
   }

   return;
} // end of cct_to_xbe_1
// -----------------------------------------------------------------------------
void cct_to_xbe_all(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct) {

   int i_xbeu,i_xbel,i_xbeu_vr,n_vr1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     n_vr1 = xbe_lib[i_xbel].n_vr;
     for (int i=0; i < n_vr1; i++) {
       i_xbeu_vr = xbe_usr[i_xbeu].vr[i];
       xbe_usr[i_xbeu].val_vr[i] = cct.val_xvr[i_xbeu_vr];
     }
   }
   return;
} // end of cct_to_xbe_all
// -----------------------------------------------------------------------------
void xbe_to_cct_op_1(
   const int i_xbeu,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct) {

   int i_xbel,i_xbeu_vr,n_ipvr1,n_vr1;

   i_xbel = xbe_usr[i_xbeu].index_xbel;
   n_vr1 = xbe_lib[i_xbel].n_vr;
   n_ipvr1 = xbe_lib[i_xbel].n_ipvr;

   for (int i = n_ipvr1; i < n_vr1; i++) {
     i_xbeu_vr = xbe_usr[i_xbeu].vr[i];
     cct.val_xvr[i_xbeu_vr] = xbe_usr[i_xbeu].val_vr[i];
   }

   return;
} // end of xbe_to_cct_op_1
// -----------------------------------------------------------------------------
void xbe_to_cct_1(
   const int i_xbeu,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct) {

   int i_xbel,i_xbeu_vr,n_vr1;

   i_xbel = xbe_usr[i_xbeu].index_xbel;
   n_vr1 = xbe_lib[i_xbel].n_vr;

   for (int i=0; i < n_vr1; i++) {
     i_xbeu_vr = xbe_usr[i_xbeu].vr[i];
     cct.val_xvr[i_xbeu_vr] = xbe_usr[i_xbeu].val_vr[i];
   }

   return;
} // end of xbe_to_cct_1
// -----------------------------------------------------------------------------
void xbe_to_cct_all(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct) {

   int i_xbeu,i_xbel,i_xbeu_vr,n_vr1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     n_vr1 = xbe_lib[i_xbel].n_vr;

     for (int i=0; i < n_vr1; i++) {
       i_xbeu_vr = xbe_usr[i_xbeu].vr[i];
       cct.val_xvr[i_xbeu_vr] = xbe_usr[i_xbeu].val_vr[i];
     }
   }
   return;
} // end of xbe_to_cct_all
// -----------------------------------------------------------------------------
void get_xbe_1(
   const int i_xbeu,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global) {

   int i_xbel;

   i_xbel = xbe_usr[i_xbeu].index_xbel;

// assign xbe_usr.val_vr from cct.val_xvr
   cct_to_xbe_1(i_xbeu,xbe_lib,xbe_usr,cct);

   get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbeu]);

// assign cct.val_xvr from xbe_usr.val_vr only for the output nodes
   xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);

   return;
} // end of get_xbe_1
// -----------------------------------------------------------------------------
void one_time_parms_x(
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;

   global.flags[global.i_one_time_parms] = true;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbel]);
   }

   global.flags[global.i_one_time_parms] = false;

   return;
} // end of one_time_parms_x
// -----------------------------------------------------------------------------
void save_history_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;

   global.flags[global.i_save_history] = true;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_savehist) {
       get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbel]);
     }
   }
   global.flags[global.i_save_history] = false;

   return;
} // end of save_history_x
// -----------------------------------------------------------------------------
void init_sol_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   SysMat &smat,
   Global &global) {

   if (slv.flag_read_solution) {
     read_solution(xbe_lib,xbe_usr,slv,cct);
   } else if (slv.flag_init_solution) {
     assign_const_1<bool>(global.flags,false);
     clear_solvec_x(xbe_lib,xbe_usr,cct,smat);
     clear_solvec_x_1(xbe_lib,xbe_usr,cct);
     global.time_given_x = global.time_begin;
     xbe_init_guess(xbe_lib,xbe_usr,xbe_jac,cct,global);
   } else if (slv.flag_prev_solution) {
     if (!slv.flag_prev_solution_exists) {
       cout << "init_sol_x: previous solution does not exist?" << endl;
       cout << "   check initial_sol statement. Halting..." << endl;
       exit(1);
     }
   } else {
     clear_solvec_x(xbe_lib,xbe_usr,cct,smat);
   }

   return;
} // end of init_sol_x
// -----------------------------------------------------------------------------
void xbe_init_guess(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;

   assign_const_1<bool>(global.flags,false);
   global.flags[global.i_init_guess] = true;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       get_xbe_1(i_xbeu,xbe_lib,xbe_usr,xbe_jac,cct,global);
     }
   }
   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_evaluate) {
       if (xbe_lib[i_xbel].flag_source) {
         get_xbe_1(i_xbeu,xbe_lib,xbe_usr,xbe_jac,cct,global);
       }
     }
   }

   for (int i_pass=0; i_pass < cct.x_n_pass; i_pass++) {
     for (int i=0; i < cct.x_pass_n_beu[i_pass]; i++) {
       i_xbeu = cct.x_pass_beu[i_pass][i];
       get_xbe_1(i_xbeu,xbe_lib,xbe_usr,xbe_jac,cct,global);
     }
   }

   global.flags[global.i_init_guess] = false;

   return;
} // end of xbe_init_guess
// -----------------------------------------------------------------------------
void read_solution(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SolveBlocks &slv,
   Circuit &cct) {

   fstream inf;
   std::string s1;
   int i_xbeu,i_xbel;

   inf.open(slv.infile_sol,ios::in|ios::binary);

   if (!inf) {
     cout << "read_solution: file " << slv.infile_sol
       << " does not exist. Halting..." << endl; exit(1);
   }

   s1 = next_string_1(inf);
   if (s1 != "cct.val_vr:") {
     cout << "expected cct.val_xvr: in file " << slv.infile_sol << endl;
     cout << "  but found " << s1 << ". Halting..." << endl; exit(1);
   }
   read_vec_double_1(inf,cct.val_xvr);

   s1 = next_string_1(inf);
   if (s1 != "xbe_usr.val_aux:") {
     cout << "expected xbe_usr.val_aux: in file " << slv.infile_sol << endl;
     cout << "  but found " << s1 << ". Halting..." << endl; exit(1);
   }
   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     for (int i=0; i < xbe_lib[i_xbel].n_aux; i++) {
       xbe_usr[i_xbeu].val_aux[i] = next_double_1(inf);
     }
   }

   inf.close();

   return;
} // end of read_solution
// -----------------------------------------------------------------------------
void write_iter_x(
   SolveBlocks &slv) {

   if (slv.iter_trns_x == 0) {
     cout << "Transient simulation starts..." << endl;
     if (slv.flag_write_time_x) {
       cout << "i=" << slv.iter_trns_x
            << " time=" << slv.time_present_x << endl;
     } else {
       cout << "i=" << slv.iter_trns_x << endl;
     }
   } else {
     if (slv.write_iter_n1_x == slv.write_iter_n_x) {
       if (slv.flag_write_time_x) {
         cout << "i=" << slv.iter_trns_x
              << ", time=" << slv.time_present_x << endl;
       } else {
         cout << "i=" << slv.iter_trns_x << endl;
       }
       slv.write_iter_n1_x = 0;
     }
   }
   slv.write_iter_n1_x++;

   return;
} // end of write_iter_x
// -----------------------------------------------------------------------------
void write_solution(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SolveBlocks &slv,
   Circuit &cct) {

   int width_r;
   int i_xbeu,i_xbel;
   int i0;

   i0 = slv.index_file_solution;
   width_r = slv.outf_sol_word_width_real;

   slv.f_output[i0] << "cct.val_xvr:" << endl;
   print_vec_double_1(slv.f_output[i0],cct.val_xvr,width_r);

   slv.f_output[i0] << "xbe_usr.val_aux:" << endl;
   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     for (int i=0; i < xbe_lib[i_xbel].n_aux; i++) {
       slv.f_output[i0] << setw(width_r)
         << xbe_usr[i_xbeu].val_aux[i] << endl;
     }
   }

   return;
} // end of write_solution
// -----------------------------------------------------------------------------
void write_trns(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   int n_var,i_ov;
   bool flag_write;

   if (cct_file.n_ov_xbe > 0) {
     xbe_ov_prm(xbe_lib,xbe_usr,xbe_jac,cct,cct_file,global);
   }

   for (int i_file=0; i_file < slv.n_outfile; i_file++) {
     if (slv.flag_solution[i_file]) continue;

     flag_write = false;

     if (slv.time_write <= slv.out_tend[i_file]) {
       if (slv.flag_out_delt_fixed[i_file]) {
         if (slv.time_write >= slv.out_tnext[i_file]) {
           flag_write = true;
           slv.out_tnext[i_file] = slv.time_write + slv.out_delt[i_file];
         }
       } else {
         if (slv.time_write >= slv.out_tstart[i_file]) {
           flag_write = true;
         }
       }
     }

     if (flag_write) {
       slv.total_lines[i_file]++;
       if (slv.total_lines[i_file] > slv.limit_lines[i_file]) {
         cout << "write_trns: limit_lines is exceeded." << endl;
         cout << "  index_solve = " << slv.index_solve
              << ", i_file = " << i_file << endl;
         cout << "  Halting..." << endl; exit (1);
       }
       n_var = slv.out_nvar[i_file];

       for (int i_var=0; i_var < n_var; i_var++) {
         i_ov = slv.out_var[i_file][i_var];
         assign_ov_prm(i_ov,i_file,i_var,xbe_usr,
           cct,slv,cct_file,global);
       }
       slv.f_output[i_file]
         << setw(slv.outf_real_word_width)
         << slv.time_write << " ";
       for (int i_var=0; i_var < n_var; i_var++) {
         slv.f_output[i_file]
           << setw(slv.outf_real_word_width)
           << slv.outvar_temp[i_file][i_var] << " ";
       }
       slv.f_output[i_file] << endl;
     }
   }
   return;
} // end of write_trns
// -----------------------------------------------------------------------------
void write_dc_startup(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   int n_var,i_ov;

   if (cct_file.n_ov_xbe > 0) {
     xbe_ov_prm(xbe_lib,xbe_usr,xbe_jac,cct,cct_file,global);
   }

   for (int i_file=0; i_file < slv.n_outfile; i_file++) {
     if (slv.flag_solution[i_file]) continue;
     n_var = slv.out_nvar[i_file];

     for (int i_var=0; i_var < n_var; i_var++) {
       i_ov = slv.out_var[i_file][i_var];

       assign_ov_prm(i_ov,i_file,i_var,xbe_usr,
         cct,slv,cct_file,global);
     }
     for (int i_var=0; i_var < n_var; i_var++) {
       slv.f_output[i_file]
         << setw(slv.outf_real_word_width)
         << slv.outvar_temp[i_file][i_var] << " ";
     }
     slv.f_output[i_file] << endl;
   }
   return;
} // end of write_dc_startup
// -----------------------------------------------------------------------------
void solve_trns_x_exp(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global) {

   assign_const_1<bool>(global.flags,false);

   init_sol_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,smat,global);

   if (cct.flag_alg_loop) {
     smat.mat_startup_1_x(xbe_lib,xbe_usr,cct,cct_file);

     smat.w_x.allocate_1(smat.m_x.n_nz,smat.m_x.n_row,smat.m_x.n_col);
     smat.mo_x.allocate_1(smat.m_x.n_row);
   }

   if (slv.x_algo_feuler) {
     solve_trns_x_feuler(xbe_lib,xbe_usr,xbe_jac,
       slv,cct,cct_file,smat,global);
   } else if (slv.x_algo_rk4) {
     solve_trns_x_rk4(xbe_lib,xbe_usr,xbe_jac,
       slv,cct,cct_file,smat,global);
   } else if (slv.x_algo_rkf45) {
     solve_trns_x_rkf45(xbe_lib,xbe_usr,xbe_jac,
       slv,cct,cct_file,smat,global);
   } else if (slv.x_algo_bs23) {
     solve_trns_x_bs23(xbe_lib,xbe_usr,xbe_jac,
       slv,cct,cct_file,smat,global);
   } else if (slv.x_algo_meuler) {
     solve_trns_x_meuler(xbe_lib,xbe_usr,xbe_jac,
       slv,cct,cct_file,smat,global);
   } else if (slv.x_algo_heun) {
     solve_trns_x_heun(xbe_lib,xbe_usr,xbe_jac,
       slv,cct,cct_file,smat,global);
   } else {
     cout << "solve_trns_x_exp:" << endl;
     cout << "  wrong trns method. Halting..." << endl;
     exit(1);
   }

// clear flags once again
   assign_const_1<bool>(global.flags,false);

   return;
} // end of solve_trns_x_exp
// -----------------------------------------------------------------------------
void solve_trns_x_feuler(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global) {

   double time_next_1;
   bool flag_tend_reached;

   cout << "solve_trns_x_feuler starts..." << endl;
   solve_trns_x_common(xbe_lib,xbe_usr,xbe_jac,
     slv,cct,cct_file,global);

   flag_tend_reached = false;
   slv.iter_trns_x = -1;
   global.iter_trns_x = slv.iter_trns_x;

   while (!flag_tend_reached) {
     slv.iter_trns_x++; global.iter_trns_x = slv.iter_trns_x;
     write_iter_x(slv);

     if (slv.iter_trns_x != 0) {
       if (cct.flag_modulo_x) {
         xbe_modulo(cct,xbe_usr);
       }
     }

     slv.time_next_x  = slv.time_present_x + slv.delt_x;
     global.time_given_x = slv.time_next_x;

//   evaluate/integrate and update:

     if (cct.flag_alg_loop) {
       xbe_feuler_al(xbe_lib,xbe_usr,xbe_jac,
         cct,slv,smat,global);
     } else {
       xbe_feuler(xbe_lib,xbe_usr,xbe_jac,cct,slv,global);
     }

//   xelements using reset:
     if (cct.flag_reset_x) {
       xbe_reset_1(cct.flag_alg_loop,xbe_lib,xbe_usr,xbe_jac,cct,smat,global);
     }

     slv.time_present_x = slv.time_next_x;
     slv.time_write     = slv.time_next_x;

     write_trns(xbe_lib,xbe_usr,xbe_jac,
       cct,slv,cct_file,global);

//   find the new time step if necessary

     slv.delt_x = slv.delt0_x;
     slv.delt_x = max(slv.delt_x,slv.delt_min_x);

     if (cct.flag_limit_tstep_x) {
       get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
     }
     time_next_1 = slv.time_present_x + slv.delt_x;
     if (time_next_1 >= global.time_end) {
       slv.delt_x = global.time_end - slv.time_present_x + slv.delt_small;
       slv.delt_x = max(slv.delt_x,slv.delt_min_x);
     }
     if (cct.flag_save_history_x) {
       save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);
     }
     if (slv.flag_limit_iter_trns) {
       if (slv.iter_trns_x > slv.itmax_trns) {
         cout << "solve_trns_x_feuler: itmax_trns exceeded." << endl;
         break;
       }
     }
     if ((slv.time_present_x+slv.delt_small) >= global.time_end) {
       flag_tend_reached = true;
     }
   }
   cout << "solve_trns_x_feuler ends" << endl;
   return;
} // end of solve_trns_x_feuler
// -----------------------------------------------------------------------------
void solve_trns_x_rk4(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global) {

   double time_next_1;
   bool flag_tend_reached;

   solve_trns_x_common(xbe_lib,xbe_usr,xbe_jac,
     slv,cct,cct_file,global);

   cout << "solve_trns_x_rk4 starting..." << endl;
   cout << "solve_trns_x_rk4: slv.itmax_trns = " << slv.itmax_trns << endl;

   flag_tend_reached = false;
   slv.iter_trns_x = -1;
   global.iter_trns_x = slv.iter_trns_x;

   while (!flag_tend_reached) {
     slv.iter_trns_x++; global.iter_trns_x = slv.iter_trns_x;
     write_iter_x(slv);

     if (slv.iter_trns_x != 0) {
       if (cct.flag_modulo_x) {
         xbe_modulo(cct,xbe_usr);
       }
     }
     slv.time_next_x  = slv.time_present_x + slv.delt_x;
//   cout << "solve_trns_x_rk4: i0 = " << i0 << endl;

//   evaluate/integrate and update:

     if (cct.flag_alg_loop) {
       xbe_rk4_al(xbe_lib,xbe_usr,xbe_jac,
         cct,slv,smat,global);
     } else {
       xbe_rk4(xbe_lib,xbe_usr,xbe_jac,cct,slv,global);
     }

     if (cct.flag_reset_x) {
       xbe_reset_1(false,xbe_lib,xbe_usr,xbe_jac,cct,smat,global);
     }

     slv.time_present_x = slv.time_next_x;
     slv.time_write     = slv.time_next_x;

     write_trns(xbe_lib,xbe_usr,xbe_jac,
       cct,slv,cct_file,global);

//   find the new time step if necessary

     slv.delt_x = slv.delt0_x;
     slv.delt_x = max(slv.delt_x,slv.delt_min_x);

     if (cct.flag_limit_tstep_x) {
       get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
     }
     time_next_1 = slv.time_present_x + slv.delt_x;
     if (time_next_1 >= global.time_end) {
       slv.delt_x = global.time_end - slv.time_present_x + slv.delt_small;
       slv.delt_x = max(slv.delt_x,slv.delt_min_x);
     }

     if (cct.flag_save_history_x) {
       save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);
     }
     if (slv.flag_limit_iter_trns) {
       if (slv.iter_trns_x > slv.itmax_trns) {
         cout << "solve_trns_x_rk4: itmax_trns exceeded." << endl;
         break;
       }
     }
     if ((slv.time_present_x+slv.delt_small) >= global.time_end) {
       flag_tend_reached = true;
     }
   }
   cout << "solve_trns_x_rk4 ends" << endl;

   return;
} // end of solve_trns_x_rk4
// -----------------------------------------------------------------------------
void solve_trns_x_rkf45(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global) {

   double time_next_1,t_n,delta;
   bool flag_tend_reached;

   cout << "solve_trns_x_rkf45 starts..." << endl;

   solve_trns_x_common(xbe_lib,xbe_usr,xbe_jac,
     slv,cct,cct_file,global);

   slv.rkf45_n_accept = 0;
   slv.rkf45_n_reject = 0;

   xbeu_copy_1(xbe_lib,xbe_usr,cct);

   flag_tend_reached = false;
   slv.iter_trns_x = -1;
   global.iter_trns_x = slv.iter_trns_x;

   while (!flag_tend_reached) {
     slv.iter_trns_x++; global.iter_trns_x = slv.iter_trns_x;
//   cout << "solve_trns_x_rkf45: slv.iter_trns_x = " << slv.iter_trns_x << endl;
     write_iter_x(slv);

     if (slv.iter_trns_x != 0) {
       if (cct.flag_modulo_x) {
         xbe_modulo(cct,xbe_usr);
       }
     }
     if (cct.flag_alg_loop) {
       xbe_rkf45_al(global.I_RKF45_LTE ,xbe_lib,xbe_usr,xbe_jac,
         cct,slv,smat,global);
       xbe_rkf45_al(global.I_RKF45_NORM,xbe_lib,xbe_usr,xbe_jac,
         cct,slv,smat,global);
     } else {
       xbe_rkf45(global.I_RKF45_LTE ,xbe_lib,xbe_usr,xbe_jac,cct,slv,global);
       xbe_rkf45(global.I_RKF45_NORM,xbe_lib,xbe_usr,xbe_jac,cct,slv,global);
     }

//   If LTE is not within tolerance, change delt_x;
//   If LTE is within tolerance, accept the time step.

     if (slv.rkf45_norm2 <= slv.rkf45_tolr) {
       slv.rkf45_n_accept++;
     } else {
       slv.rkf45_n_reject++;
       if (slv.delt_x == slv.delt_min_x) {
         cout << "solve_trns_x_rkf45:" << endl;
         cout << "  tolerance not met even with" << endl;
         cout << "  the smallest time step." << endl;
         cout << "  rkf45_norm2=" << slv.rkf45_norm2
           << ",  rkf45_tolr=" << slv.rkf45_tolr << endl;
         cout << "  delt_x=" << slv.delt_x
           << ",  delt_min_x=" << slv.delt_min_x << endl;
         cout << "  Halting..." << endl;
         exit(1);
       }
     }

     if (slv.rkf45_norm2 <= slv.rkf45_tolr) {
       t_n = slv.time_present_x + slv.delt_x;
       slv.time_present_x = t_n;
       slv.time_write     = t_n;

//     compute solution with the RK4 part of RKF45:

       if (cct.flag_alg_loop) {
         xbe_rkf45_al(global.I_RKF45_SVEC,xbe_lib,xbe_usr,xbe_jac,
           cct,slv,smat,global);
       } else {
         xbe_rkf45(global.I_RKF45_SVEC,xbe_lib,xbe_usr,xbe_jac,cct,slv,global);
       }

//     If the solution is accepted, treat intgrtr_reset kind of elements.

       if (cct.flag_reset_x) {
         xbe_reset_1(false,xbe_lib,xbe_usr,xbe_jac,cct,smat,global);
       }
       write_trns(xbe_lib,xbe_usr,xbe_jac,
         cct,slv,cct_file,global);
       xbeu_copy_1(xbe_lib,xbe_usr,cct);
     }

//   compute the next time step (whether or not LTE is within tolerance)

//   this is used to control LTE, rather than LTE/h
     delta = 0.87*pow((slv.rkf45_tolr/slv.rkf45_norm2),0.2);

     if (delta < slv.rkf45_fctr_min) {
       slv.delt_x = slv.rkf45_fctr_min*slv.delt_x;
     } else if (delta >= slv.rkf45_fctr_max) {
       slv.delt_x = slv.rkf45_fctr_max*slv.delt_x;
     } else {
       slv.delt_x= delta*slv.delt_x;
     }

//   limit the time step:

     slv.delt_x = min(slv.delt_x,slv.delt_max_x);
     slv.delt_x = max(slv.delt_x,slv.delt_min_x);

//   check if the time step should be limited by elements:

     if (slv.delt_x > slv.delt_min_x) {
       if (cct.flag_limit_tstep_x) {
         get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
       }
     }
     time_next_1 = slv.time_present_x + slv.delt_x;
     if (time_next_1 >= global.time_end) {
       slv.delt_x = global.time_end - slv.time_present_x + slv.delt_small;
       slv.delt_x = max(slv.delt_x,slv.delt_min_x);
     }
     if (cct.flag_save_history_x) {
       save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);
     }
     if (slv.flag_limit_iter_trns) {
       if (slv.iter_trns_x > slv.itmax_trns) {
         cout << "solve_trns_x_rkf45: itmax_trns exceeded." << endl;
         break;
       }
     }
//   cout << "solve_trns_x_rkf45:"
//     << " slv.time_present_x = " << slv.time_present_x
//     << " slv.delt_small = " << slv.delt_small
//     << " slv.delt_max_x = " << slv.delt_max_x
//     << " global.time_end = " << global.time_end
//     << endl;
     if ((slv.time_present_x+slv.delt_small) >= global.time_end) {
       flag_tend_reached = true;
     }
   }
   cout << "solve_trns_x_rkf45: rkf45_n_accept="
     << slv.rkf45_n_accept << endl;
   cout << "solve_trns_x_rkf45: rkf45_n_reject="
     << slv.rkf45_n_reject << endl;

   cout << "solve_trns_x_rkf45 ends" << endl;

   return;
} // end of solve_trns_x_rkf45
// -----------------------------------------------------------------------------
void solve_trns_x_bs23(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global) {

   double t_n,time_next_1,delta;
   bool flag_tend_reached;

   solve_trns_x_common(xbe_lib,xbe_usr,xbe_jac,
     slv,cct,cct_file,global);

   slv.bs23_n_accept = 0;
   slv.bs23_n_reject = 0;

   xbeu_copy_1(xbe_lib,xbe_usr,cct);

   flag_tend_reached = false;
   slv.iter_trns_x = -1;
   global.iter_trns_x = slv.iter_trns_x;

   while (!flag_tend_reached) {
     slv.iter_trns_x++; global.iter_trns_x = slv.iter_trns_x;
     write_iter_x(slv);

     if (slv.iter_trns_x != 0) {
       if (cct.flag_modulo_x) {
         xbe_modulo(cct,xbe_usr);
       }
     }
//   compute local error estimate (denote it by slv.bs23_norm2)

     if (cct.flag_alg_loop) {
       xbe_bs23_al(global.I_BS23_LTE ,xbe_lib,xbe_usr,xbe_jac,
         cct,slv,smat,global);
       xbe_bs23_al(global.I_BS23_NORM,xbe_lib,xbe_usr,xbe_jac,
         cct,slv,smat,global);
     } else {
       xbe_bs23(global.I_BS23_LTE ,xbe_lib,xbe_usr,xbe_jac,cct,slv,global);
       xbe_bs23(global.I_BS23_NORM,xbe_lib,xbe_usr,xbe_jac,cct,slv,global);
     }

//   If LTE is not within tolerance, change delt_x;
//   If LTE is within tolerance, accept the time step.

     if (slv.bs23_norm2 <= slv.bs23_tolr) {
       slv.bs23_n_accept++;
     } else {
       slv.bs23_n_reject++;
       if (slv.delt_x == slv.delt_min_x) {
         cout << "solve_trns_x_bs23:" << endl;
         cout << "  tolerance not met even with" << endl;
         cout << "  the smallest time step." << endl;
         cout << "  bs23_norm2=" << slv.bs23_norm2
           << ",  bs23_tolr=" << slv.bs23_tolr << endl;
         cout << "  delt_x=" << slv.delt_x
           << ",  delt_min_x=" << slv.delt_min_x << endl;
         cout << "  Halting..." << endl;
         exit(1);
       }
     }

     if (slv.bs23_norm2 <= slv.bs23_tolr) {
       t_n = slv.time_present_x + slv.delt_x;
       slv.time_write     = t_n;
       slv.time_present_x = t_n;

//     compute solution with the 2nd order part of bs23:

       if (cct.flag_alg_loop) {
         xbe_bs23_al(global.I_BS23_SVEC,xbe_lib,xbe_usr,xbe_jac,
           cct,slv,smat,global);
       } else {
         xbe_bs23(global.I_BS23_SVEC,xbe_lib,xbe_usr,xbe_jac,cct,slv,global);
       }

       if (cct.flag_reset_x) {
         xbe_reset_1(false,xbe_lib,xbe_usr,xbe_jac,cct,smat,global);
       }

       write_trns(xbe_lib,xbe_usr,xbe_jac,
         cct,slv,cct_file,global);

       xbeu_copy_1(xbe_lib,xbe_usr,cct);
     }
//   compute the next time step (whether or not LTE is within tolerance)
//   (1/2)**0.5 = 0.7071
//   (1/2)**0.333 = 0.7937

//   this is used to control LTE, rather than LTE/h
     delta = 0.7937*pow((slv.bs23_tolr/slv.bs23_norm2),0.3333333);
  
     if (delta < slv.bs23_fctr_min) {
       slv.delt_x = slv.bs23_fctr_min*slv.delt_x;
     } else if (delta >= slv.bs23_fctr_max) {
       slv.delt_x = slv.bs23_fctr_max*slv.delt_x;
     } else {
       slv.delt_x= delta*slv.delt_x;
     }

//   limit the time step:

     slv.delt_x = min(slv.delt_x,slv.delt_max_x);
     slv.delt_x = max(slv.delt_x,slv.delt_min_x);

//   check if the time step should be limited by elements:

     if (slv.delt_x > slv.delt_min_x) {
       if (cct.flag_limit_tstep_x) {
         get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
       }
     }
     time_next_1 = slv.time_present_x + slv.delt_x;
     if (time_next_1 >= global.time_end) {
       slv.delt_x = global.time_end - slv.time_present_x + slv.delt_small;
       slv.delt_x = max(slv.delt_x,slv.delt_min_x);
     }
     if (cct.flag_save_history_x) {
       save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);
     }
     if (slv.flag_limit_iter_trns) {
       if (slv.iter_trns_x > slv.itmax_trns) {
         cout << "solve_trns_x_bs23: itmax_trns exceeded." << endl;
         break;
       }
     }
     if ((slv.time_present_x+slv.delt_small) >= global.time_end) {
       flag_tend_reached = true;
     }
   }
   cout << "solve_trns_x_bs23: bs23_n_accept="
     << slv.bs23_n_accept << endl;
   cout << "solve_trns_x_bs23: bs23_n_reject="
     << slv.bs23_n_reject << endl;

   cout << "solve_trns_x_bs23 ends" << endl;

   return;
} // end of solve_trns_x_bs23
// -----------------------------------------------------------------------------
void solve_trns_x_meuler(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global) {

   double t_n,time_next_1;
   bool flag_tend_reached;

   cout << "solve_trns_x_meuler starts..." << endl;

   solve_trns_x_common(xbe_lib,xbe_usr,xbe_jac,
     slv,cct,cct_file,global);

   flag_tend_reached = false;
   slv.iter_trns_x = -1;
   global.iter_trns_x = slv.iter_trns_x;

   while (!flag_tend_reached) {
     slv.iter_trns_x++; global.iter_trns_x = slv.iter_trns_x;
     write_iter_x(slv);

     if (slv.iter_trns_x != 0) {
       if (cct.flag_modulo_x) {
         xbe_modulo(cct,xbe_usr);
       }
     }
     t_n = slv.time_present_x + slv.delt_x;

//   evaluate/integrate and update:

     if (cct.flag_alg_loop) {
       xbe_meuler_al(xbe_lib,xbe_usr,xbe_jac,
         cct,slv,smat,global);
     } else {
       xbe_meuler(xbe_lib,xbe_usr,xbe_jac,cct,slv,global);
     }

     if (cct.flag_reset_x) {
       xbe_reset_1(false,xbe_lib,xbe_usr,xbe_jac,cct,smat,global);
     }

     slv.time_present_x = t_n;
     slv.time_write     = t_n;

     write_trns(xbe_lib,xbe_usr,xbe_jac,
       cct,slv,cct_file,global);

//   find the new time step if necessary

     slv.delt_x = slv.delt0_x;
     slv.delt_x = max(slv.delt_x,slv.delt_min_x);

     if (cct.flag_limit_tstep_x) {
       get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
     }
     time_next_1 = slv.time_present_x + slv.delt_x;
     if (time_next_1 >= global.time_end) {
       slv.delt_x = global.time_end - slv.time_present_x + slv.delt_small;
       slv.delt_x = max(slv.delt_x,slv.delt_min_x);
     }
     if (cct.flag_save_history_x) {
       save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);
     }
     if (slv.flag_limit_iter_trns) {
       if (slv.iter_trns_x > slv.itmax_trns) {
         cout << "solve_trns_x_meuler: itmax_trns exceeded." << endl;
         break;
       }
     }
     if ((slv.time_present_x+slv.delt_small) >= global.time_end) {
       flag_tend_reached = true;
     }
   }
   cout << "solve_trns_x_meuler ends" << endl;

   return;
} // end of solve_trns_x_meuler
// -----------------------------------------------------------------------------
void solve_trns_x_heun(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global) {

   double t_n,time_next_1;
   bool flag_tend_reached;

   cout << "solve_trns_x_heun starts..." << endl;

   solve_trns_x_common(xbe_lib,xbe_usr,xbe_jac,
     slv,cct,cct_file,global);

   flag_tend_reached = false;
   slv.iter_trns_x = -1;
   global.iter_trns_x = slv.iter_trns_x;

   while (!flag_tend_reached) {
     slv.iter_trns_x++; global.iter_trns_x = slv.iter_trns_x;
     write_iter_x(slv);

     if (slv.iter_trns_x != 0) {
       if (cct.flag_modulo_x) {
         xbe_modulo(cct,xbe_usr);
       }
     }
     t_n  = slv.time_present_x + slv.delt_x;

//   evaluate/integrate and update:
     if (cct.flag_alg_loop) {
       xbe_heun_al(xbe_lib,xbe_usr,xbe_jac,
         cct,slv,smat,global);
     } else {
       xbe_heun(xbe_lib,xbe_usr,xbe_jac,cct,slv,global);
     }

     if (cct.flag_reset_x) {
       xbe_reset_1(false,xbe_lib,xbe_usr,xbe_jac,cct,smat,global);
     }

     slv.time_present_x = t_n;
     slv.time_write     = t_n;

     write_trns(xbe_lib,xbe_usr,xbe_jac,
       cct,slv,cct_file,global);

//   find the new time step if necessary

     slv.delt_x = slv.delt0_x;
     slv.delt_x = max(slv.delt_x,slv.delt_min_x);

     if (cct.flag_limit_tstep_x) {
       get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
     }
     time_next_1 = slv.time_present_x + slv.delt_x;
     if (time_next_1 >= global.time_end) {
       slv.delt_x = global.time_end - slv.time_present_x + slv.delt_small;
       slv.delt_x = max(slv.delt_x,slv.delt_min_x);
     }
     if (cct.flag_save_history_x) {
       save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);
     }
     if (slv.flag_limit_iter_trns) {
       if (slv.iter_trns_x > slv.itmax_trns) {
         cout << "solve_trns_x_heun: itmax_trns exceeded." << endl;
         break;
       }
     }
     if ((slv.time_present_x+slv.delt_small) >= global.time_end) {
       flag_tend_reached = true;
     }
   }
   cout << "solve_trns_x_heun ends" << endl;

   return;
} // end of solve_trns_x_heun
// -----------------------------------------------------------------------------
void xbe_ov_prm(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   CctFile &cct_file,
   Global &global) {

   int i_xbeu,i_xbel;

   global.flags[global.i_outvar] = true;

   for (int i=0; i < cct_file.n_ov; i++) {
     if (cct_file.ov_flag[i] == global.I_OV_XBE) {
       i_xbeu = cct_file.ov1[i];
       i_xbel = xbe_usr[i_xbeu].index_xbel;

//     assign xbe_usr.val_vr from cct.val_xvr
       cct_to_xbe_1(i_xbeu,xbe_lib,xbe_usr,cct);

       get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbeu]);
     }
   }
   global.flags[global.i_outvar] = false;

   return;
} // end of xbe_ov_prm
// -----------------------------------------------------------------------------
void assign_ov_prm(
   const int i_ov,
   const int i_file,
   const int i_var,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   int i_xbeu,i_rprm,i_xvr;

   if (cct_file.ov_flag[i_ov] == global.I_OV_XBE) {
     i_xbeu = cct_file.ov1[i_ov];
     i_rprm = cct_file.ov2[i_ov];
     slv.outvar_temp[i_file][i_var] = xbe_usr[i_xbeu].outprm[i_rprm];
   } else if (cct_file.ov_flag[i_ov] == global.I_OV_XVR) {
     i_xvr = cct_file.ov1[i_ov];
     slv.outvar_temp[i_file][i_var] = cct.val_xvr[i_xvr];
   }

   return;
} // end of assign_ov_prm
// -----------------------------------------------------------------------------
void xbe_find_nextbreak(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;

   global.flags[global.i_trns] = true;
   global.flags[global.i_next_time] = true;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     global.time_nextbreak_x = xbe_usr[i_xbeu].next_break;

     if (xbe_lib[i_xbel].flag_lmttstep) {
       get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbeu]);
       xbe_usr[i_xbeu].next_break = global.time_nextbreak_x;
     }
   }
   global.flags[global.i_trns] = false;
   global.flags[global.i_next_time] = false;

   return;
} // end of xbe_find_nextbreak
// -----------------------------------------------------------------------------
void get_tnext_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;
   double time_next_1,t_next,delt1;

   global.flags[global.i_trns] = true;
   global.flags[global.i_next_time] = true;

   global.time_given_x = slv.time_present_x;
   time_next_1 = slv.time_present_x + slv.delt_x;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     if ((time_next_1 >= xbe_usr[i_xbeu].next_break) ||
         (xbe_lib[i_xbel].flag_savehist)) {

       global.time_nextbreak_x = global.time_end;
       cct_to_xbe_1(i_xbeu,xbe_lib,xbe_usr,cct);

       get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbeu]);
       xbe_usr[i_xbeu].next_break = global.time_nextbreak_x;
     }
   }
   global.flags[global.i_trns] = false;
   global.flags[global.i_next_time] = false;

   if (slv.delt_x != slv.delt_min_x) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       if (xbe_lib[i_xbel].flag_lmttstep) {
         t_next = xbe_usr[i_xbeu].next_break;
         if (t_next <= time_next_1) {
           time_next_1 = t_next;
         }
       }
     }
     delt1 = time_next_1 - slv.time_present_x;
     slv.delt_x = max(delt1,slv.delt_min_x);

     slv.delt_x = min(slv.delt_x,slv.delt_max_x);
     slv.delt_x = max(slv.delt_x,slv.delt_min_x);
   }

   return;
} // end of get_tnext_x
// -----------------------------------------------------------------------------
void xbe_modulo(
   Circuit &cct,
   vector<XbeUsr> &xbe_usr) {

   int i_xbeu;
   int i_xbeu_vr;
   int n1,i_xbeu_1,i_vr_1;
   double x1,x2,delx,x_old,x_new;

   for (int i=0; i < cct.nttl_xbeu_modulo; i++) {
     i_xbeu = cct.xbeu_modulo_map[i];

     x1 = xbe_usr[i_xbeu].rprm[0];
     x2 = xbe_usr[i_xbeu].rprm[1];
     delx = x2-x1;

//   Expect the modulo type xbe to have only one var.

     i_xbeu_vr = xbe_usr[i_xbeu].vr[0];
     x_old = cct.val_xvr[i_xbeu_vr];

     if (x_old > x1) {
       x_new = x1 + fmod((x_old-x1),delx);
     } else {
       x_new = x2 - fmod((x1-x_old),delx);
     }

     if (x_old != x_new) {
       cct.val_xvr[i_xbeu_vr] = x_new;
       n1 = cct.n_map_vr[i_xbeu_vr];

       for (int i1=0; i1 < n1; i1++) {
         i_xbeu_1 = cct.map_vr_1[i_xbeu_vr][i1];
         i_vr_1   = cct.map_vr_2[i_xbeu_vr][i1];
         xbe_usr[i_xbeu_1].val_vr[i_vr_1] = x_new;
       }
     }
   }
   return;
} // end of xbe_modulo
// -----------------------------------------------------------------------------
void xbe_modulo_implicit(
   Circuit &cct,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SysMat &smat) {

   int i_xbeu;
   int i_xbeu_vr,i_svec;
   int n1,i_xbeu_1,i_vr_1;
   double x1,x2,delx,x_old,x_new;

   for (int i=0; i < cct.nttl_xbeu_modulo; i++) {
     i_xbeu = cct.xbeu_modulo_map[i];

     x1 = xbe_usr[i_xbeu].rprm[0];
     x2 = xbe_usr[i_xbeu].rprm[1];
     delx = x2-x1;

//   Expect the modulo type xbe to have only one var.

     i_xbeu_vr = xbe_usr[i_xbeu].vr[0];
     x_old = cct.val_xvr[i_xbeu_vr];

     if (x_old > x1) {
       x_new = x1 + fmod((x_old-x1),delx);
     } else {
       x_new = x2 - fmod((x1-x_old),delx);
     }

     if (x_old != x_new) {
       cct.val_xvr[i_xbeu_vr] = x_new;

       i_svec = cct.map_xbeuvr_to_svec[i_xbeu_vr];
       smat.svec_old_1_x[i_svec] = x_new;

       n1 = cct.n_map_vr[i_xbeu_vr];

       for (int i1=0; i1 < n1; i1++) {
         i_xbeu_1 = cct.map_vr_1[i_xbeu_vr][i1];
         i_vr_1   = cct.map_vr_2[i_xbeu_vr][i1];
         xbe_usr[i_xbeu_1].val_vr[i_vr_1] = x_new;
       }
     }
   }

   return;
} // end of xbe_modulo_implicit
// -----------------------------------------------------------------------------
void xbe_reset_1(
   const int flag_implicit,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SysMat &smat,
   Global &global) {

   int i_xbeu,i_xbel;

   global.flags[global.i_trns] = true;
   global.flags[global.i_reset_x] = true;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     if (xbe_lib[i_xbel].flag_reset) {
       if (!flag_implicit) {
         cct_to_xbe_1(i_xbeu,xbe_lib,xbe_usr,cct);
       }
       get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbeu]);

       if (flag_implicit) {
         xbe_to_cct_all(xbe_lib,xbe_usr,cct);
         form_solvec_x(xbe_lib,xbe_usr,smat,cct);
       } else {
         xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
       }
     }
   }
   global.flags[global.i_trns] = false;
   global.flags[global.i_reset_x] = false;

   return;
} // end of xbe_reset_1
// -----------------------------------------------------------------------------
void xbe_feuler(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   int i_xbeu,i_xbel,n_f1;
   int var_flag,var_number;

   xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;

       for (int i_func=0; i_func < n_f1; i_func++) {
         var_flag   = xbe_lib[i_xbel].ddt_varflag  [i_func];
         var_number = xbe_lib[i_xbel].ddt_varnumber[i_func];

         if (var_flag == global.I_XVR) {
           xbe_usr[i_xbeu].val_vr[var_number] +=
             slv.delt_x*xbe_usr[i_xbeu].f[i_func];
         } else if (var_flag == global.I_XAUX) {
           xbe_usr[i_xbeu].val_aux[var_number] +=
             slv.delt_x*xbe_usr[i_xbeu].f[i_func];
         } else {
           cout << "xbe_feuler: var_flag = " << var_flag
             << " does not make sense. Halting..." << endl;
           exit(1);
         }
       }
       xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
     }
   }
   xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

   return;
} // end of xbe_feuler
// -----------------------------------------------------------------------------
void xbe_feuler_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global) {

   int i_xbeu,i_xbel,n_f1;
   int var_flag,var_number;

// evaluate ddt elements
   xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;

       for (int i_func=0; i_func < n_f1; i_func++) {
         var_flag   = xbe_lib[i_xbel].ddt_varflag  [i_func];
         var_number = xbe_lib[i_xbel].ddt_varnumber[i_func];

         if (var_flag == global.I_XVR) {
           xbe_usr[i_xbeu].val_vr[var_number] +=
             slv.delt_x*xbe_usr[i_xbeu].f[i_func];
           xbe_usr[i_xbeu].val_vr_u[var_number] =
           xbe_usr[i_xbeu].val_vr  [var_number];
         } else if (var_flag == global.I_XAUX) {
           xbe_usr[i_xbeu].val_aux[var_number] +=
             slv.delt_x*xbe_usr[i_xbeu].f[i_func];
           xbe_usr[i_xbeu].val_aux_u[var_number] =
           xbe_usr[i_xbeu].val_aux  [var_number];
         } else {
           cout << "xbe_feuler_al: var_flag = " << var_flag
             << " does not make sense. Halting..." << endl;
           exit(1);
         }
       }
       xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
     }
   }

   if (cct.flag_linear_x) {
     solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
   } else {
     solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,global);
   }
   return;
} // end of xbe_feuler_al
// -----------------------------------------------------------------------------
void xbe_rk4(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   double h;

   h = slv.delt_x;

// copy to old
   xbeu_copy_1(xbe_lib,xbe_usr,cct);

// RK4: stage 1

   global.time_given_x = slv.time_present_x;
   if (flag_nan(slv.time_present_x)) {
     cout << "xbe_rk4: slv.time_present_x is NAN. Halting..." << endl;
     exit(1);
   }

   xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(0,xbe_lib,xbe_usr,cct);

   update_rk4(1,(0.5*h),xbe_lib,xbe_usr,cct,global);

// RK4: stage 2
   global.time_given_x = slv.time_present_x + 0.5*h;

   xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(1,xbe_lib,xbe_usr,cct);

   update_rk4(2,(0.5*h),xbe_lib,xbe_usr,cct,global);

// RK4: stage 3

   xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(2,xbe_lib,xbe_usr,cct);

   update_rk4(3,h,xbe_lib,xbe_usr,cct,global);

   global.time_given_x = slv.time_present_x + h;

   xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(3,xbe_lib,xbe_usr,cct);

   update_rk4(4,(h/6.0),xbe_lib,xbe_usr,cct,global);

// final update of non-integrator variables

// RK4: stage 4
   xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

   return;
} // end of xbe_rk4
// -----------------------------------------------------------------------------
void xbe_rk4_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global) {

   double h;

   h = slv.delt_x;

// copy to old
   xbeu_copy_1(xbe_lib,xbe_usr,cct);

// RK4: stage 1

   global.time_given_x = slv.time_present_x;
   if (flag_nan(slv.time_present_x)) {
     cout << "xbe_rk4_al: slv.time_present_x is NAN. Halting..." << endl;
     exit(1);
   }
   xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(0,xbe_lib,xbe_usr,cct);

   update_rk4_al(1,(0.5*h),xbe_lib,xbe_usr,cct,global);

// RK4: stage 2
   global.time_given_x = slv.time_present_x + 0.5*h;

   if (cct.flag_linear_x) {
     solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
   } else {
     solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,global);
   }

   xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(1,xbe_lib,xbe_usr,cct);

   update_rk4_al(2,(0.5*h),xbe_lib,xbe_usr,cct,global);

// RK4: stage 3

   if (cct.flag_linear_x) {
     solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
   } else {
     solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,global);
   }

   xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(2,xbe_lib,xbe_usr,cct);

   update_rk4_al(3,h,xbe_lib,xbe_usr,cct,global);

   global.time_given_x = slv.time_present_x + h;

   if (cct.flag_linear_x) {
     solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
   } else {
     solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,global);
   }

   xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(3,xbe_lib,xbe_usr,cct);

   update_rk4_al(4,(h/6.0),xbe_lib,xbe_usr,cct,global);

// final update of non-integrator variables

// RK4: stage 4
   if (cct.flag_linear_x) {
     solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
   } else {
     solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,global);
   }

   return;
} // end of xbe_rk4_al
// -----------------------------------------------------------------------------
void xbe_rkf45(
   const int flag,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   double h;
   h = slv.delt_x;

   if (flag == global.I_RKF45_LTE) {
//   RKF45: stage 1

     global.time_given_x = slv.time_present_x;
     if (flag_nan(slv.time_present_x)) {
       cout << "xbe_rkf45: slv.time_present_x is NAN. Halting..." << endl;
       exit(1);
     }
     xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(0,xbe_lib,xbe_usr,cct);

     update_rkf45(1,h,xbe_lib,xbe_usr,cct,slv,global);

//   RKF45: stage 2
     global.time_given_x = slv.time_present_x + slv.rkf45_a1*h;

     xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(1,xbe_lib,xbe_usr,cct);

     update_rkf45(2,h,xbe_lib,xbe_usr,cct,slv,global);

//   RKF45: stage 3
     global.time_given_x = slv.time_present_x + slv.rkf45_a2*h;

     xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(2,xbe_lib,xbe_usr,cct);

     update_rkf45(3,h,xbe_lib,xbe_usr,cct,slv,global);

//   RKF45: stage 4
     global.time_given_x = slv.time_present_x + slv.rkf45_a3*h;

     xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(3,xbe_lib,xbe_usr,cct);

     update_rkf45(4,h,xbe_lib,xbe_usr,cct,slv,global);

//   RKF45: stage 5
     global.time_given_x = slv.time_present_x + h;

     xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(4,xbe_lib,xbe_usr,cct);

     update_rkf45(5,h,xbe_lib,xbe_usr,cct,slv,global);

//   RKF45: stage 5A
     global.time_given_x = slv.time_present_x + slv.rkf45_a5*h;

     xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(5,xbe_lib,xbe_usr,cct);
   } else if (flag == global.I_RKF45_SVEC) {
//   RKF45: stage 6
//   this stage computes the updated solvec (using the 4th order
//   part of RKF45) and is required only if the current step
//   is accepted.

     update_rkf45(6,h,xbe_lib,xbe_usr,cct,slv,global);

//   final update of non-integrator variables
//   Note that time_present_x has already been updated to (t + h).

     global.time_given_x = slv.time_present_x;

     xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);
   } else if (flag == global.I_RKF45_NORM) {
//   compute norm_2 (only for the integrator variables)

     slv.rkf45_norm2 = 0.0;
     update_rkf45(7,h,xbe_lib,xbe_usr,cct,slv,global);
     slv.rkf45_norm2 = sqrt(slv.rkf45_norm2);
   }
   return;
} // end of xbe_rkf45
// -----------------------------------------------------------------------------
void xbe_rkf45_al(
   const int flag,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global) {

   double h;
   h = slv.delt_x;

   if (flag == global.I_RKF45_LTE) {
//   RKF45: stage 1

     global.time_given_x = slv.time_present_x;
     if (flag_nan(slv.time_present_x)) {
       cout << "xbe_rkf45_al: slv.time_present_x is NAN. Halting..." << endl;
       exit(1);
     }
     xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(0,xbe_lib,xbe_usr,cct);

     update_rkf45_al(1,h,xbe_lib,xbe_usr,cct,slv,global);

//   RKF45: stage 2
     global.time_given_x = slv.time_present_x + slv.rkf45_a1*h;

     if (cct.flag_linear_x) {
       solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);
     }

     xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(1,xbe_lib,xbe_usr,cct);

     update_rkf45_al(2,h,xbe_lib,xbe_usr,cct,slv,global);

//   RKF45: stage 3
     global.time_given_x = slv.time_present_x + slv.rkf45_a2*h;

     if (cct.flag_linear_x) {
       solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);
     }

     xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(2,xbe_lib,xbe_usr,cct);

     update_rkf45_al(3,h,xbe_lib,xbe_usr,cct,slv,global);

//   RKF45: stage 4
     global.time_given_x = slv.time_present_x + slv.rkf45_a3*h;

     if (cct.flag_linear_x) {
       solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);
     }

     xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(3,xbe_lib,xbe_usr,cct);

     update_rkf45_al(4,h,xbe_lib,xbe_usr,cct,slv,global);

//   RKF45: stage 5
     global.time_given_x = slv.time_present_x + h;

     if (cct.flag_linear_x) {
       solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);
     }

     xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(4,xbe_lib,xbe_usr,cct);

     update_rkf45_al(5,h,xbe_lib,xbe_usr,cct,slv,global);

//   RKF45: stage 5A
     global.time_given_x = slv.time_present_x + slv.rkf45_a5*h;

     if (cct.flag_linear_x) {
       solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);
     }

     xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(5,xbe_lib,xbe_usr,cct);
   } else if (flag == global.I_RKF45_SVEC) {
//   RKF45: stage 6
//   this stage computes the updated solvec (using the 4th order
//   part of RKF45) and is required only if the current step
//   is accepted.

     update_rkf45_al(6,h,xbe_lib,xbe_usr,cct,slv,global);

//   final update of non-integrator variables
//   Note that time_present_x has already been updated to (t + h).

     global.time_given_x = slv.time_present_x;

     if (cct.flag_linear_x) {
       solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);
     }
   } else if (flag == global.I_RKF45_NORM) {
//   compute norm_2 (only for the integrator variables)

     slv.rkf45_norm2 = 0.0;
     update_rkf45_al(7,h,xbe_lib,xbe_usr,cct,slv,global);
     slv.rkf45_norm2 = sqrt(slv.rkf45_norm2);
   }
   return;
} // end of xbe_rkf45_al
// -----------------------------------------------------------------------------
void xbe_bs23(
   const int flag,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   double h;

   h = slv.delt_x;

   if (flag == global.I_BS23_LTE) {
//   bs23: stage 1

     global.time_given_x = slv.time_present_x;
     if (flag_nan(slv.time_present_x)) {
       cout << "xbe_bs23: slv.time_present_x is NAN. Halting..." << endl;
       exit(1);
     }
     xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(0,xbe_lib,xbe_usr,cct);
     update_bs23(1,h,xbe_lib,xbe_usr,cct,slv,global);

//   bs23: stage 2
     global.time_given_x = slv.time_present_x + slv.bs23_a1*h;

     xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(1,xbe_lib,xbe_usr,cct);
     update_bs23(2,h,xbe_lib,xbe_usr,cct,slv,global);

//   bs23: stage 3
     global.time_given_x = slv.time_present_x + slv.bs23_a2*h;

     xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(2,xbe_lib,xbe_usr,cct);
     update_bs23(3,h,xbe_lib,xbe_usr,cct,slv,global);

//   bs23: stage 4
     global.time_given_x = slv.time_present_x + h;

     xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(3,xbe_lib,xbe_usr,cct);
   } else if (flag == global.I_BS23_SVEC) {
//   bs23: stage 5
//   this stage computes the updated solvec (using the 2nd order
//   part of bs23) and is required only if the current step
//   is accepted.

     update_bs23(4,h,xbe_lib,xbe_usr,cct,slv,global);

//   final update of non-integrator variables
//   Note that time_present_x has already been updated to (t + h).
     global.time_given_x = slv.time_present_x;

     xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);
   } else if (flag == global.I_BS23_NORM) {
//   compute norm_2 (only for the integrator variables)

     slv.bs23_norm2 = 0.0;
     update_bs23(5,h,xbe_lib,xbe_usr,cct,slv,global);
     slv.bs23_norm2 = sqrt(slv.bs23_norm2);
   }
   return;
} // end of xbe_bs23
// -----------------------------------------------------------------------------
void xbe_bs23_al(
   const int flag,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global) {

   double h;

   h = slv.delt_x;

   if (flag == global.I_BS23_LTE) {
//   bs23: stage 1

     global.time_given_x = slv.time_present_x;
     if (flag_nan(slv.time_present_x)) {
       cout << "xbe_bs23_al: slv.time_present_x is NAN. Halting..." << endl;
       exit(1);
     }
     xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(0,xbe_lib,xbe_usr,cct);
     update_bs23_al(1,h,xbe_lib,xbe_usr,cct,slv,global);

//   bs23: stage 2
     global.time_given_x = slv.time_present_x + slv.bs23_a1*h;

     if (cct.flag_linear_x) {
       solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);
     }

     xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(1,xbe_lib,xbe_usr,cct);
     update_bs23_al(2,h,xbe_lib,xbe_usr,cct,slv,global);

//   bs23: stage 3
     global.time_given_x = slv.time_present_x + slv.bs23_a2*h;

     if (cct.flag_linear_x) {
       solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);
     }

     xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(2,xbe_lib,xbe_usr,cct);
     update_bs23_al(3,h,xbe_lib,xbe_usr,cct,slv,global);

//   bs23: stage 4
     global.time_given_x = slv.time_present_x + h;

     if (cct.flag_linear_x) {
       solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);
     }

     xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

     xbeu_copy_2(3,xbe_lib,xbe_usr,cct);
   } else if (flag == global.I_BS23_SVEC) {
//   bs23: stage 5
//   this stage computes the updated solvec (using the 2nd order
//   part of bs23) and is required only if the current step
//   is accepted.

     update_bs23_al(4,h,xbe_lib,xbe_usr,cct,slv,global);

//   final update of non-integrator variables
     global.time_given_x = slv.time_present_x;

     if (cct.flag_linear_x) {
       solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);
     }
   } else if (flag == global.I_BS23_NORM) {
//   compute norm_2 (only for the integrator variables)

     slv.bs23_norm2 = 0.0;
     update_bs23_al(5,h,xbe_lib,xbe_usr,cct,slv,global);
     slv.bs23_norm2 = sqrt(slv.bs23_norm2);
   }
   return;
} // end of xbe_bs23_al
// -----------------------------------------------------------------------------
void xbe_meuler(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   double h;

   h = slv.delt_x;

   xbeu_copy_1(xbe_lib,xbe_usr,cct);

// meuler: stage 1

   global.time_given_x = slv.time_present_x;
   if (flag_nan(slv.time_present_x)) {
     cout << "xbe_rk4_meuler: slv.time_present_x is NAN. Halting..." << endl;
     exit(1);
   }
   xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(0,xbe_lib,xbe_usr,cct);
   update_meuler(1,h,xbe_lib,xbe_usr,cct,global);

// meuler: stage 2
   global.time_given_x = slv.time_present_x + h;

   xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(1,xbe_lib,xbe_usr,cct);
   update_meuler(2,(0.5*h),xbe_lib,xbe_usr,cct,global);

// final update of non-integrator variables

   xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

   return;
} // end of xbe_meuler
// -----------------------------------------------------------------------------
void xbe_meuler_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global) {

   double h;

   h = slv.delt_x;

   xbeu_copy_1(xbe_lib,xbe_usr,cct);

// meuler: stage 1

   global.time_given_x = slv.time_present_x;
   if (flag_nan(slv.time_present_x)) {
     cout << "xbe_meuler_al: slv.time_present_x is NAN. Halting..." << endl;
     exit(1);
   }
   xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(0,xbe_lib,xbe_usr,cct);
   update_meuler_al(1,h,xbe_lib,xbe_usr,cct,global);

// meuler: stage 2
   global.time_given_x = slv.time_present_x + h;

   if (cct.flag_linear_x) {
     solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
   } else {
     solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,global);
   }

   xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);
   xbeu_copy_2(1,xbe_lib,xbe_usr,cct);
   update_meuler_al(2,(0.5*h),xbe_lib,xbe_usr,cct,global);

// final update of non-integrator variables

   if (cct.flag_linear_x) {
     solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
   } else {
     solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,global);
   }

   return;
} // end of xbe_meuler_al
// -----------------------------------------------------------------------------
void xbe_heun(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   double h;

   h = slv.delt_x;

   xbeu_copy_1(xbe_lib,xbe_usr,cct);

// heun: stage 1

   global.time_given_x = slv.time_present_x;
   if (flag_nan(slv.time_present_x)) {
     cout << "xbe_heun: slv.time_present_x is NAN. Halting..." << endl;
     exit(1);
   }
   xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(0,xbe_lib,xbe_usr,cct);

   update_heun(1,h,xbe_lib,xbe_usr,cct,slv,global);

// heun: stage 2
   global.time_given_x = slv.time_present_x + (slv.heun_a1*h);

   xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbe_evaluate(global.I_INTEGRATE,xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(1,xbe_lib,xbe_usr,cct);
   update_heun(2,h,xbe_lib,xbe_usr,cct,slv,global);

// final update of non-integrator variables

   global.time_given_x = slv.time_present_x + h;

   xbe_evaluate_2(xbe_lib,xbe_usr,xbe_jac,cct,global);

   return;
} // end of xbe_heun
// -----------------------------------------------------------------------------
void xbe_heun_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global) {

   double h;

   h = slv.delt_x;
   xbeu_copy_1(xbe_lib,xbe_usr,cct);

// heun: stage 1

   global.time_given_x = slv.time_present_x;
   if (flag_nan(slv.time_present_x)) {
     cout << "xbe_heun_al: slv.time_present_x is NAN. Halting..." << endl;
     exit(1);
   }
   xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(0,xbe_lib,xbe_usr,cct);

   update_heun_al(1,h,xbe_lib,xbe_usr,cct,slv,global);

// heun: stage 2
   global.time_given_x = slv.time_present_x + (slv.heun_a1*h);

   if (cct.flag_linear_x) {
     solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
   } else {
     solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,global);
   }

   xbe_evaluate_al(xbe_lib,xbe_usr,xbe_jac,cct,global);

   xbeu_copy_2(1,xbe_lib,xbe_usr,cct);
   update_heun_al(2,h,xbe_lib,xbe_usr,cct,slv,global);

// final update of non-integrator variables

   global.time_given_x = slv.time_present_x + h;

   if (cct.flag_linear_x) {
     solve_trns_linear_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
   } else {
     solve_trns_nonlinear_x_al(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,global);
   }

   return;
} // end of xbe_heun_al
// -----------------------------------------------------------------------------
void xbe_evaluate(
   const int flag,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;

   global.flags[global.i_trns] = true;
   global.flags[global.i_explicit] = true;

   if (flag == global.I_INTEGRATE) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       if (xbe_lib[i_xbel].flag_integrate) {
         xbe_evaluate_1(flag,i_xbeu,xbe_lib,xbe_usr,xbe_jac,cct,global);
       }
     }
   } else if (flag == global.I_DELAY) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       if (xbe_lib[i_xbel].flag_delay) {
         xbe_evaluate_1(flag,i_xbeu,xbe_lib,xbe_usr,xbe_jac,cct,global);
       }
     }
   } else if (flag == global.I_EVAL_SRC) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       if (xbe_lib[i_xbel].flag_evaluate) {
         if (xbe_lib[i_xbel].flag_source) {
           xbe_evaluate_1(flag,i_xbeu,xbe_lib,xbe_usr,xbe_jac,cct,global);
         }
       }
     }
   } else if (flag == global.I_EVAL_NONSRC) {
     for (int i_pass=0; i_pass < cct.x_n_pass; i_pass++) {
       for (int i=0; i < cct.x_pass_n_beu[i_pass]; i++) {
         i_xbeu = cct.x_pass_beu[i_pass][i];
         xbe_evaluate_1(flag,i_xbeu,xbe_lib,xbe_usr,xbe_jac,cct,global);
       }
     }
   }
   global.flags[global.i_trns] = false;
   global.flags[global.i_explicit] = false;

   return;
} // end of xbe_evaluate
// -----------------------------------------------------------------------------
void xbe_evaluate_1(
   const int flag,
   const int i_xbeu,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global) {

   int i_xbel;

   i_xbel = xbe_usr[i_xbeu].index_xbel;

   cct_to_xbe_1(i_xbeu,xbe_lib,xbe_usr,cct);

   get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbeu]);

   if (flag != global.I_INTEGRATE) {
     xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
   }
   return;
} // end of xbe_evaluate_1
// -----------------------------------------------------------------------------
void xbe_evaluate_2(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global) {

   xbe_evaluate(global.I_DELAY      ,xbe_lib,xbe_usr,xbe_jac,cct,global);
   xbe_evaluate(global.I_EVAL_SRC   ,xbe_lib,xbe_usr,xbe_jac,cct,global);
   xbe_evaluate(global.I_EVAL_NONSRC,xbe_lib,xbe_usr,xbe_jac,cct,global);

   return;
} // end of xbe_evaluate_2
// -----------------------------------------------------------------------------
void xbe_evaluate_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;

   global.flags[global.i_trns] = true;
   global.flags[global.i_explicit] = true;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       cct_to_xbe_1(i_xbeu,xbe_lib,xbe_usr,cct);
       get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbeu]);
     }
   }
   global.flags[global.i_trns] = false;
   global.flags[global.i_explicit] = false;

   return;
} // end of xbe_evaluate
// -----------------------------------------------------------------------------
void solve_trns_x_common(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   Global &global) {

   for (int i_file=0;i_file < slv.n_outfile;i_file++) {
     if (slv.flag_out_delt_fixed[i_file]) {
       slv.out_tnext[i_file] = slv.out_tstart[i_file];
     }
   }
   slv.time_present_x = global.time_begin;
   slv.time_write     = global.time_begin;

   write_trns(xbe_lib,xbe_usr,xbe_jac,
     cct,slv,cct_file,global);

   slv.delt_x = slv.delt0_x;

   slv.write_iter_n1_x = 0;
   slv.iter_trns_x = -1; 
   global.iter_trns_x = -1;

   if (cct.flag_limit_tstep_x) {
     x_assign_nextbreak_1(xbe_usr,cct,global);
     xbe_find_nextbreak(xbe_lib,xbe_usr,xbe_jac,cct,global);

     get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
   }
   return;
} // end of solve_trns_x_common
// -----------------------------------------------------------------------------
void update_rk4(
   const int stage,
   const double h1,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;
   int var_flag,var_number;
   int n_f1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;

       for (int i_func=0; i_func < n_f1; i_func++) {
         var_flag   = xbe_lib[i_xbel].ddt_varflag  [i_func];
         var_number = xbe_lib[i_xbel].ddt_varnumber[i_func];

         if (var_flag == global.I_XVR) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_vr[var_number] =
             xbe_usr[i_xbeu].val_vr_0[var_number] + h1*xbe_usr[i_xbeu].f0[i_func];
           } if (stage == 2) {
             xbe_usr[i_xbeu].val_vr[var_number] =
             xbe_usr[i_xbeu].val_vr_0[var_number] + h1*xbe_usr[i_xbeu].f1[i_func];
           } if (stage == 3) {
             xbe_usr[i_xbeu].val_vr[var_number] =
             xbe_usr[i_xbeu].val_vr_0[var_number] + h1*xbe_usr[i_xbeu].f2[i_func];
           } if (stage == 4) {
             xbe_usr[i_xbeu].val_vr[var_number] =
             xbe_usr[i_xbeu].val_vr_0[var_number] +
             h1*(xbe_usr[i_xbeu].f0[i_func] +
               2*(xbe_usr[i_xbeu].f1[i_func] + xbe_usr[i_xbeu].f2[i_func]) +
               xbe_usr[i_xbeu].f3[i_func]);
           }
         } else if (var_flag == global.I_XAUX) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_aux[var_number] =
             xbe_usr[i_xbeu].val_aux_0[var_number] + h1*xbe_usr[i_xbeu].f0[i_func];
           } if (stage == 2) {
             xbe_usr[i_xbeu].val_aux[var_number] =
             xbe_usr[i_xbeu].val_aux_0[var_number] + h1*xbe_usr[i_xbeu].f1[i_func];
           } if (stage == 3) {
             xbe_usr[i_xbeu].val_aux[var_number] =
             xbe_usr[i_xbeu].val_aux_0[var_number] + h1*xbe_usr[i_xbeu].f2[i_func];
           } if (stage == 4) {
             xbe_usr[i_xbeu].val_aux[var_number] =
             xbe_usr[i_xbeu].val_aux_0[var_number] +
             h1*(xbe_usr[i_xbeu].f0[i_func] +
               2*(xbe_usr[i_xbeu].f1[i_func] + xbe_usr[i_xbeu].f2[i_func]) +
               xbe_usr[i_xbeu].f3[i_func]);
           }
         } else {
           cout << "update_rk4: var_flag=" << var_flag
             << " does not make sense. Halting..." << endl;
           exit(1);
         }
       }
       xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
     }
   }

   return;
} // end of update_rk4
// -----------------------------------------------------------------------------
void update_rk4_al(
   const int stage,
   const double h1,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;
   int var_flag,var_number;
   int n_f1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;

       for (int i_func=0; i_func < n_f1; i_func++) {
         var_flag   = xbe_lib[i_xbel].ddt_varflag  [i_func];
         var_number = xbe_lib[i_xbel].ddt_varnumber[i_func];

         if (var_flag == global.I_XVR) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_vr[var_number] =
             xbe_usr[i_xbeu].val_vr_0[var_number] + h1*xbe_usr[i_xbeu].f0[i_func];
           } if (stage == 2) {
             xbe_usr[i_xbeu].val_vr[var_number] =
             xbe_usr[i_xbeu].val_vr_0[var_number] + h1*xbe_usr[i_xbeu].f1[i_func];
           } if (stage == 3) {
             xbe_usr[i_xbeu].val_vr[var_number] =
             xbe_usr[i_xbeu].val_vr_0[var_number] + h1*xbe_usr[i_xbeu].f2[i_func];
           } if (stage == 4) {
             xbe_usr[i_xbeu].val_vr[var_number] =
             xbe_usr[i_xbeu].val_vr_0[var_number] +
             h1*(xbe_usr[i_xbeu].f0[i_func] +
               2*(xbe_usr[i_xbeu].f1[i_func] + xbe_usr[i_xbeu].f2[i_func]) +
               xbe_usr[i_xbeu].f3[i_func]);
           }
           xbe_usr[i_xbeu].val_vr_u[var_number] =
           xbe_usr[i_xbeu].val_vr  [var_number];
         } else if (var_flag == global.I_XAUX) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_aux[var_number] =
             xbe_usr[i_xbeu].val_aux_0[var_number] + h1*xbe_usr[i_xbeu].f0[i_func];
           } if (stage == 2) {
             xbe_usr[i_xbeu].val_aux[var_number] =
             xbe_usr[i_xbeu].val_aux_0[var_number] + h1*xbe_usr[i_xbeu].f1[i_func];
           } if (stage == 3) {
             xbe_usr[i_xbeu].val_aux[var_number] =
             xbe_usr[i_xbeu].val_aux_0[var_number] + h1*xbe_usr[i_xbeu].f2[i_func];
           } if (stage == 4) {
             xbe_usr[i_xbeu].val_aux[var_number] =
             xbe_usr[i_xbeu].val_aux_0[var_number] +
             h1*(xbe_usr[i_xbeu].f0[i_func] +
               2*(xbe_usr[i_xbeu].f1[i_func] + xbe_usr[i_xbeu].f2[i_func]) +
               xbe_usr[i_xbeu].f3[i_func]);
           }
           xbe_usr[i_xbeu].val_aux_u[var_number] =
           xbe_usr[i_xbeu].val_aux  [var_number];
         } else {
           cout << "update_rk4_al: var_flag=" << var_flag
             << " does not make sense. Halting..." << endl;
           exit(1);
         }
       }
       xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
     }
   }

   return;
} // end of update_rk4_al
// -----------------------------------------------------------------------------
void update_rkf45(
   const int stage,
   const double h,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   double a;
   int i_xbeu,i_xbel;
   int var_flag,var_number;
   int n_f1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;

       for (int i_func=0; i_func < n_f1; i_func++) {
         if (stage == 7) {
           a = slv.rkf45_e0*xbe_usr[i_xbeu].f0[i_func]
             + slv.rkf45_e2*xbe_usr[i_xbeu].f2[i_func]
             + slv.rkf45_e3*xbe_usr[i_xbeu].f3[i_func]
             + slv.rkf45_e4*xbe_usr[i_xbeu].f4[i_func]
             + slv.rkf45_e5*xbe_usr[i_xbeu].f5[i_func];
           a= a*h;
           slv.rkf45_norm2 += a*a;
           continue;
         }

         var_flag   = xbe_lib[i_xbel].ddt_varflag  [i_func];
         var_number = xbe_lib[i_xbel].ddt_varnumber[i_func];

         if (var_flag == global.I_XVR) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*(slv.rkf45_b10*xbe_usr[i_xbeu].f0[i_func]);
           } else if (stage == 2) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.rkf45_b20*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b21*xbe_usr[i_xbeu].f1[i_func]);
           } else if (stage == 3) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.rkf45_b30*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b31*xbe_usr[i_xbeu].f1[i_func]
                 + slv.rkf45_b32*xbe_usr[i_xbeu].f2[i_func]);
           } else if (stage == 4) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.rkf45_b40*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b41*xbe_usr[i_xbeu].f1[i_func]
                 + slv.rkf45_b42*xbe_usr[i_xbeu].f2[i_func]
                 + slv.rkf45_b43*xbe_usr[i_xbeu].f3[i_func]);
           } else if (stage == 5) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.rkf45_b50*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b51*xbe_usr[i_xbeu].f1[i_func]
                 + slv.rkf45_b52*xbe_usr[i_xbeu].f2[i_func]
                 + slv.rkf45_b53*xbe_usr[i_xbeu].f3[i_func]
                 + slv.rkf45_b54*xbe_usr[i_xbeu].f4[i_func]);
           } else if (stage == 6) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.rkf45_4_g0*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_4_g2*xbe_usr[i_xbeu].f2[i_func]
                 + slv.rkf45_4_g3*xbe_usr[i_xbeu].f3[i_func]
                 + slv.rkf45_4_g4*xbe_usr[i_xbeu].f4[i_func]);
           }
         } else if (var_flag == global.I_XAUX) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*(slv.rkf45_b10*xbe_usr[i_xbeu].f0[i_func]);
           } else if (stage == 2) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.rkf45_b20*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b21*xbe_usr[i_xbeu].f1[i_func]);
           } else if (stage == 3) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.rkf45_b30*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b31*xbe_usr[i_xbeu].f1[i_func]
                 + slv.rkf45_b32*xbe_usr[i_xbeu].f2[i_func]);
           } else if (stage == 4) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.rkf45_b40*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b41*xbe_usr[i_xbeu].f1[i_func]
                 + slv.rkf45_b42*xbe_usr[i_xbeu].f2[i_func]
                 + slv.rkf45_b43*xbe_usr[i_xbeu].f3[i_func]);
           } else if (stage == 5) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.rkf45_b50*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b51*xbe_usr[i_xbeu].f1[i_func]
                 + slv.rkf45_b52*xbe_usr[i_xbeu].f2[i_func]
                 + slv.rkf45_b53*xbe_usr[i_xbeu].f3[i_func]
                 + slv.rkf45_b54*xbe_usr[i_xbeu].f4[i_func]);
           } else if (stage == 6) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.rkf45_4_g0*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_4_g2*xbe_usr[i_xbeu].f2[i_func]
                 + slv.rkf45_4_g3*xbe_usr[i_xbeu].f3[i_func]
                 + slv.rkf45_4_g4*xbe_usr[i_xbeu].f4[i_func]);
           }
         } else {
           cout << "update_rkf45: var_flag=" << var_flag
             << " does not make sense. Halting..." << endl;
           exit(1);
         }
       }
       xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
     }
   }
   return;
} // end of update_rkf45
// -----------------------------------------------------------------------------
void update_rkf45_al(
   const int stage,
   const double h,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   double a;
   int i_xbeu,i_xbel;
   int var_flag,var_number;
   int n_f1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;

       for (int i_func=0; i_func < n_f1; i_func++) {
         if (stage == 7) {
           a = slv.rkf45_e0*xbe_usr[i_xbeu].f0[i_func]
             + slv.rkf45_e2*xbe_usr[i_xbeu].f2[i_func]
             + slv.rkf45_e3*xbe_usr[i_xbeu].f3[i_func]
             + slv.rkf45_e4*xbe_usr[i_xbeu].f4[i_func]
             + slv.rkf45_e5*xbe_usr[i_xbeu].f5[i_func];
           a= a*h;
           slv.rkf45_norm2 += a*a;
           continue;
         }

         var_flag   = xbe_lib[i_xbel].ddt_varflag  [i_func];
         var_number = xbe_lib[i_xbel].ddt_varnumber[i_func];

         if (var_flag == global.I_XVR) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*(slv.rkf45_b10*xbe_usr[i_xbeu].f0[i_func]);
           } else if (stage == 2) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.rkf45_b20*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b21*xbe_usr[i_xbeu].f1[i_func]);
           } else if (stage == 3) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.rkf45_b30*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b31*xbe_usr[i_xbeu].f1[i_func]
                 + slv.rkf45_b32*xbe_usr[i_xbeu].f2[i_func]);
           } else if (stage == 4) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.rkf45_b40*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b41*xbe_usr[i_xbeu].f1[i_func]
                 + slv.rkf45_b42*xbe_usr[i_xbeu].f2[i_func]
                 + slv.rkf45_b43*xbe_usr[i_xbeu].f3[i_func]);
           } else if (stage == 5) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.rkf45_b50*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b51*xbe_usr[i_xbeu].f1[i_func]
                 + slv.rkf45_b52*xbe_usr[i_xbeu].f2[i_func]
                 + slv.rkf45_b53*xbe_usr[i_xbeu].f3[i_func]
                 + slv.rkf45_b54*xbe_usr[i_xbeu].f4[i_func]);
           } else if (stage == 6) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.rkf45_4_g0*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_4_g2*xbe_usr[i_xbeu].f2[i_func]
                 + slv.rkf45_4_g3*xbe_usr[i_xbeu].f3[i_func]
                 + slv.rkf45_4_g4*xbe_usr[i_xbeu].f4[i_func]);
           }
           xbe_usr[i_xbeu].val_vr_u[var_number] =
           xbe_usr[i_xbeu].val_vr  [var_number];
         } else if (var_flag == global.I_XAUX) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*(slv.rkf45_b10*xbe_usr[i_xbeu].f0[i_func]);
           } else if (stage == 2) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.rkf45_b20*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b21*xbe_usr[i_xbeu].f1[i_func]);
           } else if (stage == 3) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.rkf45_b30*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b31*xbe_usr[i_xbeu].f1[i_func]
                 + slv.rkf45_b32*xbe_usr[i_xbeu].f2[i_func]);
           } else if (stage == 4) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.rkf45_b40*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b41*xbe_usr[i_xbeu].f1[i_func]
                 + slv.rkf45_b42*xbe_usr[i_xbeu].f2[i_func]
                 + slv.rkf45_b43*xbe_usr[i_xbeu].f3[i_func]);
           } else if (stage == 5) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.rkf45_b50*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_b51*xbe_usr[i_xbeu].f1[i_func]
                 + slv.rkf45_b52*xbe_usr[i_xbeu].f2[i_func]
                 + slv.rkf45_b53*xbe_usr[i_xbeu].f3[i_func]
                 + slv.rkf45_b54*xbe_usr[i_xbeu].f4[i_func]);
           } else if (stage == 6) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.rkf45_4_g0*xbe_usr[i_xbeu].f0[i_func]
                 + slv.rkf45_4_g2*xbe_usr[i_xbeu].f2[i_func]
                 + slv.rkf45_4_g3*xbe_usr[i_xbeu].f3[i_func]
                 + slv.rkf45_4_g4*xbe_usr[i_xbeu].f4[i_func]);
           }
           xbe_usr[i_xbeu].val_aux_u[var_number] =
           xbe_usr[i_xbeu].val_aux  [var_number];
         } else {
           cout << "update_rkf45: var_flag=" << var_flag
             << " does not make sense. Halting..." << endl;
           exit(1);
         }
       }
       xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
     }
   }
   return;
} // end of update_rkf45_al
// -----------------------------------------------------------------------------
void update_bs23(
   const int stage,
   const double h,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   double a;
   int i_xbeu,i_xbel;
   int var_flag,var_number;
   int n_f1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;

       for (int i_func=0; i_func < n_f1; i_func++) {
         if (stage == 5) {
           a = slv.bs23_e0*xbe_usr[i_xbeu].f0[i_func]
             + slv.bs23_e1*xbe_usr[i_xbeu].f1[i_func]
             + slv.bs23_e2*xbe_usr[i_xbeu].f2[i_func]
             + slv.bs23_e3*xbe_usr[i_xbeu].f3[i_func];
           slv.bs23_norm2 += a*a;
           continue;
         }
         var_flag   = xbe_lib[i_xbel].ddt_varflag  [i_func];
         var_number = xbe_lib[i_xbel].ddt_varnumber[i_func];

         if (var_flag == global.I_XVR) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*(slv.bs23_b10*xbe_usr[i_xbeu].f0[i_func]);
           } else if (stage == 2) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.bs23_b21*xbe_usr[i_xbeu].f1[i_func]);
           } else if (stage == 3) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.bs23_b30*xbe_usr[i_xbeu].f0[i_func]
                 + slv.bs23_b31*xbe_usr[i_xbeu].f1[i_func]
                 + slv.bs23_b32*xbe_usr[i_xbeu].f2[i_func]);
           } else if (stage == 4) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.bs23_2_g0*xbe_usr[i_xbeu].f0[i_func]
                 + slv.bs23_2_g1*xbe_usr[i_xbeu].f1[i_func]
                 + slv.bs23_2_g2*xbe_usr[i_xbeu].f2[i_func]
                 + slv.bs23_2_g3*xbe_usr[i_xbeu].f3[i_func]);
           }
         } else if (var_flag == global.I_XAUX) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*(slv.bs23_b10*xbe_usr[i_xbeu].f0[i_func]);
           } else if (stage == 2) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.bs23_b21*xbe_usr[i_xbeu].f1[i_func]);
           } else if (stage == 3) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.bs23_b30*xbe_usr[i_xbeu].f0[i_func]
                 + slv.bs23_b31*xbe_usr[i_xbeu].f1[i_func]
                 + slv.bs23_b32*xbe_usr[i_xbeu].f2[i_func]);
           } else if (stage == 4) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.bs23_2_g0*xbe_usr[i_xbeu].f0[i_func]
                 + slv.bs23_2_g1*xbe_usr[i_xbeu].f1[i_func]
                 + slv.bs23_2_g2*xbe_usr[i_xbeu].f2[i_func]
                 + slv.bs23_2_g3*xbe_usr[i_xbeu].f3[i_func]);
           }
         } else {
           cout << "update_bs23: var_flag=" << var_flag
             << " does not make sense. Halting..." << endl;
           exit(1);
         }
       }
       xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
     }
   }
   return;
} // end of update_bs23
// -----------------------------------------------------------------------------
void update_bs23_al(
   const int stage,
   const double h,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   double a;
   int i_xbeu,i_xbel;
   int var_flag,var_number;
   int n_f1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;

       for (int i_func=0; i_func < n_f1; i_func++) {
         if (stage == 5) {
           a = slv.bs23_e0*xbe_usr[i_xbeu].f0[i_func]
             + slv.bs23_e1*xbe_usr[i_xbeu].f1[i_func]
             + slv.bs23_e2*xbe_usr[i_xbeu].f2[i_func]
             + slv.bs23_e3*xbe_usr[i_xbeu].f3[i_func];
           slv.bs23_norm2 += a*a;
           continue;
         }
         var_flag   = xbe_lib[i_xbel].ddt_varflag  [i_func];
         var_number = xbe_lib[i_xbel].ddt_varnumber[i_func];

         if (var_flag == global.I_XVR) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*(slv.bs23_b10*xbe_usr[i_xbeu].f0[i_func]);
           } else if (stage == 2) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.bs23_b21*xbe_usr[i_xbeu].f1[i_func]);
           } else if (stage == 3) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.bs23_b30*xbe_usr[i_xbeu].f0[i_func]
                 + slv.bs23_b31*xbe_usr[i_xbeu].f1[i_func]
                 + slv.bs23_b32*xbe_usr[i_xbeu].f2[i_func]);
           } else if (stage == 4) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.bs23_2_g0*xbe_usr[i_xbeu].f0[i_func]
                 + slv.bs23_2_g1*xbe_usr[i_xbeu].f1[i_func]
                 + slv.bs23_2_g2*xbe_usr[i_xbeu].f2[i_func]
                 + slv.bs23_2_g3*xbe_usr[i_xbeu].f3[i_func]);
           }
           xbe_usr[i_xbeu].val_vr_u[var_number] =
           xbe_usr[i_xbeu].val_vr  [var_number];
         } else if (var_flag == global.I_XAUX) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*(slv.bs23_b10*xbe_usr[i_xbeu].f0[i_func]);
           } else if (stage == 2) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.bs23_b21*xbe_usr[i_xbeu].f1[i_func]);
           } else if (stage == 3) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.bs23_b30*xbe_usr[i_xbeu].f0[i_func]
                 + slv.bs23_b31*xbe_usr[i_xbeu].f1[i_func]
                 + slv.bs23_b32*xbe_usr[i_xbeu].f2[i_func]);
           } else if (stage == 4) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.bs23_2_g0*xbe_usr[i_xbeu].f0[i_func]
                 + slv.bs23_2_g1*xbe_usr[i_xbeu].f1[i_func]
                 + slv.bs23_2_g2*xbe_usr[i_xbeu].f2[i_func]
                 + slv.bs23_2_g3*xbe_usr[i_xbeu].f3[i_func]);
           }
           xbe_usr[i_xbeu].val_aux_u[var_number] =
           xbe_usr[i_xbeu].val_aux  [var_number];
         } else {
           cout << "update_bs23: var_flag=" << var_flag
             << " does not make sense. Halting..." << endl;
           exit(1);
         }
       }
       xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
     }
   }
   return;
} // end of update_bs23_al
// -----------------------------------------------------------------------------
void update_meuler(
   const int stage,
   const double h1,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;
   int var_flag,var_number;
   int n_f1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;

       for (int i_func=0; i_func < n_f1; i_func++) {
         var_flag   = xbe_lib[i_xbel].ddt_varflag  [i_func];
         var_number = xbe_lib[i_xbel].ddt_varnumber[i_func];

         if (var_flag == global.I_XVR) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h1*xbe_usr[i_xbeu].f0[i_func];
           } if (stage == 2) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h1*( xbe_usr[i_xbeu].f0[i_func] 
                  + xbe_usr[i_xbeu].f1[i_func]);
           }
         } else if (var_flag == global.I_XAUX) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h1*xbe_usr[i_xbeu].f0[i_func];
           } if (stage == 2) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h1*( xbe_usr[i_xbeu].f0[i_func] 
                  + xbe_usr[i_xbeu].f1[i_func]);
           }
         } else {
           cout << "update_meuler: var_flag=" << var_flag
             << " does not make sense. Halting..." << endl;
           exit (1);
         }
       }
       xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
     }
   }
   return;
} // end of update_meuler
// -----------------------------------------------------------------------------
void update_meuler_al(
   const int stage,
   const double h1,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;
   int var_flag,var_number;
   int n_f1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;

       for (int i_func=0; i_func < n_f1; i_func++) {
         var_flag   = xbe_lib[i_xbel].ddt_varflag  [i_func];
         var_number = xbe_lib[i_xbel].ddt_varnumber[i_func];

         if (var_flag == global.I_XVR) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h1*xbe_usr[i_xbeu].f0[i_func];
           } if (stage == 2) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h1*( xbe_usr[i_xbeu].f0[i_func] 
                  + xbe_usr[i_xbeu].f1[i_func]);
           }
           xbe_usr[i_xbeu].val_vr_u[var_number] =
           xbe_usr[i_xbeu].val_vr  [var_number];
         } else if (var_flag == global.I_XAUX) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h1*xbe_usr[i_xbeu].f0[i_func];
           } if (stage == 2) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h1*( xbe_usr[i_xbeu].f0[i_func] 
                  + xbe_usr[i_xbeu].f1[i_func]);
           }
           xbe_usr[i_xbeu].val_aux_u[var_number] =
           xbe_usr[i_xbeu].val_aux  [var_number];
         } else {
           cout << "update_meuler: var_flag=" << var_flag
             << " does not make sense. Halting..." << endl;
           exit (1);
         }
       }
       xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
     }
   }
   return;
} // end of update_meuler_al
// -----------------------------------------------------------------------------
void update_heun(
   const int stage,
   const double h,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   int i_xbeu,i_xbel;
   int var_flag,var_number;
   int n_f1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;

       for (int i_func=0; i_func < n_f1; i_func++) {
         var_flag   = xbe_lib[i_xbel].ddt_varflag  [i_func];
         var_number = xbe_lib[i_xbel].ddt_varnumber[i_func];

         if (var_flag == global.I_XVR) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + slv.heun_b10*h*xbe_usr[i_xbeu].f0[i_func];
           } if (stage == 2) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.heun_g0*xbe_usr[i_xbeu].f0[i_func] 
                 + slv.heun_g1*xbe_usr[i_xbeu].f1[i_func]);
           }
         } else if (var_flag == global.I_XAUX) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + slv.heun_b10*h*xbe_usr[i_xbeu].f0[i_func];
           } if (stage == 2) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.heun_g0*xbe_usr[i_xbeu].f0[i_func] 
                 + slv.heun_g1*xbe_usr[i_xbeu].f1[i_func]);
           }
         } else {
           cout << "update_heun: var_flag=" << var_flag
             << " does not make sense. Halting..." << endl;
           exit (1);
         }
       }
       xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
     }
   }
   return;
} // end of update_heun
// -----------------------------------------------------------------------------
void update_heun_al(
   const int stage,
   const double h,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   int i_xbeu,i_xbel;
   int var_flag,var_number;
   int n_f1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;

       for (int i_func=0; i_func < n_f1; i_func++) {
         var_flag   = xbe_lib[i_xbel].ddt_varflag  [i_func];
         var_number = xbe_lib[i_xbel].ddt_varnumber[i_func];

         if (var_flag == global.I_XVR) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + slv.heun_b10*h*xbe_usr[i_xbeu].f0[i_func];
           } if (stage == 2) {
             xbe_usr[i_xbeu].val_vr[var_number]
             = xbe_usr[i_xbeu].val_vr_0[var_number]
             + h*( slv.heun_g0*xbe_usr[i_xbeu].f0[i_func] 
                 + slv.heun_g1*xbe_usr[i_xbeu].f1[i_func]);
           }
           xbe_usr[i_xbeu].val_vr_u[var_number] =
           xbe_usr[i_xbeu].val_vr  [var_number];
         } else if (var_flag == global.I_XAUX) {
           if (stage == 1) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + slv.heun_b10*h*xbe_usr[i_xbeu].f0[i_func];
           } if (stage == 2) {
             xbe_usr[i_xbeu].val_aux[var_number]
             = xbe_usr[i_xbeu].val_aux_0[var_number]
             + h*( slv.heun_g0*xbe_usr[i_xbeu].f0[i_func] 
                 + slv.heun_g1*xbe_usr[i_xbeu].f1[i_func]);
           }
           xbe_usr[i_xbeu].val_aux_u[var_number] =
           xbe_usr[i_xbeu].val_aux  [var_number];
         } else {
           cout << "update_heun: var_flag=" << var_flag
             << " does not make sense. Halting..." << endl;
           exit (1);
         }
       }
       xbe_to_cct_op_1(i_xbeu,xbe_lib,xbe_usr,cct);
     }
   }
   return;
} // end of update_heun_al
// -----------------------------------------------------------------------------
void xbeu_copy_1(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct) {

   int i_xbeu,i_xbel;
   int n_vr1,n_aux1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     n_vr1 = xbe_lib[i_xbel].n_vr;
     for (int i=0; i < n_vr1; i++) {
       xbe_usr[i_xbeu].val_vr_0[i] = xbe_usr[i_xbeu].val_vr[i];
     }
     n_aux1 = xbe_lib[i_xbel].n_aux;
     for (int i=0; i < n_aux1; i++) {
       xbe_usr[i_xbeu].val_aux_0[i] = xbe_usr[i_xbeu].val_aux[i];
     }
   }
   return;
} //end of xbe_copy_1
// -----------------------------------------------------------------------------
void xbeu_copy_2(
   const int flag,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct) {

// copy f to f0 or f1 ...

   int i_xbeu,i_xbel;
   int n_func1;

   if (flag == 0) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       n_func1 = xbe_lib[i_xbel].n_f;
       for (int i=0; i < n_func1; i++) {
         xbe_usr[i_xbeu].f0[i] = xbe_usr[i_xbeu].f[i];
       }
     }
   } else if (flag == 1) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       n_func1 = xbe_lib[i_xbel].n_f;
       for (int i=0; i < n_func1; i++) {
         xbe_usr[i_xbeu].f1[i] = xbe_usr[i_xbeu].f[i];
       }
     }
   } else if (flag == 2) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       n_func1 = xbe_lib[i_xbel].n_f;
       for (int i=0; i < n_func1; i++) {
         xbe_usr[i_xbeu].f2[i] = xbe_usr[i_xbeu].f[i];
       }
     }
   } else if (flag == 3) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       n_func1 = xbe_lib[i_xbel].n_f;
       for (int i=0; i < n_func1; i++) {
         xbe_usr[i_xbeu].f3[i] = xbe_usr[i_xbeu].f[i];
       }
     }
   } else if (flag == 4) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       n_func1 = xbe_lib[i_xbel].n_f;
       for (int i=0; i < n_func1; i++) {
         xbe_usr[i_xbeu].f4[i] = xbe_usr[i_xbeu].f[i];
       }
     }
   } else if (flag == 5) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       n_func1 = xbe_lib[i_xbel].n_f;
       for (int i=0; i < n_func1; i++) {
         xbe_usr[i_xbeu].f5[i] = xbe_usr[i_xbeu].f[i];
       }
     }
   }
   return;
} //end of xbe_copy_2
// -----------------------------------------------------------------------------
void copy_func_to_old_x(
   const int flag_1,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;
   int n_g1;

   if (flag_1 == global.I_COPY_0_TO_1) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       n_g1 = xbe_lib[i_xbel].n_g;
       for (int i=0; i < n_g1; i++) {
         xbe_usr[i_xbeu].g_old_1[i] = xbe_usr[i_xbeu].g[i];
       }
     }
   } else if (flag_1 == global.I_COPY_1_TO_2) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       n_g1 = xbe_lib[i_xbel].n_g;
       for (int i=0; i < n_g1; i++) {
         xbe_usr[i_xbeu].g_old_2[i] = xbe_usr[i_xbeu].g_old_1[i];
       }
     }
   } else if (flag_1 == global.I_COPY_2_TO_0) {
     for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
       i_xbel = xbe_usr[i_xbeu].index_xbel;
       n_g1 = xbe_lib[i_xbel].n_g;
       for (int i=0; i < n_g1; i++) {
         xbe_usr[i_xbeu].g[i] = xbe_usr[i_xbeu].g_old_2[i];
       }
     }
   }
   return;
} //end of copy_func_to_old_x
// -----------------------------------------------------------------------------
void solve_startup(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   solve_startup_x_exp(xbe_lib,xbe_usr,xbe_jac,
     cct,slv,cct_file,global);

   return;
} // end of solve_startup
// -----------------------------------------------------------------------------
void solve_startup_x_exp(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   int i_xbeu,i_xbel;

   cout << "solve_startup_x_exp starting..." << endl;

   assign_const_1<bool>(global.flags,false);
   global.time_given_x = 0.0;

   global.flags[global.i_startup] = true;
   global.flags[global.i_explicit] = true;

// step 1: assign start-up values to integrate and source type elements

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_integrate) {
       get_xbe_1(i_xbeu,xbe_lib,xbe_usr,xbe_jac,cct,global);
     }
   }
   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     if (xbe_lib[i_xbel].flag_evaluate) {
       if (xbe_lib[i_xbel].flag_source) {
         get_xbe_1(i_xbeu,xbe_lib,xbe_usr,xbe_jac,cct,global);
       }
     }
   }
// step 2: process the remaining elements

   for (int i_pass=0; i_pass < cct.x_n_pass; i_pass++) {
     for (int i=0; i < cct.x_pass_n_beu[i_pass]; i++) {
       i_xbeu = cct.x_pass_beu[i_pass][i];
       get_xbe_1(i_xbeu,xbe_lib,xbe_usr,xbe_jac,cct,global);
     }
   }
   global.flags[global.i_startup] = false;
   global.flags[global.i_explicit] = false;

   write_dc_startup(xbe_lib,xbe_usr,xbe_jac,
     cct,slv,cct_file,global);

   cout << "solve_startup_x_exp over" << endl;

   return;
} // end of solve_startup_x_exp
// -----------------------------------------------------------------------------
void solve_startup_x_imp(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   init_sol_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,smat,global);
   smat.mat_startup_1_x(xbe_lib,xbe_usr,cct,cct_file);
   form_solvec_x(xbe_lib,xbe_usr,smat,cct);

   assign_const_1<bool>(global.flags,false);
   global.time_given_x = 0.0;

   if (cct.flag_linear_x) {
     solve_startup_linear_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
   } else {
     solve_startup_nonlinear_x(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,global);
   }
   write_dc_startup(xbe_lib,xbe_usr,xbe_jac,
     cct,slv,cct_file,global);

   return;
} // end of solve_startup_x_imp
// -----------------------------------------------------------------------------
void solve_startup_linear_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   form_jac_rhs_startup_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);

   form_solvec_x(xbe_lib,xbe_usr,smat,cct);

   negative_double_1(smat.m_x.n_row,smat.rhs_m_x);

   smat.w_x.allocate_1(smat.m_x.n_nz,smat.m_x.n_row,smat.m_x.n_col);
   smat.mo_x.allocate_1(smat.m_x.n_row);

   solve_jac_1_x(smat,slv,global);

   copy_array_1<double>(smat.m_x.n_row,smat.svec_orig_x,smat.delsvec_x);

   add_arrays_1<double>(smat.m_x.n_row,smat.delsvec_x,smat.svec_x);
   dcmp_solvec_x(xbe_lib,xbe_usr,smat,cct);

   return;
} // end of solve_startup_linear_x
// -----------------------------------------------------------------------------
void solve_startup_nonlinear_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   int i_newt;
   bool flag_nan_1;

   check_convergence_count_x(slv);
   slv.get_dmp();

   smat.w_x.allocate_1(smat.m_x.n_nz,smat.m_x.n_row,smat.m_x.n_col);
   smat.mo_x.allocate_1(smat.m_x.n_row);

   for (i_newt=0; i_newt < slv.x_nr_itermax_a; i_newt++) {
     cout << "solve_startup_nonlinear_x: i_newt = " << i_newt << endl;
     slv.iter_newton = i_newt;
     form_jac_rhs_startup_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);

     negative_double_1(smat.m_x.n_row,smat.rhs_m_x);

     if (i_newt == 0) {
       solve_jac_1_x(smat,slv,global);
     } else {
       solve_jac_2_x(smat,slv,global);
     }
     copy_array_1<double>(smat.m_x.n_row,smat.svec_orig_x,smat.delsvec_x);

     if ((slv.x_nr_flag_dmp_a) && (i_newt <= slv.x_nr_dmp_itermax_a)) {
       mult_array_1<double>(smat.m_x.n_row,smat.delsvec_x,slv.x_nr_dmp_k_a);
     }
     add_arrays_1<double>(smat.m_x.n_row,smat.delsvec_x,smat.svec_x);
     check_array_for_nan_2(smat.m_x.n_row,smat.svec_x,flag_nan_1);
     if (flag_nan_1) {
       cout << "solve_startup_nonlinear_x: svec_x has a NaN" << endl;
       slv.flag_nr_converged = false;
       slv.flag_nr_norm_large = true;
       break;
     } else {
       dcmp_solvec_x(xbe_lib,xbe_usr,smat,cct);
       check_convergence_x(smat,slv);
       if (slv.flag_nr_converged) {
         cout << "solve_startup_nonlinear_x: NR converged." << endl;
         break;
       }
     }
   }
   return;
} // end of solve_startup_nonlinear_x
// -----------------------------------------------------------------------------
void solve_trns(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   if (cct.flag_x_explicit) {
     solve_trns_x_exp(xbe_lib,xbe_usr,xbe_jac,
       slv,cct,cct_file,smat,global);
   } else {
     solve_trns_x_imp(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,cct_file,global);
   }

   return;
} // end of solve_trns
// -----------------------------------------------------------------------------
void solve_trns_x_imp(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   init_sol_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,smat,global);

   form_solvec_x(xbe_lib,xbe_usr,smat,cct);
   form_map_xbeuvr_1(smat,cct);

   smat.mat_trns_1_x(xbe_lib,xbe_usr,cct,global,cct_file);
   smat.w_x.allocate_1(smat.m_x.n_nz,smat.m_x.n_row,smat.m_x.n_col);
   smat.mo_x.allocate_1(smat.m_x.n_row);

   slv.x_algo_be0 = slv.x_algo_be || slv.x_algo_be_auto;
   slv.x_algo_trz0 = slv.x_algo_trz || slv.x_algo_trz_auto;
   slv.x_algo_auto = slv.x_algo_be_auto || slv.x_algo_trz_auto;

   for (int i_file=0; i_file < slv.n_outfile; i_file++) {
     if (slv.flag_out_delt_fixed[i_file]) {
       slv.out_tnext[i_file] = slv.out_tstart[i_file];
     }
   }

   if (slv.x_algo_be) {
     solve_trns_x_be(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,cct_file,global);
   } else if (slv.x_algo_trz) {
     solve_trns_x_trz(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,cct_file,global);
   } else if (slv.x_algo_be_auto) {
     solve_trns_x_be_auto(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,cct_file,global);
   } else if (slv.x_algo_trz_auto) {
     solve_trns_x_trz_auto(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,cct_file,global);
   } else if (slv.x_algo_trbdf2) {
     solve_trns_x_trbdf2(xbe_lib,xbe_usr,xbe_jac,
       smat,cct,slv,cct_file,global);
   }

   return;
} // end of solve_trns_x_imp
// -----------------------------------------------------------------------------
void solve_trns_x_be(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   double time_next_1;
   bool flag_tend_reached;

   flag_tend_reached = false;

   slv.write_iter_n1_x = 0;
   slv.iter_trns_x = -1;

   slv.time_present_x = global.time_begin;
   slv.time_write     = global.time_begin;

   write_trns(xbe_lib,xbe_usr,xbe_jac,
     cct,slv,cct_file,global);

   slv.delt_x = slv.delt0_x;

   if (cct.flag_limit_tstep_x) {
     x_assign_nextbreak_1(xbe_usr,cct,global);
     xbe_find_nextbreak(xbe_lib,xbe_usr,xbe_jac,cct,global);

     get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
   }
   copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);

   while (!flag_tend_reached) {
     slv.iter_trns_x++;
     global.iter_trns_x = slv.iter_trns_x;
     write_iter_x(slv);

     if (slv.iter_trns_x != 0) {
       if (cct.flag_modulo_x) {
         xbe_modulo_implicit(cct,xbe_lib,xbe_usr,smat);
       }
     }
     global.flags[global.i_trns_first_x] = (slv.iter_trns_x == 0);

     slv.time_next_x  = slv.time_present_x + slv.delt_x;
     global.time_given_x = slv.time_next_x;
     if (flag_nan(slv.time_next_x)) {
       cout << "solve_trns_x_be: slv.time_next_x is NAN. Halting..." << endl;
       exit(1);
     }

     slv.trns_constants_2_x();
     form_solvec_x(xbe_lib,xbe_usr,smat,cct);

     if (cct.flag_linear_x) {
       solve_trns_linear_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);

       if (!slv.flag_nr_converged) {
         if (slv.flag_write_solution) {
           write_solution(xbe_lib,xbe_usr,slv,cct);
         }
         cout << "solve_trns_x_be: N-R iterations did not converge." << endl;
         cout << "  iter_trns_x =" << slv.iter_trns_x << endl;
         cout << "  time ="
           << scientific << setprecision(6)
           << slv.time_present_x << endl;
         cout << "  Halting..." << endl;
         exit (1);
       }
     }
     slv.time_present_x = slv.time_next_x;
     slv.time_write     = slv.time_next_x;

     write_trns(xbe_lib,xbe_usr,xbe_jac,
       cct,slv,cct_file,global);

     copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);

     slv.delt_x = slv.delt0_x;
     slv.delt_x = max(slv.delt_x,slv.delt_min_x);

     if (cct.flag_limit_tstep_x) {
       get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
     }
     time_next_1 = slv.time_present_x + slv.delt_x;
     if (time_next_1 >= global.time_end) {
       slv.delt_x = global.time_end - slv.time_present_x + slv.delt_small;
       slv.delt_x = max(slv.delt_x,slv.delt_min_x);
     }
     if (cct.flag_save_history_x) {
       save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);
     }
     if (slv.flag_limit_iter_trns) {
       if (slv.iter_trns_x > slv.itmax_trns) {
         cout << "solve_trns_x_be: itmax_trns exceeded." << endl;
         break;
       }
     }
     if ((slv.time_present_x+slv.delt_small) >= global.time_end) {
       flag_tend_reached = true;
     }
   }
   return;
} // end of solve_trns_x_be
// -----------------------------------------------------------------------------
void solve_trns_x_trz(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   double time_next_1;
   bool flag_tend_reached;

   flag_tend_reached = false;

   slv.write_iter_n1_x = 0;
   slv.iter_trns_x = -1;

   slv.time_present_x = global.time_begin;
   slv.time_write     = global.time_begin;

   write_trns(xbe_lib,xbe_usr,xbe_jac,
     cct,slv,cct_file,global);

   slv.delt_x = slv.delt0_x;

   if (cct.flag_limit_tstep_x) {
     x_assign_nextbreak_1(xbe_usr,cct,global);
     xbe_find_nextbreak(xbe_lib,xbe_usr,xbe_jac,cct,global);

     get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
   }
   copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);

   slv.time_next_x = slv.time_present_x;
   global.time_given_x = slv.time_next_x;
   if (flag_nan(slv.time_next_x)) {
     cout << "solve_trns_x_trz: slv.time_next_x is NAN. Halting..." << endl;
     exit(1);
   }
   slv.trns_constants_2_x();
   find_functions_trns_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);

   while (!flag_tend_reached) {
     slv.iter_trns_x++;
     global.iter_trns_x = slv.iter_trns_x;
     write_iter_x(slv);

     if (slv.iter_trns_x != 0) {
       if (cct.flag_modulo_x) {
         xbe_modulo_implicit(cct,xbe_lib,xbe_usr,smat);
       }
     }
     global.flags[global.i_trns_first_x] = (slv.iter_trns_x == 0);

     slv.time_next_x  = slv.time_present_x + slv.delt_x;
     global.time_given_x = slv.time_next_x;

     form_solvec_x(xbe_lib,xbe_usr,smat,cct);

     if (cct.flag_linear_x) {
       find_functions_trns_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);
     }
     copy_func_to_old_x(global.I_COPY_0_TO_1,xbe_lib,xbe_usr,cct,global);
     slv.trns_constants_2_x();

     if (cct.flag_linear_x) {
       solve_trns_linear_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);

       if (!slv.flag_nr_converged) {
         if (slv.flag_write_solution) {
           write_solution(xbe_lib,xbe_usr,slv,cct);
         }
         cout << "solve_trns_x_trz: N-R iterations did not converge." << endl;
         cout << "  iter_trns_x =" << slv.iter_trns_x << endl;
         cout << "  time ="
           << scientific << setprecision(6)
           << slv.time_present_x << endl;
         cout << "  Halting..." << endl;
         exit (1);
       }
     }
     slv.time_present_x = slv.time_next_x;
     slv.time_write     = slv.time_next_x;

     write_trns(xbe_lib,xbe_usr,xbe_jac,
       cct,slv,cct_file,global);

     copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);

//   find the new time step if necessary

     slv.delt_x = slv.delt0_x;
     slv.delt_x = max(slv.delt_x,slv.delt_min_x);

     if (cct.flag_limit_tstep_x) {
       get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
     }
     time_next_1 = slv.time_present_x + slv.delt_x;
     if (time_next_1 >= global.time_end) {
       slv.delt_x = global.time_end - slv.time_present_x + slv.delt_small;
       slv.delt_x = max(slv.delt_x,slv.delt_min_x);
     }
     if (cct.flag_save_history_x) {
       save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);
     }

     if (slv.flag_limit_iter_trns) {
       if (slv.iter_trns_x > slv.itmax_trns) {
         cout << "solve_trns_e_trz: itmax_trns exceeded." << endl;
         break;
       }
     }
     if ((slv.time_present_x+slv.delt_small) >= global.time_end) {
       flag_tend_reached = true;
     }
   }
   return;
} // end of solve_trns_x_trz
// -----------------------------------------------------------------------------
void solve_trns_x_be_auto(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   double time_next_1;
   bool flag_tend_reached,flag_repeat_step;
   int iter_stepred;

   flag_tend_reached = false;
   flag_repeat_step = false;

   slv.write_iter_n1_x = 0;
   slv.iter_trns_x = -1;

   slv.time_present_x = global.time_begin;
   slv.time_write     = global.time_begin;

   write_trns(xbe_lib,xbe_usr,xbe_jac,
     cct,slv,cct_file,global);

   slv.delt_x = slv.delt0_x;

   if (cct.flag_limit_tstep_x) {
     x_assign_nextbreak_1(xbe_usr,cct,global);
     xbe_find_nextbreak(xbe_lib,xbe_usr,xbe_jac,cct,global);

     get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
   }
   copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);

   while (!flag_tend_reached) {
     if (!flag_repeat_step) {
       slv.iter_trns_x++;
       global.iter_trns_x = slv.iter_trns_x;
       write_iter_x(slv);

       if (slv.iter_trns_x != 0) {
         if (cct.flag_modulo_x) {
           xbe_modulo_implicit(cct,xbe_lib,xbe_usr,smat);
         }
       }
       global.flags[global.i_trns_first_x] = (slv.iter_trns_x == 0);
       iter_stepred = 0;
     }
     flag_repeat_step = false;

     slv.time_next_x  = slv.time_present_x + slv.delt_x;
     global.time_given_x = slv.time_next_x;
     if (flag_nan(slv.time_next_x)) {
       cout << "solve_trns_x_be_auto: slv.time_next_x is NAN. Halting..." << endl;
       exit(1);
     }

     slv.trns_constants_2_x();

     form_solvec_x(xbe_lib,xbe_usr,smat,cct);

     if (cct.flag_linear_x) {
       solve_trns_linear_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);

       if (!slv.flag_nr_converged) {
         if (slv.delt_x == slv.delt_min_x) {
           cout << "solve_trns_x_be_auto: no convergence even with" << endl;
           cout << "  the smallest time step." << endl;
           cout << "iter_trns_x=" << slv.iter_trns_x
                << ", time =" << slv.time_present_x << endl;
           cout << "  Halting..." << endl;
           exit (1);
         } else {
           iter_stepred++;
           if (iter_stepred > slv.itmax_stepred) {
             cout << "solve_trns_x_be_auto: iter_stepred has exceeded" << endl;
             cout << "  itmax_stepred." << endl;
             cout << "  iter_stepred=" << iter_stepred << endl;
             cout << "  itmax_stepred=" << slv.itmax_stepred << endl;
             cout << "  Halting..." << endl;
             exit (1);
           }
           slv.delt_x = slv.factor_stepdec*slv.delt_x;
           slv.delt_x = max(slv.delt_x,slv.delt_min_x);

           copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);
           flag_repeat_step = true;
         }
       } else {
         if (slv.delt_x != slv.delt_max_x) {
           slv.delt_x = slv.factor_stepinc*slv.delt_x;
           slv.delt_x = min(slv.delt_x,slv.delt_max_x);
         }
         flag_repeat_step = false;
       }
     }
     if (!flag_repeat_step) {
       slv.time_present_x = slv.time_next_x;
       slv.time_write     = slv.time_next_x;

       write_trns(xbe_lib,xbe_usr,xbe_jac,
         cct,slv,cct_file,global);

       copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);

       slv.delt_x = max(slv.delt_x,slv.delt_min_x);

       if (cct.flag_limit_tstep_x) {
         get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
       }
       time_next_1 = slv.time_present_x + slv.delt_x;
       if (time_next_1 >= global.time_end) {
         slv.delt_x = global.time_end - slv.time_present_x + slv.delt_small;
         slv.delt_x = max(slv.delt_x,slv.delt_min_x);
       }
       if (cct.flag_save_history_x) {
         save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);
       }
       if (slv.flag_limit_iter_trns) {
         if (slv.iter_trns_x > slv.itmax_trns) {
           cout << "solve_trns_e_be_auto: itmax_trns exceeded." << endl;
           break;
         }
       }
       if ((slv.time_present_x+slv.delt_small) >= global.time_end) {
         flag_tend_reached = true;
       }
     }
   }
   return;
} // end of solve_trns_x_be_auto
// -----------------------------------------------------------------------------
void solve_trns_x_trz_auto(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   double time_next_1;
   bool flag_tend_reached,flag_repeat_step;
   int iter_stepred;

   flag_tend_reached = false;
   flag_repeat_step = false;

   slv.write_iter_n1_x = 0;
   slv.iter_trns_x = -1;

   slv.time_present_x = global.time_begin;
   slv.time_write     = global.time_begin;

   write_trns(xbe_lib,xbe_usr,xbe_jac,
     cct,slv,cct_file,global);

   slv.delt_x = slv.delt0_x;

   if (cct.flag_limit_tstep_x) {
     x_assign_nextbreak_1(xbe_usr,cct,global);
     xbe_find_nextbreak(xbe_lib,xbe_usr,xbe_jac,cct,global);

     get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
   }
   copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);

   slv.time_next_x = slv.time_present_x;
   global.time_given_x = slv.time_next_x;
   if (flag_nan(slv.time_next_x)) {
     cout << "solve_trns_x_trz_auto: slv.time_next_x is NAN. Halting..." << endl;
     exit(1);
   }
   slv.trns_constants_2_x();
   find_functions_trns_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);

   while (!flag_tend_reached) {
     if (!flag_repeat_step) {
       slv.iter_trns_x++;
       global.iter_trns_x = slv.iter_trns_x;
       write_iter_x(slv);

       if (slv.iter_trns_x != 0) {
         if (cct.flag_modulo_x) {
           xbe_modulo_implicit(cct,xbe_lib,xbe_usr,smat);
         }
       }
       global.flags[global.i_trns_first_x] = (slv.iter_trns_x == 0);
       iter_stepred = 0;
     }
     flag_repeat_step = false;

     slv.time_next_x  = slv.time_present_x + slv.delt_x;
     global.time_given_x = slv.time_next_x;

     form_solvec_x(xbe_lib,xbe_usr,smat,cct);

     if (cct.flag_linear_x) {
       find_functions_trns_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);
     }
     copy_func_to_old_x(global.I_COPY_0_TO_1,xbe_lib,xbe_usr,cct,global);
     slv.trns_constants_2_x();

     if (cct.flag_linear_x) {
       solve_trns_linear_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);
     } else {
       solve_trns_nonlinear_x(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);

       if (!slv.flag_nr_converged) {
         if (slv.delt_x == slv.delt_min_x) {
           cout << "solve_trns_x_trzauto: no convergence even with" << endl;
           cout << "  the smallest time step." << endl;
           cout << "iter_trns_x=" << slv.iter_trns_x
                << ", time =" << slv.time_present_x << endl;
           cout << "  Halting..." << endl;
           exit (1);
         } else {
           iter_stepred++;
           if (iter_stepred > slv.itmax_stepred) {
             cout << "solve_trns_x_trzauto: iter_stepred has exceeded" << endl;
             cout << "  itmax_stepred." << endl;
             cout << "  iter_stepred=" << iter_stepred << endl;
             cout << "  itmax_stepred=" << slv.itmax_stepred << endl;
             cout << "  Halting..." << endl;
             exit (1);
           }
           slv.delt_x = slv.factor_stepdec*slv.delt_x;
           slv.delt_x = max(slv.delt_x,slv.delt_min_x);

           copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);

           flag_repeat_step = true;
         }
       } else {
         if (slv.delt_x != slv.delt_max_x) {
           slv.delt_x = slv.factor_stepinc*slv.delt_x;
           slv.delt_x = min(slv.delt_x,slv.delt_max_x);
         }
         flag_repeat_step = false;
       }
     }
     if (!flag_repeat_step) {
       slv.time_present_x = slv.time_next_x;
       slv.time_write     = slv.time_next_x;

       write_trns(xbe_lib,xbe_usr,xbe_jac,
         cct,slv,cct_file,global);

       copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);

       slv.delt_x = max(slv.delt_x,slv.delt_min_x);

       if (cct.flag_limit_tstep_x) {
         get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
       }
       time_next_1 = slv.time_present_x + slv.delt_x;
       if (time_next_1 >= global.time_end) {
         slv.delt_x = global.time_end - slv.time_present_x + slv.delt_small;
         slv.delt_x = max(slv.delt_x,slv.delt_min_x);
       }
       if (cct.flag_save_history_x) {
         save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);
       }
       if (slv.flag_limit_iter_trns) {
         if (slv.iter_trns_x > slv.itmax_trns) {
           cout << "solve_trns_e_trz_auto: itmax_trns exceeded." << endl;
           break;
         }
       }
       if ((slv.time_present_x+slv.delt_small) >= global.time_end) {
         flag_tend_reached = true;
       }
     }
   }
   return;
} // end of solve_trns_x_trz_auto
// -----------------------------------------------------------------------------
void solve_trns_x_trbdf2(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global) {

   double delt_tmp,t_present_tmp,time_next_1;
   bool flag_tend_reached,flag_repeat_step;
   int iter_trbdf2,iter_trz,iter_bdf2,n_reject;

   slv.x_algo_trz0 = true;

   flag_tend_reached = false;
   flag_repeat_step = false;

   slv.write_iter_n1_x = 0;
   slv.iter_trns_x = -1;

   slv.time_present_x = global.time_begin;
   slv.time_write     = global.time_begin;

   write_trns(xbe_lib,xbe_usr,xbe_jac,
     cct,slv,cct_file,global);

   slv.delt_x = slv.delt0_x;

   if (cct.flag_limit_tstep_x) {
     x_assign_nextbreak_1(xbe_usr,cct,global);
     xbe_find_nextbreak(xbe_lib,xbe_usr,xbe_jac,cct,global);

     get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);
   }
   copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);

   slv.trns_constants_2_x();
   find_functions_trns_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);

   while (!flag_tend_reached) {
     if (!flag_repeat_step) {
       slv.iter_trns_x++;
       global.iter_trns_x = slv.iter_trns_x;
       write_iter_x(slv);

       if (slv.iter_trns_x != 0) {
         if (cct.flag_modulo_x) {
           xbe_modulo_implicit(cct,xbe_lib,xbe_usr,smat);
         }
       }
       global.flags[global.i_trns_first_x] = (slv.iter_trns_x == 0);

       iter_trbdf2 = 0;
       iter_trz = 0;
       iter_bdf2 = 0;
       n_reject = 0;
     }
     flag_repeat_step = false;

     iter_trbdf2++;
     if (iter_trbdf2 > slv.itmax_trbdf2) {
       cout << "solve_trns_x_trbdf2: iter_trbdf2 has exceeded itmax_trbdf2" << endl;
       cout << "  iter_trbdf2=" << iter_trbdf2 << endl;
       cout << "  itmax_trbdf2=" << slv.itmax_trbdf2 << endl;
       cout << "  n_reject=" << n_reject << endl;
       cout << "  time_given=" << global.time_given_x << endl;
       cout << "  Halting..." << endl; exit (1);
     }
     delt_tmp = slv.delt_x;
     t_present_tmp = slv.time_present_x;

     slv.delt_x = slv.bank_gamma*delt_tmp;
     slv.time_next_x  = slv.time_present_x + slv.delt_x;
     global.time_given_x = slv.time_next_x;
     if (flag_nan(slv.time_next_x)) {
       cout << "solve_trns_x_trbdf2: slv.time_next_x is NAN. Halting..." << endl;
       exit(1);
     }

//   solve trz part:

     slv.x_algo_trz0 = true;
     slv.x_algo_bdf2 = false;
     iter_trz++;

     copy_func_to_old_x(global.I_COPY_0_TO_1,xbe_lib,xbe_usr,cct,global);
     slv.trns_constants_2_x();

     form_solvec_x(xbe_lib,xbe_usr,smat,cct);

     if (cct.flag_linear_x) {

       solve_trns_linear_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);

       find_functions_trns_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);
     } else {
       solve_trns_nonlinear_x(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,global);

       if (!slv.flag_nr_converged) {
         if (slv.delt_x == slv.delt_min_x) {
           cout << "solve_trns_x_trbdf2: no convergence even with" << endl;
           cout << "  the smallest time step." << endl;
           cout << "  iter_trz=" << iter_trz << endl;
           cout << "  iter_trbdf2=" << iter_trbdf2 << endl;
           cout << "  n_reject=" << n_reject << endl;
           cout << "  iter_trns_x="   << slv.iter_trns_x << endl;
           cout << "  Halting..." << endl;
           exit (1);
         }

         slv.delt_x = slv.factor_stepdec*delt_tmp;

         slv.delt_x = max(slv.delt_x,slv.delt_min_x);

         copy_array_1<double>(smat.n_solvec_x,smat.svec_old_1_x,smat.svec_x);
         flag_repeat_step = true;
       }
     }

     if (!flag_repeat_step) {
       copy_array_1<double>(smat.n_solvec_x,smat.svec_old_1_x,smat.svec_old_2_x);
       copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);

       copy_func_to_old_x(global.I_COPY_1_TO_2,xbe_lib,xbe_usr,cct,global);
       copy_func_to_old_x(global.I_COPY_0_TO_1,xbe_lib,xbe_usr,cct,global);

//     solve bdf2 part:

       slv.x_algo_trz0 = false;
       slv.x_algo_bdf2 = true;
       iter_bdf2++;

       slv.time_next_x  = t_present_tmp + delt_tmp;
       global.time_given_x = slv.time_next_x;
       slv.delt_x = delt_tmp;
       slv.trns_constants_2_x();

       form_solvec_x(xbe_lib,xbe_usr,smat,cct);

       if (cct.flag_linear_x) {

         solve_trns_linear_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,slv,global);

         find_functions_trns_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);
       } else {
         solve_trns_nonlinear_x(xbe_lib,xbe_usr,xbe_jac,
           smat,cct,slv,global);

         if (!slv.flag_nr_converged) {
           if (slv.delt_x == slv.delt_min_x) {
             cout << "solve_trns_x_trbdf2: no convergence even with" << endl;
             cout << "  the smallest time step." << endl;
             cout << "  iter_trz=" << iter_trz << endl;
             cout << "  iter_trbdf2=" << iter_trbdf2 << endl;
             cout << "  n_reject=" << n_reject << endl;
             cout << "  iter_trns_x=" << slv.iter_trns_x << endl;
             cout << "  Halting..." << endl;
             exit (1);
           }
           slv.delt_x = slv.factor_stepdec*delt_tmp;
           slv.delt_x = max(slv.delt_x,slv.delt_min_x);

           copy_array_1<double>(smat.n_solvec_x,smat.svec_old_2_x,smat.svec_old_1_x);
           copy_array_1<double>(smat.n_solvec_x,smat.svec_old_2_x,smat.svec_x);

           flag_repeat_step = true;
         }
       }
     }
     if (!flag_repeat_step) {
//     check whether the solution should be accepted

       trzbdf2_1_x(xbe_usr,smat,slv);

       slv.delt_x = slv.delt_new_x;

       slv.delt_x = max(slv.delt_x,slv.delt_min_x);
       slv.delt_x = min(slv.delt_x,slv.delt_max_x);

       if (!slv.flag_accept_sol) {
         n_reject++;

         copy_array_1<double>(smat.n_solvec_x,smat.svec_old_2_x,smat.svec_old_1_x);
         copy_array_1<double>(smat.n_solvec_x,smat.svec_old_2_x,smat.svec_x);

         copy_func_to_old_x(global.I_COPY_2_TO_0,xbe_lib,xbe_usr,cct,global);
         dcmp_solvec_x(xbe_lib,xbe_usr,smat,cct);

         flag_repeat_step = true;
       }
     }
     if (!flag_repeat_step) {
       slv.time_present_x = slv.time_next_x;
       slv.time_write     = slv.time_next_x;

       write_trns(xbe_lib,xbe_usr,xbe_jac,
         cct,slv,cct_file,global);

       copy_array_1<double>(smat.n_solvec_x,smat.svec_x,smat.svec_old_1_x);

       slv.delt_x = max(slv.delt_x,slv.delt_min_x);

       if (cct.flag_limit_tstep_x) {

         get_tnext_x(xbe_lib,xbe_usr,xbe_jac,slv,cct,global);

       }
       time_next_1 = slv.time_present_x + slv.delt_x;
       if (time_next_1 >= global.time_end) {
         slv.delt_x = global.time_end - slv.time_present_x + slv.delt_small;
         slv.delt_x = max(slv.delt_x,slv.delt_min_x);
       }
       if (cct.flag_save_history_x) {
         save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);
       }
       if (slv.flag_limit_iter_trns) {
         if (slv.iter_trns_x > slv.itmax_trns) {
           cout << "solve_trns_x_trbdf2: itmax_trns exceeded." << endl;
           break;
         }
       }
       if ((slv.time_present_x+slv.delt_small) >= global.time_end) {
         flag_tend_reached = true;
       }
     }
   }
   return;
} // end of solve_trns_x_trbdf2
// -----------------------------------------------------------------------------
void solve_trns_linear_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   bool flag_write;

   form_jac_rhs_trns_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);
   add_trns_terms_x(xbe_usr,smat,slv);

   negative_double_1(smat.m_x.n_row,smat.rhs_m_x);

   if ((slv.iter_trns_x == 0) && (!slv.x_algo_bdf2)) {
     solve_jac_1_x(smat,slv,global);
   } else {
     solve_jac_2_x(smat,slv,global);
   }

   copy_array_1<double>(smat.m_x.n_row,smat.svec_orig_x,smat.delsvec_x);

   add_arrays_1<double>(smat.m_x.n_row,smat.delsvec_x,smat.svec_x);
   dcmp_solvec_x(xbe_lib,xbe_usr,smat,cct);

   xbe_reset_1(true,xbe_lib,xbe_usr,xbe_jac,cct,smat,global);
// save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);

   return;
} // end of solve_trns_linear_x
// -----------------------------------------------------------------------------
void solve_trns_nonlinear_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   bool flag_nan_1,flag_write;
   int i_newt;

   check_convergence_count_x(slv);
   slv.get_dmp();

   for (i_newt=0; i_newt < slv.x_nr_itermax_a; i_newt++) {
     slv.iter_newton = i_newt;

     form_jac_rhs_trns_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);
     add_trns_terms_x(xbe_usr,smat,slv);

     negative_double_1(smat.m_x.n_row,smat.rhs_m_x);

     if (slv.iter_newton == 0) {
       if ((slv.iter_trns_x == 0) && (!slv.x_algo_bdf2)) {
         solve_jac_1_x(smat,slv,global);
       } else {
         solve_jac_2_x(smat,slv,global);
       }
     } else {
       solve_jac_2_x(smat,slv,global);
     }

     copy_array_1<double>(smat.m_x.n_row,smat.svec_orig_x,smat.delsvec_x);
     if ((slv.x_nr_flag_dmp_a) && (i_newt <= slv.x_nr_dmp_itermax_a)) {
       mult_array_1<double>(smat.m_x.n_row,smat.delsvec_x,slv.x_nr_dmp_k_a);
     }
     add_arrays_1<double>(smat.m_x.n_row,smat.delsvec_x,smat.svec_x);
     check_array_for_nan_2(smat.m_x.n_row,smat.svec_x,flag_nan_1);
     if (flag_nan_1) {
       cout << "solve_trns_nonlinear_x: svec_x has a NaN" << endl;
       slv.flag_nr_converged = false;
       slv.flag_nr_norm_large = true;
       break;
     } else {
       dcmp_solvec_x(xbe_lib,xbe_usr,smat,cct);
       check_convergence_x(smat,slv);
       if (slv.flag_nr_converged) {
         break;
       }
     }
   }
// save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);

   return;
} // end of solve_trns_nonlinear_x
// -----------------------------------------------------------------------------
void solve_trns_linear_x_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

   bool flag_write;

// solve algebraic loop equations in only xbe's/explicit trns case.

   form_jac_rhs_trns_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);

   form_solvec_x(xbe_lib,xbe_usr,smat,cct);

   negative_double_1(smat.m_x.n_row,smat.rhs_m_x);

   if (slv.iter_trns_x == 0) {
     solve_jac_1_x(smat,slv,global);
   } else {
     solve_jac_2_x(smat,slv,global);
   }

   copy_array_1<double>(smat.m_x.n_row,smat.svec_orig_x,smat.delsvec_x);

   add_arrays_1<double>(smat.m_x.n_row,smat.delsvec_x,smat.svec_x);
   dcmp_solvec_x(xbe_lib,xbe_usr,smat,cct);

   xbe_reset_1(true,xbe_lib,xbe_usr,xbe_jac,cct,smat,global);
// save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);

   return;
} // end of solve_trns_linear_x_al
// -----------------------------------------------------------------------------
void solve_trns_nonlinear_x_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global) {

// solve algebraic loop equations in only xbe's/explicit trns case.

   bool flag_nan_1,flag_write;
   int i_newt;

   check_convergence_count_x(slv);
   slv.get_dmp();

   for (i_newt=0; i_newt < slv.x_nr_itermax_a; i_newt++) {
     slv.iter_newton = i_newt;

     form_jac_rhs_trns_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);

     negative_double_1(smat.m_x.n_row,smat.rhs_m_x);

     if (slv.iter_newton == 0) {
       if (slv.iter_trns_x == 0) {
         solve_jac_1_x(smat,slv,global);
       } else {
         solve_jac_2_x(smat,slv,global);
       }
     } else {
       solve_jac_2_x(smat,slv,global);
     }

     copy_array_1<double>(smat.m_x.n_row,smat.svec_orig_x,smat.delsvec_x);
     if ((slv.x_nr_flag_dmp_a) && (i_newt <= slv.x_nr_dmp_itermax_a)) {
       mult_array_1<double>(smat.m_x.n_row,smat.delsvec_x,slv.x_nr_dmp_k_a);
     }
     add_arrays_1<double>(smat.m_x.n_row,smat.delsvec_x,smat.svec_x);
     check_array_for_nan_2(smat.m_x.n_row,smat.svec_x,flag_nan_1);
     if (flag_nan_1) {
       cout << "solve_trns_nonlinear_x_al: svec_x has a NaN" << endl;
       slv.flag_nr_converged = false;
       slv.flag_nr_norm_large = true;
       break;
     } else {
       dcmp_solvec_x(xbe_lib,xbe_usr,smat,cct);
       check_convergence_x(smat,slv);
       if (slv.flag_nr_converged) {
         break;
       }
     }
   }
   if (!slv.flag_nr_converged) {
     cout << "solve_trns_nonlinear_x_al: NR did not converge. Halting..."
       << endl; exit(1);
   }
   xbe_reset_1(true,xbe_lib,xbe_usr,xbe_jac,cct,smat,global);
// save_history_x(xbe_lib,xbe_usr,xbe_jac,cct,global);

   return;
} // end of solve_trns_nonlinear_x_al
// -----------------------------------------------------------------------------
void solve_jac_1_x(
   SysMat &smat,
   SolveBlocks &slv,
   Global &global) {

   solve_jac(
     global.I_NO_RHS_MTOW,
     global.I_MAT_MTOW,
     global.I_GAUSS1A,
     global.I_NO_GAUSS2,
     global.I_NO_LU_BSUB,
     slv.flag_debug_gauss1,slv.flag_debug_gauss2,
     slv.gauss_epsln,slv.zero_piv,
     smat.m_x,smat.w_x,smat.mo_x,smat.rhs_m_x,
     smat.rhs_w_x,smat.svec_w_x,smat.svec_orig_x,global);
   solve_jac(
     global.I_RHS_MTOW,
     global.I_NO_MAT_MTOW,
     global.I_NO_GAUSS1,
     global.I_GAUSS2,
     global.I_LU_BSUB,
     slv.flag_debug_gauss1,slv.flag_debug_gauss2,
     slv.gauss_epsln,slv.zero_piv,
     smat.m_x,smat.w_x,smat.mo_x,smat.rhs_m_x,
     smat.rhs_w_x,smat.svec_w_x,smat.svec_orig_x,global);

   return;
} // end of solve_jac_1_x
// -----------------------------------------------------------------------------
void solve_jac_2_x(
   SysMat &smat,
   SolveBlocks &slv,
   Global &global) {

   solve_jac(
     global.I_RHS_MTOW,
     global.I_NO_MAT_MTOW,
     global.I_NO_GAUSS1,
     global.I_GAUSS2,
     global.I_LU_BSUB,
     slv.flag_debug_gauss1,slv.flag_debug_gauss2,
     slv.gauss_epsln,slv.zero_piv,
     smat.m_x,smat.w_x,smat.mo_x,smat.rhs_m_x,
     smat.rhs_w_x,smat.svec_w_x,smat.svec_orig_x,global);

   return;
} // end of solve_jac_2_x
// -----------------------------------------------------------------------------
void solve_jac(
   const int flag_copy_rhs,
   const int flag_copy_m,
   const int flag_call_gauss1,
   const int flag_call_gauss2,
   const int flag_gauss2,
   const bool flag_debug_gauss1,
   const bool flag_debug_gauss2,
   const double gauss_epsln,
   const double zero_piv,
   KnuthMat &m,
   KnuthMat &w,
   MatOp &mo,
   double *rhs_m,
   double *rhs_w,
   double *svec,
   double *svec_orig,
   Global &global) {

   bool flag_error,flag_error_1;

   if (flag_copy_rhs == global.I_RHS_MTOW) {
     copy_array_1<double>(m.n_row,rhs_m,rhs_w);
   }
   if (flag_copy_m == global.I_MAT_MTOW) {
     knuth_copy_1(m,w);
   }
   if (flag_call_gauss1 == global.I_GAUSS1) {
     gauss_1(flag_debug_gauss1,m,w,mo);
   } else if (flag_call_gauss1 == global.I_GAUSS1A) {
     gauss_1a(flag_debug_gauss1,gauss_epsln,m,w,mo);
   }

   if (flag_call_gauss2 == global.I_GAUSS2) {
     gauss_2(flag_gauss2,flag_debug_gauss2,zero_piv,flag_error,global,
       m,w,mo,rhs_w,svec,svec_orig);
     if (flag_error) {
       cout << "solve_jac: error in gauss_2." << endl;
       cout << "  Calling gauss_1a (to reorder)..." << endl;

//     need to allocate w and mo again:
//     Note: memory allocation for mo will not change, but we still need to
//       call delete_1/allocate_1 here because allocate_1 will also initialise
//       values.

       if (w.flag_delete) w.delete_1();
       if (mo.flag_delete) mo.delete_1(m.n_row);

       w.allocate_1(m.n_nz,m.n_row,m.n_col);
       mo.allocate_1(m.n_row);

       knuth_copy_1(m,w);

       gauss_1a(flag_debug_gauss1,gauss_epsln,m,w,mo);
       copy_array_1<double>(m.n_row,rhs_m,rhs_w);
       gauss_2(global.I_LU_BSUB,flag_debug_gauss2,zero_piv,flag_error_1,global,
         m,w,mo,rhs_w,svec,svec_orig);

       if (flag_error_1) {
         cout << "solve_jac: error in gauss_2 again!" << endl;
         cout << "  (even after reordering)" << endl;
         cout << "  Halting..." << endl;
         exit (1);
       }
     }
   }
   return;
} //end of solve_jac
// -----------------------------------------------------------------------------
void form_jac_rhs_startup_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global) {

   cct_to_xbe_all(xbe_lib,xbe_usr,cct);
   mat_startup_2_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);

   return;
} // end of form_jac_rhs_startup_x
// -----------------------------------------------------------------------------
void form_jac_rhs_trns_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global) {

   cct_to_xbe_all(xbe_lib,xbe_usr,cct);
   mat_trns_2_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);

   return;
} // end of form_jac_rhs_trns_x
// -----------------------------------------------------------------------------
void form_jac_rhs_trns_x_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global) {

   cct_to_xbe_all(xbe_lib,xbe_usr,cct);
   mat_trns_2_x_al(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);

   return;
} // end of form_jac_rhs_trns_x_al
// -----------------------------------------------------------------------------
void find_functions_trns_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global) {

   cct_to_xbe_all(xbe_lib,xbe_usr,cct);
   mat_trns_2a_x(xbe_lib,xbe_usr,xbe_jac,smat,cct,global);

   return;
} // end of find_functions_trns_x
// -----------------------------------------------------------------------------
void form_solvec_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SysMat &smat,
   Circuit &cct) {

   int i_svec;
   int i_xbeu,i_xbel,i_xbeu_vr;
   int n_aux1;

   i_svec = smat.offs_xvr;
   for (i_xbeu_vr=0; i_xbeu_vr < cct.n_xbeu_vr; i_xbeu_vr++) {
//   cout << "form_solvec_x: cct.val_xvr[i_xbeu_vr]: "
//      << cct.val_xvr[i_xbeu_vr] << endl;
     smat.svec_x[i_svec] = cct.val_xvr[i_xbeu_vr];
     i_svec++;
   }
   i_svec = smat.offs_xaux;
   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     n_aux1 = xbe_lib[i_xbel].n_aux;

     for (int i=0; i < n_aux1; i++) {
       smat.svec_x[i_svec] = xbe_usr[i_xbeu].val_aux[i];
       i_svec++;
     }
   }
   return;
} // end of form_solvec_x
// -----------------------------------------------------------------------------
void form_map_xbeuvr_1(
   SysMat &smat,
   Circuit &cct) {

   int i_svec;
   int i_xbeu,i_xbel,i_xbeu_vr;

   i_svec = smat.offs_xvr;
   for (i_xbeu_vr=0; i_xbeu_vr < cct.n_xbeu_vr; i_xbeu_vr++) {
     cct.map_xbeuvr_to_svec[i_xbeu_vr] = i_svec;
     i_svec++;
   }
   return;
} // end of form_solvec_x
// -----------------------------------------------------------------------------
void dcmp_solvec_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SysMat &smat,
   Circuit &cct) {

   int i_svec;
   int i_xbeu,i_xbel,i_xbeu_vr;
   int n_aux1;

   i_svec = smat.offs_xvr;
   for (i_xbeu_vr=0; i_xbeu_vr < cct.n_xbeu_vr; i_xbeu_vr++) {
     cct.val_xvr[i_xbeu_vr] = smat.svec_x[i_svec];
     i_svec++;
   }
   cct_to_xbe_all(xbe_lib,xbe_usr,cct);

   i_svec = smat.offs_xaux;
   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     n_aux1 = xbe_lib[i_xbel].n_aux;
     for (int i=0; i < n_aux1; i++) {
       xbe_usr[i_xbeu].val_aux[i] = smat.svec_x[i_svec];
       i_svec++;
     }
   }
   return;
} // end of dcmp_solvec_x
// -----------------------------------------------------------------------------
void mat_startup_2_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;

   global.flags[global.i_startup] = true;
   global.flags[global.i_implicit] = true;
   global.flags[global.i_function] = true;
   global.flags[global.i_jacobian] = true;

   xbe_init_jac_startup_x(xbe_lib,xbe_usr,cct,smat);

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbel]);

     mat_startup_3_x(i_xbeu,i_xbel,xbe_lib,xbe_usr,xbe_jac,smat,global);
   }
   global.flags[global.i_startup] = false;
   global.flags[global.i_implicit] = false;
   global.flags[global.i_function] = false;
   global.flags[global.i_jacobian] = false;

   return;
} // end of mat_startup_2_x
// -----------------------------------------------------------------------------
void mat_startup_3_x(
   const int i_xbeu,
   const int i_xbel,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Global &global) {

   int i_f,i_g,n_f1,n_g1,n_gvar1,k,row0;
   int var_flag,var_number;
   double val;

   if (xbe_lib[i_xbel].flag_integrate) {
     n_f1 = xbe_lib[i_xbel].n_f;

     for (i_f=0; i_f < n_f1; i_f++) {
       row0 = smat.xbe_f_to_row[i_xbeu][i_f];
       if (row0 == -1) {
         cout << "mat_startup_3_x: row0 = -1 is not expected. Halting.." << endl;
         exit(1);
       } else {
         smat.rhs_m_x[row0] = xbe_usr[i_xbeu].h[i_f];
       }
     }
   } else {
     n_g1 = xbe_lib[i_xbel].n_g;

     for (i_g=0; i_g < n_g1; i_g++) {
       n_gvar1 = xbe_lib[i_xbel].n_gvar[i_g];
       for (int i_gvar=0; i_gvar < n_gvar1; i_gvar++) {
         k = smat.map_gvar_to_xbe[i_xbeu][i_g][i_gvar];
         if (k == -1) {
           cout << "mat_startup_3_x: k = -1 is not expected. Halting.." << endl;
           exit(1);
         } else {
           var_flag   = xbe_lib[i_xbel].gvar_flag [i_g][i_gvar];
           var_number = xbe_lib[i_xbel].gvar_index[i_g][i_gvar];

           if (var_flag == global.I_XVR) {
             val = xbe_jac[i_xbel].dgdvr[i_g][var_number];
           } else if (var_flag == global.I_XAUX) {
             val = xbe_jac[i_xbel].dgdaux[i_g][var_number];
           } else {
             cout << "mat_startup_3_x: incorrect value of var_flag." << endl;
             cout << "  i_g = " << i_g << ", i_gavr = " << i_gvar << endl;
             cout << "  xbe_lib is " << xbe_lib[i_xbel].name << ". Halting.." << endl;
             exit(1);
           }
           smat.m_x.val[k] = val;
         }
       }
       row0 = smat.xbe_g_to_row[i_xbeu][i_g];
       if (row0 == -1) {
         cout << "mat_startup_3_x: row0 = -1 is not expected. Halting.." << endl;
         exit(1);
       } else {
         cout << "mat_startup_3_x: row0 = " << row0 << endl;
         smat.rhs_m_x[row0] = xbe_usr[i_xbeu].g[i_g];
       }
     }
   }

   return;
} // end of mat_startup_3_x
// -----------------------------------------------------------------------------
void mat_trns_2_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;

   xbe_init_jac_trns_x(xbe_lib,xbe_usr,cct,smat);

   global.flags[global.i_trns] = true;
   global.flags[global.i_function] = true;
   global.flags[global.i_jacobian] = true;
   global.flags[global.i_implicit] = true;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbel]);
     mat_trns_3_x(i_xbeu,i_xbel,xbe_lib,xbe_usr,xbe_jac,smat,global);
   }
   global.flags[global.i_trns] = false;
   global.flags[global.i_function] = false;
   global.flags[global.i_jacobian] = false;
   global.flags[global.i_implicit] = false;

   return;
} // end of mat_trns_2_x
// -----------------------------------------------------------------------------
void mat_trns_2a_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;

   global.flags[global.i_trns] = true;
   global.flags[global.i_function] = true;
   global.flags[global.i_implicit] = true;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbel]);
     mat_trns_3a_x(i_xbeu,i_xbel,xbe_lib,xbe_usr,smat);
   }
   global.flags[global.i_trns] = false;
   global.flags[global.i_function] = false;
   global.flags[global.i_implicit] = false;

   return;
} // end of mat_trns_2a_x
// -----------------------------------------------------------------------------
void mat_trns_2_x_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global) {

   int i_xbeu,i_xbel;

   global.flags[global.i_trns] = true;

   xbe_init_jac_startup_x(xbe_lib,xbe_usr,cct,smat);

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     if (xbe_lib[i_xbel].flag_integrate) {
       global.flags[global.i_explicit] = true;
       global.flags[global.i_alg_loop] = true;
     } else {
       global.flags[global.i_implicit] = true;
       global.flags[global.i_function] = true;
       global.flags[global.i_jacobian] = true;
     }

     get_xbe(i_xbel,global,xbe_usr[i_xbeu],xbe_jac[i_xbel]);

     if (xbe_lib[i_xbel].flag_integrate) {
       global.flags[global.i_explicit] = false;
       global.flags[global.i_alg_loop] = false;
     } else {
       global.flags[global.i_implicit] = false;
       global.flags[global.i_function] = false;
       global.flags[global.i_jacobian] = false;
     }

     mat_startup_3_x(i_xbeu,i_xbel,xbe_lib,xbe_usr,xbe_jac,smat,global);
   }
   global.flags[global.i_trns] = false;

   return;
} // end of mat_trns_2_x_al
// -----------------------------------------------------------------------------
void mat_trns_3_x(
   const int i_xbeu,
   const int i_xbel,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Global &global) {

   int i_g,n_g1,n_gvar1,k,row0;
   int var_flag,var_number;
   double val;

   n_g1 = xbe_lib[i_xbel].n_g;

   for (i_g=0; i_g < n_g1; i_g++) {
     n_gvar1 = xbe_lib[i_xbel].n_gvar[i_g];
     for (int i_gvar=0; i_gvar < n_gvar1; i_gvar++) {
       k = smat.map_gvar_to_xbe[i_xbeu][i_g][i_gvar];
       if (k == -1) {
         cout << "mat_trns_3_x: k = -1 is not expected. Halting.." << endl;
         exit(1);
       } else {
         var_flag   = xbe_lib[i_xbel].gvar_flag [i_g][i_gvar];
         var_number = xbe_lib[i_xbel].gvar_index[i_g][i_gvar];

         if (var_flag == global.I_XVR) {
           val = xbe_jac[i_xbel].dgdvr[i_g][var_number];
         } else if (var_flag == global.I_XAUX) {
           val = xbe_jac[i_xbel].dgdaux[i_g][var_number];
         } else {
           cout << "mat_trns_3_x: incorrect value of var_flag." << endl;
           cout << "  i_g = " << i_g << ", i_gavr = " << i_gvar << endl;
           cout << "  xbe_lib is " << xbe_lib[i_xbel].name << ". Halting.." << endl;
           exit(1);
         }
         smat.m_x.val[k] = val;
       }
     }
     if (xbe_lib[i_xbel].flag_integrate) {
       row0 = smat.xbe_f_to_row[i_xbeu][i_g];
       if (row0 == -1) {
         cout << "mat_trns_3_x: row0 = -1 is not expected. Halting.." << endl;
         exit(1);
       } else {
         smat.rhs_m_x[row0] = xbe_usr[i_xbeu].g[i_g];
       }
     } else {
       row0 = smat.xbe_g_to_row[i_xbeu][i_g];
       if (row0 == -1) {
         cout << "mat_trns_3_x: row0 = -1 is not expected. Halting.." << endl;
         exit(1);
       } else {
         smat.rhs_m_x[row0] = xbe_usr[i_xbeu].g[i_g];
       }
     }
   }
   return;
} // end of mat_trns_3_x
// -----------------------------------------------------------------------------
void mat_trns_3a_x(
   const int i_xbeu,
   const int i_xbel,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SysMat &smat) {

   int i_g,n_g1,row0;

   n_g1 = xbe_lib[i_xbel].n_g;

   for (i_g=0; i_g < n_g1; i_g++) {
     if (xbe_lib[i_xbel].flag_integrate) {
       row0 = smat.xbe_f_to_row[i_xbeu][i_g];
       if (row0 == -1) {
         cout << "mat_trns_3a_x: row0 = -1 is not expected. Halting.." << endl;
         exit(1);
       } else {
         smat.rhs_m_x[row0] = xbe_usr[i_xbeu].g[i_g];
       }
     } else {
       row0 = smat.xbe_g_to_row[i_xbeu][i_g];
       if (row0 == -1) {
         cout << "mat_trns_3a_x: row0 = -1 is not expected. Halting.." << endl;
         exit(1);
       } else {
         smat.rhs_m_x[row0] = xbe_usr[i_xbeu].g[i_g];
       }
     }
   }
   return;
} // end of mat_trns_3a_x
// -----------------------------------------------------------------------------
void add_trns_terms_x(
   vector<XbeUsr> &xbe_usr,
   SysMat &smat,
   SolveBlocks &slv) {

   int i_xbeu,i_f;
   double x_0,x_1,x_2,f_1;
   int i_rhs,i_var,pntr;

   if (slv.x_algo_be0) {
     for (i_rhs=0; i_rhs < smat.m_x.n_row; i_rhs++) {
       if (smat.xbe_rhs_ddt_flag[i_rhs]) {
         i_var = smat.xbe_rhs_ddt_varnumber[i_rhs];
         pntr = smat.xbe_rhs_ddt_pntr[i_rhs];
         x_0 = smat.svec_x      [i_var];
         x_1 = smat.svec_old_1_x[i_var];
         smat.rhs_m_x[i_rhs] +=
           -slv.beuler_1_x*(x_0-x_1);

         smat.m_x.val[pntr] += -slv.beuler_1_x;

       }
     }
   } else if (slv.x_algo_trz0) {
     for (i_rhs=0; i_rhs < smat.m_x.n_row; i_rhs++) {
       if (smat.xbe_rhs_ddt_flag[i_rhs]) {
         i_var = smat.xbe_rhs_ddt_varnumber[i_rhs];
         pntr = smat.xbe_rhs_ddt_pntr[i_rhs];
         x_0 = smat.svec_x      [i_var];
         x_1 = smat.svec_old_1_x[i_var];

         i_xbeu = smat.xbe_rhs_ddt_i_xbeu[i_rhs];
         i_f = smat.xbe_rhs_ddt_i_f[i_rhs];
         f_1 = xbe_usr[i_xbeu].g_old_1[i_f];

         smat.rhs_m_x[i_rhs] +=
           -slv.trz_1_x*(x_0-x_1) + f_1;
         smat.m_x.val[pntr] += -slv.trz_1_x;
       }
     }
   } else if (slv.x_algo_bdf2) {
     for (i_rhs=0; i_rhs < smat.m_x.n_row; i_rhs++) {
       if (smat.xbe_rhs_ddt_flag[i_rhs]) {
         i_var = smat.xbe_rhs_ddt_varnumber[i_rhs];
         pntr = smat.xbe_rhs_ddt_pntr[i_rhs];
         x_0 = smat.svec_x      [i_var];
         x_1 = smat.svec_old_1_x[i_var];
         x_2 = smat.svec_old_2_x[i_var];

         smat.rhs_m_x[i_rhs] +=
           -(slv.bdf2_1_x*x_0 - slv.bdf2_2_x*x_1 + slv.bdf2_3_x*x_2);
         smat.m_x.val[pntr] += -slv.bdf2_1_x;
       }
     }
   }
// cout << "add_trns_terms_x ending" << endl;

   return;
} // end of add_trns_terms_x
// -----------------------------------------------------------------------------
void xbe_init_jac_startup_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SysMat &smat) {

   int i_xbeu,i_xbel,k,i_g,n_g1,n_gvar1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     if (!xbe_lib[i_xbel].flag_integrate) {
       n_g1 = xbe_lib[i_xbel].n_g;
       for (i_g=0; i_g < n_g1; i_g++) {
         n_gvar1 = xbe_lib[i_xbel].n_gvar[i_g];
         for (int i_gvar=0; i_gvar < n_gvar1; i_gvar++) {
           k = smat.map_gvar_to_xbe[i_xbeu][i_g][i_gvar];
           smat.m_x.val[k] = 0.0;
         }
       }
     }
   }
   return;
} // end of xbe_init_jac_startup_x
// -----------------------------------------------------------------------------
void xbe_init_jac_trns_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SysMat &smat) {

   int i_xbeu,i_xbel,i_f,row0,k,i_g,n_f1,n_g1,n_gvar1;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;
       for (i_f=0; i_f < n_f1; i_f++) {
         row0 = smat.xbe_f_to_row[i_xbeu][i_f];
         k = smat.xbe_rhs_ddt_pntr[row0];
         smat.m_x.val[k] = 0.0;
       }
     }
     n_g1 = xbe_lib[i_xbel].n_g;
     for (i_g=0; i_g < n_g1; i_g++) {
       n_gvar1 = xbe_lib[i_xbel].n_gvar[i_g];
       for (int i_gvar=0; i_gvar < n_gvar1; i_gvar++) {
         k = smat.map_gvar_to_xbe[i_xbeu][i_g][i_gvar];
         smat.m_x.val[k] = 0.0;
       }
     }
   }
   return;
} // end of xbe_init_jac_trns_x
// -----------------------------------------------------------------------------
void check_convergence_x(
   SysMat &smat,
   SolveBlocks &slv) {

   int precision_real = 4;
   int width = precision_real + 7;
   bool flag_1;

   cout << scientific; cout << setprecision(precision_real);

// Compute norms if check/write is specified:

   if (slv.x_nr_flag_check_rhs2 || slv.x_nr_flag_write_rhs2) {
     slv.x_nr_norm_rhs2 = norm_2(smat.n_solvec_x,smat.rhs_m_x);
   }
   if (slv.x_nr_flag_write_rhsinf) {
     slv.x_nr_norm_rhsinf = norm_inf(smat.n_solvec_x,smat.rhs_m_x);
   }

// write norms to console

   if (slv.x_nr_flag_write_rhs2) {
     cout << slv.iter_newton << " " << "x_nr_norm_rhs2 =" << setw(width)
       << slv.x_nr_norm_rhs2 << endl;
   }
   if (slv.x_nr_flag_write_rhsinf) {
     cout << slv.iter_newton << " " << "x_nr_norm_rhsinf =" << setw(width)
       << slv.x_nr_norm_rhsinf << endl;
   }

// Check convergence

   slv.flag_nr_converged = true;
   slv.flag_nr_norm_large = false;

   if (slv.x_nr_flag_check_rhs2) {
     if (slv.iter_newton > 0) {
       flag_1 = flag_nan(slv.x_nr_norm_rhs2) ||
         (slv.x_nr_norm_rhs2 > slv.nr_norm_large);
       if (flag_1) {
         cout << "check_convergence_x: x_nr_norm_rhs2 = " << setw(width)
           << slv.x_nr_norm_rhs2 << " is too large." << endl;
         slv.flag_nr_converged = false;
         slv.flag_nr_norm_large = true;
       }
     }
     if (slv.x_nr_norm_rhs2 > slv.x_nr_eps_rhs) {
       slv.flag_nr_converged = false;
     }
   }

   return;
} // end of check_convergence_x
// -----------------------------------------------------------------------------
void check_convergence_count_x(
   SolveBlocks &slv) {

// for now, we will only allow check_rhs2, but keep the following anyway
// for later use.

   int count;

   count = 0;
   if (slv.x_nr_flag_check_rhs2) count++;

   if (count == 0) {
     cout << "check_convergence_count_x: no convergence criterion?" << endl;
     cout << "  Halting..." << endl; exit(1);
   }
   return;
} // end of check_convergence_count_x
// -----------------------------------------------------------------------------
void trzbdf2_1_x(
   vector<XbeUsr> &xbe_usr,
   SysMat &smat,
   SolveBlocks &slv) {

   double norm2,norm,r,hstar,step1,step2,x,c1,c2,c3;
   double g_0,g_1,g_2;
   int i_xbeu;
   int i_f,i_rhs;

   norm2 = 0.0;
   c1 = 1.0/slv.bank_gamma;
   c3 = 1.0/(1.0-slv.bank_gamma);
   c2 = c1*c3;

   for (i_rhs=0; i_rhs < smat.m_x.n_row; i_rhs++) {
     if (smat.xbe_rhs_ddt_flag[i_rhs]) {
       i_xbeu = smat.xbe_rhs_ddt_i_xbeu[i_rhs];
       i_f = smat.xbe_rhs_ddt_i_f[i_rhs];

       g_0 = xbe_usr[i_xbeu].g      [i_f];
       g_1 = xbe_usr[i_xbeu].g_old_1[i_f];
       g_2 = xbe_usr[i_xbeu].g_old_2[i_f];

       x = c1*g_2 - c2*g_1 + c3*g_0;

       norm2 = norm2 + x*x;
       if (flag_nan(norm2)) {
         cout << "trzbdf2_1_x: norm2 is NAN. Halting..." << endl; exit(1);
       }
     }
   }

   norm = 2.0*fabs(slv.bank_c)*slv.delt_x*sqrt(norm2);

   r = norm/slv.bank_tolr;

   hstar = slv.delt_x*pow(r,-0.333333333);

   if (r >= 2.0) {
     slv.flag_accept_sol = false;
     slv.delt_new_x = 0.9*hstar;
   } else {
     slv.flag_accept_sol = true;
     step1 = slv.bank_theta1*hstar;
     step2 = 2.0*slv.delt_x;
     slv.delt_new_x = min(step1,step2);
   }

   return;
} // end of trzbdf2_1_x
// -----------------------------------------------------------------------------
void x_assign_nextbreak_1(
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global) {

   for (int i=0; i < cct.n_xbeu; i++) {
     xbe_usr[i].next_break = global.time_end;
   }
   return;
} // end of x_assign_nextbreak_1
// -----------------------------------------------------------------------------
