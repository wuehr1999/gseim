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

#ifndef SOLVE_BLOCKS_H
#define SOLVE_BLOCKS_H

#include <iostream>
#include <cstring>
#include <fstream>
#include <math.h>

#include "global.h"
#include "cctfile.h"
#include "circuit.h"
#include "utils.h"
#include "xbelib.h"
#include "xbeusr.h"
#include "solveparm.h"

using namespace std; 

class SolveBlocks {

public:
   vector<SolveParm> parms;

   bool flag_startup,flag_trns,flag_dc;

   int outf_real_precision,outf_real_word_width;
   int outf_time_precision,outf_time_word_width;
   int outf_sol_precision,outf_sol_word_width_real;

   int index_file_solution;

   int index_solve;
   int n_outfile,solve_type;
   int solve_type_previous;

   vector<fstream> f_output;
   vector<int> total_lines;

   bool flag_read_solution;
   bool flag_init_solution;
   bool flag_prev_solution;
   bool flag_prev_solution_exists;
   bool flag_write_solution;

   bool flag_write_time_e;
   bool flag_write_time_x;

   int write_iter_n_e,write_iter_n1_e;
   int write_iter_n_x,write_iter_n1_x;

   std::string infile_sol;

   vector<int> limit_lines;
   vector<int> out_nvar;;
   vector<std::string> outf_name;

   vector<bool> flag_solution;

   int flag_algo_x;

   bool x_algo_startup_exp;

   bool x_algo_feuler;
   bool x_algo_rk4;
   bool x_algo_rkf45;
   bool x_algo_bs23;
   bool x_algo_meuler;
   bool x_algo_heun;
   bool x_algo_be;
   bool x_algo_be_auto;
   bool x_algo_trz;
   bool x_algo_trz_auto;
   bool x_algo_trbdf2;

   bool x_algo_be0,x_algo_trz0,x_algo_auto,x_algo_bdf2;

   vector<bool> flag_out_delt_fixed;

   vector<double> out_delt;
   vector<double> out_tstart;
   vector<double> out_tend;
   vector<double> out_tnext;

   vector< vector<int> > out_var;
   vector< vector<double> > outvar_temp;

   bool flag_limit_iter_trns;
   int itmax_trns,iter_trns_x,iter_trns_e;
   int itmax_stepred,itmax_trbdf2;

//*skip for gseim
   int x_trns_iter_debug;
//*end skip for gseim

   int iter_newton;

   bool flag_fixed_delt_x;
   double time_present_x,time_next_x,time_write_x;
   double delt_min_x,delt_max_x;
   double delt_x,delt0_x;

   double delt_small;
   double time_write;

   double time_startup;
   double zero_piv,gauss_epsln;
   bool flag_debug_gauss1,flag_debug_gauss2;

//*skip for gseim
   bool flag_write_matrix_x; // wrt_debug_smat in Sequel
//*skip for gseim

// NR parameters:

   bool flag_nr_converged;
   bool flag_nr_norm_large;

   double nr_norm_large;

   int x_nr_itermax0;
   bool x_nr_flag_dmp0;
   int x_nr_dmp_itermax0;
   double x_nr_dmp_k0;

   bool x_nr_flag_dmp;
   int x_nr_itermax;
   int x_nr_dmp_itermax;
   double x_nr_dmp_k;

   bool x_nr_flag_dmp_a;
   int x_nr_itermax_a;
   int x_nr_dmp_itermax_a;
   double x_nr_dmp_k_a;

   int x_nr_iter_debug;

   bool x_nr_flag_check_rhs2;
   bool x_nr_flag_write_rhs2;
   bool x_nr_flag_write_rhsinf;

   double x_nr_eps_rhs;

   double x_nr_norm_rhs2;
   double x_nr_norm_rhsinf;

// RKF45 constants:
// "a1" stands for alpha(1),
// "b10" stands for beta(1,0),
// "g1"  stands for gamma(1), etc.

   double rkf45_a1;
   double rkf45_a2;
   double rkf45_a3;
   double rkf45_a4;
   double rkf45_a5;

   double rkf45_b10;

   double rkf45_b20;
   double rkf45_b21;

   double rkf45_b30;
   double rkf45_b31;
   double rkf45_b32;

   double rkf45_b40;
   double rkf45_b41;
   double rkf45_b42;
   double rkf45_b43;

   double rkf45_b50;
   double rkf45_b51;
   double rkf45_b52;
   double rkf45_b53;
   double rkf45_b54;

   double rkf45_4_g0;
   double rkf45_4_g1;
   double rkf45_4_g2;
   double rkf45_4_g3;
   double rkf45_4_g4;

// "e" stands for error:

   double rkf45_e0;
   double rkf45_e1;
   double rkf45_e2;
   double rkf45_e3;
   double rkf45_e4;
   double rkf45_e5;

   double rkf45_tolr,rkf45_norm2;
   double rkf45_fctr_min,rkf45_fctr_max;

   int rkf45_n_accept,rkf45_n_reject;

// variables for Bogacki-Shampine method:

   double bs23_a1;
   double bs23_a2;

   double bs23_b10;

   double bs23_b21;

   double bs23_b30;
   double bs23_b31;
   double bs23_b32;

   double bs23_b40;
   double bs23_b41;
   double bs23_b42;

// 3rd order gamma's

   double bs23_3_g0;
   double bs23_3_g1;
   double bs23_3_g2;

// 2nd order gamma's

   double bs23_2_g0;
   double bs23_2_g1;
   double bs23_2_g2;
   double bs23_2_g3;

// "e" stands for error:

   double bs23_e0;
   double bs23_e1;
   double bs23_e2;
   double bs23_e3;

   double bs23_tolr,bs23_norm2;
   double bs23_fctr_min,bs23_fctr_max;

   int bs23_n_accept,bs23_n_reject;

   double heun_a1,heun_b10,heun_g0,heun_g1;

   double bank_gamma,bank_c,bank_tolr,bank_theta1;

   double beuler_1_x,trz_1_x,bdf2_1_x,bdf2_2_x,bdf2_3_x;

   bool flag_accept_sol;
   double delt_new_x;

   double factor_stepinc,factor_stepdec;

// these are used to check if a specific parameter has been specified
// in the circuit file (and generate suitable error message if it is not):
   bool flag_t_start,flag_t_end;
   bool flag_delt0_x,flag_delt_min_x,flag_delt_max_x;

public:
  SolveBlocks();

  void set_values_1(
   const std::string &filename,
   Circuit &cct,
   Global &global,
   CctFile &cct_file);

  void method_default(
   Global &global);

  void trns_constants_1();

  void trns_constants_2_x();

  void open_output_files();

  void close_output_files();

  void get_dmp();

  void assign_parm(
   const std::string s2,
   const std::string s3,
   Global &global);

};
#endif
