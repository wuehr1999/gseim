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

#ifndef ROUTINES1_H
#define ROUTINES1_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <math.h>
#include <iomanip>

using namespace std;

#include "global.h"
#include "xbelib.h"
#include "xbeusr.h"
#include "xbejac.h"
#include "solveblocks.h"
#include "circuit.h"
#include "sysmat.h"
#include "matgnrl1.h"
#include "utils.h"
#include "get_yyy.h"

void clear_solvec_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SysMat &smat);

void clear_solvec_x_1(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct);

void cct_to_xbe_1(
   const int i_xbeu,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct);

void cct_to_xbe_all(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct);

void xbe_to_cct_op_1(
   const int i_xbeu,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct);

void xbe_to_cct_1(
   const int i_xbeu,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct);

void xbe_to_cct_all(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct);

void get_xbe_1(
   const int i_xbeu,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);

void one_time_parms_x(
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);

void save_history_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);

void init_sol_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   SysMat &smat,
   Global &global);

void xbe_init_guess(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);

void read_solution(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SolveBlocks &slv,
   Circuit &cct);

void write_iter_x(
   SolveBlocks &slv);

void write_solution(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SolveBlocks &slv,
   Circuit &cct);

void write_solution_1_x(
   std::string filename,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SolveBlocks &slv,
   Circuit &cct);

void write_trns(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void write_dc_startup(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void solve_trns_x_exp(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global);

void solve_trns_x_feuler(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global);

void solve_trns_x_rk4(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global);

void solve_trns_x_rkf45(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global);

void solve_trns_x_bs23(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global);

void solve_trns_x_meuler(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global);

void solve_trns_x_heun(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   SysMat &smat,
   Global &global);

void xbe_ov_prm(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   CctFile &cct_file,
   Global &global);

void assign_ov_prm(
   const int i_ov,
   const int i_file,
   const int i_var,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void xbe_find_nextbreak(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);

void get_tnext_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   Global &global);

void xbe_modulo(
   Circuit &cct,
   vector<XbeUsr> &xbe_usr);

void xbe_modulo_implicit(
   Circuit &cct,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SysMat &smat);

void xbe_reset_1(
   const int flag_implicit,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SysMat &smat,
   Global &global);

void xbe_feuler(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void xbe_feuler_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global);

void xbe_rk4(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void xbe_rk4_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global);

void xbe_rk4_exs_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global);

void xbe_rkf45(
   const int flag,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void xbe_rkf45_al(
   const int flag,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global);

void xbe_bs23(
   const int flag,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void xbe_bs23_al(
   const int flag,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global);

void xbe_meuler(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void xbe_meuler_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global);

void xbe_meuler_exs_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global);

void xbe_heun(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void xbe_heun_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   SysMat &smat,
   Global &global);

void xbe_evaluate(
   const int flag,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);

void xbe_evaluate_1(
   const int flag,
   const int i_xbeu,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);

void xbe_evaluate_2(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);

void xbe_evaluate_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);

void solve_trns_x_common(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   CctFile &cct_file,
   Global &global);

void update_rk4(
   const int stage,
   const double h1,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global);

void update_rk4_al(
   const int stage,
   const double h1,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global);

void update_rkf45(
   const int stage,
   const double h,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void update_rkf45_al(
   const int stage,
   const double h,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void update_bs23(
   const int stage,
   const double h,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void update_bs23_al(
   const int stage,
   const double h,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void update_meuler(
   const int stage,
   const double h1,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global);

void update_meuler_al(
   const int stage,
   const double h1,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global);

void update_heun(
   const int stage,
   const double h,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void update_heun_al(
   const int stage,
   const double h,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void xbeu_copy_1(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct);

void xbeu_copy_2(
   const int flag,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct);

void copy_func_to_old_x(
   const int flag_1,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global);

void solve_dc(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void solve_startup(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void solve_startup_x_exp(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void solve_startup_x_imp(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void solve_startup_linear_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void solve_startup_nonlinear_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void solve_trns(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void solve_trns_x_imp(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void solve_trns_x_be(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void solve_trns_x_trz(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void solve_trns_x_be_auto(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void solve_trns_x_trz_auto(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void solve_trns_x_trbdf2(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   CctFile &cct_file,
   Global &global);

void solve_trns_linear_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void solve_trns_nonlinear_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void solve_trns_linear_x_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void solve_trns_nonlinear_x_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void exs_x_feuler(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   SysMat &smat,
   Global &global);

void exs_x_rk4(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   SysMat &smat,
   Global &global);

void exs_x_meuler(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   SysMat &smat,
   Global &global);

void exs_x_heun(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   SysMat &smat,
   Global &global);

void exs_x_be(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void exs_x_trz(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   SolveBlocks &slv,
   Global &global);

void solve_jac_1_x(
   SysMat &smat,
   SolveBlocks &slv,
   Global &global);

void solve_jac_2_x(
   SysMat &smat,
   SolveBlocks &slv,
   Global &global);

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
   Global &global);

void form_jac_rhs_startup_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global);

void form_jac_rhs_trns_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global);

void form_jac_rhs_trns_x_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global);

void find_functions_trns_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global);

void form_solvec_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SysMat &smat,
   Circuit &cct);

void form_map_xbeuvr_1(
   SysMat &smat,
   Circuit &cct);

void dcmp_solvec_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SysMat &smat,
   Circuit &cct);

void mat_startup_2_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global);

void mat_startup_3_x(
   const int i_xbeu,
   const int i_xbel,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Global &global);

void mat_trns_2_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global);

void mat_trns_2a_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global);

void mat_trns_2_x_al(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Circuit &cct,
   Global &global);

void mat_trns_3_x(
   const int i_xbeu,
   const int i_xbel,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SysMat &smat,
   Global &global);

void mat_trns_3a_x(
   const int i_xbeu,
   const int i_xbel,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SysMat &smat);

void add_trns_terms_x(
   vector<XbeUsr> &xbe_usr,
   SysMat &smat,
   SolveBlocks &slv);

void xbe_init_jac_startup_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SysMat &smat);

void xbe_init_jac_trns_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SysMat &smat);

void check_convergence_x(
   SysMat &smat,
   SolveBlocks &slv);

void check_convergence_count_x(
   SolveBlocks &slv);

void trzbdf2_1_x(
   vector<XbeUsr> &xbe_usr,
   SysMat &smat,
   SolveBlocks &slv);

void x_assign_nextbreak_1(
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global);

#endif
