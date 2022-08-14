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

#ifndef ROUTINES2_H
#define ROUTINES2_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <math.h>
#include <iomanip>

using namespace std;

#include "global.h"
#include "ebelib.h"
#include "ebeusr.h"
#include "ebejac.h"
#include "xbelib.h"
#include "xbeusr.h"
#include "xbejac.h"
#include "solveblocks.h"
#include "circuit.h"
#include "sysmat.h"
#include "matgnrl1.h"
#include "utils.h"
#include "get_yyy.h"

void clear_solvec_e(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct,
   SysMat &smat);
void clear_solvec_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   SysMat &smat);
void clear_solvec_ex(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
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
void cct_to_ebe_nd_1(
   const int i_ebeu,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct);
void cct_to_ebe_nd_all(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct);
void cct_to_ebe_xvr_all(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
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
void assign_xvr_ebe_in(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct);
void get_xbe_1(
   const int i_xbeu,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);
void form_solvec_e(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   SysMat &smat,
   SolveBlocks &slv,
   Circuit &cct);
void form_solvec_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SysMat &smat,
   Circuit &cct);
void form_map_xbeuvr_1(
   SysMat &smat,
   Circuit &cct);
void form_solvec_ex(
   vector<EbeLib> &ebe_lib,
   vector<XbeLib> &xbe_lib,
   vector<EbeUsr> &ebe_usr,
   vector<XbeUsr> &xbe_usr,
   SysMat &smat,
   SolveBlocks &slv,
   Circuit &cct);
void dcmp_solvec_e(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   SysMat &smat,
   SolveBlocks &slv,
   Circuit &cct);
void dcmp_solvec_ssw_e(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   SysMat &smat,
   SolveBlocks &slv,
   Circuit &cct);
void dcmp_solvec_ssw_ex(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   SysMat &smat,
   SolveBlocks &slv,
   Circuit &cct);
void dcmp_solvec_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SysMat &smat,
   Circuit &cct);
void dcmp_solvec_ex(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   SysMat &smat,
   SolveBlocks &slv,
   Circuit &cct);
void one_time_parms_x(
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);
void one_time_parms_e(
   vector<EbeUsr> &ebe_usr,
   vector<EbeJac> &ebe_jac,
   Circuit &cct,
   Global &global);
void save_history_e(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   vector<EbeJac> &ebe_jac,
   Circuit &cct,
   Global &global);
void save_history_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);
void init_sol_e(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   vector<EbeJac> &ebe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   SysMat &smat,
   Global &global);
void init_sol_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   SolveBlocks &slv,
   Circuit &cct,
   SysMat &smat,
   Global &global);
void init_sol_ex(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   vector<EbeJac> &ebe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   SysMat &smat,
   Global &global);
void ebe_init_guess(
   vector<EbeUsr> &ebe_usr,
   vector<EbeJac> &ebe_jac,
   Circuit &cct,
   Global &global);
void xbe_init_guess(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   Global &global);
void read_solution(
   vector<XbeLib> &xbe_lib,
   vector<EbeLib> &ebe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<EbeUsr> &ebe_usr,
   SolveBlocks &slv,
   Circuit &cct);
void write_iter_x(
   SolveBlocks &slv);
void write_iter_e(
   SolveBlocks &slv);
void write_solution(
   vector<XbeLib> &xbe_lib,
   vector<EbeLib> &ebe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<EbeUsr> &ebe_usr,
   SolveBlocks &slv,
   Circuit &cct);
void write_solution_1_e(
   std::string filename,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   SolveBlocks &slv,
   Circuit &cct);
void write_solution_1_x(
   std::string filename,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   SolveBlocks &slv,
   Circuit &cct);
void xbe_ov_prm(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
   CctFile &cct_file,
   Global &global);
void ebe_ov_prm(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   vector<EbeJac> &ebe_jac,
   Circuit &cct,
   CctFile &cct_file,
   Global &global);
void assign_ov_prm(
   const int i_ov,
   const int i_file,
   const int i_var,
   vector<XbeUsr> &xbe_usr,
   vector<EbeUsr> &ebe_usr,
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
void ebe_find_nextbreak(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   vector<EbeJac> &ebe_jac,
   Circuit &cct,
   Global &global);
void get_tnext_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   Global &global);
void get_tnext_e(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   vector<EbeJac> &ebe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   Global &global);
void get_tnext_ex(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   vector<EbeJac> &ebe_jac,
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
void xbe_reset_1_exc(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   vector<EbeJac> &ebe_jac,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   SolveBlocks &slv,
   Circuit &cct,
   SysMat &smat,
   Global &global);
void xbe_time_parms(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   vector<XbeJac> &xbe_jac,
   Circuit &cct,
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
void ebeu_copy_stv_1(
   const int &flag_1,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct,
   Global &global);
void copy_func_to_old_e(
   const int flag_1,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct,
   Global &global);
void copy_cur_nd_nr_1(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct);
void copy_func_to_old_x(
   const int flag_1,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global);
void copy_func_to_old_ex(
   const int flag_1,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global);
void ebe_form_arrays_1(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct);
void ebe_form_arrays_ssw_e(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   SysMat &smat,
   Circuit &cct);
void ebe_form_arrays_ssw_ex(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   SysMat &smat,
   Circuit &cct);
void ebe_update_stv(
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   vector<EbeJac> &ebe_jac,
   Circuit &cct,
   Global &global);
void e_assign_nextbreak_1(
   vector<EbeUsr> &ebe_usr,
   Circuit &cct,
   Global &global);
void x_assign_nextbreak_1(
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global);
void ex_assign_nextbreak_1(
   vector<EbeUsr> &ebe_usr,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global);
void check_convergence_e(
   SysMat &smat,
   SolveBlocks &slv,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct);
void check_convergence_x(
   SysMat &smat,
   SolveBlocks &slv);
void check_convergence_ex(
   SysMat &smat,
   SolveBlocks &slv,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct);
void compute_norm_spice_e(
   SysMat &smat,
   SolveBlocks &slv,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct);
void compute_norm_spice_ex(
   SysMat &smat,
   SolveBlocks &slv,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct);
void check_spice_e(
   SysMat &smat,
   SolveBlocks &slv,
   vector<EbeLib> &ebe_lib,
   vector<EbeUsr> &ebe_usr,
   Circuit &cct,
   bool &flag_converged,
   bool &flag_norm_large);
void check_delx_all_e(
   SysMat &smat,
   SolveBlocks &slv);
void check_delx_all_x(
   SysMat &smat,
   SolveBlocks &slv);
void check_delx_all_ex(
   SysMat &smat,
   SolveBlocks &slv);
void check_convergence_count_e(
   SolveBlocks &slv);
void check_convergence_count_x(
   SolveBlocks &slv);
void check_convergence_count_ex(
   SolveBlocks &slv);

#endif