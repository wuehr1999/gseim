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

#ifndef SYSMAT_H
#define SYSMAT_H

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <vector>
#include "utils.h"
#include "solveblocks.h"
#include "xbeusr.h"
#include "knuthmat.h"
#include "matgnrl1.h"
#include "matop.h"

using namespace std;

class SysMat {

  public:

   int n_solvec_x;
   int n_solvec_x0;
   int n_rhs_x0;

   int offs_xvr,offs_xaux;

   int offs_xvr_in,offs_xvr_out;
   int offs_row_xvr_in;

   int n_xvr,n_xaux;

   vector<int> offs;
   vector<int> svec_flag_x;

   KnuthMat m_x,w_x;

   MatOp mo_x;

   bool flag_alloc,flag_delete;

   double *delsvec_x,*svec_x,*svec_w_x,*svec_orig_x;

   double *svec_old_1_x,*svec_old_2_x;

   double *rhs_m_x,*rhs_w_x;

   int ***map_gvar_to_xbe;

   int **xbe_f_to_row;
   int **xbe_g_to_row;

   bool *xbe_rhs_ddt_flag;
   int *xbe_rhs_ddt_varnumber;
   int *xbe_rhs_ddt_pntr;
   int *xbe_rhs_ddt_varflag;
   int *xbe_rhs_ddt_i_xbeu;
   int *xbe_rhs_ddt_i_f;

   int *map_jac_xbe_to_m;
   int *map_rhs_xbe_to_m;

   int *eqn_xbe_be_index;
   int *eqn_xbe_be_eqn;
   int *xbe_ddt_varnumber;

  public:

  SysMat();

  void allocate_1(
   Global &global,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct);

  void set_values_1(
   Global &global,
   SolveBlocks &slv,
   Circuit &cct);

  void mat_startup_1_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   CctFile &cct_file);

  void mat_trns_1_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global,
   CctFile &cct_file);

  void delete_1();

};
#endif
