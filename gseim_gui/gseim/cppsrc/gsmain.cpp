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

#include <iostream>
#include <cstring>
#include <fstream>
#include <sys/time.h>
#include <sys/resource.h>

#include "global.h"
#include "cctfile.h"
#include "utils.h"
#include "xbelib.h"
#include "xbeusr.h"
#include "circuit.h"
#include "solveblocks.h"
#include "routines1.h"

using namespace std;

int main(int argc, char** argv) {

   std::string file_xbe_list;
   std::string file_circuit;
   int n_xbel;

   std::string homedir;

   Global global;
   CctFile cct_file;
   Circuit cct;
   SolveBlocks slv;
   SysMat smat;

   vector<XbeLib> xbe_lib;
   vector<XbeUsr> xbe_usr;

   std::ofstream outf;

   homedir = getenv("HOME");
   file_xbe_list = homedir + "/gseim_gui/gseim/cppsrc/xbe.aux";

   if (argc == 1) {
     file_circuit = "qqtest.in";
   } else {
     file_circuit = argv[1];
   }
   cout << "main: file_circuit: " << file_circuit << endl;

   cct_file = CctFile(file_circuit);

   xbe_lib = get_lib_elements<XbeLib>(global,file_xbe_list);

   n_xbel = xbe_lib.size();

   vector<XbeJac> xbe_jac;
   xbe_jac.resize(n_xbel);
   for (int i=0; i < n_xbel; i++) {
     xbe_jac[i].initialise(xbe_lib,i);
   }

   cct = Circuit(file_circuit,xbe_lib,xbe_usr,global,cct_file);

// compute one-time parameters:

   if (cct.flag_x) one_time_parms_x(xbe_usr,xbe_jac,cct,global);

   slv.solve_type_previous = global.I_CLEAR;
   slv.flag_prev_solution_exists = false;

   cct.xbe_map_vr(xbe_lib,xbe_usr);

   cct.check_save_history(xbe_lib,xbe_usr);

   cct.check_reset_x(xbe_lib,xbe_usr);
   cct.check_modulo_x(xbe_lib,xbe_usr);

   cct.check_limit_tstep(xbe_lib,xbe_usr);

   cct.assign_flag_linear_x(xbe_lib,xbe_usr);

// fixed size allocation (not dependent on solve blocks):
   smat.allocate_1(global,xbe_lib,xbe_usr,cct);

   for (int i_solve=0; i_solve < cct_file.n_solve; i_solve++) {
     slv.index_solve = i_solve;
     assign_const_1<bool>(global.flags,false);
     cout << "main: i_solve = " << i_solve << endl;
     slv.set_values_1(file_circuit,cct,global,cct_file);

     slv.open_output_files();

//   computation of offs, n_solvec:
     smat.set_values_1(global,slv,cct);

     if (slv.flag_startup) {
       cout << "main: calling solve_startup" << endl;
       solve_startup(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,cct_file,global);
     }
     if (slv.flag_trns) {
       cout << "main: calling solve_trns" << endl;
       solve_trns(xbe_lib,xbe_usr,xbe_jac,
         smat,cct,slv,cct_file,global);
     }
     if (slv.flag_write_solution) {
       write_solution(xbe_lib,xbe_usr,slv,cct);
     }
     slv.close_output_files();

     slv.flag_prev_solution_exists = true;
     slv.solve_type_previous = slv.solve_type;

     smat.delete_1();
   }

   cout << "GSEIM: Program completed." << endl;

   return 0;
}
