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

#include "sysmat.h"

SysMat::SysMat() {
   flag_alloc = true; flag_delete = false;
   return;
}
// -----------------------------------------------------------------------------
void SysMat::allocate_1(
   Global &global,
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct) {

// allocate variables whose size will not change (so they will not need to
// be deleted)

   int n, i_xbeu,i_xbel;
   int n_f1,n_g1,n_h1;
   int n_fvar1,n_gvar1,n_hvar1;
   int i_f,i_g,i_h;

   n = 0;
   n = max(n,global.I_XVR);
   n = max(n,global.I_XAUX);
   n = max(n,global.I_EAUXS);

   offs.resize(n+1);

   map_gvar_to_xbe = new int**[cct.n_xbeu];
   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     n_g1 = xbe_lib[i_xbel].n_g;
     map_gvar_to_xbe[i_xbeu] = new int*[n_g1];

     for (i_g=0; i_g < n_g1; i_g++) {
       n_gvar1 = xbe_lib[i_xbel].n_gvar[i_g];
       map_gvar_to_xbe[i_xbeu][i_g] = new int[n_gvar1];
     }
   }

   xbe_f_to_row = new int*[cct.n_xbeu];
   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     n_f1 = xbe_lib[i_xbel].n_f;
     xbe_f_to_row[i_xbeu] = new int[n_f1];
   }

   xbe_g_to_row = new int*[cct.n_xbeu];
   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     n_g1 = xbe_lib[i_xbel].n_g;
     xbe_g_to_row[i_xbeu] = new int[n_g1];
   }

   rhs_m_x = NULL;
   rhs_w_x = NULL;

   return;
} // end of SysMat::allocate_1
// -----------------------------------------------------------------------------
void SysMat::set_values_1(
   Global &global,
   SolveBlocks &slv,
   Circuit &cct) {

// allocate/assign variable values which will depend on the solve section.
// initialise (probably not required; never mind)
   offs_xvr   = 0;
   offs_xaux  = 0;

// added to handle exs:

   n_solvec_x  = 0;

   n_rhs_x0 = 0;

   n_solvec_x0  = cct.n_xbeu_vr + cct.x_n_ttlaux;

   svec_flag_x.clear();

   n_xvr = cct.n_xbeu_vr;
   n_xaux = cct.x_n_ttlaux;

   if (cct.flag_x_matrix) {
     n_solvec_x = n_solvec_x0;
   
     if (slv.flag_trns) {
       n_rhs_x0 = n_solvec_x0;
     }
   
     offs_xvr = 0;
     offs_xaux = offs_xvr + cct.n_xbeu_vr;
   
     for (int i=offs_xvr; i < offs_xaux; i++) {
       svec_flag_x.push_back(global.I_XVR);
     }
     for (int i=offs_xaux; i < n_solvec_x; i++) {
       svec_flag_x.push_back(global.I_XAUX);
     }
   }

   offs[global.I_XVR  ] = offs_xvr;
   offs[global.I_XAUX ] = offs_xaux;

   cout << "SysMat::set_values_1:" << endl;
   cout << "  n_solvec_x = " << n_solvec_x << endl;

   if (flag_alloc) {
     svec_x  = new double[n_solvec_x ];

     svec_w_x  = new double[n_solvec_x ];

     svec_orig_x  = new double[n_solvec_x ];

     svec_old_1_x  = new double[n_solvec_x ];
     svec_old_2_x  = new double[n_solvec_x ];

     delsvec_x  = new double[n_solvec_x ];

     xbe_rhs_ddt_flag      = new bool[n_rhs_x0];
     xbe_rhs_ddt_varnumber = new int [n_rhs_x0];
     xbe_rhs_ddt_pntr      = new int [n_rhs_x0];
     xbe_rhs_ddt_varflag   = new int [n_rhs_x0];
     xbe_rhs_ddt_i_xbeu    = new int [n_rhs_x0];
     xbe_rhs_ddt_i_f       = new int [n_rhs_x0];

     flag_alloc = false; flag_delete = true;
   } else {
     cout << "SysMat::set_values_1: trying to allocate" << endl;
     cout << "  with flag_alloc = false? Halting..." << endl; exit(1);
   }

   return;
} // end of SysMat::set_values_1
// -----------------------------------------------------------------------------
void SysMat::mat_startup_1_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   CctFile &cct_file) {

   int i_xbel,i_xbeu,n_f1,n_g1,n_gvar1;
   int n_col1,n_row1;
   int var_flag,var_number;
   int n_rows,row0,col0;
   int i_f,i_g;

   n_row1 = cct.x_n_ttlg;
   n_col1 = n_solvec_x;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     n_f1 = xbe_lib[i_xbel].n_f;
     n_g1 = xbe_lib[i_xbel].n_g;

     for (i_f=0; i_f < n_f1; i_f++) {
       xbe_f_to_row[i_xbeu][i_f] = -1;
     }
     for (i_g=0; i_g < n_g1; i_g++) {
       n_gvar1 = xbe_lib[i_xbel].n_gvar[i_g];
       for (int i_gvar=0; i_gvar < n_gvar1; i_gvar++) {
         map_gvar_to_xbe[i_xbeu][i_g][i_gvar] = -1;
       }
       xbe_g_to_row[i_xbeu][i_g] = -1;
     }
   }

   m_x.allocate_1(0,n_row1,n_col1);
   knuth_zero_1(m_x);
   n_rows = 0;

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;

     if (xbe_lib[i_xbel].flag_integrate) {
       n_f1 = xbe_lib[i_xbel].n_f;
       for (i_f=0; i_f < n_f1; i_f++) {
         var_number = xbe_usr[i_xbeu].fvar[i_f];
         var_flag = xbe_lib[i_xbel].ddt_varflag[i_f];

         col0 = offs[var_flag] + var_number;
         row0 = n_rows;
         xbe_f_to_row[i_xbeu][i_f] = row0;
         knuth_addentry(m_x,row0,col0,1.0);
         n_rows++;
       }
     } else {
       n_g1 = xbe_lib[i_xbel].n_g;
       for (i_g=0; i_g < n_g1; i_g++) {
         n_gvar1 = xbe_lib[i_xbel].n_gvar[i_g];
         for (int i_gvar=0; i_gvar < n_gvar1; i_gvar++) {
           var_number = xbe_usr[i_xbeu].gvar[i_g][i_gvar];
           var_flag = xbe_lib[i_xbel].gvar_flag[i_g][i_gvar];

           col0 = offs[var_flag] + var_number;
           row0 = n_rows;
           xbe_g_to_row[i_xbeu][i_g] = row0;
           map_gvar_to_xbe[i_xbeu][i_g][i_gvar] = m_x.n_nz;
           knuth_addentry(m_x,row0,col0,0.0);
         }
         n_rows++;
       }
     }
   }

   check_knuth(m_x,(char*)"mat_startup_1_x: matrix error");

   map_rhs_xbe_to_m = new int[m_x.n_row];
   map_jac_xbe_to_m = new int[m_x.n_nz];

   rhs_m_x = new double[m_x.n_row];
   rhs_w_x = new double[m_x.n_row];

   return;
} //end of SysMat::mat_startup_1_x
// -----------------------------------------------------------------------------
void SysMat::mat_trns_1_x(
   vector<XbeLib> &xbe_lib,
   vector<XbeUsr> &xbe_usr,
   Circuit &cct,
   Global &global,
   CctFile &cct_file) {

   int i_xbel,i_xbeu,n_f1,n_g1,n_gvar1;
   int n_col1,n_row1;
   int var_flag,var_number;
   int n_rows,row0,col0;
   int *flag_col,*nnz_dummy;
   int flag_col_1;
   int i_f,i_g;

   n_row1 = cct.x_n_ttlg;

   n_col1 = n_solvec_x;

   cout << "SysMat::mat_trns_1_x: n_row1 = " << n_row1
     << " n_col1 = " << n_col1 << endl;

   flag_col  = new int[n_col1];
   nnz_dummy = new int[n_col1];

   for (int i=0; i < n_col1; i++) {
     flag_col [i] = -1;
     nnz_dummy[i] = -1;
   }

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     n_f1 = xbe_lib[i_xbel].n_f;
     for (i_f=0; i_f < n_f1; i_f++) {
       xbe_f_to_row[i_xbeu][i_f] = -1;
     }
     n_g1 = xbe_lib[i_xbel].n_g;
     for (i_g=0; i_g < n_g1; i_g++) {
       xbe_g_to_row[i_xbeu][i_g] = -1;

       n_gvar1 = xbe_lib[i_xbel].n_gvar[i_g];
       for (int i_gvar=0; i_gvar < n_gvar1; i_gvar++) {
         map_gvar_to_xbe[i_xbeu][i_g][i_gvar] = -1;
       }
     }
   }

   m_x.allocate_1(0,n_row1,n_col1);
   knuth_zero_1(m_x);
   n_rows = 0;

   assign_array_1<bool>(xbe_rhs_ddt_flag,     n_rhs_x0,false);
   assign_array_1<int> (xbe_rhs_ddt_varnumber,n_rhs_x0,-1);
   assign_array_1<int> (xbe_rhs_ddt_varflag,  n_rhs_x0,-1);
   assign_array_1<int> (xbe_rhs_ddt_pntr,     n_rhs_x0,-1);
   assign_array_1<int> (xbe_rhs_ddt_i_xbeu,   n_rhs_x0,-1);
   assign_array_1<int> (xbe_rhs_ddt_i_f,      n_rhs_x0,-1);

   for (i_xbeu=0; i_xbeu < cct.n_xbeu; i_xbeu++) {
     i_xbel = xbe_usr[i_xbeu].index_xbel;
     n_g1 = xbe_lib[i_xbel].n_g;

     if (xbe_lib[i_xbel].flag_integrate) {
       for (i_f=0; i_f < n_g1; i_f++) {
         var_number = xbe_usr[i_xbeu].fvar[i_f];
         var_flag = xbe_lib[i_xbel].ddt_varflag[i_f];

         col0 = offs[var_flag] + var_number;
         row0 = n_rows;

         flag_col_1 = i_f + col0;
         flag_col[col0] = flag_col_1;

         xbe_rhs_ddt_flag[row0] = true;
         xbe_rhs_ddt_i_xbeu[row0] = i_xbeu;
         xbe_rhs_ddt_i_f[row0] = i_f;
         xbe_rhs_ddt_varflag[row0] = var_flag;
         xbe_rhs_ddt_varnumber[row0] = col0;

         xbe_f_to_row[i_xbeu][i_f] = row0;
         nnz_dummy[col0] = m_x.n_nz;
         xbe_rhs_ddt_pntr[row0] = m_x.n_nz;

         knuth_addentry(m_x,row0,col0,0.0);

//       the non_ddt terms will come from g

         n_gvar1 = xbe_lib[i_xbel].n_gvar[i_f];
         for (int i_gvar=0; i_gvar < n_gvar1; i_gvar++) {
           var_number = xbe_usr[i_xbeu].gvar[i_f][i_gvar];
           var_flag = xbe_lib[i_xbel].gvar_flag[i_f][i_gvar];
           col0 = offs[var_flag] + var_number;

           flag_col_1 = i_f + col0;

           if (flag_col[col0] != flag_col_1) {

             map_gvar_to_xbe[i_xbeu][i_f][i_gvar] = m_x.n_nz;
             flag_col[col0] = flag_col_1;
             nnz_dummy[col0] = m_x.n_nz;

             knuth_addentry(m_x,row0,col0,0.0);

           } else {
             map_gvar_to_xbe[i_xbeu][i_f][i_gvar] = nnz_dummy[col0];
           }
         }
         n_rows++;
       }
     } else {
       for (i_g=0; i_g < n_g1; i_g++) {
         n_gvar1 = xbe_lib[i_xbel].n_gvar[i_g];
         for (int i_gvar=0; i_gvar < n_gvar1; i_gvar++) {
           var_number = xbe_usr[i_xbeu].gvar[i_g][i_gvar];
           var_flag = xbe_lib[i_xbel].gvar_flag[i_g][i_gvar];
           col0 = offs[var_flag] + var_number;
           row0 = n_rows;
           xbe_g_to_row[i_xbeu][i_g] = row0;

           map_gvar_to_xbe[i_xbeu][i_g][i_gvar] = m_x.n_nz;
           knuth_addentry(m_x,row0,col0,0.0);
         }
         n_rows++;
       }
     }
   }
   check_knuth(m_x,(char*)"mat_trns_1_x: matrix error");

   map_rhs_xbe_to_m = new int[m_x.n_row];
   map_jac_xbe_to_m = new int[m_x.n_nz];

   rhs_m_x = new double[m_x.n_row];
   rhs_w_x = new double[m_x.n_row];

   delete[] flag_col;
   delete[] nnz_dummy;

   return;
} //end of SysMat::mat_trns_1_x
// -----------------------------------------------------------------------------
void SysMat::delete_1() {

   if (rhs_m_x != NULL) {
     delete[] rhs_m_x; rhs_m_x = NULL;
   }
   if (rhs_w_x != NULL) {
     delete[] rhs_w_x; rhs_w_x = NULL;
   }

   if (flag_delete) {
     delete[] svec_x;

     delete[] svec_w_x;

     delete[] svec_orig_x;

     delete[] svec_old_1_x;
     delete[] svec_old_2_x;

     delete[] delsvec_x;

     flag_alloc = true; flag_delete = false;
   } else {
     cout << "SysMat::delete_1: trying to delete" << endl;
     cout << "  with flag_delete = false? Halting..." << endl; exit(1);
   }

// delete arrays in Knuth matrices and matop.

   if (m_x.flag_delete) m_x.delete_1();
   if (w_x.flag_delete) w_x.delete_1();

   return;
} // end of SysMat::delete_1
