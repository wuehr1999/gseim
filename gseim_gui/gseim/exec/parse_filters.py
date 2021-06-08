"""
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
"""

# Parsing of filter_1 elements.
# Read the circuit file and rewrite after parsing filter_1 elements.
#

import sys
import os
import numpy as np
import gutils_gseim as gu
from ctypes import cdll

class PICDummy(object):
    def __init__(self):
        self.obj = lib.PICDummy_new()
  
    def filter_compute_coef(self):
        lib.PICDummy_filter_compute_coef(self.obj)

if len(sys.argv) != 3:
    print('parse_filters.py: need 2 arguments. Halting...')
    sys.exit()

cct_file = sys.argv[1]
dir_exec = sys.argv[2]

# load the library
lib = cdll.LoadLibrary(dir_exec + 'libfilter.so')
  
temp_file_1 = 'cct_temp_1.in'
temp_file_2 = 'cct_temp_2.in'

cmd = 'cp ' + cct_file + ' ' + temp_file_1
os.system(cmd)
cmd = 'cp ' + cct_file + ' ' + temp_file_2
os.system(cmd)

nmax = 80
indent1 = '   '
indent2 = '+     '

flag_break = False

dummy = PICDummy()
n_filters = 0

f_1 = open(os.path.expanduser(temp_file_1), 'r')

while True:
#   Check if there is any filter_1
    line = f_1.readline()
    l = line.split()
    if l[0] == 'end_cf':
        f_1.close()
        break
    elif l[0] == 'xelement':
        if l[1].split('=')[-1].startswith('filter_1'):
            flag_break = False
            n_filters += 1
            f_1.close()
            f_1 = open(os.path.expanduser(temp_file_1), 'r')
            f_2 = open(os.path.expanduser(temp_file_2), 'w')
            while True:
                if flag_break:
                    break
                pos = f_1.tell()
                line = f_1.readline()
                l = line.split()
                if l[0] == 'xelement':
                    if l[1].split('=')[-1].startswith('filter_1'):
                        element_name = l[1].split('=')[-1]

                        f_1.seek(pos)
                        pos, l0 = gu.next_line_1(f_1, pos)

                        x_node = l0[l0.index('x') + 1]
                        y_node = l0[l0.index('y') + 1]

#                       keep values of a0,a1,.. as strings (in l2); no need to convert to float

                        l1 = ['a0', 'a1', 'a2', 'a3', 'a4', 'a5',
                              'b0', 'b1', 'b2', 'b3', 'b4', 'b5',
                              'scale_coef', 'f0']
                        d = {}

                        for k1 in l1:
                            d[k1] = l0[l0.index(k1) + 1] if k1 in l0 else '0'

                        flag_scale = False

                        if d['scale_coef'] != '0':
                            k_scale = 2.0*np.pi*float(d['f0'])
                            flag_scale = True

                        f_temp = open('temp_coef_1.dat', 'w')
                        for k in l1:
                            f_temp.write(d[k] + '  -- '+ k + '\n')
                        f_temp.close()

#                       call C++ program to compute partial fractions

                        dummy.filter_compute_coef()

#                       read partial fraction info from file

                        f_temp = open('temp_coef_2.dat', 'r')

                        line = f_temp.readline()
                        l = line.split()
                        n_roots = int(l[0])
                        gain = float(l[1])
                        a0p = float(l[2])

                        print('n_roots:', n_roots, 'gain:', gain, 'a0p:', a0p)

                        flag_real = []
                        root_power = []
                        alpha = []
                        beta = []
                        a = []
                        b = []

                        for i in range(n_roots):
                            line = f_temp.readline()
                            l = line.split()
                            flag_real.append(int(l[0]))
                            power_dummy = int(l[1])
                            root_power.append(power_dummy)

                            if flag_scale:
                                alpha_scale = k_scale
                                if power_dummy == 1:
                                    a_scale = k_scale
                                elif power_dummy == 2:
                                    a_scale = k_scale*k_scale
                            else:
                                alpha_scale = 1.0
                                a_scale = 1.0

                            line = f_temp.readline()
                            l = line.split()
                            alpha_dummy = float(l[0])
                            beta_dummy = float(l[1])
                            alpha.append(alpha_scale*alpha_dummy)
                            beta.append(alpha_scale*beta_dummy)

                            line = f_temp.readline()
                            l = line.split()
                            a_dummy = float(l[0])
                            b_dummy = float(l[1])
                            a.append(a_scale*a_dummy)
                            b.append(a_scale*b_dummy)

                        f_temp.close()

                        for i in range(n_roots):
                            print(
                              'flag_real:', flag_real[i],
                              'root_power:', root_power[i],
                              'alpha:', alpha[i],
                              'beta:', beta[i],
                              'a:', a[i],
                              'b:', b[i],
                            )

                        for i in range(n_roots):
                            a[i] *= gain
                            b[i] *= gain

                        n_elements_ttl = 0

                        for i in range(n_roots):
                            if a[i] != 0.0 or b[i] != 0.0:
                                n_elements_ttl += 1
                        if a0p != 0.0:
                            n_elements_ttl += 1

                        n_elements = 0

                        for i in range(n_roots):
                            if a[i] != 0.0 or b[i] != 0.0:
                                n_elements += 1
                                if n_elements_ttl > 1:
                                    y_nd = x_node + "_f_" + str(n_filters) + "_" + str(n_elements)
                                else:
                                    y_nd = y_node

                                if flag_real[i] == 1:
                                    l_line = ['xelement']
                                    l_line.append('type=pole_real_order_' + str(root_power[i]))
                                    l_line.append('a=' + ('%14.7e' % (a[i])).replace(" ", ""))
                                    l_line.append('alpha=' + ('%14.7e' % (alpha[i])).replace(" ", ""))
                                else:
                                    l_line = ['xelement']
                                    l_line.append('type=pole_complex_order_' + str(root_power[i]))
                                    l_line.append('a=' + ('%14.7e' % (a[i])).replace(" ", ""))
                                    l_line.append('b=' + ('%14.7e' % (b[i])).replace(" ", ""))
                                    l_line.append('alpha=' + ('%14.7e' % (alpha[i])).replace(" ", ""))
                                    l_line.append('beta=' + ('%14.7e' % (beta[i])).replace(" ", ""))

                                l_line.append('x=' + x_node)
                                l_line.append('y=' + y_nd)
                                if n_elements == 1:
                                    y_st = float(l0[l0.index("y_st") + 1] if "y_st" in l0 else '0')
                                    l_line.append('y_st=' + ('%14.7e' % (y_st)).replace(" ", ""))
                                gu.format_string_4a(f_2, nmax, indent1, indent2, l_line)

                        if a0p != 0.0:
                            n_elements += 1
                            l_line = ['xelement type=multscl']
                            y_nd = x_node + "_f_" + str(n_filters) + "_" + str(n_elements)
                            l_line.append('x=' + x_node)
                            l_line.append('y=' + y_nd)
                            l_line.append('k=' + ('%14.7e' % (gain*a0p)).replace(" ", ""))
                            gu.format_string_4a(f_2, nmax, indent1, indent2, l_line)

                        if n_elements > 1:
                            l_line = ['xelement']
                            l_line.append('type=sum_' + str(n_elements))
                            for i in range(n_elements):
                                x_nd = x_node + "_f_" + str(n_filters) + "_" + str(i+1)
                                l_line.append('x' + str(i+1) + '=' + str(x_nd))
                            l_line.append('y=' + y_node)
                            gu.format_string_4a(f_2, nmax, indent1, indent2, l_line)

                        f_1.seek(pos)
#                       now look for outvar statements with element_name
                        while True:
                            pos = f_1.tell()
                            line = f_1.readline()
                            l = line.split()
                            if l[0] == 'outvar:':
                                if l[1].split('_of_')[-1] == element_name:
                                    if l[1].split('=')[-1].split('_of_')[0] == 'x':
                                        ov_name = l[1].split('=')[0]
                                        f_2.write(indent1 + 'outvar: ' +ov_name +
                                           '=xvar_of_' + x_node + '\n')
                                    elif l[1].split('=')[-1].split('_of_')[0] == 'y':
                                        ov_name = l[1].split('=')[0]
                                        f_2.write(indent1 + 'outvar: ' +ov_name +
                                           '=xvar_of_' + y_node + '\n')
                                else:
                                    f_2.write(line)
                            else:
                                f_2.write(line)
                                if line.split()[0] == 'end_cf':
                                    f_1.close()
                                    f_2.close()
                                    flag_break = True
                                    break
                    else:
                        f_2.write(line)
                else:
                    f_2.write(line)

            f_1.close()
            f_2.close()
            cmd = 'cp ' + temp_file_2 + ' ' + temp_file_1
            os.system(cmd)
            f_1 = open(os.path.expanduser(temp_file_1), 'r')

# At this stage, temp_file_1 is the parsed version of cct_file.
# temp_file_2 can be discarded.

