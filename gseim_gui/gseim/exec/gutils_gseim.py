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

import sys
import os

def next_line(f, pos_in, string_list):

    f.seek(pos_in)
    for i in range(500):
        pos = f.tell()
        line = f.readline()

        s2 = line.replace("\n","")
        s3 = s2.replace("="," ")
        s1 = s3.split()
        if (len(s1) == 0):
            continue
        if (s1[0] == '#'):
            continue
        if (len(string_list) == 0):
            for k in s1:
                string_list.append(k)
        else:
            if (s1[0] == '+'):
                for j in range(1,len(s1)):
                    string_list.append(s1[j])
            else:
                break
    return pos

def next_line_1(f, pos_in):

    string_list = []
    f.seek(pos_in)
    for i in range(500):
        pos = f.tell()
        line = f.readline()
        s2 = line.replace("\n","")
        s3 = s2.replace("="," ")
        s1 = s3.split()
        if (len(s1) == 0):
            continue
        if (s1[0] == '#'):
            continue
        if (len(string_list) == 0):
            for k in s1:
                string_list.append(k)
        else:
            if (s1[0] == '+'):
                for j in range(1,len(s1)):
                    string_list.append(s1[j])
            else:
                break
    return pos, string_list

def format_string_1(nmax, prefix, indent1, indent2, slist_in, slist_out):

    i_word = 0
    i_line = 0
    new_line = True
    n = len(slist_in)
    flag_done = False

    while True:
        if (flag_done):
            break
        if (new_line):
            if (i_line == 0):
                string1 = indent1 + prefix
            else:
                string1 = indent1 + indent2
            new_line = False
        for i in range(i_word, n):
            if (max(len(indent1 + prefix),
                    len(indent1 + indent2))
                + len(slist_in[i]) + 1) > nmax:
                print('format_string_1: error. Halting...')
                sys.exit()
            n1 = len(string1) + len(slist_in[i]) + 1
            if (n1 > nmax):
                new_line = True
                i_line += 1
                i_word = i
                slist_out.append(string1)
                break
            else:
                if (i == n-1):
                    string1 += slist_in[i] + ';'
                    slist_out.append(string1)
                    flag_done = True
                    break
                else:
                    string1 += slist_in[i] + ','

def format_string_1a(f, nmax, prefix, indent1, indent2, slist_in, slist_out):
    format_string_1(nmax, prefix, indent1, indent2, slist_in, slist_out)
    for i in slist_out:
        f.write(i + '\n')

def format_string_2(prefix, offset, indent1, slist_in, slist_out):
    n = len(slist_in)
    for i in range(n):
        string1 = indent1 + prefix + slist_in[i] + ' = ' + str(i + offset) + ';'
        slist_out.append(string1)

def format_string_2a(f, offset, prefix, indent1, slist_in, slist_out):
    format_string_2(prefix, offset, indent1, slist_in, slist_out)
    for i in slist_out:
        f.write(i + '\n')

def format_string_3(prefix, indent1, n, slist_out):
    for i in range(n):
        string1 = indent1 + prefix + str(i+1) + ' = ' + str(i) + ';'
        slist_out.append(string1)

def format_string_3a(f, prefix, indent1, n, slist_out):
    format_string_3(prefix, indent1, n, slist_out)
    for i in slist_out:
        f.write(i + '\n')

def format_string_4(nmax, indent1, indent2, slist_in, slist_out):

    i_word = 0
    i_line = 0
    new_line = True
    n = len(slist_in)
    flag_done = False

    while True:
        if (flag_done):
            break
        if (new_line):
            if (i_line == 0):
                string1 = indent1
            else:
                string1 = indent2
            new_line = False
        for i in range(i_word, n):
            if (max(len(indent1), len(indent2)) + len(slist_in[i]) + 1) > nmax:
                print('format_string_4: error. Halting...')
                sys.exit()
            n1 = len(string1) + len(slist_in[i]) + 1
            if (n1 > nmax):
                new_line = True
                i_line += 1
                i_word = i
                slist_out.append(string1)
                break
            else:
                if (i == n-1):
                    string1 += slist_in[i]
                    slist_out.append(string1)
                    flag_done = True
                    break
                else:
                    string1 += slist_in[i] + ' '

def format_string_4a(f, nmax, indent1, indent2, slist_in):
    slist_out = []
    format_string_4(nmax, indent1, indent2, slist_in, slist_out)
    for i in slist_out:
        f.write(i.rstrip() + '\n')

def extract_strings_1(filename, keyword, string_list):

    f = open(os.path.expanduser(filename), 'r')
    pos_in = f.tell()
    flag_found = False

    nmaxlines = 10000

    for i in range(nmaxlines):
        s1 = []
        pos_in = next_line(f, pos_in, s1)
        if (s1[0] == keyword):
            flag_found = True
            for j in range(1,len(s1),2):
                string_list.append(s1[j])
            break

    if (not flag_found):
        print('extract_strings_1: ',keyword,' not found')
        sys.exit()

    f.close()

def extract_strings_2(filename, keyword, string_list):

    f = open(os.path.expanduser(filename), 'r')
    pos_in = f.tell()
    flag_found = False

    nmaxlines = 10000

    for i in range(nmaxlines):
        s1 = []
        pos_in = next_line(f, pos_in, s1)
        if (s1[0] == keyword):
            flag_found = True
            for j in range(1,len(s1),1):
                string_list.append(s1[j])
            break

    if (not flag_found):
        print('extract_strings_2: ',keyword,' not found')
        sys.exit()

    f.close()

def extract_dict_1(filename, keyword, d1):

    f = open(os.path.expanduser(filename), 'r')
    pos_in = f.tell()
    flag_found = False

    nmaxlines = 10000

    for i in range(nmaxlines):
        s1 = []
        pos_in = next_line(f, pos_in, s1)
        if (s1[0] == keyword):
            flag_found = True
            for j in range(1,len(s1),2):
                d1[s1[j]] = s1[j+1]
            break

    if (not flag_found):
        print('extract_dict_1: ',keyword,' not found')
        sys.exit()

    f.close()

def write_elements_2(el_list,auxdir,filename,path,ext):

    outfile = auxdir + '/' + filename
    f_el = open(os.path.expanduser(outfile), 'w')
    
    for i in el_list:
        f_el.write(path + i + ext + '\n')
    f_el.close()

def extract_lines_1(filename, keyword1, keyword2, string_list):

    flag_keyword1 = False
    flag_keyword2 = False
    f = open(os.path.expanduser(filename), 'r')

    for i in range(10000):
        line = f.readline()
        if (line == ""):
            print('extract_lines_1: unexpected eof reached. Halting...')
            sys.exit()
        s2 = line.replace("\n","")
        s1 = s2.split()
        if (len(s1) > 0):
            if (s1[0] == keyword2):
                flag_keyword2 = True
                break
        if (flag_keyword1):
            string_list.append(line)
            continue
        if (len(s1) > 0):
            if (s1[0] == keyword1):
                flag_keyword1 = True

    f.close()

    if (not flag_keyword1):
        print('extract_lines_1: ',keyword1,' not found. Halting...')
        sys.exit()
    if (not flag_keyword2):
        print('extract_lines_1: ',keyword2,' not found. Halting...')
        sys.exit()

def extract_lines_1a(f_out, filename, keyword1, keyword2):
    s1 = []
    extract_lines_1(filename, keyword1, keyword2, s1)
    for i in s1:
        f_out.write(i)

def extract_int_1(filename, keyword):

    f = open(os.path.expanduser(filename), 'r')
    pos_in = f.tell()
    flag_found = False

    for i in range(10000):
        s1 = []
        pos_in = next_line(f, pos_in, s1)
        if (s1[0] == keyword):
            flag_found = True
            break

    if (not flag_found):
        print('extract_int_1: ',keyword,' not found')
        sys.exit()
    else:
        f.close()
        return int(s1[1])

def get_parms_1(filename, d_parms):

    f = open(os.path.expanduser(filename), 'r')
    pos_in = f.tell()
    pos_in, l_strings = next_line_1(f, pos_in)
    if l_strings[0] != 'begin_file':
        print('get_parms: expected begin_file. Halting...')
        sys.exit()
    while True:
        pos_in, l_strings = next_line_1(f, pos_in)
        if l_strings[0] == 'end_file':
            f.close()
            break
        if l_strings[0] != 'begin_parm':
            print('get_parms: expected begin_parm. Halting...')
            sys.exit()
        flag_force_write = True if 'force_write' in l_strings else False

        pos_in, l_strings = next_line_1(f, pos_in)
        if l_strings[0] != 'keyword:':
            print('get_parms: expected keyword. Halting...')
            sys.exit()
        parm_name = l_strings[1]

        pos_in, l_strings = next_line_1(f, pos_in)
        if l_strings[0] != 'options:':
            print('get_parms: expected options. Halting...')
            sys.exit()
        l_options = l_strings[1:]

        pos_in, l_strings = next_line_1(f, pos_in)
        if l_strings[0] != 'default:':
            print('get_parms: expected default. Halting...')
            sys.exit()
        parm_default = l_strings[1]

        pos_in, l_strings = next_line_1(f, pos_in)
        if l_strings[0] != 'end_parm':
            print('get_parms: expected end_parm. Halting...')
            sys.exit()

        d_parms[parm_name] = {'options': l_options, 'default': parm_default}
        d_parms[parm_name]['flag_force_write'] = flag_force_write

