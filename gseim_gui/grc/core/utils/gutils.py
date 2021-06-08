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
import fcntl

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

def get_parms(filename, d_parms):

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

def find_nth(s_in, c, n):
#   https://stackoverflow.com/questions/30833409/python-deleting-the-first-2-lines-of-a-string

    start = s_in.find(c)
    while start >= 0 and n > 1:
        start = s_in.find(c, start+len(c))
        n -= 1
    return start

def is_running(pid):
    if os.path.isdir('/proc/{}'.format(pid)):
        return True
    return False

def get_sec(time_str):
#   https://stackoverflow.com/questions/6402812/how-to-convert-an-hmmss-time-string-to-seconds-in-python/6402934
    """Get Seconds from time."""
    h, m, s = time_str.split(':')
    return int(h)*3600 + int(m)*60 + int(s)

def non_block_read(output):
    ''' even in a thread, a normal read with block until the buffer is full '''
    fd = output.fileno()
    fl = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)

    op = output.read()

    if op == None:
        return ''
    return op.decode('utf-8')

