"""
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
"""


# begin_coord

# 1 wdg -> 2 wdg xfmr

#delta = 3
delta = Constants.CANVAS_GRID_SIZE

#port_sep_x = 3
#port_sep_y = 1
#port_block_x = 3
#port_block_y = 12

port_sep_x = 6
port_sep_y = 4
port_block_x = 6
port_block_y = 24

delxb2 = delta*(port_block_x + port_sep_x)
delyb2 = delta*(port_block_y + port_sep_y)

delx = 2*delxb2
dely = 2*delyb2

k1x = 0.25

dy = dely/14

d2 = 2.0*delta*port_sep_x - dy
d3 = k1x*d2
d1 = d2 + dy

# vertical lines (in the centre)

x1 = delxb2 - 0.5*d3
y1 = dy

x2 = x1
y2 = dely - dy

x3 = delxb2 + 0.5*d3
y3 = y1

x4 = x3
y4 = y2

# vertical lines (at the end of coils)

# left winding

x5 = delxb2 - d1/2
y5 = 0

x6 = x5
y6 = dy*4

x7 = x5
y7 = dely - dy*4

x8 = x5
y8 = dely

# top right winding

x9 = delxb2 + d1/2
y9 = 0

x10 = x9
y10 = dy

x11 = x9
y11 = y10 + 4*dy

x12 = x9
y12 = y11 + dy

# bottom right winding

x13 = x9
y13 = y12 + 2*dy

x14 = x9
y14 = y13 + dy

x15 = x9
y15 = y14 + 4*dy

x16 = x9
y16 = y15 + dy

# horizontal contact lines

x17 = x9
y17 = y12

x18 = delx
y18 = y12

x19 = x9
y19 = y13

x20 = delx
y20 = y13

# semi circles
# A: move_to locations, B: centres

A_left = []
B_left = []
A_top_right = []
B_top_right = []
A_bottom_right = []
B_bottom_right = []

y0 = 4*dy
for i in range(6):
    A_left.append((x5, y0))
    B_left.append((x5, y0 + dy/2))
    y0 += dy

y0 = dy
for i in range(4):
    A_top_right.append((x9, y0))
    B_top_right.append((x9, y0 + dy/2))
    y0 += dy

y0 = 9*dy
for i in range(4):
    A_bottom_right.append((x9, y0))
    B_bottom_right.append((x9, y0 + dy/2))
    y0 += dy

c_ = []
c_.append((x1, y1))  # 0
c_.append((x2, y2))  # 1
c_.append((x3, y3))  # 2
c_.append((x4, y4))  # 3
c_.append((x5, y5))  # 4
c_.append((x6, y6))  # 5
c_.append((x7, y7))  # 6
c_.append((x8, y8))  # 7
c_.append((x9, y9))  # 8
c_.append((x10, y10))  # 9
c_.append((x11, y11))  # 10
c_.append((x12, y12))  # 11
c_.append((x13, y13))  # 12
c_.append((x14, y14))  # 13
c_.append((x15, y15))  # 14
c_.append((x16, y16))  # 15
c_.append((x17, y17))  # 16
c_.append((x18, y18))  # 17
c_.append((x19, y19))  # 18
c_.append((x20, y20))  # 19

for i in range(6):
    c_.append(A_left[i])
    c_.append(B_left[i])

for i in range(4):
    c_.append(A_top_right[i])
    c_.append(B_top_right[i])

for i in range(4):
    c_.append(A_bottom_right[i])
    c_.append(B_bottom_right[i])

radius = dy/2

a_ = [-90, 90, -90, -270]
t_ = [radius]

# end_coord

# begin_draw

cr.move_to(*c_[0])
cr.line_to(*c_[1])

cr.move_to(*c_[2])
cr.line_to(*c_[3])

cr.move_to(*c_[4])
cr.line_to(*c_[5])

cr.move_to(*c_[6])
cr.line_to(*c_[7])

cr.move_to(*c_[8])
cr.line_to(*c_[9])

cr.move_to(*c_[10])
cr.line_to(*c_[11])

cr.move_to(*c_[12])
cr.line_to(*c_[13])

cr.move_to(*c_[14])
cr.line_to(*c_[15])

cr.move_to(*c_[16])
cr.line_to(*c_[17])

cr.move_to(*c_[18])
cr.line_to(*c_[19])

radius = t_[0]
a1 = math.pi/180
alpha1 = a1*a_[0]
alpha2 = a1*a_[1]
alpha3 = a1*a_[2]
alpha4 = a1*a_[3]

k = 20

for i in range(6):
    cr.move_to(*c_[k])
    k += 1
    cr.arc(*c_[k], radius, alpha1, alpha2)
    k += 1

for i in range(4):
    cr.move_to(*c_[k])
    k += 1
    cr.arc_negative(*c_[k], radius, alpha3, alpha4)
    k += 1

for i in range(4):
    cr.move_to(*c_[k])
    k += 1
    cr.arc_negative(*c_[k], radius, alpha3, alpha4)
    k += 1

cr.stroke()

# end_draw
