
GSEIM (General-purpose Simulator with Explicit and Implicit
Methods) is meant for solving a set of ordinary differential
equations (ODEs) arising from a system (circuit) schematic
diagram provided by the user. In terms of basic functionality,
it is similar to
[Matlab-Simulink](https://in.mathworks.com/products/simulink.html)
and
[Scicos](http://www.scicos.org/).
GSEIM includes a schematic entry GUI (python) adapted from the
[GNU Radio](https://www.gnuradio.org//) package, a C++ solver,
a Qt-python plotting GUI,
and a few additional python pieces.

# Installation instructions (for Linux only)

- Make sure you have ``` python3 ``` and ``` g++ ``` installed.
- Download the ```gseim_gui``` folder in your home directory.
- ```cd ~/gseim_gui/gseim/cppsrc```
- ```make -f mkgseim```
- ```make -f mkfilter_so```
- ```mv libfilter.so ~/gseim_gui/gseim/exec/.```
- ```cd ~/gseim_gui```
- ```python3 ~/gseim_gui/grc/scripts/run_gseim```
- The GSEIM GUI should show up. Change the window to full size;
  quit the GUI and run ```run_gseim``` again.
- Follow instructions in the **Getting started** page of the
  [GSEIM manual](https://gseim.github.io).

# Documentation

The GSEIM documentation can be accessed
online [(https://gseim.github.io)](https://gseim.github.io).

# License

The use and redistribution of GSEIM is governed by GPLv3.

**Contributions** are sought in the following areas:

- Schematic entry GUI
  - allow images within element symbols
  - allow rectilinear wires, with features like this
    [Simulink example](https://www.youtube.com/watch?v=iOmqgewj5XI)
  - allow text on canvas with user-specified font, symbols;
    allow moving text around in the canvas with the mouse.
  - bus facility, as in
    [this example](https://www.youtube.com/watch?v=LWnZXSuRt0I)
  - allow panning with the mouse (in addition to using scroll bars)
- Plotting GUI
  - include post-processing features such as average, rms, Fourier spectrum
  - dynamic plotting: plot should get updated as simulation proceeds

In addition, suggestions for improving the package, examples,
documentation are welcome.

**Future Plans:**

- enhance the element library
- include additional examples
- incorporate *electrical* type elements
- develop a teaching/RnD product for power electronics
- extension to other areas
