A library for calculating static stability metrics of walking systems 
written by Simon Wilshin, swilshin@gmail.com, Mar 2016.

This work was completed with support from the Royal Veterinary College, 
London (www.rvc.ac.uk).


Full documentation is provided by epydoc, a script for producing this 
documentation on linux systems is provided (epydoc.sh). Demo programs compute 
the foot placements of a spider and demonstrate the ability to construct 
the convex hull of the feet in contact with the ground 
(staticstab/examples/demo.py), to calculate the longitudinal 
stability[1] of walks close to trot (staticstab/examples/mfmargins.py), to 
calculate the stability difference between a quadruped using a forward and 
backward walking gait, and to produce a video of a X-RHEX[2] type robot 
walking forward and backward with convex hull.

[1] McGhee, Robert B., and Frank, Andrew A. "On the stability properties of 
quadruped creeping gaits." Mathematical Biosciences 3 (1968): 331-351.
[2] Kevin C. Galloway, G. C. Haynes, B. Deniz Ilhan, Aaron M. Johnson, Ryan 
Knopf, Goran Lynch, Benjamin Plotnick, Mackenzie White and D. E. Koditschek
X-RHex: A Highly Mobile Hexapedal Robot for Sensorimotor Tasks, University 
of Pennsylvania Technical Report (2010)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
