# mpLBP-for-triangulations 
MATLAB code which implements the mpLBP for triangulations

# File list:
- (folder) docs
  - citeme.bib
- mpLBP.m
- mpLBPp.m
- LICENCE.md
- README.md (this file)

# Description of the package:

This is a MATLAB implementation of the mpLBP as presented in the paper "mpLBP: A point-based representation for surface pattern description". You can find the complete reference in the docs folder. If you use this code, or its variations, in your work, PLEASE cite us.

This version works only on triangulations. The only function you need to run the code is mpLBPp.m. No additional toolboxes are required. The code is not fully optimized. I'm working on a better-structured version of the code, eventually, I will update this archive with that version. This code was tested on Matlab 2021b but it should work in previous versions as well. You can also try to use it in Octave but I never tested it. If you have the Parallel Computing Toolbox, you can use the mpLBPp.m function, reducing the running time (otherwise, the two functions are identical).

Both functions are fully commented on. If you have any problems and/or questions, please contact me at elia.moscoso@ge.imati.cnr.it. 

Copyright (c) Moscoso~Thompson Elia, 2020

# FAQ:

Q: I'm using your code on mymesh.off, and I get an error. Why?  
A: Usually, this is due to the mesh geometry. Check for duplicated points, intersecting triangles and other triangulation irregularities.


Q: Which radius should I use? multP? n_rad?  
A: A quick read of the paper mentioned above can help you in this. About the radius, try to set it such that a sphere of that radius contains the feature that is repeated in the pattern. This is not a general rule, but practically speaking it is the most efficient one. About the other parameter, the most efficient settings we found across multiple datasets are P=4, nrad=7.


Q: I'm looking at your code, the structure at line xx is odd/doesn't make sense. Are you hiding something?  
A: This is because the code I'm giving you is a simplified/cut version of the advanced code I'm currently working on. Right now, I have no time for developing a new and simpler version from scratches. Also, I'm "rushing" this release and I'm currently out of time to do something smoother.


Q: I'm using your code on this dataset you have used on your paper but my results are different!  
A: Of course, other factors must be considered, like how the descriptor is computed, which implementation of the distance measure between descriptors is used, and so on. If you need help with this, contact me. I'll answer when I can.
