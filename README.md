# NewBESO
The BESO79 and BESO94 codes implement the Bi-directional Evolutionary Structural Optimization (BESO) method for topology optimization of structures.
•	BESO79 is designed for 2D problems, making it suitable for optimizing planar structures such as beams, frames, or plates.
•	BESO94 extends this functionality to 3D problems, allowing optimization of volumetric structures like bridges, shells, or trusses.
# Methodology
The BESO method iteratively adjusts material distribution within a predefined design domain to achieve an optimized topology under given constraints and loads. Material is added or removed based on a sensitivity analysis of the structural response, aiming to minimize compliance (maximize stiffness) or achieve other design objectives.
# BESO Strategies:
-	Mesh: Structured rectangular/hexahedral elements
-	Design variable: Element density (0/1)
-	Filter: Sensitivity filter
-	Projection: No need
-	Sensitivity stabilization: Weighted average with the previous iteration
-	Optimizer: Modified OC
-	Step length control: Evolutionary ratio
-	Post-processing: No need
# Special thanks
[1] X. Huang, Y.M. Xie, Evolutionary topology optimization of continuum structures: Methods and applications, Wiley, Chichester, 2010.

[2] E. Andreassen, A. Clausen, M. Schevenels, B.S. Lazarov, O. Sigmund, Efficient topology optimization in MATLAB using 88 lines of code, Structural and Multidisciplinary Optimization, 43 (2011) 1-16.

[3] F. Ferrari, O. Sigmund, A new generation 99 line Matlab code for compliance topology optimization and its extension to 3D, Structural and Multidisciplinary Optimization, 62 (2020) 2211-2228.

[4] Z.H. Zuo, Y.M. Xie, A simple and compact Python code for complex 3D topology optimization, Advances in Engineering Software, 85 (2015) 1-11.
