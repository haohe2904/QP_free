This package contains implementations of the algorithms described in the paper:
Chunming Tang, Hao He, Wen Huang, Jinbao Jian, and Ruyun Li. "A globally and superlinearly convergent QO-free method for nonlinear optimization on Riemannian manifolds."

The repository includes five subfolders:

1. solvers – main solver routines.

2. Clients – intermediate scripts that link the problem formulation to the solver.

3. fix_rank_function – basic utilities for fixed-rank manifolds.

4. connectors – LAPACK interface.

5. supportingfunctions – additional supporting functions.

In addition to these subfolders, five demonstration scripts are provided in the root directory:
demo_fixedrank_acc.m, demo_fixedrank_rub.m, demo_fixedrank_plot.m, demo_Oblique_acc.m, and demo_Oblique_rub.m.
These illustrate the usage of all solvers implemented in this package.

Prerequisites:

1. Ensure that the MEX environment for C is properly configured in MATLAB.

2. The Manopt toolbox (https://www.manopt.org) must be installed and set up correctly.

Run Instructions:

1. Run importconstrained to add the required paths.

2. Navigate to the connectors subfolder and execute generate_mex.m to generate the necessary MEX files.

3. To test the fixed-rank matrix completion problem, run: demo_fixedrank_acc.m \ demo_fixedrank_rub.m \ demo_fixedrank_plot.m.

4. To test the nonnegative PCA problem, run : demo_Oblique_acc.m \ demo_Oblique_rub.m.

For any questions or issues, please contact: mahehao@mail.scut.edu.cn
