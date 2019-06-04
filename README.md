# linear_algebra
some usefull stuff for computations in linear algebra with C++

1. Matrix determinant calculator with Gaussian elimination. Output example:

        Input matrix size n: 5
        Input elements of 5x5 matrix:
        2 3 4 1 8
        3 9 7 1 2
        3 9 4 7 1
        5 3 5 9 1
        4 0 8 1 2

        Det(A) = -2247.000

        Matrix after Gaussian elimination:
            5.000    3.000    5.000    9.000    0.000
            0.000    7.200    4.000   -4.400    0.000
            0.000    0.000    5.333   -7.667    0.000
            0.000    0.000    0.000    1.688    0.000
            0.000    0.000    0.000    0.000    6.935

2. System of linear equations solver. Output example:
                
        Input equations number: 3
        Input variables number: 3

        Example of input 3x4 augmented matrix (A|b):
        a1 a2 a3 b1
        a4 a5 a6 b2
        a7 a8 a9 b3

        Your input:
        4 2 1 1
        7 8 9 1
        9 1 3 2

        The system has a single unique solution <= rank(A|b) = rank(A) = n = 3

        Solution:
        X1 =     0.260870
        X2 =     0.043478
        X3 =    -0.130435

        Cramer's rule: Xi = det_i / det(A)

        Det(A) = 115.000000

        Cramer's determinants:
        det #1 =       30.000
        det #2 =        5.000
        det #3 =      -15.000

        (A|b) after Gauss elimination:
           9.000   1.000   0.000   2.391
           0.000   7.222   0.000   0.314
           0.000   0.000  -1.769   0.231
