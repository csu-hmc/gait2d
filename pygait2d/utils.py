#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sympy import srepr, Matrix, sympify


def save_sympy_matrix(matrix, filename):
    """Writes a matrix to file in the SymPy representation (srepr)."""
    num_rows, num_cols = matrix.shape
    with open(filename, 'w') as f:
        f.write(str(num_rows) + "\n")
        f.write(str(num_cols) + "\n")
        for expr in matrix:
            f.write(srepr(expr) + "\n")


def load_sympy_matrix(filename):
    """Loads a matrix from file created with save_sympy_matrix."""
    exprs = []
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                num_rows = int(line.strip())
            elif i == 1:
                num_cols = int(line.strip())
            else:
                exprs.append(sympify(line.strip()))
    return Matrix(exprs).reshape(num_rows, num_cols)
