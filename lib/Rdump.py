#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 12:50:59 2017

@author: meysamhashemi
"""

from six import PY2, string_types
from collections import Callable, Sequence
import numpy as np
from numbers import Number

import inspect
import io
import itertools
import logging
import math
import os
import random
import re
import sys
import tempfile
import time
import warnings




def is_legal_stan_vname(name):
    stan_kw1 = ('for', 'in', 'while', 'repeat', 'until', 'if', 'then', 'else',
                'true', 'false')
    stan_kw2 = ('int', 'real', 'vector', 'simplex', 'ordered', 'positive_ordered',
                'row_vector', 'matrix', 'corr_matrix', 'cov_matrix', 'lower', 'upper')
    stan_kw3 = ('model', 'data', 'parameters', 'quantities', 'transformed', 'generated')
    cpp_kw = ("alignas", "alignof", "and", "and_eq", "asm", "auto", "bitand", "bitor", "bool",
              "break", "case", "catch", "char", "char16_t", "char32_t", "class", "compl",
              "const", "constexpr", "const_cast", "continue", "decltype", "default", "delete",
              "do", "double", "dynamic_cast", "else", "enum", "explicit", "export", "extern",
              "false", "float", "for", "friend", "goto", "if", "inline", "int", "long", "mutable",
              "namespace", "new", "noexcept", "not", "not_eq", "nullptr", "operator", "or", "or_eq",
              "private", "protected", "public", "register", "reinterpret_cast", "return",
              "short", "signed", "sizeof", "static", "static_assert", "static_cast", "struct",
              "switch", "template", "this", "thread_local", "throw", "true", "try", "typedef",
              "typeid", "typename", "union", "unsigned", "using", "virtual", "void", "volatile",
              "wchar_t", "while", "xor", "xor_eq")
    illegal = stan_kw1 + stan_kw2 + stan_kw3 + cpp_kw
    if re.findall(r'(\.|^[0-9]|__$)', name):
        return False
    return not name in illegal


def _dict_to_rdump(data):
    parts = []
    for name, value in data.items():
        if isinstance(value, (Sequence, Number, np.number, np.ndarray, int, bool, float)) \
           and not isinstance(value, string_types):
            value = np.asarray(value)
        else:
            raise ValueError("Variable {} is not a number and cannot be dumped.".format(name))

        if value.dtype == np.bool:
            value = value.astype(int)

        if value.ndim == 0:
            s = '{} <- {}\n'.format(name, str(value))
        elif value.ndim == 1:
            s = '{} <-\nc({})\n'.format(name, ', '.join(str(v) for v in value))
        elif value.ndim > 1:
            tmpl = '{} <-\nstructure(c({}), .Dim = c({}))\n'
            # transpose value as R uses column-major
            # 'F' = Fortran, column-major
            s = tmpl.format(name,
                            ', '.join(str(v) for v in value.flatten(order='F')),
                            ', '.join(str(v) for v in value.shape))
        parts.append(s)
    return ''.join(parts)


def stan_rdump(data, filename):
    """
    Dump a dictionary with model data into a file using the R dump format that
    Stan supports.
    Parameters
    ----------
    data : dict
    filename : str
    """
    for name in data:
        if not is_legal_stan_vname(name):
            raise ValueError("Variable name {} is not allowed in Stan".format(name))
    with open(filename, 'w') as f:
         f.write(_dict_to_rdump(data))
