import ast

import numpy as np


def wrap_longitude(longitudes, base=0.0):
    """
    Convert points in longitude to range from base (default zero) to base+360
    """
    if isinstance(longitudes, list):
        longitudes = np.array(longitudes)
        intype = 'list'
    elif isinstance(longitudes, tuple):
        longitudes = np.array(longitudes)
        intype = 'tuple'
    else:
        intype = None
    wrapped_longs = ((longitudes - base + 720.) % 360.) + base

    if intype == 'list':
        wrapped_longs = wrapped_longs.tolist()
    if intype == 'tuple':
        wrapped_longs = tuple(wrapped_longs.tolist())

    return wrapped_longs


def safe_eval(expr, allow_calls=None, allow_methods=None, **kwargs):
    """
    Overly safe evaluation of arbitrary python code.

    Allows the same types as ast.literal_eval as well as:
        * Local namespace references (must be defined using keywords)
        * Object attributes
        * Binary, unary, comparison, boolean and arithmetic operators
        * Comprehension statements (e.g. list)
        * Subscripts and indexing
        * Explicitly allowed calls and methods
    """

    _whitelist = (ast.Expression, ast.Load,
                  ast.Name, ast.Attribute, ast.Subscript, ast.Index,
                  ast.boolop, ast.unaryop, ast.operator, ast.cmpop, ast.comprehension,
                  ast.BinOp,
                  ast.List, ast.Tuple, ast.Dict,
                  ast.Str, ast.Num)

    if allow_calls and (not hasattr(allow_calls, '__iter__')):
        allow_calls = [allow_calls]
    if allow_methods and (not hasattr(allow_methods, '__iter__')):
        allow_methods = [allow_methods]

    args_wl = []

    def check_wlist(n):
        return isinstance(n, _whitelist)

    def check_call(n):
        return allow_calls \
            and isinstance(n, ast.Call) and isinstance(n.func, ast.Name) \
            and (allow_calls == '*' or n.func.id in allow_calls)

    def check_method(n):
        return allow_methods \
            and isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) \
            and (allow_methods == '*' or n.func.attr in allow_methods)

    def check_args(n):
        return n in args_wl

    def check_node(n):
        return check_wlist(n) or check_call(n) or check_method(n) or check_args(n)

    tree = ast.parse(expr.strip(), mode='eval')
    nodes = [node for node in ast.walk(tree)]

    for node in nodes:
        if check_call(node) or check_method(node):
            args_wl += node.keywords

    unsafe = [type(node).__name__ for node in ast.walk(tree) if not check_node(node)]

    if unsafe:
        msg = 'The string to be evaluated ("{}") has one or more unpermitted expression types:\n    {}'
        raise TypeError(msg.format(expr, unsafe))

    try:
        result = eval(compile(tree, '<AST>', mode='eval'), kwargs)
    except NameError:
        msg = 'Expression "{}" is not fully described by keys {}'
        raise ValueError(msg.format(expr, kwargs.keys()))

    return result
