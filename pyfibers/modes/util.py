from functools import wraps
from scipy.special import hankel1, hankel2
from inspect import ismethod
import numpy as np

def hankel(j, nu, z):
    if j==1:
        return hankel1(nu, z)
    elif j==2:
        return hankel2(nu, z)
    raise ValueError("j must be either 1 or 2")

def hankelp(j, nu, z):
    """
    The derivative of hankel(j, nu, z) with respect to z
    """
    return hankel(j, nu-1, z) - nu*hankel(j, nu, z)/z


def vectorize_method(func):
    """
    A wrapper for numpy's vectorize which works for methods without explicitly passing the object
    :return:
    """
    @wraps(func)
    def inner(obj, *iargs, **ikwargs):

        def functionalized_method(*args, **kwargs):
            return func(obj, *args, **kwargs)

        return np.vectorize(functionalized_method, otypes=[np.dtype('c16')])(*iargs, **ikwargs)
        #return np.vectorize(functionalized_method)(*iargs, **ikwargs)
    return inner

