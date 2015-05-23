from scipy.special import hankel1, hankel2
import inspect

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


