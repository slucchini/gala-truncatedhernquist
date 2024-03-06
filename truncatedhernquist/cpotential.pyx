# coding: utf-8
# cython: boundscheck=False
# cython: nonecheck=False
# cython: cdivision=True
# cython: wraparound=False

# Third-party
import numpy as np
cimport numpy as np
np.import_array()

# Project
from gala.potential.potential.cpotential cimport CPotentialWrapper
from gala.potential.potential.cpotential cimport energyfunc, gradientfunc, densityfunc, hessianfunc
from gala.potential import CPotentialBase, PotentialParameter

cdef extern from "src/potential.h":
    double truncated_hernquist_energy(double t, double *pars, double *q, int n_dim) nogil
    void truncated_hernquist_gradient(double t, double *pars, double *q, int n_dim, double *grad) nogil
    double truncated_hernquist_density(double t, double *pars, double *q, int n_dim) nogil

__all__ = ['TruncatedHernquistPotential']


cdef class TruncatedHernquistWrapper(CPotentialWrapper):

    def __init__(self, G, parameters, q0, R):
        self.init([G] + list(parameters),
                  np.ascontiguousarray(q0),
                  np.ascontiguousarray(R))
        self.cpotential.value[0] = <energyfunc>(truncated_hernquist_energy)
        self.cpotential.gradient[0] = <gradientfunc>(truncated_hernquist_gradient)
        self.cpotential.density[0] = <densityfunc>(truncated_hernquist_density)


class TruncatedHernquistPotential(CPotentialBase):
    r"""
    Parameters
    ----------
    m : :class:`~astropy.units.Quantity`, numeric [mass]
        Particle mass.
    units : `~gala.units.UnitSystem` (optional)
        Set of non-reducable units that specify (at minimum) the
        length, mass, time, and angle units.

    """
    m = PotentialParameter("m",physical_type='mass')
    c = PotentialParameter("c",physical_type='length')
    rmax = PotentialParameter('rmax',physical_type='length')

    Wrapper = TruncatedHernquistWrapper
