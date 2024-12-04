# Copyright (c) 2009-2024 The Regents of the University of Michigan.
# Part of HOOMD-blue, released under the BSD 3-Clause License.

"""Example pair potential."""

# Import the C++ module.
from hoomd.capillary_soft_shell import _capillary_soft_shell

# Impot the hoomd Python package and other necessary components.
from hoomd.md import pair
from hoomd.data.parameterdicts import TypeParameterDict
from hoomd.data.typeparam import TypeParameter


class CapillaryInteraction(pair.Pair):
    """Example pair potential."""

    # set static class data
    _ext_module = _capillary_soft_shell
    _cpp_class_name = "CapillaryInteraction"
    _accepted_modes = ("none", "shift", "xplor")

    def __init__(self, nlist, default_r_cut=None, default_r_on=0., mode='none'):
        super().__init__(nlist, default_r_cut, default_r_on, mode)
        params = TypeParameter('params', 'particle_types',
                               TypeParameterDict(q=float, len_keys=2))
        self._add_typeparam(params)
class SoftShell(pair.Pair):
    """Example pair potential."""

    # set static class data
    _ext_module = _capillary_soft_shell
    _cpp_class_name = "SoftShell"
    _accepted_modes = ("none", "shift", "xplor")

    def __init__(self, nlist, default_r_cut=None, default_r_on=0., mode='none'):
        super().__init__(nlist, default_r_cut, default_r_on, mode)
        params = TypeParameter('params', 'particle_types',
                               TypeParameterDict(f=float,r_cor=float,len_keys=2))
        self._add_typeparam(params)
