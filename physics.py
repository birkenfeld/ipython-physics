# -*- coding: utf-8 -*-
"""
IPython (0.11) extension for physical quantity input.
See README.txt for usage examples.

Author: Georg Brandl <georg@python.org>.
This file has been placed in the public domain.
"""

import re
import sys
from math import pi
from IPython.core import ipapi

from Scientific.Physics.PhysicalQuantities import PhysicalQuantity, _addUnit

name = r'([_a-zA-Z]\w*)'
number = r'(-?[\d0-9.eE]+)'
unit = r'([a-zA-Z1][a-zA-Z0-9/*^-]*)'
quantity = number + r'\s*' + unit

inline_unit_re = re.compile(r'\((%s)\)' % quantity)
slash_conv_re = re.compile(r'^(.*?)//\s*%s$' % unit)
trailing_conv_re = re.compile(r'\s*//\s*%s$' % unit)
nice_assign_re = re.compile(r'^%s\s*=\s*(%s)$' % (name, quantity))
quantity_re = re.compile(quantity)
subst_re = re.compile(r'\?' + name)

def replace_inline(match):
    return 'Q(\'' + match.group(1).replace('^', '**') + '\')'
def replace_slash(match):
    expr = match.group(1)
    unit = str(match.group(2))  # PhysicalQuantity doesn't like Unicode strings
    if quantity_re.match(expr):
        return 'Q(\'' + expr + '\').inUnitsOf(%r)' % unit
    elif not expr:
        expr = '_'
    return '(' + expr + ').inUnitsOf(%r)' % unit
def replace_conv(match):
    return 'Q(\'' + match.group(1).replace('^', '**') + '\').inUnitsOf(%r)' % \
        str(match.group(4))
def replace_assign(match):
    return '%s = Q(\'%s\')' % (match.group(1), match.group(2).replace('^', '**'))

class QTransformer(object):
    # XXX: inheriting from PrefilterTransformer as documented gives TypeErrors,
    # but apparently is not needed after all
    priority = 99
    enabled = True
    def transform(self, line, continue_prompt):
        line = inline_unit_re.sub(replace_inline, line)
        if not continue_prompt:
            line = slash_conv_re.sub(replace_slash, line)
            line = nice_assign_re.sub(replace_assign, line)
        return line

def Q(v):
    try: return PhysicalQuantity(v)
    except NameError: raise ValueError('invalid unit in %r' % v)

def tbl_magic(shell, arg):
    """tbl <expr>: Evaluate <expr> for a range of parameters, given
    as "?name" in the expr.
    """
    unit = None
    match = trailing_conv_re.search(arg)
    if match:
        arg = arg[:match.start()]
        unit = match.group(1)
    substs = sorted(set(subst_re.findall(arg)))
    if not substs:
        raise ValueError('no substitutions in expr')
    while 1:
        expr = arg
        for subst in substs:
            try:
                val = raw_input('%s = ' % subst)
            except EOFError:
                sys.stdout.write('\n')
                return
            if not val:
                return
            if quantity_re.match(val):
                val = '(' + val + ')'
            expr = expr.replace('?' + subst, val)
        if unit:
            expr = '(' + expr + ').inUnitsOf("' + unit + '")'
        shell.run_cell(expr, False)

# monkey-patch a little
global_precision = [8]
PhysicalQuantity.__str__ = \
    lambda self: '%.*g %s' % (global_precision[0], self.value,
                              self.unit.name().replace('**', '^'))
PhysicalQuantity.__repr__ = PhysicalQuantity.__str__
PhysicalQuantity.__truediv__ = PhysicalQuantity.__div__
PhysicalQuantity.__rtruediv__ = PhysicalQuantity.__rdiv__
PhysicalQuantity.base = property(lambda self: self.inBaseUnits())
PhysicalQuantity.units = PhysicalQuantity.inUnitsOf

ip = ipapi.get()
# set up simplified quantity input
ip.user_ns['Q'] = Q
ip.prefilter_manager.register_transformer(QTransformer())
# setter for custom precision
ip.user_ns['setprec'] = lambda p: global_precision.__setitem__(0, p)
# quick evaluator
ip.define_magic('tbl', tbl_magic)

# active true float division
exec ip.compile('from __future__ import division', '<input>', 'single') \
    in ip.user_ns

# add some well used constants
ip.user_ns['c0'] = Q('299792458. m/s')
ip.user_ns['mu0'] = Q('4.e-7 pi*N/A**2').base
ip.user_ns['eps0'] = Q('1 1/mu0/c**2').base
ip.user_ns['Grav'] = Q('6.67259e-11 m**3/kg/s**2')
ip.user_ns['hpl'] = Q('6.62606957e-34 J*s')
ip.user_ns['hbar'] = ip.user_ns['hpl']/(2*pi)
ip.user_ns['e0'] = Q('1.60217733e-19 C')
ip.user_ns['me'] = Q('9.1093897e-31 kg')
ip.user_ns['mp'] = Q('1.6726231e-27 kg')
ip.user_ns['mn'] = Q('1.6749274e-27 kg')
ip.user_ns['NA'] = Q('6.0221367e23 1/mol')
ip.user_ns['kb'] = Q('1.380658e-23 J/K')

# essential units :)
_addUnit('furlong', '201.168*m', 'furlongs')
_addUnit('fortnight', '1209600*s', '14 days')

print
print 'Unit calculation and physics extensions activated.'
