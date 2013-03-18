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

import numpy as np

# allow uncertain values if the "uncertainties" package is available
try:
    from uncertainties import ufloat, Variable, AffineScalarFunc
    import uncertainties.umath as unp
    uncertain = (Variable, AffineScalarFunc)
    def valuetype((v, u)):
        if isinstance(v, uncertain):
            return v
        return ufloat((v, u))
except ImportError:
    uncertain = ()
    valuetype = lambda (v, u): v
    unp = np


class UnitError(ValueError):
    pass

# Adapted from ScientificPython:
# Written by Konrad Hinsen <hinsen@cnrs-orleans.fr>
# with contributions from Greg Ward
# last revision: 2007-5-25

class NumberDict(dict):
    """Dictionary storing numerical values.

    An instance of this class acts like an array of number with generalized
    (non-integer) indices. A value of zero is assumed for undefined
    entries. NumberDict instances support addition, and subtraction with other
    NumberDict instances, and multiplication and division by scalars.
    """

    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            return 0

    def __coerce__(self, other):
        if type(other) == type({}):
            other = NumberDict(other)
        return self, other

    def __add__(self, other):
        sum_dict = NumberDict()
        for key in self.keys():
            sum_dict[key] = self[key]
        for key in other.keys():
            sum_dict[key] = sum_dict[key] + other[key]
        return sum_dict

    def __sub__(self, other):
        sum_dict = NumberDict()
        for key in self.keys():
            sum_dict[key] = self[key]
        for key in other.keys():
            sum_dict[key] = sum_dict[key] - other[key]
        return sum_dict

    def __mul__(self, other):
        new = NumberDict()
        for key in self.keys():
            new[key] = other*self[key]
        return new
    __rmul__ = __mul__

    def __div__(self, other):
        new = NumberDict()
        for key in self.keys():
            new[key] = self[key]/other
        return new


# Type checks

def isPhysicalUnit(x):
    return hasattr(x, 'factor') and hasattr(x, 'powers')

def isPhysicalQuantity(x):
    return hasattr(x, 'value') and hasattr(x, 'unit')


class PhysicalUnit(object):
    """Physical unit.

    A physical unit is defined by a name (possibly composite), a scaling factor,
    and the exponentials of each of the SI base units that enter into it. Units
    can be multiplied, divided, and raised to integer powers.
    """

    def __init__(self, names, factor, powers, offset=0):
        """
        @param names: a dictionary mapping each name component to its
                      associated integer power (e.g. C{{'m': 1, 's': -1}})
                      for M{m/s}). As a shorthand, a string may be passed
                      which is assigned an implicit power 1.
        @param factor: a scaling factor
        @param powers: the integer powers for each of the nine base units
        @param offset: an additive offset to the base unit (used only for
                       temperatures)
        """
        if isinstance(names, basestring):
            self.names = NumberDict()
            self.names[names] = 1
        else:
            self.names = names
        self.factor = factor
        self.offset = offset
        self.powers = powers

    def set_name(self, name):
        self.names = NumberDict()
        self.names[name] = 1

    def name(self):
        num = ''
        denom = ''
        for unit in self.names.keys():
            power = self.names[unit]
            if power < 0:
                denom = denom + '/' + unit
                if power < -1:
                    denom = denom + '**' + str(-power)
            elif power > 0:
                num = num + '*' + unit
                if power > 1:
                    num = num + '**' + str(power)
        if len(num) == 0:
            num = '1'
        else:
            num = num[1:]
        return num + denom

    @property
    def is_dimensionless(self):
        return not reduce(lambda a, b: a or b, self.powers)

    @property
    def is_angle(self):
        return self.powers[7] == 1 and \
               reduce(lambda a, b: a + b, self.powers) == 1

    def __str__(self):
        return self.name()

    def __repr__(self):
        return '<PhysicalUnit ' + self.name() + '>'

    def __cmp__(self, other):
        if self.powers != other.powers:
            raise UnitError('Incompatible units')
        return cmp(self.factor, other.factor)

    def __mul__(self, other):
        if self.offset != 0 or (isPhysicalUnit(other) and other.offset != 0):
            raise UnitError('Cannot multiply units with non-zero offset')
        if isPhysicalUnit(other):
            return PhysicalUnit(self.names + other.names,
                                self.factor * other.factor,
                                map(lambda a,b: a+b, self.powers, other.powers))
        else:
            return PhysicalUnit(self.names + {str(other): 1},
                                self.factor * other, self.powers,
                                self.offset * other)

    __rmul__ = __mul__

    def __div__(self, other):
        if self.offset != 0 or (isPhysicalUnit(other) and other.offset != 0):
            raise UnitError('Cannot divide units with non-zero offset')
        if isPhysicalUnit(other):
            return PhysicalUnit(self.names - other.names,
                                self.factor / other.factor,
                                map(lambda a, b: a-b, self.powers, other.powers))
        else:
            return PhysicalUnit(self.names+{str(other): -1},
                                self.factor/other, self.powers)

    def __rdiv__(self, other):
        if self.offset != 0 or (isPhysicalUnit(other) and other.offset != 0):
            raise UnitError('Cannot divide units with non-zero offset')
        if isPhysicalUnit(other):
            return PhysicalUnit(other.names - self.names,
                                other.factor/self.factor,
                                map(lambda a,b: a-b, other.powers, self.powers))
        else:
            return PhysicalUnit({str(other): 1} - self.names,
                                other / self.factor,
                                map(lambda x: -x, self.powers))

    def __pow__(self, other):
        if self.offset != 0:
            raise UnitError('Cannot exponentiate units with non-zero offset')
        if isinstance(other, int):
            return PhysicalUnit(other*self.names, pow(self.factor, other),
                                map(lambda x,p=other: x*p, self.powers))
        if isinstance(other, float):
            inv_exp = 1./other
            rounded = int(np.floor(inv_exp + 0.5))
            if abs(inv_exp-rounded) < 1.e-10:
                if reduce(lambda a, b: a and b,
                          map(lambda x, e=rounded: x%e == 0, self.powers)):
                    f = pow(self.factor, other)
                    p = map(lambda x,p=rounded: x/p, self.powers)
                    if reduce(lambda a, b: a and b,
                              map(lambda x, e=rounded: x%e == 0,
                                  self.names.values())):
                        names = self.names/rounded
                    else:
                        names = NumberDict()
                        if f != 1.:
                            names[str(f)] = 1
                        for i in range(len(p)):
                            names[_base_names[i]] = p[i]
                    return PhysicalUnit(names, f, p)
                else:
                    raise UnitError('Illegal exponent %f' % other)
        raise UnitError('Only integer and inverse integer exponents allowed')

    def conversion_factor_to(self, other):
        """Return conversion factor to another unit."""
        if self.powers != other.powers:
            raise UnitError('Incompatible units')
        if self.offset != other.offset and self.factor != other.factor:
            raise UnitError(('Unit conversion (%s to %s) cannot be expressed ' +
                             'as a simple multiplicative factor') %
                             (self.name(), other.name()))
        return self.factor/other.factor

    def conversion_tuple_to(self, other):
        """Return conversion factor and offset to another unit."""
        if self.powers != other.powers:
            raise UnitError('Incompatible units')

        # let (s1,d1) be the conversion tuple from 'self' to base units
        #   (ie. (x+d1)*s1 converts a value x from 'self' to base units,
        #   and (x/s1)-d1 converts x from base to 'self' units)
        # and (s2,d2) be the conversion tuple from 'other' to base units
        # then we want to compute the conversion tuple (S,D) from
        #   'self' to 'other' such that (x+D)*S converts x from 'self'
        #   units to 'other' units
        # the formula to convert x from 'self' to 'other' units via the
        #   base units is (by definition of the conversion tuples):
        #     ( ((x+d1)*s1) / s2 ) - d2
        #   = ( (x+d1) * s1/s2) - d2
        #   = ( (x+d1) * s1/s2 ) - (d2*s2/s1) * s1/s2
        #   = ( (x+d1) - (d1*s2/s1) ) * s1/s2
        #   = (x + d1 - d2*s2/s1) * s1/s2
        # thus, D = d1 - d2*s2/s1 and S = s1/s2
        factor = self.factor / other.factor
        offset = self.offset - (other.offset * other.factor / self.factor)
        return (factor, offset)


# Helper functions

def _findUnit(unit):
    if isinstance(unit, basestring):
        name = unit.strip().replace('^', '**').replace('µ', 'mu').replace('°', 'deg')
        try:
            unit = eval(name, _unit_table)
        except NameError:
            raise UnitError('Invalid or unknown unit in %r' % unit)
        for cruft in ['__builtins__', '__args__']:
            try: del _unit_table[cruft]
            except: pass
    if not isPhysicalUnit(unit):
        raise UnitError(str(unit) + ' is not a unit')
    return unit


def _convertValue(value, src_unit, target_unit):
    (factor, offset) = src_unit.conversion_tuple_to(target_unit)
    return (value + offset) * factor


class PhysicalQuantity(object):
    """Physical quantity with units.

    PhysicalQuantity instances allow addition, subtraction, multiplication, and
    division with each other as well as multiplication, division, and
    exponentiation with numbers.  Addition and subtraction check that the units
    of the two operands are compatible and return the result in the units of the
    first operand. A limited set of mathematical functions (from numpy) is
    applicable as well.
    """

    global_precision = 8

    _number = re.compile(r'([+-]?[0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?)'
                         r'(?:\s+\+\/-\s+([+-]?[0-9]+(?:\.[0-9]*)?(?:[eE][+-]?[0-9]+)?))?')

    def __init__(self, value, unit=None, stdev=None):
        """There are two constructor calling patterns:

        1. PhysicalQuantity(value, unit), where value is any number and unit is
           a string defining the unit

        2. PhysicalQuantity(value_with_unit), where value_with_unit is a string
           that contains both the value and the unit, i.e. '1.5 m/s'. This form
           is provided for more convenient interactive use.
        """
        if unit is not None:
            self.value = valuetype((value, stdev or 0))
            self.unit = _findUnit(unit)
        else:
            s = value.strip()
            match = self._number.match(s)
            if match is None:
                raise UnitError('No number found in %r' % value)
            self.value = valuetype((float(match.group(1)),
                                    float(match.group(2) or 0)))
            self.unit = _findUnit(s[match.end(0):])

    def __str__(self):
        prec = self.global_precision
        unit = self.unit.name().replace('**', '^')
        if isinstance(self.value, uncertain):
            stdev = self.value.std_dev()
            if stdev:
                return '%.*g +/- %.*g %s' % (prec, self.value.nominal_value,
                                             prec, stdev, unit)
            return '%.*g %s' % (prec, self.value.nominal_value, unit)
        return '%.*g %s' % (prec, self.value, unit)

    def __repr__(self):
        return self.__str__()

    def _sum(self, other, sign1, sign2):
        if not isPhysicalQuantity(other):
            raise UnitError('Incompatible types')
        new_value = sign1 * self.value + \
            sign2 * other.value * other.unit.conversion_factor_to(self.unit)
        return self.__class__(new_value, self.unit)

    def __add__(self, other):
        return self._sum(other, 1, 1)

    __radd__ = __add__

    def __sub__(self, other):
        return self._sum(other, 1, -1)

    def __rsub__(self, other):
        return self._sum(other, -1, 1)

    def __cmp__(self, other):
        diff = self._sum(other, 1, -1)
        return cmp(diff.value, 0)

    def __mul__(self, other):
        if not isPhysicalQuantity(other):
            return self.__class__(self.value * other, self.unit)
        value = self.value * other.value
        unit = self.unit * other.unit
        if unit.is_dimensionless:
            return value * unit.factor
        else:
            return self.__class__(value, unit)

    __rmul__ = __mul__

    def __div__(self, other):
        if not isPhysicalQuantity(other):
            return self.__class__(self.value / other, self.unit)
        value = self.value / other.value
        unit = self.unit / other.unit
        if unit.is_dimensionless:
            return value * unit.factor
        else:
            return self.__class__(value, unit)

    def __rdiv__(self, other):
        if not isPhysicalQuantity(other):
            return self.__class__(other / self.value, pow(self.unit, -1))
        value = other.value / self.value
        unit = other.unit / self.unit
        if unit.is_dimensionless:
            return value * unit.factor
        else:
            return self.__class__(value, unit)

    __truediv__ = __div__
    __rtruediv__ = __rdiv__

    def __pow__(self, other):
        if isPhysicalQuantity(other):
            raise UnitError('Exponents must be dimensionless')
        return self.__class__(pow(self.value, other), pow(self.unit, other))

    def __rpow__(self, other):
        raise UnitError('Exponents must be dimensionless')

    def __abs__(self):
        return self.__class__(abs(self.value), self.unit)

    def __pos__(self):
        return self

    def __neg__(self):
        return self.__class__(-self.value, self.unit)

    def __nonzero__(self):
        return self.value != 0

    def __format__(self, *args, **kw):
        return "{1:{0}} {2}".format(args[0],self.value, self.unit)
    
    def convert(self, unit):
        """Change the unit and adjust the value such that the combination is
        equivalent to the original one. The new unit must be compatible with the
        previous unit of the object.
        """
        unit = _findUnit(unit)
        self.value = _convertValue(self.value, self.unit, unit)
        self.unit = unit

    def _round(self, x):
        if np.greater(x, 0.):
            return np.floor(x)
        else:
            return np.ceil(x)

    def to(self, *units):
        """Express the quantity in different units. If one unit is specified, a
        new PhysicalQuantity object is returned that expresses the quantity in
        that unit. If several units are specified, the return value is a tuple
        of PhysicalObject instances with with one element per unit such that the
        sum of all quantities in the tuple equals the the original quantity and
        all the values except for the last one are integers. This is used to
        convert to irregular unit systems like hour/minute/second.
        """
        units = map(_findUnit, units)
        if len(units) == 1:
            unit = units[0]
            value = _convertValue(self.value, self.unit, unit)
            return self.__class__(value, unit)
        else:
            units.sort()
            result = []
            value = self.value
            unit = self.unit
            for i in range(len(units)-1,-1,-1):
                value = value*unit.conversion_factor_to(units[i])
                if i == 0:
                    rounded = value
                else:
                    rounded = self._round(value)
                result.append(self.__class__(rounded, units[i]))
                value = value - rounded
                unit = units[i]
            return tuple(result)

    @staticmethod
    def any_to(qty, unit):
        if not isPhysicalQuantity(qty):
            qty = PhysicalQuantity(qty, 'rad')
        return qty.to(unit)

    @property
    def base(self):
        """Returns the same quantity converted to base units."""
        new_value = self.value * self.unit.factor
        num = ''
        denom = ''
        for i in xrange(9):
            unit = _base_names[i]
            power = self.unit.powers[i]
            if power < 0:
                denom += '/' + unit
                if power < -1:
                    denom += '**' + str(-power)
            elif power > 0:
                num += '*' + unit
                if power > 1:
                    num += '**' + str(power)
        if len(num) == 0:
            num = '1'
        else:
            num = num[1:]
        return self.__class__(new_value, num + denom)

    @property
    def cgs(self):
        """Returns the same quantity converted to cgs units."""
        new_value = self.value * self.unit.factor
        num = ''
        denom = ''
        for i in xrange(9):

            unit_name = _base_names[i]
            cgs_name = _cgs_names[i]
            power = self.unit.powers[i]

            conversion_factor = Q('1 '+unit_name).to(cgs_name).value
            new_value *= conversion_factor**power

            if power < 0:
                denom += '/' + cgs_name
                if power < -1:
                    denom += '**' + str(-power)
            elif power > 0:
                num += '*' + cgs_name
                if power > 1:
                    num += '**' + str(power)
        if len(num) == 0:
            num = '1'
        else:
            num = num[1:]

        return self.__class__(new_value, num + denom)

    # implementations of special functions, used by numpy ufuncs

    def sqrt(self):
        return pow(self, 0.5)

    def sin(self):
        if self.unit.is_angle:
            return unp.sin(self.value *
                           self.unit.conversion_factor_to(_unit_table['rad']))
        else:
            raise UnitError('Argument of sin must be an angle')

    def cos(self):
        if self.unit.is_angle:
            return unp.cos(self.value *
                           self.unit.conversion_factor_to(_unit_table['rad']))
        else:
            raise UnitError('Argument of cos must be an angle')

    def tan(self):
        if self.unit.is_angle:
            return unp.tan(self.value *
                           self.unit.conversion_factor_to(_unit_table['rad']))
        else:
            raise UnitError('Argument of tan must be an angle')


Q = PhysicalQuantity


# SI unit definitions

_base_names = ['m', 'kg', 's', 'A', 'K', 'mol', 'cd', 'rad', 'sr']

_base_units = [
    ('m',   PhysicalUnit('m',   1.,    [1,0,0,0,0,0,0,0,0])),
    ('g',   PhysicalUnit('g',   0.001, [0,1,0,0,0,0,0,0,0])),
    ('s',   PhysicalUnit('s',   1.,    [0,0,1,0,0,0,0,0,0])),
    ('A',   PhysicalUnit('A',   1.,    [0,0,0,1,0,0,0,0,0])),
    ('K',   PhysicalUnit('K',   1.,    [0,0,0,0,1,0,0,0,0])),
    ('mol', PhysicalUnit('mol', 1.,    [0,0,0,0,0,1,0,0,0])),
    ('cd',  PhysicalUnit('cd',  1.,    [0,0,0,0,0,0,1,0,0])),
    ('rad', PhysicalUnit('rad', 1.,    [0,0,0,0,0,0,0,1,0])),
    ('sr',  PhysicalUnit('sr',  1.,    [0,0,0,0,0,0,0,0,1])),
]

_cgs_names = ['cm', 'g', 's', 'abA', 'K', 'mol', 'cd', 'rad', 'sr']


_prefixes = [
    ('Y',  1.e24), ('Z',  1.e21), ('E',  1.e18), ('P',  1.e15), ('T',  1.e12),
    ('G',  1.e9),  ('M',  1.e6),  ('k',  1.e3),  ('h',  1.e2),  ('da', 1.e1),
    ('d',  1.e-1), ('c',  1.e-2), ('m',  1.e-3), ('mu', 1.e-6), ('n',  1.e-9),
    ('p',  1.e-12), ('f',  1.e-15), ('a',  1.e-18), ('z',  1.e-21),
    ('y',  1.e-24),
]

_unit_table = {}

for unit in _base_units:
    _unit_table[unit[0]] = unit[1]

def _addUnit(name, unit, comment=''):
    if _unit_table.has_key(name):
        raise KeyError('Unit ' + name + ' already defined')
    if type(unit) == type(''):
        unit = eval(unit, _unit_table)
        for cruft in ['__builtins__', '__args__']:
            try: del _unit_table[cruft]
            except: pass
    unit.set_name(name)
    _unit_table[name] = unit

def _addPrefixed(unit):
    _prefixed_names = []
    for prefix in _prefixes:
        name = prefix[0] + unit
        _addUnit(name, prefix[1]*_unit_table[unit])
        _prefixed_names.append(name)


# SI derived units; these automatically get prefixes
_unit_table['kg'] = PhysicalUnit('kg',   1., [0,1,0,0,0,0,0,0,0])

_addUnit('Hz', '1/s', 'Hertz')
_addUnit('N', 'm*kg/s**2', 'Newton')
_addUnit('Pa', 'N/m**2', 'Pascal')
_addUnit('J', 'N*m', 'Joule')
_addUnit('W', 'J/s', 'Watt')
_addUnit('C', 's*A', 'Coulomb')
_addUnit('V', 'W/A', 'Volt')
_addUnit('F', 'C/V', 'Farad')
_addUnit('ohm', 'V/A', 'Ohm')
_addUnit('S', 'A/V', 'Siemens')
_addUnit('Wb', 'V*s', 'Weber')
_addUnit('T', 'Wb/m**2', 'Tesla')
_addUnit('H', 'Wb/A', 'Henry')
_addUnit('lm', 'cd*sr', 'Lumen')
_addUnit('lx', 'lm/m**2', 'Lux')
_addUnit('Bq', '1/s', 'Becquerel')
_addUnit('Gy', 'J/kg', 'Gray')
_addUnit('Sv', 'J/kg', 'Sievert')
_addUnit('kat', 'mol/s', 'Katal')

_addUnit('abA', '10*A', 'Abampere')

del _unit_table['kg']

for unit in _unit_table.keys():
    _addPrefixed(unit)

# Fundamental constants, as far as needed to define other units
_unit_table['pi'] = np.pi
_addUnit('c0', '299792458.*m/s', 'speed of light')
_addUnit('mu0', '4.e-7*pi*N/A**2', 'permeability of vacuum')
_addUnit('eps0', '1/mu0/c0**2', 'permittivity of vacuum')
_addUnit('hplanck', '6.62606957e-34*J*s', 'Planck constant')
_addUnit('hbar', 'hplanck/(2*pi)', 'Planck constant / 2pi')
_addUnit('e0', '1.602176565e-19*C', 'elementary charge')
_addUnit('me', '9.10938291e-31*kg', 'electron mass')
_addUnit('kb', '1.3806488e-23*J/K', 'Boltzmann constant')

# Time units
_addUnit('min', '60*s', 'minute')
_addUnit('h', '60*min', 'hour')
_addUnit('d', '24*h', 'day')
_addUnit('wk', '7*d', 'week')
_addUnit('yr', '365.25*d', 'year')
_addPrefixed('yr')
_addUnit('fortnight', '1209600*s', '14 days')

# Length units
_addUnit('inch', '2.54*cm', 'inch')
_addUnit('ft', '12*inch', 'foot')
_addUnit('yd', '3*ft', 'yard')
_addUnit('mi', '5280.*ft', '(British) mile')
_addUnit('nmi', '1852.*m', 'Nautical mile')
_addUnit('Ang', '1.e-10*m', 'Angstrom')
_addUnit('AA', '1.e-10*m', 'Angstrom')
_addUnit('lyr', 'c0*yr', 'light year')
_addUnit('Bohr', '4*pi*eps0*hbar**2/me/e0**2', 'Bohr radius')
_addUnit('furlong', '201.168*m', 'furlongs')
_addUnit('au', '149597870691*m', 'astronomical unit')

# Area units
_addUnit('ha', '10000*m**2', 'hectare')
_addUnit('acres', 'mi**2/640', 'acre')
_addUnit('b', '1.e-28*m', 'barn')

# Volume units
_addUnit('l', 'dm**3', 'liter')
_addUnit('dl', '0.1*l', 'deci liter')
_addUnit('cl', '0.01*l', 'centi liter')
_addUnit('ml', '0.001*l', 'milli liter')
_addUnit('mul', '0.000001*l', 'micro liter')
_addUnit('tsp', '4.92892159375*ml', 'teaspoon')
_addUnit('tbsp', '3*tsp', 'tablespoon')
_addUnit('floz', '2*tbsp', 'fluid ounce')
_addUnit('cup', '8*floz', 'cup')
_addUnit('pt', '16*floz', 'pint')
_addUnit('qt', '2*pt', 'quart')
_addUnit('galUS', '4*qt', 'US gallon')
_addUnit('galUK', '4.54609*l', 'British gallon')

# Mass units
_addUnit('t', '1000*kg', 'Metric ton')
_addUnit('amu', '1.660538921e-27*kg', 'atomic mass units')
_addUnit('Da', '1*amu', 'Dalton')
_addUnit('oz', '28.349523125*g', 'ounce')
_addUnit('lb', '16*oz', 'pound')
_addUnit('ton', '2000*lb', 'US ton')

# Force units
_addUnit('dyn', '1.e-5*N', 'dyne (cgs unit)')

# Energy units
_addUnit('erg', '1.e-7*J', 'erg (cgs unit)')
_addUnit('eV', 'e0*V', 'electron volt')
_addUnit('Hartree', 'me*e0**4/16/pi**2/eps0**2/hbar**2', 'Wavenumbers/inverse cm')
_addUnit('Ken', 'kb*K', 'Kelvin as energy unit')
_addUnit('cal', '4.184*J', 'thermochemical calorie')
_addUnit('kcal', '1000*cal', 'thermochemical kilocalorie')
_addUnit('cali', '4.1868*J', 'international calorie')
_addUnit('kcali', '1000*cali', 'international kilocalorie')
_addUnit('Btu', '1055.05585262*J', 'British thermal unit')

_addPrefixed('eV')

# Electromagnetic units
_addUnit('G', '1e-4*T', 'Gauss')
_addUnit('Oe', '79.5774715*A/m', 'Oersted')

_addPrefixed('G')
_addPrefixed('Oe')

# Power units
_addUnit('hp', '745.7*W', 'horsepower')

# Pressure units
_addUnit('bar', '1.e5*Pa', 'bar (cgs unit)')
_addUnit('mbar', '1.e2*Pa', 'millibar')
_addUnit('kbar', '1.e8*Pa', 'kilobar')
_addUnit('atm', '101325.*Pa', 'standard atmosphere')
_addUnit('torr', 'atm/760', 'torr = mm of mercury')
_addUnit('psi', '6894.75729317*Pa', 'pounds per square inch')

# Angle units
_addUnit('deg', 'pi*rad/180', 'degrees')
_addUnit('arcmin', 'pi*rad/180/60', 'minutes of arc')
_addUnit('arcsec', 'pi*rad/180/3600', 'seconds of arc')
_unit_table['cycles'] = 2*np.pi

# Temperature units -- can't use the 'eval' trick that _addUnit provides
# for degC and degF because you can't add units
kelvin = _findUnit('K')
_addUnit('degR', '(5./9.)*K', 'degrees Rankine')
_addUnit('degC', PhysicalUnit(None, 1.0, kelvin.powers, 273.15),
         'degrees Celcius')
_addUnit('degF', PhysicalUnit(None, 5./9., kelvin.powers, 459.67),
         'degree Fahrenheit')
del kelvin

# Radiation-related units
_addUnit('Ci', '3.7e10*Bq', 'Curie')
_addUnit('rem', '0.01*Sv', 'Rem')

_addPrefixed('Ci')
_addPrefixed('rem')

# Astronomical units
_addUnit('Msol', '1.98892e30*kg', 'solar mass')
_addUnit('Lsol', '3.839e26*W', 'solar luminosity')
_addUnit('pc', '3.08568025e16*m')
_addPrefixed('pc')

# Important physical constants
_constants = [
    ('pi', np.pi),
    ('e', np.e),
    ('c0', Q('299792458. m/s')),
    ('mu0', Q('4.e-7 pi*N/A**2').base),
    ('eps0', Q('1 1/mu0/c0**2').base),
    ('Grav', Q('6.67384e-11 m**3/kg/s**2')),
    ('hpl', Q('6.62606957e-34 J*s')),
    ('hbar', Q('6.62606957e-34 J*s')/(2*pi)),
    ('e0', Q('1.602176565e-19 C')),
    ('me', Q('9.10938291e-31 kg')),
    ('mp', Q('1.672621777e-27 kg')),
    ('mn', Q('1.674927351e-27 kg')),
    ('NA', Q('6.02214129e23 1/mol')),
    ('kb', Q('1.3806488e-23 J/K')),
    ('g0', Q('9.80665 m/s**2')),
    ('R', Q('8.3144621 J/mol/K')),
    ('alpha', 7.2973525698e-3),
    ('Ry', Q('10973731.568539 1/m')),
    ('mu_n', Q('-0.96623647e-26 J/T')),
    ('gamma', Q('183.247179 MHz/T')),
    ('h0', 0.704),  # WMAP-7 + BAO constraint
    ('sigmaT', Q('6.652453e-29 m**2')),
]

name = r'([_a-zA-Z]\w*)'
number = r'(-?[\d0-9.eE-]+)'
unit = r'([a-zA-Z1°µ][a-zA-Z0-9°µ/*^-]*)'
quantity = number + r'(?:\s+\+\/-\s+' + number + ')?' + r'\s+' + unit

inline_unit_re = re.compile(r'\((%s)\)' % quantity)
slash_conv_re = re.compile(r'^(.*?)//\s*%s$' % unit)
slash_last_re = re.compile(r'^()\(/, %s\)$' % unit)
trailing_conv_re = re.compile(r'\s*//\s*%s$' % unit)
nice_assign_re = re.compile(r'^%s\s*=\s*(%s)$' % (name, quantity))
quantity_re = re.compile(quantity)
subst_re = re.compile(r'\?' + name)

def replace_inline(match):
    """Replace an inline unit expression, e.g. ``(1 m)``, by valid Python code
    using a Quantity call.
    """
    return '(Quantity(\'' + match.group(1) + '\'))'

def replace_slash(match):
    """Replace a double-slash unit conversion, e.g. ``c // km/s``, by valid
    Python code using a Quantity call.
    """
    expr = match.group(1)
    unit = str(match.group(2))  # PhysicalQuantity doesn't like Unicode strings
    if quantity_re.match(expr):
        expr = 'Quantity(\'' + expr + '\')'
    elif not expr:
        expr = '_'
    else:
        expr = '(' + expr + ')'
    if unit == 'base':
        return '(' + expr + ').base'
    if unit == 'cgs':
        return '(' + expr + ').cgs'
    else:
        return 'Quantity.any_to(%s, %r)' % (expr, unit)

def replace_assign(match):
    """Replace a pretty assignment, e.g. ``B = 1 T``, by valid Python code using
    a Quantity call.
    """
    return '%s = Quantity(\'%s\')' % (match.group(1), match.group(2))


class QTransformer(object):
    """IPython command line transformer that recognizes and replaces unit
    expressions.
    """
    # XXX: inheriting from PrefilterTransformer as documented gives TypeErrors,
    # but apparently is not needed after all
    priority = 99
    enabled = True
    def transform(self, line, continue_prompt):
        line = inline_unit_re.sub(replace_inline, line)
        if not continue_prompt:
            line = slash_conv_re.sub(replace_slash, line)
            line = nice_assign_re.sub(replace_assign, line)
            # lines that look like ``(/, unit)`` have been ``// unit`` but
            # already preprocessed by IPython, let's recognize them
            line = slash_last_re.sub(replace_slash, line)
        return line


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
            expr = 'Quantity.any_to((' + expr + '), %r)' % unit
        if hasattr(shell, 'shell'):  # later IPython versions
            shell.shell.run_cell(expr, False)
        else:
            shell.run_cell(expr, False)


q_transformer = QTransformer()


def load_ipython_extension(ip):
    # set up simplified quantity input
    ip.user_ns['Q'] = Q
    ip.user_ns['Quantity'] = Q
    ip.prefilter_manager.register_transformer(q_transformer)
    # setter for custom precision
    ip.user_ns['setprec'] = \
        lambda p: setattr(PhysicalQuantity, 'global_precision', p)
    # quick evaluator
    ip.define_magic('tbl', tbl_magic)

    # active true float division
    exec ip.compile('from __future__ import division', '<input>', 'single') \
        in ip.user_ns

    # add constants of nature
    for const, value in _constants:
        ip.user_ns[const] = value

    print 'Unit calculation and physics extensions activated.'

def unload_ipython_extension(ip):
    ip.prefilter_manager.unregister_transformer(q_transformer)
