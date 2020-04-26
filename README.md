ipython-physics
---------------

This is an extension for IPython 0.11+ that at the moment mainly enables easy
input of physical quantities (i.e. numbers with units).  The implementation
is adapted from Konrad Hinsen's "Scientific.Physics.PhysicalQuantities" module.

If the "uncertainties" module (http://github.com/newville/uncertainties/) is
installed, you can also enter values with standard errors wherever you can enter
quantities, see the examples below.

You can have a look at a short notebook tutorial put together by Hans Fangohr at
http://www.southampton.ac.uk/~fangohr/blog/physical-quantities-numerical-value-with-units-in-python.html.

Quick installation advice:

  Place `physics.py` in any directory on PYTHONPATH, or add its directory to
  `sys.path` in your IPython `ipython_config.py` file.  Then you can load it
  either in the config file or on the command line as described here:
  http://ipython.org/ipython-doc/dev/config/extensions/index.html

Quick usage examples:

```
  In:  1 m // cm                        # convert between units
  Out: 100 cm                           # (syntax inspired by Mathematica)

  In:  (1 m)/(1 s)                      # sugar for inline quantity input
  Out: 1 m/s                            # in arbitrary expressions

  In:  Quantity('1 m')/Quantity('1 s')  # this is the desugared form
  Out: 1 m/s

  In:  // furlong/fortnight             # convert units in last result
  Out: 6012.8848 furlong/fortnight

  In:  alpha = 90 deg                   # more sugar for assignment: no
                                        # parentheses needed

  In:  sin(alpha)                       # angle units work with NumPy
  Out: 1.0                              # trigonometric functions

  In:  m = 80 +/- 5 kg                  # calculating with uncertainties
  In:  v = 130 +/- 10 m/s               # (needs the "uncertainties" module)
  In:  0.5 * m * v**2 // kJ
  Out: 676 +/- 112.25445 kJ

  In:  %tbl sqrt(?x**2 + ?y**2) // cm   # quickly tabulate a formula:
  x = 1 m                               # provide some values
  y = 2 m
  Out: 223.6068 cm                      # and get the result
  x = 3 m                               # ... this continues as long as you
  y = 4 m                               # enter new values
  Out: 500 cm

  In:  c0                               # important physical constants
  Out: 2.9979246e+08 m/s
  In:  setprec(4)                       # set the display precision
  In:  c0
  Out: 2.998e+08 m/s
```

The predefined constants are:

```
  pi
  e
  c0    -- vacuum speed of light
  mu0   -- magnetic constant
  eps0  -- electric constant
  Grav  -- Newton's constant
  hpl   -- Planck's constant
  hbar  -- Planck's constant / 2pi
  e0    -- elementary charge
  me    -- electron mass
  mp    -- proton mass
  mn    -- neutron mass
  NA    -- Avogadro's number
  kb    -- Boltzmann constant
  g0    -- Standard earth gravity
  R     -- Universal gas constant
  alpha -- fine structure constant
  Ry    -- Rydberg constant
  mu_n  -- Magnetic moment of the neutron
  gamma -- Gyromagnetic ratio of the neutron
  h0    -- dimensionless Hubble parameter
```

Please let me know if anything is missing.
