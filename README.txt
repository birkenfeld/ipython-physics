ipython-physics
---------------

This is an extension for IPython 0.11 that at the moment mainly enables easy
input of physical quantities (i.e. numbers with units).  It requires the
"ScientificPython" (not SciPy) package by Konrad Hinsen.

Quick usage examples:

  In:  1 m // cm                        # convert between units
  Out: 100 cm                           # (syntax inspired by Mathematica)

  In:  (1 m)/(1 s)                      # sugar for inline quantity input
  Out: 1 m/s                            # in arbitrary expressions

  In:  Q('1 m')/Q('1 s')                # this is the desugared form
  Out: 1 m/s

  In:  // furlong/fortnight             # convert units in last result
  Out: 6012.8848 furlong/fortnight

  In:  alpha = 90 deg                   # more sugar for assignment: no
                                        # parentheses needed

  In:  sin(alpha)                       # angle units work with NumPy
  Out: 1.0                              # trigonometric functions

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

The predefined physical constants are:

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

Please let me know if anything is missing.
