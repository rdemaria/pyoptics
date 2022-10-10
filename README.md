# pyoptics
Classes to work with mad and LHC optics

## install
```
pip install git+https://github.com/rdemaria/pyoptics.git
```

Content:
- generic table manager with specific twiss/track/survey post-prossinging [src](https://github.com/rdemaria/pyoptics/blob/master/pyoptics/optics.py)
  - row filtering on regexp and column expression, new column on expressions
  - plots with lattices (automatic reload on changes)
  - reload on change
  - interpolation with s
  - table pretty printer
  - twiss/survey:
  - interpolation with s, cycle w/wo reset cumulative quantity, range extraction, append
  - twiss:
    - calculate beta max inside quadrupole, integral k*beta, normalized dispersion, transfer matrix in between points, 2D normalization matrix, partial integrals to compute Q', Q'
    - add errors from error tables
    - compute driving terms kernels
    - compute detuning factor
    - compute footprint based on detuning equations
  - survey:
     - redefine origin to any point (2D only)
     - compute rotation matrix
  - tfs reader/basic writer [src](https://github.com/rdemaria/pyoptics/blob/master/pyoptics/tfsdata/tfsdata.py)
- beam envelope class [src](https://github.com/rdemaria/pyoptics/blob/master/pyoptics/aperture.py)
     - draw 2d section of beam-envelopes and apertures x-y, s-x, s-y around reference orbit or lab frame
     - compute aperture margin
     - plot beam-beam separations
- madx language parser [src](https://github.com/rdemaria/pyoptics/tree/master/pyoptics/madlang/madlang.py)
     - parse MAD-X definitions (no support for macros) and build data structure including expressions
     - compute dependency matrix
 - harmonic fit [src](https://github.com/rdemaria/pyoptics/blob/master/pyoptics/harmonic_fit.py)
     - various tune fitting routines on turn-by-turn data
     - optics transition tables (table of strength) [src](https://github.com/rdemaria/pyoptics/blob/master/pyoptics/irplots.py)
     - plot, polynomial fit with consraints on deritives, convert fit to madx expressions
     - jacobian matching method in python [src](https://github.com/rdemaria/pyoptics/blob/master/pyoptics/jacobian.py)
 - LHC circuit model [src](https://github.com/rdemaria/pyoptics/blob/master/pyoptics/lhccircuit.py)
     - compute max ktoI and ItoK, max V, I', I'' using data from LSA or from IREF/VREF fit
     - compute minimum beam-process time base on I', I'' max
 - sdds reader [src](https://github.com/rdemaria/pyoptics/blob/master/pyoptics/sddsdata.py)
 - yasp reader [src](https://github.com/rdemaria/pyoptics/blob/master/pyoptics/yaspdata.py)
 - lot of obselete files ;-)
      
