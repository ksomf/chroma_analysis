#!/bin/env python

import numpy as np
from hl import *

dot   = lambda (a1,a2,a3), (b1,b2,b3): a1*b1 + a2*b2 + a3*b3
omega = lambda q, p: dot(q,p) / float(dot(q,q))
lattice_mom_to_discrete_mom = lambda momenta: np.sum( (np.asarray(momenta)*(2*np.pi/nx))**2, 0 )
dispersion_relation = lambda momenta: np.sqrt( (mass_boot**2) + lattice_mom_to_discrete_mom(momenta) )
nt = 64
nx = 32

for q in [ (2,1,0), (3,1,0), (3,2,0), (4,2,0), (4,3,0) ]:
  mf = LoadJSONFile( 'q{}{}{}'.format(*q)+'_v3/fit/nucleon_effmass_pp+0+0+0.json' )
  mass_boot = np.asarray( mf['eff']['boot'] )
  
  for base_momenta in [ (0,1,), (1,0,), (0,2,), (2,0,) ]:
    points = []
    _swaptuple = lambda (a,b,): (-a,-b,)
    def _modtuple( (a,b,c), z ):
      return (a,b,z)
    momenta = [ base_momenta + (0,), _swaptuple(base_momenta) + (0,) ]
    momenta = momenta + [ _modtuple( m, z ) for m in momenta for z in [ -1, 1 ] ]
    #print momenta
    format_string = lambda q: 'q{}{}{}'.format(*q)+'_v3/fit/nucleon_simple_even_quark2_pp{:+}{:+}{:+}_g3_cos_op_q' + '{:+}{:+}{:+}'.format(*q) + '_0.0125single_mom.json'
    for pp in momenta:
      jf = LoadJSONFile( format_string(q).format( *pp ) )
      boot = np.asarray(jf['eff']['boot']) * 2 * dispersion_relation( pp )
      w = omega(q,pp)
      points.append( (w,pp[-1],boot) )
    
    vals = map( Thd, points )
    F1_strict = np.mean( vals[0:2], 0 )
    F1_mixed  = np.mean( vals[2:], 0 )
    F2 = F1_strict - F1_mixed
    F1m, F1s = np.mean( F1_strict ), np.std( F1_strict, ddof=1 )
    F2m = np.mean( F2 )
    F2s = np.std( F2, ddof=1 )
    print q, base_momenta, points[0][0], ErrVal2( F1m, F1s ), ErrVal2( F2m, F2s )
