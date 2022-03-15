import os 

def cor_path( cor, key ):
  return os.path.join(cor['path'],cor[key])

parameter_map = {
  'umd'   : 'quark1_minus_quark2'
, 'upd'   : 'quark1_plus_quark2'
, 'u'     : 'quark1'
, 'd'     : 'quark2'
, 'q1mq2' : 'quark1_minus_quark2'
, 'q1pq2' : 'quark1_plus_quark2'
, 'q1'    : 'quark1'
, 'q2'    : 'quark2'
, 'o'     : 'spin_odd'
, 'so'    : 'simple_odd'
, 'se'    : 'simple_even'
, 'ph'    : 'spin_partial_hilbert'
, 'e'     : 'spin_even'
, 'sr'    : 'simple_ratio'
, 'sioo'  : 'spin_independent_odd_odd_ratio'
, 'sipoo' : 'spin_independent_phase_odd_odd_ratio'
, 'ppoo'  : 'polarised_phase_odd_odd_ratio'
, 'oo'    : 'spin_odd_odd'
, 'oe'    : 'odd_even_ratio'
, 'eo'    : 'even_odd_ratio'
, 'ee'    : 'even_even_ratio'
, 'm'     : 'mass'
, 'cor'   : 'correlator'
, 'rw'    : 'reweight'
, 'lams'  : 'lambdas'
, '2d3'   : 'twodimthreepoint'
, '2d3r'  : 'twodimthreepointr'
}

def ChromaPrinter( (key, val) ):
  if key == 'dim':
    return '{:}x{:}'.format( val[0], val[-1] )
  elif key == 'kappas':
    return 'kappa' + 'x'.join( map(str,val) )
  elif key == 'qsq':
    return 'qsq' + str(val)
  elif key == 'plowersq':
    return 'plowersq' + str(val)
  elif key == 'puppersq':
    return 'puppersq' + str(val)
  elif key == 'label':
    return val
  else:
    print 'unrecognised key: ', key
    exit(1)
