import sys
import json
import os.path
import base64
import stat

from functools import partial, reduce
from itertools import chain, compress, izip_longest, groupby
from ast       import literal_eval
from scipy.special import factorial
from matplotlib.path import Path

import numpy as np
import scipy.linalg as la
import scipy.optimize as op

##-- BASIC FUNCTIONS --##
Fst  = lambda xs: xs[0]
Snd  = lambda xs: xs[1]
Thd  = lambda xs: xs[2]
Lst  = lambda xs: xs[-1]
Not  = lambda x: not x
Neg  = lambda x: -x
Id   = lambda x: x
Get  = lambda key: lambda x: x[key]
Mult = lambda x,y: x*y
MethodToFunction = lambda method, *args: lambda x: x[method](*args)
Plus = lambda x,y: x+y
Minus = lambda x,y: x-y

##-- STRING FUNCTIONS --##
splitfirst = lambda s, t: (s.split(t)[0], t.join(s.split(t)[1:]))

##-- LIST FUNCTIONS --##
flatten        = lambda xs: list(chain(*xs))
IsSubset       = lambda ys, xs: all( [x in ys for x in xs] )
OneIsSubset    = lambda ys, xs: IsSubset(ys, xs) if len(ys) >= len(xs) else IsSubset(xs, ys)
product        = lambda xs: reduce( Mult, xs )
sumsquares     = lambda xs: sum([x**2 for x in xs])
AllEqual       = lambda xs: len(set(map(str,xs))) <= 1
unique         = lambda xs: list(set(xs))

##-- DICTIONARY FUNCTIONS --##
sub_object      = lambda except_keys, obj: { k:v for k,v in obj.items() if k not in except_keys }
EqualExcept     = lambda except_keys, o1, o2: sub_object(except_keys,o1) == sub_object(except_keys,o2)
sub_dict        = lambda keys, obj: { k:obj[k] for k in keys }
SubsetEqual     = lambda keys, o1, o2: sub_dict( keys, o1 ) == sub_dict( keys, o2 )
ValuesAreSubset = lambda subset_key , o1, o2: OneIsSubset( o1[subset_key], o2[subset_key] )
extend          = lambda objects: { k:v for object in objects for k,v in object.items() }

##-- FILTER FUNCTIONS --##
cmap           = lambda f, pred, xs: [ f(x) if p else x for p,x in zip(pred,xs) ]
common_items   = lambda xs, key=Id: list( set.intersection( *[ set(key(x)) for x in xs ] ) )
FilterMethod   = lambda xs, method, args=None: filter( lambda x: method(x,args), xs )
_SplitFilter   = lambda selector, xs: ( list(compress(xs,selector)), list(compress(xs,map(Not,selector))) )
SplitFilter    = lambda func, xs: _SplitFilter( map( func, xs ), xs )
overlap        = lambda xs, ys: [ x for x in xs if x in ys ]
exclude        = lambda xs, ys: [ x for x in xs if x not in ys ]

def Partition( transative_comparison, xs, aux_match=None, aux_match_comparison=None ):
  res = []
  while len(xs):
    #print map(partial(transative_comparison,xs[0]),xs)
    part,xs = SplitFilter( partial(transative_comparison, xs[0]), xs )
    if aux_match and aux_match_comparison:
      for y in aux_match:
        if aux_match_comparison( y, part[0] ):
          part.append(y)
    res.append( part )
  return res

def RePartition( transative_comparison, xss ):
  return flatten( map( partial(Partition,transative_comparison), xss ) )

def ArgvToArgsAndFiles( arg_map, *file_extensions ):
  proto_files = filter( lambda s: '=' not in s, sys.argv[1:] )
  files = [ FilterMethod( proto_files, str.endswith, args='.'+file_extension ) for file_extension in file_extensions ]
  args  = [ arg for arg in sys.argv[1:] if arg not in flatten( files )  ]
  args = { k[1:]:v.replace( '~', ' ' ) for k,v in map( lambda s: splitfirst(s,'=') if '=' in s else (s,''), args ) }
  args = { k:arg_map[v] if v in arg_map else v for k,v in args.items() }
  args = dict([ (k,v) if v != '' else (k,True) for k,v in args.items() ])
  args = dict([ (arg_map[k],v) if k in arg_map else (k,v) for k,v in args.items() ])
  print sys.argv[0], args
  if len(files):
    return tuple( [args] + files )
  else:
    return args
  

sng_report = lambda a,b:   '{:40} {:<11}'       .format( a, b )
dbl_report = lambda a,b,c: '{:40} {:<11} -> {:}'.format( a, b, c )
running_rep = lambda digits=2: '\r{:40} {:'+str(digits)+'}/{:'+str(digits)+'}  '
def progmap( msg, f, xs, *args, **kwargs ):
  length = len(xs)
  digits = len(str(length))
  fmt = running_rep(digits)
  res = []
  for i,x in enumerate(xs):
    sys.stdout.write( fmt.format( msg, i, length ) )
    sys.stdout.flush()
    res.append( f(x,*args,**kwargs) )
  #sys.stdout.write( '\r' )
  sys.stdout.write( fmt.format( msg, length, length ) )
  sys.stdout.write( '\n' )
  sys.stdout.flush()
  return res

def progpoolmap( msg, pool, f, xs ):
  length = len(xs)
  digits = len(str(length))
  fmt    = running_rep(digits)
  for i,_ in enumerate(pool.imap_unordered( f, xs )):
    sys.stdout.write( fmt.format( msg, i, length ) )
    sys.stdout.flush()
  sys.stdout.write( fmt.format( msg, length, length ) )
  sys.stdout.write( '\n' )
  sys.stdout.flush()

def progprogmap( msg, f, xs ):
  length = len(xs)
  digits = len(str(length))
  fmt_gen = lambda d1, d2: '\r{:40} {:'+str(d1)+'}/{:'+str(d1)+'}  {:'+str(d2)+'}/{:'+str(d2)+'}'
  fmt    = fmt_gen(digits,digits)
  res = []
  def _print( msg, i, l, j, len2 ):
    fmt = fmt_gen( digits, len(str(len2)) )
    sys.stdout.write( fmt.format( msg, i, l, j, len2 ) )
    sys.stdout.flush()
  _print( msg, 0, length, '?', '?' )
  for i,x in enumerate(xs):
    res.append( f(lambda j,l: _print( msg, i, length, j, l), x) )
  sys.stdout.write( '\r' )
  return res

def nmap( *args ):
  functions = args[:-1]
  xs = args[-1]
  composed_function = reduce( lambda f, g: lambda x: f(g(x)), functions, lambda x: x )
  return map( composed_function, xs )
  #return reduce( lambda f, xs: map(f,xs), args[:-1:-1], args[-1] )

def get_path_filename( file ):
  return os.path.split(file)

def change_extension( file, ext ):
  return '.'.join( file.split('.')[:-1] + [ext] )

##-- NUMPY STUFF --##
Chi2 = lambda data_mean, data_err, fit, dof=1, **extra: np.sum( ( (np.asarray(data_mean) - np.asarray(fit)) / np.asarray(data_err) )[:]**2 ) / dof
def Chi2Cov( data_mean, cov, fit, dof=1 ):
  mat = fit-data_mean
  try:
    invcov = np.linalg.inv( cov )
    return np.dot( np.dot( np.transpose(mat), invcov ), mat ) / dof
  except:
    return None

def cov_object( boot, dict ):
  if type(dict) == np.ndarray and len(boot.shape) == 2:
    dict['cov'] = np.cov(boot, rowvar=False) 
  return dict

def mean_var( xss, axis=1 ):
  return np.mean( xss, axis ), np.std( xss, axis, ddof=1 )
def mystd( xss, axis, ddof=1 ):
  if  'dtype' in dir(xss) and xss.dtype == np.complex128:
    return np.std(xss.real, axis, ddof=ddof) + 1.0j*np.std(xss.imag, axis, ddof=ddof)
  else:
    return np.std(xss,axis,ddof=ddof)
boot_object = lambda boot, **extra: extend( [ { 'boot':np.asarray(boot), 'std':mystd(boot,0,ddof=extra.pop('ddof',1.0)), 'mean':np.mean(boot,0) }, cov_object(boot,extra) ] )

def FitPoly1D( xs, ys, ss, coefs ):
  def _poly_func(x, *cs):
    return np.sum([ c*xs**cp for c,cp in zip(cs,coefs) ],0)
  return op.curve_fit( _poly_func, xs, ys, sigma=ss, absolute_sigma=True, p0=np.ones(len(coefs)) )
    
def FitPoly1D( xs, ys, ss, coefs ):
  A = np.array( [ x**i for x in xs for i in coefs ] ).reshape( len( ys ), len( coefs ) )
  A = np.transpose( np.multiply( np.transpose( A ), np.reciprocal( ss ) ) )
  B = np.divide( ys, ss )
  return la.lstsq( A, B )
def BootFitPoly1D( xs, yss, ss, coefs, **kwargs ):
  factor    = kwargs.pop( 'factor', 1.0 )
  boot_func = kwargs.pop( 'boot_func', boot_object )
  data_type = kwargs.pop( 'data_type', 'real' )
  xs  = np.asarray( xs )
  yss = np.asarray( yss )
  ss  = np.asarray( ss )
  if data_type == 'real':
    yss = yss.real
    ss  = ss.real
  elif data_type == 'imag':
    yss = yss.imag
    ss  = ss.imag
  else:
    print 'unknown data_type', data_type
  result = []
  for i in range(yss.shape[0]):
    result.append( FitPoly1D( xs, yss[i,:], ss, coefs )[0] / factor )
  result = boot_func( result, coefs=coefs, type='poly', **kwargs )
  return result

def FitExponential1D( xs, ys, ss ):
  result,residual,rank,singular = FitPoly1D( xs, np.log(np.abs(ys)), np.log(np.abs(ss)), [ 0, 1 ] )
  result[ 0 ] = np.exp( result[ 0 ] )
  return result

def BootFitExponential1D( xs, yss, ss, **kwargs ):
  xs = np.asarray(xs)
  factor = kwargs.pop( 'factor', 1.0 )
  result = []
  data_mean = np.mean( yss, 0 )
  data_err  = np.std(  yss, 0 )
  for i in range(yss.shape[0]):
    result.append( FitExponential1D( xs, yss[i,:], ss ) )
    result[i][1] /= factor
  booter = kwargs.pop( 'boot_func', boot_object )
  result = booter( result, **kwargs )
  cov    = np.cov( yss, rowvar=False )
  fit = result['mean'][0] * np.exp( result['mean'][1] * xs )
  result['chi2dof'] = Chi2Cov( data_mean, cov, fit, dof=(data_mean.shape[0]-3) )
  return result

cosh_func = lambda x, a, b: a * np.cosh( b * x )
def FitCosh1D( xs, ys, ss ):
  res, covariance_matrix = op.curve_fit( cosh_func, xs, ys, sigma=ss, absolute_sigma=True )
  return res
def BootFitCosh1D( xs, yss, ss, **kwargs ):
  xs = np.asarray(xs)
  offset = kwargs.get( 'offset', 0.0 )
  result = []
  data_mean = np.mean( yss, 0 )
  data_err  = np.std(  yss, 0 )
  for i in range(yss.shape[0]):
    result.append( FitCosh1D( xs - offset, yss[i,:], ss ) )
  booter = kwargs.pop( 'boot_func', boot_object )
  result = booter( result, **kwargs )
  cov = np.cov( yss, rowvar=False )
  fit = result['mean'][0] * np.cosh( result['mean'][1] * ( xs - offset ) )
  result['chi2dof'] = Chi2Cov( data_mean, cov, fit, dof=(data_mean.shape[0]-3) )
  return result

dipole_func = 'lambda x, a, b: a / ( 1 + x/b**2 )**2'
normed_dipole_func = 'lambda x, b: 1.0 / ( 1.0 + x/b**2 )**2'
poly_deg1_imp = lambda x,a2,b2: a2 + b2 * x
dipole_func_imp = eval(dipole_func)
dipole_jac  = lambda x, a, b: np.asarray([ 1 / ( 1 + x/b**2 )**2, -2.0 / b**2 / ( 1 + x/b**2 )**3 ])

def FitCurve1D( xs, ys, ss, curve, **kwargs ):
  #print xs, ys, ss
  #print xs.shape, ys.shape, ss.shape
  res, covariance_matrix = op.curve_fit( curve, xs, ys, sigma=ss, absolute_sigma=True, maxfev=10000000 )
  return res
def BootFitCurve1DOverBootX( xss, yss, ss, curve, **kwargs ):
  print xss.shape, yss.shape, ss.shape
  evalled_curve = eval(curve)
  jac = kwargs.pop('jac',None)
  data_mean = np.mean( yss, 0 )
  data_err  = np.std(  yss, 0 )
  result = progmap( 'Performing boot curve fit', lambda i: FitCurve1D( xss[i,:], yss[i,:], ss, evalled_curve, jac=kwargs.get('jac',None) ), range(yss.shape[0]) )
  booter = kwargs.pop( 'boot_func', boot_object )
  result = booter( result, type='function', function=curve, **kwargs )
  print sng_report('Performed boot curve fit', result['boot'].shape)
  cov    = np.cov( yss, rowvar=False )
  #print np.mean(xss,0).shape
  def _configure_axes( xs ):
    if len(np.asarray(xs).shape) > 1:
      return np.swapaxes(xs,0,1)
    else:
      return xs
  fit    = np.asarray([ evalled_curve( x, *result['mean'].tolist() ) for x in _configure_axes(np.mean(xss,0).tolist())  ])
  result['chi2dof'] = Chi2Cov( data_mean, cov, fit, dof=(data_mean.shape[0]-len(result['mean'].tolist())-1) )
  return result
def BootFitCurve1D( xs, yss, ss, curve, **kwargs ):
  #print xs.shape, yss.shape, np.asarray([xs for i in range(yss.shape[0])]).shape
  return BootFitCurve1DOverBootX( np.asarray([xs for i in range(yss.shape[0])]), yss, ss, curve, **kwargs )

def DecomposeLinearSystem( A, bm, bs ):
  A, bm, bs = np.asarray( A ), np.asarray( bm ), np.asarray( bs )
  f = lambda xs: (np.tensordot( A, xs, (1,0) ) - bm) / bs
  #print A.shape, bm.shape, bs.shape
  res = op.least_squares( f, np.zeros(A.shape[-1]) )
  return res['x']

def BootDecomposeLinearSystem( A, bb ):
  A, bb = np.asarray( A ), np.swapaxes(np.asarray( bb ),0,1)
  bs = np.std( bb, ddof=1 )
  nb = bb.shape[0]
  res = np.zeros( (nb,A.shape[-1]) )
  for i in range(nb):
    res[i,:] = DecomposeLinearSystem( A, bb[i,:], bs )
  return res

class NumpyJsonEncoder(json.JSONEncoder):
  def default(self,obj):
    if isinstance(obj, np.ndarray):
      #return obj.tolist()
      return dict(__ndarray__=base64.b64encode(obj.data),dtype=str(obj.dtype),shape=obj.shape)
    return json.JSONEncoder.default(self,obj)

def NumpyJSONHook(dct):
  if isinstance(dct,dict) and '__ndarray__' in dct:
    data = base64.b64decode( dct['__ndarray__'] )
    return np.frombuffer( data, dct['dtype'] ).reshape( dct['shape'] )
  return dct

def LoadJSONFile( filename ):
  #print filename
  file = open( filename, 'r' )
  res  = json.load( file, object_hook=NumpyJSONHook )
  file.close()
  return res

def WriteJSONFile( filename, d ):
  json.dump( d, open( filename, 'w+' ), cls=NumpyJsonEncoder )

def ErrVal( mean, std, sig, psign=False ):
  if type(mean) is np.ndarray and len(mean.shape) == 1 and type(std) is np.ndarray and len(mean.shape) == 1:
    if mean.shape[0] == 1 and std.shape[0] == 1:
      mean = Fst(mean)
      std = Fst(std)
    elif mean.shape[0] == std.shape[0]:
      return [ ErrVal( m, s, sig, psign ) for m,s in zip(mean.tolist(),std.tolist()) ]
  if mean == 0.0 and std == 0.0:
    return '0(0)'
  if std == 0.0:
    return '{:}'.format(mean)

  mean_mag = int(np.floor( np.log10( np.abs(mean) ) )) #Get the number of digits
  std_mag  = int(np.floor( np.log10( np.abs(std)  ) ))

  sig = max(sig,std_mag+1) #If the number is to the left of the decimal increase sig figs
  digits_after_decimal = max( 0, sig - std_mag - 1 )
  fmt_str = '{:{sign_fmt}.{digits}f}'
  if sig:
    fmt_str = fmt_str + '({:.0f})'
  if psign:
    sign_fmt = '+'
  else:
    sign_fmt = ''
  return fmt_str.format( mean, std * 10**(-std_mag+sig-1), sign_fmt=sign_fmt, digits=digits_after_decimal )
def ErrVal2( mean, std ):
  return ErrVal( mean, std, 2 )
def ErrVal2D( d ):
  return ErrVal2( d['mean'], d['std'] )

def unit( xs ):
  return xs / np.linalg.norm( xs )
def GenVal( mean, std, num ):
  return mean + std*np.random.randn(num)

def MinMax( *args, **kwargs ):
  scale = kwargs.pop('scale',1.0)-1.0
  args = filter( lambda d: np.asarray(d).shape != (0,), args )
  minv, maxv = min([ np.min(d) for d in args ]), max([ np.max(d) for d in args ])
  if minv != minv:
    minv = 0
  if maxv != maxv:
    maxv = 0
  extension = (maxv-minv)*scale
  return minv-extension, maxv+extension


def PointsInside( fourpoints, points, radius=0.1 ):
  verts = fourpoints + [ fourpoints[0] ]
  codes  = [Path.MOVETO] + [Path.LINETO]*(len(verts)-2) + [Path.CLOSEPOLY]
  path   = Path( verts, codes )
  return points[ path.contains_points( points, radius=radius ) ].tolist()

def afterslash( str, tok ):
  s = str.split('/')
  s[-1] = tok + s[-1]
  return '/'.join( s )

def threevec( xs ):
  return ''.join(map( lambda x: '{0:+}'.format(x), xs ))

def RollBackwards( lattice, source, offset=0 ):
  return np.roll( np.roll( np.roll( np.roll( lattice, -source[0], 3+offset ), -source[1], 2+offset ), -source[2], 1+offset ), -source[3], 0+offset )

def chain_to_dict( d, ks, val ):
  kf = ks.pop(-1)
  for k in ks:
    if k not in d:
      d[k] = {}
    d = d[k]
  d[kf] = val

def valify( d ):
  if 'val' in d:
    return d['val']
  else:
    return { k:valify(v) for k,v in d.items() }

_boot_object_errval = lambda boot_obj, sig: extend([ boot_obj, {'val':ErrVal( boot_obj['mean'], boot_obj['std'], sig )} ])
boot_object_errval = lambda boot, sig, **extra: _boot_object_errval( boot_object( boot, **extra ), sig )
boot_object_errval_2 = lambda boot, **extra: boot_object_errval( boot, 2, **extra )

def outerlist( xss ):
  return [ xss[i] for i in range(xss.shape[0]) ]

def splitlist( xs, num ):
  return list(izip_longest(*[iter(xs)]*num, fillvalue=None))

def nullspace(A, atol=1e-13, rtol=0):
    """Compute an approximate basis for the nullspace of A.

    The algorithm used by this function is based on the singular value
    decomposition of `A`.

    Parameters
    ----------
    A : ndarray
        A should be at most 2-D.  A 1-D array with length k will be treated
        as a 2-D with shape (1, k)
    atol : float
        The absolute tolerance for a zero singular value.  Singular values
        smaller than `atol` are considered to be zero.
    rtol : float
        The relative tolerance.  Singular values less than rtol*smax are
        considered to be zero, where smax is the largest singular value.

    If both `atol` and `rtol` are positive, the combined tolerance is the
    maximum of the two; that is::
        tol = max(atol, rtol * smax)
    Singular values smaller than `tol` are considered to be zero.

    Return value
    ------------
    ns : ndarray
        If `A` is an array with shape (m, k), then `ns` will be an array
        with shape (k, n), where n is the estimated dimension of the
        nullspace of `A`.  The columns of `ns` are a basis for the
        nullspace; each element in numpy.dot(A, ns) will be approximately
        zero.
    """

    A = np.atleast_2d(A)
    u, s, vh = np.linalg.svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns

def GetDecomposableElements( xss ):
  xss = np.asarray(xss)
  nspace = nullspace( xss )
  num_combs = xss.shape[0]
  num_terms, failing_combinations = nspace.shape
  unuseable_terms = [ j for j in range(num_terms) if np.sum(np.any(nspace[j,:] != 0)) ]
  useable_combinations = map(Fst,filter(Snd,enumerate([ np.all(xss[j,unuseable_terms]==0) for j in range(num_combs) ])))
  return useable_combinations, [ j for j in range(num_terms) if j not in unuseable_terms ]

#def EqualExcept( except_keys, o1, o2 ):
#  s1 = sub_object(except_keys,o1) 
#  s2 = sub_object(except_keys,o2)
#  print s1
#  print s2
#  return s1 == s2

def create_executable_script( f, s ):
  open( f, 'w+' ).write( s )
  os.chmod( f, os.stat(f).st_mode | stat.S_IEXEC )

def Collect( keys, differ_by_keys ): #Collect the keys that differe only in specified indicies, and then sorts by them
  groupby_key = range(len(keys[0]))
  for k in sorted(differ_by_keys)[::-1]:
    groupby_key.pop(k)
  group_key  = lambda xs: tuple([xs[i] for i in groupby_key])
  sort_key = lambda xs: group_key(xs) + tuple([xs[i] for i in differ_by_keys])
  #print groupby_key, group_key(keys[0]), keys[0]
  sorted_keys = sorted(keys,key=sort_key)
  res = groupby(sorted_keys,key=group_key)
  res = [ (k,[a for a in v]) for k,v in res ]
  #print map(Fst,res)
  return res
