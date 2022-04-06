from hl import *

I = np.diag([ 1.0, 1.0, 1.0, 1.0 ])
trace    = lambda mat: np.sum(np.diag(mat))
tracemul = lambda xs: trace(reduce( np.dot, xs ))

minkowski_four_dot = lambda four_vec1, four_vec2: reduce(Minus, [ a*b for a,b in zip(four_vec1,four_vec2)])
minkowski_slash    = lambda four_vec: minkowski_four_dot( four_vec, minkowski_gammas )
minkowski_sigmavec = lambda mu, vec: 0.5j * (np.dot(minkowski_gammas[mu],minkowski_slash(vec)) - np.dot(minkowski_slash(vec),minkowski_gammas[mu]))
minkowski_dispersion_relation = lambda four_mom: np.asarray( [np.sqrt(np.sum(np.asarray(four_mom)**2))]+np.asarray(four_mom)[1:].tolist() )
def getAltMinkowskiKinematics( pp, p ):
  P, q = (pp+p)/2.0, (pp-p)
  q2   = minkowski_four_dot( q, q )
  return P, q, q2

euclidean_four_dot = lambda four_vec1, four_vec2: reduce(Plus, [ a*b for a,b in zip(four_vec1,four_vec2)])
euclidean_slash    = lambda four_vec: euclidean_four_dot( four_vec, euclidean_gammas )
euclidean_sigmavec = lambda mu, vec: 0.5j * (np.dot(euclidean_gammas[mu],euclidean_slash(vec)) - np.dot(euclidean_slash(vec),euclidean_gammas[mu]))
euclidean_dispersion_relation = lambda four_mom: np.asarray( [1.0j*np.sqrt(np.sum(np.asarray(four_mom)**2))]+np.asarray(four_mom)[1:].tolist() )
def getAltEuclideanKinematics( pp, p ):
  P, q = (pp+p)/2.0, (pp-p)
  q2   = euclidean_four_dot( q, q )
  return P, q, q2

minkowski_gamma0 = np.diag( [ 1.0+0.0j, 1.0, -1.0, -1.0 ] )
minkowski_gamma1 = -np.fliplr(minkowski_gamma0)
minkowski_gamma2 = -np.fliplr( np.diag( [ -1.0j, 1.0j, 1.0j, -1.0j ] ) )
minkowski_gamma3 = -np.asarray( [ [  0.0j, 0.0, 1.0,  0.0 ]
                                , [  0.0 , 0.0, 0.0, -1.0 ]
                                , [ -1.0 , 0.0, 0.0,  0.0 ]
                                , [  0.0 , 1.0, 0.0,  0.0 ] ] )
minkowski_gamma5 = np.abs(minkowski_gamma3)
minkowski_gammas = [ minkowski_gamma0, minkowski_gamma1, minkowski_gamma2, minkowski_gamma3, minkowski_gamma0 ]
minkowski_spin_up_projector        = (I+1.0j*np.dot(minkowski_gamma3,minkowski_gamma5))
minkowski_spin_down_projector      = (I-1.0j*np.dot(minkowski_gamma3,minkowski_gamma5))
minkowski_pol_projector            = -np.dot(minkowski_gamma3,minkowski_gamma5)
minkowski_pos_parity_projector     = (I+minkowski_gamma0)/2.0
minkowski_pos_spin_up_projector    = np.dot( minkowski_pos_parity_projector, minkowski_spin_up_projector   )
minkowski_pos_spin_down_projector  = np.dot( minkowski_pos_parity_projector, minkowski_spin_down_projector )
minkowski_pos_parity_pol_projector = np.dot( minkowski_pos_parity_projector, minkowski_pol_projector       )


euclidean_gamma1 = np.fliplr(np.diag( [ 1.0j, 1.0j, -1.0j, -1.0j ] ))
euclidean_gamma2 = np.fliplr(np.diag( [ 1.0, -1.0, -1.0, 1.0 ] ))
euclidean_gamma3 = np.asarray( [ [  0.0j, 0.0 , 1.0j,  0.0  ]
                               , [  0.0 , 0.0 , 0.0 , -1.0j ]
                               , [ -1.0j, 0.0 , 0.0 ,  0.0  ]
                               , [  0.0 , 1.0j, 0.0 ,  0.0  ] ] )
euclidean_gamma4 = np.diag([ 1.0, 1.0, -1.0, -1.0 ])
euclidean_gamma5 = np.abs(euclidean_gamma3)
euclidean_gammas = [ euclidean_gamma4, euclidean_gamma1, euclidean_gamma2, euclidean_gamma3, euclidean_gamma4 ]
euclidean_spin_up_projector         = (I+1.0j*np.dot(euclidean_gamma3,euclidean_gamma5))
euclidean_spin_down_projector       = (I-1.0j*np.dot(euclidean_gamma3,euclidean_gamma5))
euclidean_polx_projector            = 1.0j*np.dot(euclidean_gamma5,euclidean_gamma1)
euclidean_poly_projector            = 1.0j*np.dot(euclidean_gamma5,euclidean_gamma2)
euclidean_polz_projector            = 1.0j*np.dot(euclidean_gamma5,euclidean_gamma3)
euclidean_pol_projector             = euclidean_polz_projector
euclidean_pos_parity_projector      = (I+euclidean_gamma4)/2.0
euclidean_pos_spin_up_projector     = np.dot( euclidean_pos_parity_projector, euclidean_spin_up_projector   )
euclidean_pos_spin_down_projector   = np.dot( euclidean_pos_parity_projector, euclidean_spin_down_projector )
euclidean_pos_parity_polx_projector  = np.dot( euclidean_pos_parity_projector, euclidean_polx_projector       )
euclidean_pos_parity_poly_projector  = np.dot( euclidean_pos_parity_projector, euclidean_poly_projector       )
euclidean_pos_parity_polz_projector  = np.dot( euclidean_pos_parity_projector, euclidean_polz_projector       )
euclidean_pos_parity_pol_projector   = euclidean_pos_parity_polz_projector
euclidean_pos_parity_pol_tuple = np.asarray([ euclidean_pos_parity_projector, euclidean_pos_parity_polx_projector, euclidean_pos_parity_poly_projector, euclidean_pos_parity_polz_projector ])

def GluonTensorFF( mass, mu, nu, pp, p, P, q, sync_prop, src_prop, TwoPointGamma, Gamma, FFs ):
  #if any( map( lambda x:x==0, [2*np.abs(pp[0]), 2*np.abs(p[0]), tracemul([TwoPointGamma,sync_prop]), tracemul([TwoPointGamma,src_prop])] ) ):
  #print 2*np.abs(pp[0]), 2*np.abs(p[0]), tracemul([TwoPointGamma,sync_prop]), tracemul([TwoPointGamma,src_prop])
  fraction_factor =  ( I / 2.0 ) / np.sqrt( 2*np.abs(pp[0]) * 2*np.abs(p[0]) * tracemul([TwoPointGamma,sync_prop]) * tracemul([TwoPointGamma,src_prop]) )
  res = np.asarray(map(lambda inner: tracemul([Gamma,sync_prop,inner,src_prop,fraction_factor]), FFs))
  return res
def MinkowskiGluonTensorFF( mass, mu, nu, pp, p, TwoPointGamma=I, Gamma=I ):
  P, q, q2 = getAltMinkowskiKinematics( pp, p )
  sync_prop       = minkowski_slash(pp) + mass*I
  src_prop        = minkowski_slash(p)  + mass*I
  FFs  = [ minkowski_gammas[mu]*P[nu] + minkowski_gammas[nu]*P[mu]
         , 1.0j*(minkowski_sigmavec(mu,q)*P[nu] + minkowski_sigmavec(nu,q)*P[mu])/(2*mass)
         , (q[mu]*q[nu] + q[nu]*q[mu]) / mass * I ]
  return GluonTensorFF( mass, mu, nu, pp, p, P, q, sync_prop, src_prop, TwoPointGamma, Gamma, FFs ), q2
def EuclideanGluonTensorFF( mass, mu, nu, pp, p, TwoPointGamma=I, Gamma=I ):
  P, q, q2 = getAltEuclideanKinematics( pp, p )
  sync_prop       = -1.0j*euclidean_slash(pp) + mass*I
  src_prop        = -1.0j*euclidean_slash(p)  + mass*I
  FFs  = [ -1.0j*(euclidean_gammas[mu]*P[nu] + euclidean_gammas[nu]*P[mu])
         , -1.0j*(euclidean_sigmavec(mu,q)*P[nu] + euclidean_sigmavec(nu,q)*P[mu])/(2*mass)
         ,  -(q[mu]*q[nu] + q[nu]*q[mu]) / mass * I ]
  return GluonTensorFF( mass, mu, nu, pp, p, P, q, sync_prop, src_prop, TwoPointGamma, Gamma, FFs ), q2
def MinkowskiGluonIrrep1( mass, pp, p, TwoPointGamma=I, Gamma=I ):
  temporal, q2 = MinkowskiGluonTensorFF( mass, 0, 0, pp, p, TwoPointGamma, Gamma )
  spatial1, _  = MinkowskiGluonTensorFF( mass, 1, 1, pp, p, TwoPointGamma, Gamma )
  spatial2, _  = MinkowskiGluonTensorFF( mass, 2, 2, pp, p, TwoPointGamma, Gamma )
  spatial3, _  = MinkowskiGluonTensorFF( mass, 3, 3, pp, p, TwoPointGamma, Gamma )
  return (temporal + (spatial1+spatial2+spatial3)/3.0).tolist(), q2
def EuclideanGluonIrrep1( mass, pp, p, TwoPointGamma=I, Gamma=I ):
  temporal, q2 = EuclideanGluonTensorFF( mass, 0, 0, pp, p, TwoPointGamma, Gamma )
  spatial1, _  = EuclideanGluonTensorFF( mass, 1, 1, pp, p, TwoPointGamma, Gamma )
  spatial2, _  = EuclideanGluonTensorFF( mass, 2, 2, pp, p, TwoPointGamma, Gamma )
  spatial3, _  = EuclideanGluonTensorFF( mass, 3, 3, pp, p, TwoPointGamma, Gamma )
  return (-temporal + (spatial1+spatial2+spatial3)/3.0).tolist(), q2
def MinkowskiGluonIrrep2( mass, i, pp, p, TwoPointGamma=I, Gamma=I ):
  res,q2= MinkowskiGluonTensorFF( mass, 0, i, pp, p, TwoPointGamma, Gamma )
  return res.tolist(), q2
def EuclideanGluonIrrep2( mass, i, pp, p, TwoPointGamma=I, Gamma=I ):
  res,q2= EuclideanGluonTensorFF( mass, 0, i, pp, p, TwoPointGamma, Gamma )
  return res.tolist(), q2

def kinematics_transform( obj, return_q_mom2=False ): 
   pp3 = np.asarray(obj['momentum'])
   q3  = np.asarray(obj['q'])
   p3  = pp3 - q3
   sync_mom2 = int(np.sum(pp3**2))
   src_mom2  = int(np.sum(p3**2))
   q3_mom2   = int(np.sum(q3**2))
   if return_q_mom2:
     return pp3, q3, p3, sync_mom2, src_mom2, q3_mom2
   else:
     return pp3, q3, p3, sync_mom2, src_mom2

# deprecated because for cleanliness we use dispersion relation to determine this
#def four_kinematics_transform( sync_mass, pp3, src_mass, p3, nx, ny, nz ):
#  mom_conversion_factor = 2*np.pi*np.reciprocal(map(np.float,[nx,ny,nz]))
#  pp = np.asarray([1.0j*sync_mass]+(pp3*mom_conversion_factor).tolist())
#  p  = np.asarray([1.0j*src_mass] +(p3 *mom_conversion_factor).tolist())
#  q  = pp - p
#  return pp, q, p

def ff_dispersion_relation( mass, mom3, nx, ny, nz ):
  mom_conversion_factor = 2*np.pi*np.reciprocal(map(np.float,[nx,ny,nz]))
  euclidean_mom = euclidean_dispersion_relation( [mass]+(mom3*mom_conversion_factor).tolist() )
  return euclidean_mom
