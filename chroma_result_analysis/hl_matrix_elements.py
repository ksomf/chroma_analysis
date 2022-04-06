from hl import *
from itertools import groupby

I  = np.diag( [1.,1.,1.,1.] )
g0 = np.diag( [1.+0j,1.,-1.,-1.] )
g1 = np.fliplr( g0 )
g2 = np.fliplr(np.diag( [-1.j, 1.j, 1.j, -1.j ] ))
g3 = np.asarray( [[0.j,0,1.,0],[0,0,0,-1.],[-1.,0,0,0],[0,1.,0,0]] )
g5 = np.asarray( [[0.j,0,1.,0],[0,0,0,1.],[1.,0,0,0],[0,1.,0,0]] )
gamma = [ g0, g1, g2, g3 ]

#def KinematicsForConstQ2( range ):
def slash( four_vec ):
  return four_vec[0]*g0 - four_vec[1]*g1 - four_vec[2]*g2 - four_vec[3]*g3

def sigmavec( mu, vec ):
  return np.dot(gamma[mu],slash(vec)) - np.dot(slash(vec),gamma[mu])

def dispersion_relation( four_vec ):
  res = np.zeros( four_vec.shape )
  res[0] = np.sqrt(np.sum(four_vec**2))
  res[1:] = four_vec[1:]
  return res

def trace( mat ):
  return np.sum(np.diag(mat))

def TensorFF( kinematics ):
  m     = kinematics['m']
  mu    = kinematics['mu']
  nu    = kinematics['nu']
  pp    = dispersion_relation(np.asarray([m]+kinematics['pp']))
  p     = dispersion_relation(np.asarray([m]+kinematics['p']))
  Gamma = kinematics['Gamma']
  P     = (pp+p)/2.0
  q     = (pp-p)
  q2    = q[0]**2 - np.sum(q[1:4]**2)
  FFs = [ (gamma[mu]*P[nu]+gamma[nu]*P[mu])/2.0
        , (sigmavec(mu,q)*P[nu]+sigmavec(nu,q)*P[mu])/(4.0*m)
        , q[mu]*q[nu]/m ]
  def _project_and_trace( inner_ff ):
    return trace(reduce( np.dot, [Gamma, slash(pp)+m*I, inner_ff, slash(p)+m*I] ))
  return np.asarray(map( _project_and_trace, FFs ))

def GluonIrrep1( kinematics ):
  kinematics['mu'], kinematics['nu'] = 0, 0
  temporal_component = TensorFF( kinematics )
  kinematics['mu'], kinematics['nu'] = 1, 1
  spatial_component1 = TensorFF( kinematics )
  kinematics['mu'], kinematics['nu'] = 2, 2
  spatial_component2 = TensorFF( kinematics )
  kinematics['mu'], kinematics['nu'] = 3, 3
  spatial_component3 = TensorFF( kinematics )
  return temporal_component + ( spatial_component1 + spatial_component2 + spatial_component3 ) / 3.0

def GetFFCoefficients( cor ):
  op    = cor['op']
  pp    = cor['momentum']
  q     = cor['q']
  p     = (np.asarray(pp)-np.asarray(q)).tolist()
  if cor['spin'] == 'up':
    Gamma = np.dot((I+gamma[0])/2.0,(I+1.0j*np.dot(g3,g5)))
  elif cor['spin'] == 'down':
    Gamma = np.dot((I+gamma[0])/2.0,(I-1.0j*np.dot(g3,g5)))
  kinematics = { 'm':0.12334, 'pp':pp, 'p':p, 'Gamma':Gamma }

  if op == 'wilsonB2mE2':
    return str(map(str,GluonIrrep1( kinematics )))


mass=0.987
projector=np.dot((I+gamma[0])/2.0,(I+1.0j*np.dot(g3,g5)))
#projector=(I+gamma[0])/2.0
#print 1, GluonIrrep1( { 'm':mass, 'pp':[0,0,1] , 'p':[1,0,1], 'Gamma':projector } )
#print 2, GluonIrrep1( { 'm':mass, 'pp':[0,0,1] , 'p':[0,1,1], 'Gamma':projector } )
#print 3, GluonIrrep1( { 'm':mass, 'pp':[0,1,0] , 'p':[0,1,1], 'Gamma':projector } )
#print 4, GluonIrrep1( { 'm':mass, 'pp':[0,1,0] , 'p':[1,1,0], 'Gamma':projector } )
#print 5, GluonIrrep1( { 'm':mass, 'pp':[1,0,0] , 'p':[1,0,1], 'Gamma':projector } )
#print 6, GluonIrrep1( { 'm':mass, 'pp':[1,0,0] , 'p':[1,1,0], 'Gamma':projector } )
#print 7, GluonIrrep1( { 'm':mass, 'pp':[1,0,0] , 'p':[0,1,1], 'Gamma':projector } )
#print 8, GluonIrrep1( { 'm':mass, 'pp':[0,1,0] , 'p':[1,0,1], 'Gamma':projector } )
#print 9, GluonIrrep1( { 'm':mass, 'pp':[0,0,1] , 'p':[1,1,0], 'Gamma':projector } )
