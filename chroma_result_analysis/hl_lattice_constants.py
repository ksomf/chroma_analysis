hbarc = .197 #GeV nm
renorms = {
  5.50 : {
    'g3' : 0.857
  , 'ig3g5' : 0.8728
  }
  , 5.40 : {
    'g3' : 0.8327
  }
  , 5.65 : {
    'g3' : 0.8403
  }
}

spacings = {
  5.40 : 0.0818
, 5.50 : 0.074
, 5.65 : 0.068
}

def GenerateLatticeBlock( beta, kappa_l, kappa_s, dim ):
  return 'b{:1d}p{:2d}kp{}kp{}_{}x{}'.format( int(beta), int(100*(beta - int(beta))), int(kappa_l*10**6), int(kappa_s*10**6), dim[0], dim[3] )

def GenerateSimulationBlock( q, op ):
  return 'q{}{}{}_{}'.format( q[0], q[1], q[2], op )
