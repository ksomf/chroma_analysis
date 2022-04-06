import numpy as np 

from hl                   import *
from hl_naming_convention import chroma_str
from hl_plot              import _ChromaPlot, MultiPlot
from hl_load_json         import StoreFit
from hl_matrix_elements   import *

from functools import partial
from itertools import groupby

import ast

bookeeping_keys = [ 'file', 'keys', 'extraction', 'data', 'used_keys', 'used', 'key_map', 'test', 'weight', 'wgt', 'boot_file', 'three_file', 'threeboot_file', 'R_file' ]
data_keys       = [ 'raw', 'cor', 'culled', 'boot', 'culled_weight', 'cor_weight', 'threeboot', 'boot2', 'eff', 'R', 'path', 'cor_file', 'sync_cor_file', 'boot_threepoint_file', 'current_file', 'boot_sync_cor_file', 'boot_cor_file', 'boot_src_cor_file', 'src_cor_file', 'sync_boot', 'src_boot', 'threept_boot' ]
meta_keys       = bookeeping_keys + data_keys
id_keys         = [ 'hadron', 'prefix', 'weight_prefix', 'kappas', 'dim', 'mom', 'prop_modification_type' ]
lam_id_keys     = id_keys + [ 'feynman_hellmann' ]
ff_keys         = exclude( lam_id_keys, ['mom'] ) + [ 'momentum' ]


def DictDiff( d1, d2, ignore_keys ):
  if isinstance(d1,list) and isinstance(d2,list):
    return filter(Id,[ DictDiff( d1a, d2a, ignore_keys ) for d1a,d2a in zip(d1,d2) ])
  if isinstance(d1,dict) and isinstance(d2,dict):
    d1k, d2k = exclude(d1.keys(),ignore_keys), exclude(d2.keys(),ignore_keys)
    d1_differing_keys   = exclude(d1k,d2k) 
    d2_differing_keys   = exclude(d2k,d1k)
    common_keys         = overlap(d1k,d2k)
    differing_values    = [ k for k in common_keys if d1[k] != d2[k] ]
    
    res = { k:(d1[k],d2[k]) for k in differing_values }
    return res

def FeynmanEqualExcept( cor1, cor2 ):
  d = DictDiff( cor1.get('feynman_hellmann',[]), cor2.get('feynman_hellmann',[]), ['lambda'] )
  return not d and EqualExcept( meta_keys, cor1, cor2 )

def cor_id( cor, keys=id_keys, **kwargs ):
  return '_'.join( [ chroma_str(k,cor[k]) for k in keys if k in cor ] + kwargs.values() )
def cor_id_dict( cor, keys=lam_id_keys, extra_args=[], **kwargs ):
  return dict([ (k,v) for k,v in cor.items() if k in keys ] + kwargs.items() + [('file_prefix',cor_id(cor,keys))] + extra_args )


get_effmass  = lambda cor, shift=1, factor=1.0, axis=1, **kwargs: boot_object( np.log( np.abs( np.divide( cor, np.roll(cor,-shift,axis=axis) ) ) )/shift/factor, shift=shift, time=kwargs.get('absshift',0.0) + shift/2.0 + np.arange( cor.shape[-1] ), **kwargs )
phase_extraction = lambda cor, shift=1: np.divide( cor, np.roll(cor,-shift,axis=1) )
def get_plot_args(args):
  res = { k:ast.literal_eval(args[k]) for k in overlap(args.keys(),['displace','ylims']) }
  return res


def FitCorrelatorEffmass( id, cor, eff, save_path, fit_range_str, factor=1.0, allow_cosh_fit=True, save_append='', store=True ):
  fit_range = ast.literal_eval(fit_range_str)
  lower,upper = fit_range
  #print cor.keys(), eff.keys(), id.keys()
  if allow_cosh_fit and id['hadron'] == 'g5g5':
    fit = BootFitCosh1D( cor['time'][lower:upper+1]
                       , cor['boot'][:,lower:upper+1]
                       , cor['std'][lower:upper+1]
                       , offset = id['dim'][-1]/2.0, range=fit_range, factor=factor, type='cosh' )
    eff['fit'] = dict(fit)
    eff['fit']['type'] = 'eff_cosh'
  else:
    fit = BootFitExponential1D( cor['time'][lower:upper+1]
                              , cor['boot'][:,lower:upper+1]
                              , cor['std'][lower:upper+1], range=fit_range, factor=factor, type='poly', coefs=[0,1] )
    eff['fit'] = boot_object( -fit['boot'][:,1:], range=fit_range, coefs=[0], type='poly', chi2dof=fit['chi2dof'] )
  cor['fit'] = fit
  if store:
    StoreFit( eff['fit'], save_path, save_append=save_append, **id )
  return cor, eff

def FitThreePointCorrelator( id, cor, save_path, fit_range_str, flat=False, factor=1.0 ):
  lower,upper = fit_range = ast.literal_eval(fit_range_str)
  if flat:
    coefs = [ 0 ]
  else:
    coefs = [ 0, 1 ]
  fit = BootFitPoly1D( cor['time'][lower:upper+1]
                     , cor['boot'][:,lower:upper+1]
                     , cor['std'][lower:upper+1], coefs, range=fit_range, factor=factor, type='poly' )
  fit['coefs'] = coefs
  cor['fit'] = fit
  StoreFit( cor['fit'], save_path, **id )
  return cor

def FitThreePointTrapezoid( id, cor, save_path, fit_trapezoid, factor=1.0, data_type='real' ):
  nt = id['dim'][-1]
  min_tau, max_tau = map( int, MinMax( map( Fst, fit_trapezoid ) ) )
  min_t  , max_t   = map( int, MinMax( map( Snd, fit_trapezoid ) ) )
  points = np.asarray([ (j,i) for i in range(nt) for j in range(nt) ])
  contained_points = PointsInside( fit_trapezoid, points )
  correlator_boots  = []
  correlator_points = []
  correlator_stds   = []
  for it,itau in contained_points:
    correlator_boots.append(  cor['boot'][:,it,itau] )
    correlator_stds.append(   cor['std'][it,itau]    )
    correlator_points.append( it*itau )
  correlator_boots  = np.transpose( np.asarray( correlator_boots   ) )
  correlator_points = np.asarray( correlator_points )
  correlator_stds   = np.asarray( correlator_stds   )

  fit = BootFitPoly1D( correlator_points, correlator_boots, correlator_stds, [0], factor=factor, trapezoid=fit_trapezoid, data_type=data_type )
  cor['fit'] = fit
  print id
  StoreFit( cor['fit'], save_path, **id )
  return cor, [min_tau-1,max_tau+1], [min_t-1,max_t+1]

##-- FILTER/PARTITION SETS --##
_lambda_getter = lambda c: [ d['lambda'] for d in c.get('feynman_hellmann',[]) ]
_lambda_len_getter = lambda c: len(c.get('feynman_hellmann',[]))
def _hlRatioSplitByOpNum( cors, list ):
  return ( filter( lambda c: _lambda_len_getter(c)==num, cors ) for num in list )

def hlRatioFilterSpinIndependant( cors ):
  res = Partition( partial(EqualExcept, meta_keys + [ 'spin' ]), cors )
  return res

def hlRatioFilterSpinAndMomentumDirIndependant( cors ):
  res = Partition( partial(EqualExcept, meta_keys + [ 'pp', 'momentum', 'spin' ]), cors )
  return res

def hlRatioFilterLambdaSameOpSets( cors ):
  one_cors, zero_cors = _hlRatioSplitByOpNum( cors )
  res = Partition( partial(EqualExcept, meta_keys + [ 'pp', 'momentum', 'spin', 'lambdas' ]), one_cors )
  return res

def hlRatioCombiner( cors, type='tavg', data='real'):
  if   type == 'tavg':
    slice    = np.s_[:]
    avg_dims = (0,2)
  elif type == 'tfwd':
    slice    = np.s_[0]
    avg_dims = (0)
  elif type == 'tbck':
    slice    = np.s_[1]
    avg_dims = (0)
  if data == 'real':
    return np.mean( np.asarray([ c['sync_boot'][:,slice,:].real for c in cors ]), avg_dims )
  elif data == 'imag':
    return np.mean( np.asarray([ c['sync_boot'][:,slice,:].imag for c in cors ]), avg_dims )
  elif data == 'cmpl':
    return np.mean( np.asarray([ c['sync_boot'][:,slice,:]      for c in cors ]), avg_dims )

tuplise = lambda xs: tuple(xs) if isinstance(xs,list) else xs
def MassRatio( cors, args, save_path ):
  sync2_sets = hlRatioFilterSpinAndMomentumDirIndependant( cors )
  sync_sets  = hlRatioFilterSpinIndependant( cors )

  sync2_fits = {}
  for sync2_set in sync2_sets:
    #print sync2_set[0]['mom']
    #print exclude( sync2_set[0].keys(), meta_keys )
    nt = sync2_set[0]['dim'][-1]
    rat        = hlRatioCombiner( sync2_set )
    correlator = boot_object( rat, time=np.arange(nt) )
    effmass    = get_effmass( correlator['boot'] )
    id         = cor_id(sync2_set[0],keys=lam_id_keys)
    id_dict    = { k:sync2_set[0][k] for k in ['hadron','dim','kappas','mom'] }
    id_dict['fit'] = 'effmass'
    id_tuple   = tuple([ tuplise(id_dict[k])  for k in ['hadron','dim','kappas','mom'] ])
    if 'fit' in args:
      correlator, effmass = FitCorrelatorEffmass( id_dict, correlator, effmass, save_path, args['fit'] )
      sync2_fits[id_tuple] = effmass['fit']
      params = { 'xstride':12 }
      _ChromaPlot( [correlator], save=save_path+'png/cor_'+id+'.png', yscale='log', **params )
      params['xlims'] = literal_eval( args.get('xlims', '[0.0,'+str(nt)+']' ) )
      if 'ylims' in args:
        params['ylims'] = ast.literal_eval(args['ylims'])
      _ChromaPlot( [effmass], save=save_path+'png/eff_'+id+'.png', **params )

  sync_fits = {}
  for sync_set in sync_sets:
    nt = sync_set[0]['dim'][-1]
    rat        = hlRatioCombiner( sync_set )
    correlator = boot_object( rat, time=np.arange(nt) )
    effmass    = get_effmass( correlator['boot'] )
    id         = cor_id(sync_set[0],keys=lam_id_keys+['pp'])
    id_dict    = { k:sync_set[0][k] for k in ['hadron','dim','kappas','pp'] }
    id_dict['fit'] = 'effmass'
    id_tuple   = tuple([ tuplise(id_dict[k]) for k in ['hadron','dim','kappas','pp'] ])
    if 'fit' in args:
      correlator, effmass = FitCorrelatorEffmass( id_dict, correlator, effmass, save_path, args['fit'] )
      sync_fits[id_tuple] = effmass['fit']
      params = { 'xstride':12 }
      _ChromaPlot( [correlator], save=save_path+'png/cor_'+id+'.png', yscale='log', **params )
      params['xlims'] = literal_eval( args.get('xlims', '[0.0,'+str(nt)+']' ) )
      if 'ylims' in args:
        params['ylims'] = args['ylims']
      _ChromaPlot( [effmass], save=save_path+'png/eff_'+id+'.png', **params )

def SimpleRatio( cors, args, save_path ):
  nt = cors[0]['dim'][-1]
  zero_cors, other_cors  = _hlRatioSplitByOpNum( cors, [0,1] )
  print len(cors), len(other_cors), len(zero_cors)
  vals  = hlRatioCombiner( other_cors, type='tfwd', data='imag' )
  zeros = hlRatioCombiner( zero_cors , type='tfwd' )
  correlator = boot_object( vals / zeros, time=np.arange(nt) )
  effmass = get_effmass( correlator['boot'] )
  params = { 'xstride':12 }
  _ChromaPlot( [effmass], save=save_path+'png/test.png', **params )

#def SimpleOddRatio( cors, args, save_path ):
#  print 'Performing Simple Odd Ratio'
#  cors = filter( lambda cor:     cor['lambdas']!=[0.0], cors )
#  cors = filter( lambda cor: len(cor['lambdas'])==1   , cors )
#  one_cors = _hlRatioSplitByOpNum( cors, [1] )
#  if len(one_cors):
#    cor_sets = Partition( partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin']), one_cors )
#    for cor_set in cor_sets: #Should be different operators here
#      split_cor_set  = Partition( lambda cor1, cor2: abs(cor1['lambdas'][0]) == abs(cor2['lambdas'][0]), cor_set )
#      for lambda_cor_set in split_cor_set:
#        pos_lam_cors, neg_lam_cors = SplitFilter( lambda cor: cor['lambdas'][0] > 0.0, lambda_cor_set )
#        id = cor_id( pos_lam_cors[0], keys=lam_id_keys )
#        id_dict = cor_id_dict( pos_lam_cors[0], fit='simple_odd' )
#        op  = pos_lam_cors[0]['op']
#        lam = pos_lam_cors[0]['lambdas'][0]
#
#        numerator   = hlRatioCombiner( pos_lam_cors )
#        denominator = hlRatioCombiner( neg_lam_cors )
#        rat = numerator / denominator
#        #old_pos = pos_lam_cors[1]['boot'][:,1,:].real
#        #old_neg = neg_lam_cors[1]['boot'][:,1,:].real
#        #rat = old_pos / old_neg
#        correlator = boot_object( rat, time=range(rat.shape[-1]) )
#        for s in [ 1, 2, 3, 4 ]:
#          effmass    = get_effmass( correlator['boot'], factor=2*lam, shift=s )
#          if 'fit' in args:
#            correlator, effmass = FitCorrelatorEffmass( id_dict, correlator, effmass, save_path, args['fit'], factor=2*lam )
#          _ChromaPlot( [effmass], save=save_path+'png/simple_odd_ratio_effmass_'+id+'shift'+str(s)+'.png', xstride=12, xlims=ast.literal_eval(args.get('xlims','[0.0,24.0]')) )
#
#        #print np.mean( np.asarray([ np.mean( cor['boot'].real, 1 ) for cor in pos_lam_cors ]), 0 ).shape
#        #print np.asarray([ np.mean( cor['boot'].real, 1 ) for cor in neg_lam_cors ]).shape

def SimpleEvenRatio( cors, args, save_path ):
  print 'Performing Simple Even Ratio'
  #print len(cors)
  zero_cors, one_cors = _hlRatioSplitByOpNum( cors, [0,1] )
  nt = cors[0]['dim'][-1]

  extra_plot_args = get_plot_args( args ) 
  #print exclude( cors[0].keys(), data_keys+bookeeping_keys )
  #print exclude( zero_cors[0].keys(), data_keys+bookeeping_keys )
  if len(one_cors):
    zero_set, one_set = _hlRatioSplitByOpNum( cors, [0,1] )
    print [ abs(_lambda_getter(cor)[0]) for cor in one_set ]
    split_cor_set  = Partition( lambda cor1, cor2: abs(_lambda_getter(cor1)[0]) == abs(_lambda_getter(cor2)[0]), one_set )
    effmasses = []
    diagnostic_effmasses = []
    for i,lambda_cor_set in enumerate(split_cor_set):
      #print len(pos_lam_cors), len(neg_lam_cors), len(zero_cors)
      keys = exclude(lam_id_keys+['momentum'],['mom'])
      if not args.get('simplified_ratio',False):
        pos_lam_cors, neg_lam_cors  = SplitFilter( lambda cor: cor['feynman_hellmann'][0]['lambda']>0.0, lambda_cor_set )
        id = cor_id( pos_lam_cors[0], keys=keys )
        id_dict = cor_id_dict( pos_lam_cors[0], fit='simple_even', keys=keys, ff=args['ff'] )
        lam = _lambda_getter(pos_lam_cors[0])[0]
        assert(len(pos_lam_cors))
        assert(len(neg_lam_cors))
        assert(len(zero_cors))
        print len(pos_lam_cors)
        pos = hlRatioCombiner( pos_lam_cors )
        neg = hlRatioCombiner( neg_lam_cors )
        den = hlRatioCombiner( zero_cors )
        rat = (pos * neg) / (den*den)
        diagrats = []
        if not i:
          diagrats = [ (den*den,0.0) ]
        diagrats = diagrats + [ (pos*neg,lam) ]
        factor=lam**2 #/2 for even ratio, *2 for taylor 2! term
        #rat = pos / den
        #rat = pos
      else:
        id = cor_id( lambda_cor_set[0], keys=keys )
        id_dict = cor_id_dict( lambda_cor_set[0], fit='simple_even', keys=keys, ff=args['ff'] )
        lam = _lambda_getter(lambda_cor_set[0])[0]
        factor=1./2.*lam**2
        assert(len(lambda_cor_set))
        assert(len(zero_cors))
        num = hlRatioCombiner( lambda_cor_set )
        den = hlRatioCombiner( zero_cors )
        rat = num / den
        if not i:
          diagrats = [ (den,0.0) ]
        diagrats = diagrats + [ (num,lam) ]

      #print np.asarray([ c['boot'].real for c in pos_lam_cors ]).shape
      #print np.asarray([ c['boot'].real for c in zero_cors ]).shape
      correlator = boot_object( rat, time=range(nt) )
      diag_correlators = [ boot_object(r,time=range(nt),l=l) for r,l in diagrats ]
      for s in [1]: #[ 1, 2, 3, 4 ]:
        effmass = get_effmass( correlator['boot'], factor=factor, shift=s )
        diag_effmasses = [ get_effmass( c['boot'], shift=s ) for c in diag_correlators ]
        if 'fit' in args:
          correlator, effmass = FitCorrelatorEffmass( id_dict, correlator, effmass, save_path, args['fit'], factor=factor, allow_cosh_fit=False, save_append=args.get('save_append','') )
          diag_correlators, diag_effmasses = map(list,zip(*[ FitCorrelatorEffmass( id_dict, c, e, '', args['fit'], store=False ) for c,e in zip(diag_correlators,diag_effmasses) ]))
      effmasses.append(effmass)
      diagnostic_effmasses = diagnostic_effmasses + diag_effmasses
    _ChromaPlot( effmasses           , save=save_path+'png/simple_even_ratio_effmass_'     +id+args.get('save_append','')+'shift'+str(s)+'.png', xstride=12, xlims=ast.literal_eval(args.get('xlims','[0.0,24.0]')), **extra_plot_args )
    _ChromaPlot( diagnostic_effmasses, save=save_path+'png/simple_even_ratio_effmass_diag_'+id+args.get('save_append','')+'shift'+str(s)+'.png', xstride=12, xlims=[4.0,24.0], displace=0.5, legend=True )



def _SpinOddRatio( pos_lam_cors, neg_lam_cors, zero_up_cors, zero_down_cors ):
  num1 = hlRatioCombiner( pos_lam_cors  , type='tfwd' )
  num2 = hlRatioCombiner( zero_down_cors, type='tfwd' )
  num3 = hlRatioCombiner( neg_lam_cors  , type='tbck' )
  num4 = hlRatioCombiner( zero_up_cors  , type='tbck' )
  den1 = hlRatioCombiner( zero_up_cors  , type='tfwd' )
  den2 = hlRatioCombiner( neg_lam_cors  , type='tfwd' )
  den3 = hlRatioCombiner( zero_down_cors, type='tbck' )
  den4 = hlRatioCombiner( pos_lam_cors  , type='tbck' )
  rat = (num1/den1) * (num2/den2) * (num3/den3) * (num4/den4)
  return rat

def _SpinEvenRatio( pos_lam_cors, neg_lam_cors, zero_up_cors, zero_down_cors ):
  num1 = hlRatioCombiner( pos_lam_cors   )
  num2 = hlRatioCombiner( neg_lam_cors   )
  den1 = hlRatioCombiner( zero_up_cors   )
  den2 = hlRatioCombiner( zero_down_cors )
  rat = (num1/den1) * (num2/den2)
  return rat

def _SpinSimpleRatio( pos_lam_cors, neg_lam_cors, zero_up_cors, zero_down_cors ):
  set1 = hlRatioCombiner( pos_lam_cors, type='tfwd', data='imag' )
  set2 = hlRatioCombiner( zero_up_cors, type='tfwd', data='iamg' )
  rat = set1 - set2
  return rat

def extract_plot_params( params, args ):
  for k in overlap( args.keys(), [ 'ylims', 'xlabel', 'ylabel', 'fs' ] ):
    params[k] = args[k]

def _SpinDependentRatio( cors, args, save_path, ratio, fitstr, factor_pow=0, scale_factor=1 ): 
  #print [ c['type'] for c in cors ]
  #print [ c.get('lambdas',[]) for c in cors ]
  zero_cors, one_cors = _hlRatioSplitByOpNum( cors, [0,1] )
  if len(one_cors) and len(zero_cors):
    comparison            = partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin'])
    zero_match_comparison = partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin','op','type','q'])
    #print [ exclude( c.keys(), data_keys+bookeeping_keys+['lambdas','momentum','spin','op','type'] ) for c in one_cors ]
    #print [ exclude( c.keys(), data_keys+bookeeping_keys+['lambdas','momentum','spin','op','type'] ) for c in zero_cors ]
    cor_sets = Partition( comparison, one_cors, aux_match=zero_cors, aux_match_comparison=zero_match_comparison )
    for cor_set in cor_sets:
      #print len(cor_set)
      cor_set, zero_cor_set = SplitFilter( lambda cor: cor.get('lambdas',[])!=[], cor_set )
      #print len(cor_set), len(zero_cor_set)
      zero_up_cors, zero_down_cors = SplitFilter( lambda cor: cor['spin']=='up', zero_cor_set )
      split_cor_set = Partition( lambda cor1, cor2: abs(cor1['lambdas'][0]) == abs(cor2['lambdas'][0]), cor_set )
      nt = cor_set[0]['dim'][-1]
      for s in [ 1, 2, 3, 4 ]:
        res=[]
        for lambda_cor_set in split_cor_set:
          pos_lam_cors, neg_lam_cors = SplitFilter( lambda cor: cor['spin'] == 'up', lambda_cor_set )
          #print len(pos_lam_cors), len(neg_lam_cors), len(zero_up_cors), len(zero_down_cors)
          #for c in pos_lam_cors:
          #  print c['lambdas'], c['spin']
          id = cor_id( pos_lam_cors[0], keys=lam_id_keys )
          id_dict = cor_id_dict( pos_lam_cors[0], fit=fitstr )
          op  = pos_lam_cors[0]['op']
          lam = pos_lam_cors[0]['lambdas'][0]
          rat = ratio( pos_lam_cors, neg_lam_cors, zero_up_cors, zero_down_cors )
          correlator = boot_object( rat, time=np.asarray(range(rat.shape[-1])) )
          effmass_scale_factor=scale_factor*lam**factor_pow / factorial(factor_pow)
          effmass    = get_effmass( correlator['boot'], factor=effmass_scale_factor, shift=s )
          if 'fit' in args:
            correlator, effmass = FitCorrelatorEffmass( id_dict, correlator, effmass, save_path, args['fit'], factor=effmass_scale_factor )
          res.append( {'id':id,'cor':correlator,'eff':effmass,'idd':id_dict} )

        params = {
          'xlims' : ast.literal_eval(args.get('xlims','[0.0,'+str(nt/2)+']'))
        }
        extract_plot_params( params, args )

        cors = [ r['cor'] for r in res ]
        effs = [ r['eff'] for r in res ]
        _ChromaPlot( effs, save=save_path+'png/'+fitstr+'_compare_effmass_'+res[0]['id']+'shift'+str(s)+'.png', displace=True, **params )
        #_ChromaPlot( effs, save=save_path+'png/'+fitstr+'_compare_effmass_'+res[0]['id']+'shift'+str(s)+'.png', xstride=12, displace=True, **params )
        #_ChromaPlot( cors, save=save_path+'png/'+fitstr+'_compare_correlators_'+res[0]['id']+'shift'+str(s)+'.png', xstride=12, xlims=xlims, displace=True )
        for r,eff in zip(res,effs):
          id = r['id']
          _ChromaPlot( [eff], save=save_path+'png/'+fitstr+'_effmass_'+id+'shift'+str(s)+'.png', **params )
          #_ChromaPlot( [eff], save=save_path+'png/'+fitstr+'_effmass_'+id+'shift'+str(s)+'.png', xstride=12, **params )

def SpinDependentOddRatio( cors, args, save_path ):
  print 'Performing Spin Dependent Odd Ratio'
  _SpinDependentRatio( cors, args, save_path, _SpinOddRatio, 'spin_odd_ratio', factor_pow=1, scale_factor=4 )

def SpinDependentEvenRatio( cors, args, save_path ):
  print 'Performing Spin Dependent Even Ratio'
  _SpinDependentRatio( cors, args, save_path, _SpinEvenRatio, 'spin_even_ratio', factor_pow=2, scale_factor=4 )

def SpinDependentSimpleRatio( cors, args, save_path ):
  print 'Performing simple spin dependent ratio'
  _SpinDependentRatio( cors, args, save_path, _SpinSimpleRatio, 'spin_simple_ratio', factor_pow=0, scale_factor=1 )
    
def SpinDependentPartialHilbertRatio( cors, args, save_path ):
  fitstr = 'spin_partial_hilbert'
  zero_cors, one_cors = _hlRatioSplitByOpNum( cors, [0,1] )
  if len(one_cors) and len(zero_cors):
    comparison            = partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin'])
    zero_match_comparison = partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin','op','type'])
    cor_sets = Partition( comparison, one_cors, aux_match=zero_cors, aux_match_comparison=zero_match_comparison )
    for cor_set in cor_sets:
      zero_cor_set, one_cor_set = _hlRatioSplitByOpNum( cor_set, [0,1] )
      zero_up_cors, zero_down_cors = SplitFilter( lambda cor: cor['spin']=='up', zero_cor_set )
      split_cor_set = Partition( lambda cor1, cor2: abs(cor1['lambdas'][0]) == abs(cor2['lambdas'][0]), one_cor_set )
      nt = cor_set[0]['dim'][-1]
      for s in [ 1, 2, 3, 4 ]:
        res=[]
        for lambda_cor_set in split_cor_set:
          pos_lam_cors, neg_lam_cors = SplitFilter( lambda cor: cor['lambdas'][0] > 0.0, lambda_cor_set )
          id = cor_id( pos_lam_cors[0], keys=lam_id_keys )
          id_dict = cor_id_dict( pos_lam_cors[0], fit=fitstr )
          op  = pos_lam_cors[0]['op']
          lam = pos_lam_cors[0]['lambdas'][0]
          dE  = np.power( _SpinOddRatio(  pos_lam_cors, neg_lam_cors, zero_up_cors, zero_down_cors ), (1.0/4.0) )
          d2E = np.power( _SpinEvenRatio( pos_lam_cors, neg_lam_cors, zero_up_cors, zero_down_cors ), (2.0/4.0)/lam )
          ph  = np.zeros( d2E.shape )
          for i in range(nt):
            ph[:,i] = np.power( dE[:,i], np.log(dE[:,i])*(i+1) ) / d2E[:,i]
          correlator = boot_object( ph, time=range(ph.shape[-1]) )
          effmass    = get_effmass( correlator['boot'], factor=lam**1, shift=s )
          res.append( {'id':id,'cor':correlator,'eff':effmass,'idd':id_dict} )
        xlims=ast.literal_eval(args.get('xlims','[0.0,'+str(nt/2)+']'))
        effs = [ r['eff'] for r in res ]
        _ChromaPlot( effs, save=save_path+'png/'+fitstr+'_compare_effmass_'+res[0]['id']+'shift'+str(s)+'.png', xstride=12, xlims=xlims, displace=True )
        for r,eff in zip(res,effs):
          id = r['id']
          _ChromaPlot( [eff], save=save_path+'png/'+fitstr+'_effmass_'+id+'shift'+str(s)+'.png', xstride=12, xlims=xlims )

def IsCoreOpLambda( op_lams, cor ):
  feynman_hellmann_modifications = cor.get('feynman_hellmann',[])
  if len(feynman_hellmann_modifications) == len(op_lams):
    for mod in feynman_hellmann_modifications:
      if not (mod['op'],mod['lambda']) in op_lams:
        return False
    return True
  return False

def SpinIndependentOddOddRatio( cors, args, save_path ):
  print 'Spin Independent Odd Odd Ratio'
  #cors = filter( lambda cor: cor['spin'] == 'up', cors )
  zero_cors, one_cors, two_cors = _hlRatioSplitByOpNum( cors, [0,1,2] )
  extra_plot_args = get_plot_args( args )
  #print len(two_cors), len(one_cors), len(zero_cors)
  if( len(two_cors) ): #and len(one_cors) and len(zero_cors) ):
    comparison           = partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin'])
    one_match_comparison = partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin','op','type'])
    cor_sets = Partition( comparison, two_cors, aux_match=one_cors, aux_match_comparison=one_match_comparison )
    for cor_set in cor_sets:
      two_lam_cors = filter( lambda cor: len(cor['lambdas'])==2, cor_set )
      one_lam_cors = filter( lambda cor: len(cor['lambdas'])==1, cor_set )
      op_lam_list = sorted(list(set([ (o,abs(l)) for cor in cor_set for o,l in zip(cor['op'],cor['lambdas']) ])))
      pair1,pair2 = SplitFilter( lambda x: x[0] == op_lam_list[0][0], op_lam_list )
      for s in [ 1, 2, 3, 4 ]:
        effmasses = []
        for (o1,l1) in pair1:
          for (o2,l2) in pair2:
            #print map(len,[pos_pos_lam_cors,pos_neg_lam_cors,neg_pos_lam_cors,neg_neg_lam_cors])

            id = cor_id( pos_pos_lam_cors[0], keys=exclude(lam_id_keys+['momentum'],['lambdas']) )
            id_dict = cor_id_dict( pos_pos_lam_cors[0], fit='odd_odd' )
            op  = pos_pos_lam_cors[0]['op']
            lam = pos_pos_lam_cors[0]['lambdas'][0]

            num1 = hlRatioCombiner( pos_pos_lam_cors )
            num2 = hlRatioCombiner( neg_neg_lam_cors )
            den1 = hlRatioCombiner( neg_pos_lam_cors )
            den2 = hlRatioCombiner( pos_neg_lam_cors )
            rat = (num1/den1) * (num2/den2)
            correlator = boot_object( rat, time=range(rat.shape[-1]) )
            effmass = get_effmass( correlator['boot'], factor=l1*l2, shift=s )
            if 'fit' in args:
              correlator, effmass = FitCorrelatorEffmass( id_dict, correlator, effmass, save_path, args['fit'], factor=l1*l2 )
            effmasses.append(effmass)
        _ChromaPlot( effmasses, save=save_path+'png/spin_independent_odd_odd_ratio_effmass_'+id+'shift'+str(s)+'.png', xstride=12, xlims=ast.literal_eval(args.get('xlims','[0.0,24.0]')), **extra_plot_args )

def SpinIndependentPhaseOddOddRatio( cors, args, save_path ):
  print 'Spin Independent Phase Odd Odd Ratio'
  #cors = filter( lambda cor: cor['spin'] == 'up', cors )
  #print [ len(cor['lambdas']) for cor in cors ]
  #print [ cor['op'] for cor in cors ]
  #print [ cor['lambdas'] for cor in cors ]
  zero_cors, one_cors, two_cors = _hlRatioSplitByOpNum( cors, [0,1,2] )
  extra_plot_args = get_plot_args( args )
  #print len(two_cors), len(one_cors), len(zero_cors)
  if( len(two_cors) ): #and len(one_cors) and len(zero_cors) ):
    comparison           = partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin'])
    one_match_comparison = partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin','op','type'])
    cor_sets = Partition( comparison, two_cors, aux_match=one_cors, aux_match_comparison=one_match_comparison )
    for cor_set in cor_sets:
      two_lam_cors = filter( lambda cor: len(cor['feynman_hellmann'])==2, cor_set )
      one_lam_cors = filter( lambda cor: len(cor['feynman_hellmann'])==1, cor_set )
      op_lam_list = sorted(list(set( [ (mod['op'],mod['lambda']) for cor in cor_set for mod in cor.get('feynman_hellmann',[]) ] )))
      #op_lam_list = sorted(list(set([ (o,abs(l)) for cor in cor_set for o,l in zip(cor['op'],cor['lambdas']) ])))
      pair1,pair2 = SplitFilter( lambda x: x[0] == op_lam_list[0][0], op_lam_list )
      for s in [4]: #[ 1, 2, 3, 4 ]:
        effmasses = []
        for (o1,l1) in pair1:
          for (o2,l2) in pair2:
            pos_pos_lam_cors = filter( partial(IsCoreOpLambda, [(o1, l1),(o2, l2)]), two_lam_cors )
            pos_neg_lam_cors = filter( partial(IsCoreOpLambda, [(o1, l1),(o2,-l2)]), two_lam_cors )
            neg_pos_lam_cors = filter( partial(IsCoreOpLambda, [(o1,-l1),(o2, l2)]), two_lam_cors )
            neg_neg_lam_cors = filter( partial(IsCoreOpLambda, [(o1,-l1),(o2,-l2)]), two_lam_cors )
            #print map(len,[pos_pos_lam_cors,pos_neg_lam_cors,neg_pos_lam_cors,neg_neg_lam_cors])

            id = cor_id( pos_pos_lam_cors[0], keys=exclude(lam_id_keys+['momentum'],['lambdas']) )
            id_dict = cor_id_dict( pos_pos_lam_cors[0], fit='odd_odd' )
            #op  = pos_pos_lam_cors[0]['op']
            #lam = pos_pos_lam_cors[0]['lambdas'][0]

            pos_pos_set = hlRatioCombiner( pos_pos_lam_cors, type='tfwd', data='cmpl' )
            neg_neg_set = hlRatioCombiner( neg_neg_lam_cors, type='tfwd', data='cmpl' )
            neg_pos_set = hlRatioCombiner( neg_pos_lam_cors, type='tfwd', data='cmpl' )
            pos_neg_set = hlRatioCombiner( pos_neg_lam_cors, type='tfwd', data='cmpl' )

            #num1 = pos_pos_set
            #num2 = neg_neg_set
            #den1 = neg_pos_set
            #den2 = pos_neg_set
            #rat = (num1/den1) * (num2/den2)
            #time_set = np.arange(rat.shape[-1])
            #phi = boot_object( 1.0j * np.log( phase_extraction( rat, shift=s ) ), time=time_set )

            pos_pos_effmass = phase_extraction( pos_pos_set, shift=s )
            neg_neg_effmass = phase_extraction( neg_neg_set, shift=s )
            neg_pos_effmass = phase_extraction( neg_pos_set, shift=s )
            pos_neg_effmass = phase_extraction( pos_neg_set, shift=s )
            num1 = pos_pos_effmass
            num2 = neg_neg_effmass
            den1 = neg_pos_effmass
            den2 = pos_neg_effmass
            rat = (num1/den1) * (num2/den2)
            time_set = np.arange(rat.shape[-1])
            phi = boot_object( 1.0j * np.log( rat ) / (l1*l2), time=time_set )
            effmasses.append(phi)
            #correlator = boot_object( rat, time=time_set )
            #effmass = get_effmass( correlator['boot'], factor=l1*l2, shift=s )
            #if 'fit' in args:
            #  correlator, effmass = FitCorrelatorEffmass( id_dict, correlator, effmass, save_path, args['fit'], factor=l1*l2 )
            #effmasses.append(effmass)
        _ChromaPlot( effmasses, save=save_path+'png/spin_independent_phase_odd_odd_ratio_effmass_'+id+'shift'+str(s)+'.png', xstride=12, xlims=ast.literal_eval(args.get('xlims','[0.0,24.0]')), **extra_plot_args )
    
def PolarisedPhaseOddOddRatio( cors, args, save_path ):
  print 'Polarised Phase Odd Odd Ratio'
  #cors = filter( lambda cor: cor['spin'] == 'up', cors )
  #print [ len(cor['lambdas']) for cor in cors ]
  #print [ cor['op'] for cor in cors ]
  #print [ cor['lambdas'] for cor in cors ]
  zero_cors, one_cors, two_cors = _hlRatioSplitByOpNum( cors, [0,1,2] )
  extra_plot_args = get_plot_args( args )
  #print len(two_cors), len(one_cors), len(zero_cors)
  if( len(two_cors) ): #and len(one_cors) and len(zero_cors) ):
    comparison           = partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin'])
    one_match_comparison = partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin','op','type'])
    cor_sets = Partition( comparison, two_cors, aux_match=one_cors, aux_match_comparison=one_match_comparison )
    for cor_set in cor_sets:
      two_lam_cors = filter( lambda cor: len(cor['lambdas'])==2, cor_set )
      one_lam_cors = filter( lambda cor: len(cor['lambdas'])==1, cor_set )
      op_lam_list = sorted(list(set([ (o,abs(l)) for cor in cor_set for o,l in zip(cor['op'],cor['lambdas']) ])))
      pair1,pair2 = SplitFilter( lambda x: x[0] == op_lam_list[0][0], op_lam_list )
      for pol in ['polx','poly','polz']:
        polled_two_lam_cors = filter( lambda cor: cor['spin'] == pol, two_lam_cors )
        for s in [4]: #[ 1, 2, 3, 4 ]:
          effmasses = []
          for (o1,l1) in pair1:
            for (o2,l2) in pair2:
              pos_pos_lam_cors = filter( partial(IsCoreOpLambda, [(o1, l1),(o2, l2)]), polled_two_lam_cors )
              pos_neg_lam_cors = filter( partial(IsCoreOpLambda, [(o1, l1),(o2,-l2)]), polled_two_lam_cors )
              neg_pos_lam_cors = filter( partial(IsCoreOpLambda, [(o1,-l1),(o2, l2)]), polled_two_lam_cors )
              neg_neg_lam_cors = filter( partial(IsCoreOpLambda, [(o1,-l1),(o2,-l2)]), polled_two_lam_cors )
              #print map(len,[pos_pos_lam_cors,pos_neg_lam_cors,neg_pos_lam_cors,neg_neg_lam_cors])

              id = cor_id( pos_pos_lam_cors[0], keys=exclude(lam_id_keys+['momentum'],['lambdas']) )
              id_dict = cor_id_dict( pos_pos_lam_cors[0], fit='odd_odd' )
              #op  = pos_pos_lam_cors[0]['op']
              #lam = pos_pos_lam_cors[0]['lambdas'][0]

              pos_pos_set = hlRatioCombiner( pos_pos_lam_cors, type='tfwd', data='cmpl' )
              neg_neg_set = hlRatioCombiner( neg_neg_lam_cors, type='tfwd', data='cmpl' )
              neg_pos_set = hlRatioCombiner( neg_pos_lam_cors, type='tfwd', data='cmpl' )
              pos_neg_set = hlRatioCombiner( pos_neg_lam_cors, type='tfwd', data='cmpl' )

              pos_pos_effmass = phase_extraction( pos_pos_set, shift=s )
              neg_neg_effmass = phase_extraction( neg_neg_set, shift=s )
              neg_pos_effmass = phase_extraction( neg_pos_set, shift=s )
              pos_neg_effmass = phase_extraction( pos_neg_set, shift=s )
              num1 = pos_pos_effmass
              num2 = neg_neg_effmass
              den1 = neg_pos_effmass
              den2 = pos_neg_effmass
              rat = (num1/den1) * (num2/den2)
              time_set = np.arange(rat.shape[-1])
              phi = boot_object( 1.0j * np.log( rat ) / (l1*l2), time=time_set )
              effmasses.append(phi)
              #correlator = boot_object( rat, time=time_set )
              #effmass = get_effmass( correlator['boot'], factor=l1*l2, shift=s )
              #if 'fit' in args:
              #  correlator, effmass = FitCorrelatorEffmass( id_dict, correlator, effmass, save_path, args['fit'], factor=l1*l2 )
              #effmasses.append(effmass)
          _ChromaPlot( effmasses, save=save_path+'png/polarised_phase_odd_odd_ratio_effmass_'+id+'_'+pol+'_shift'+str(s)+'.png', xstride=12, xlims=ast.literal_eval(args.get('xlims','[0.0,24.0]')), **extra_plot_args )
    
def SpinDependentOddOddRatio( cors, args, save_path ):
  print 'Spin Dependent Odd Odd Ratio'
  zero_cors, one_cors, two_cors = _hlRatioSplitByOpNum( cors, [0,1,2] )
  one_cors = filter( lambda cor: cor['type']                             , one_cors )
  one_cors = filter( lambda cor: not any([ 'g5' in o for o in cor['op']]), one_cors )
  if( len(two_cors) ):
    comparison           = partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin'])
    one_match_comparison = partial(EqualExcept, data_keys+bookeeping_keys+['lambdas','momentum','spin','op','type'])
    cor_sets = Partition( comparison, two_cors, aux_match=one_cors, aux_match_comparison=one_match_comparison )
    for cor_set in cor_sets:
      two_lam_cors = filter( lambda cor: len(cor['lambdas'])==2, cor_set )
      one_lam_cors = filter( lambda cor: len(cor['lambdas'])==1, cor_set )
      op_lam_list = sorted(list(set([ (o,abs(l)) for cor in cor_set for o,l in zip(cor['op'],cor['lambdas']) ])))
      pair1,pair2 = SplitFilter( lambda x: x[0] == op_lam_list[0][0], op_lam_list )
      for (o1,l1) in pair1:
        for (o2,l2) in pair2:
          if 'g5' in o1:
            print 'error on order'
            exit(1)
          pos_pos_lam_cors = filter( partial(IsCoreOpLambda, [(o1, l1),(o2, l2)]), two_lam_cors )
          pos_neg_lam_cors = filter( partial(IsCoreOpLambda, [(o1, l1),(o2,-l2)]), two_lam_cors )
          neg_pos_lam_cors = filter( partial(IsCoreOpLambda, [(o1,-l1),(o2, l2)]), two_lam_cors )
          neg_neg_lam_cors = filter( partial(IsCoreOpLambda, [(o1,-l1),(o2,-l2)]), two_lam_cors )
          pos_lam_cors = filter( partial(IsCoreOpLambda, [(o1, l1)]), one_lam_cors )
          neg_lam_cors = filter( partial(IsCoreOpLambda, [(o1,-l1)]), one_lam_cors )
          pos_up_lam_cors, pos_down_lam_cors = SplitFilter( lambda cor: cor['spin']=='up', pos_lam_cors )
          neg_up_lam_cors, neg_down_lam_cors = SplitFilter( lambda cor: cor['spin']=='up', neg_lam_cors )

          id = cor_id( pos_pos_lam_cors[0], keys=lam_id_keys )
          id_dict = cor_id_dict( pos_pos_lam_cors[0], fit='spin_odd_odd' )
          op  = pos_pos_lam_cors[0]['op']
          lam = pos_pos_lam_cors[0]['lambdas'][0]

          cor_len_list = map(len, [ pos_pos_lam_cors, pos_neg_lam_cors, neg_pos_lam_cors, neg_neg_lam_cors
                                  , pos_up_lam_cors, pos_down_lam_cors, neg_up_lam_cors, neg_down_lam_cors ])
          if not all(cor_len_list):
            print 'error'
            print cor_len_list
            exit( 1 )

          num1 = hlRatioCombiner( pos_pos_lam_cors , type='tfwd' )
          num2 = hlRatioCombiner( neg_neg_lam_cors , type='tfwd' )
          num3 = hlRatioCombiner( neg_up_lam_cors  , type='tfwd' )
          num4 = hlRatioCombiner( pos_down_lam_cors, type='tfwd' )
          den1 = hlRatioCombiner( pos_up_lam_cors  , type='tfwd' )
          den2 = hlRatioCombiner( neg_down_lam_cors, type='tfwd' )
          den3 = hlRatioCombiner( neg_pos_lam_cors , type='tfwd' )
          den4 = hlRatioCombiner( pos_neg_lam_cors , type='tfwd' )
          rat = (num1/den1) * (num2/den2) * (num3/den3) * (num4/den4)

          #num5 = np.mean( np.asarray([ c['boot'][:,1,:].real for c in pos_up_lam_cors   ]), 0 )
          #num6 = np.mean( np.asarray([ c['boot'][:,1,:].real for c in neg_down_lam_cors ]), 0 )
          #num7 = np.mean( np.asarray([ c['boot'][:,1,:].real for c in neg_neg_lam_cors  ]), 0 )
          #num8 = np.mean( np.asarray([ c['boot'][:,1,:].real for c in pos_pos_lam_cors  ]), 0 )
          #den5 = np.mean( np.asarray([ c['boot'][:,1,:].real for c in pos_neg_lam_cors  ]), 0 )
          #den6 = np.mean( np.asarray([ c['boot'][:,1,:].real for c in neg_pos_lam_cors  ]), 0 )
          #den7 = np.mean( np.asarray([ c['boot'][:,1,:].real for c in neg_up_lam_cors   ]), 0 )
          #den8 = np.mean( np.asarray([ c['boot'][:,1,:].real for c in pos_down_lam_cors ]), 0 )
          #rat  = (num1/den1) * (num2/den2) * (num3/den3) * (num4/den4) * (num5/den5) * (num6/den6) * (num7/den7) * (num8/den8)

          correlator = boot_object( rat, time=range(rat.shape[-1]) )
          for s in [ 1, 2, 3, 4 ]:
            effmass = get_effmass( correlator['boot'], factor=l1*l2, shift=s )
            if 'fit' in args:
              correlator, effmass = FitCorrelatorEffmass( id_dict, correlator, effmass, save_path, args['fit'], factor=l1*l2 )
            _ChromaPlot( [effmass], save=save_path+'png/spin_odd_odd_ratio_effmass_'+id+'shift'+str(s)+'.png', xstride=12, xlims=ast.literal_eval(args.get('xlims','[0.0,24.0]')) )
 
def ThreePointMass( cors, args, save_path ):
  print 'Perfroming Three Point Mass Plot'
  threept_sets = Partition( partial(EqualExcept, data_keys + bookeeping_keys + [ 'momentum', 'spin' ]), cors )
  for threept_set in threept_sets:
    nt = threept_set[0]['dim'][-1]
    delta = int(args.get( 'delta', 1 ))
    twopt   = np.asarray([ cor['boot']      for cor in threept_set ]).real
    threept = np.asarray([ cor['threeboot'] for cor in threept_set ])
    id_dict     = cor_id_dict( threept_set[0], keys=exclude(lam_id_keys,['lambdas']), fit='threept_slope'     )
    id_dict_cmp = cor_id_dict( threept_set[0], keys=exclude(lam_id_keys,['lambdas']), fit='threept_slope_cmp' )
    R    = np.zeros( twopt.shape )
    for it in range(nt):
      R[:,:,:,it]    = np.sum( threept[:,:,:,it,delta:it-delta+1], -1 ) / twopt[:,:,:,it] 
    rat    = np.mean( np.sum( R[:,:,:,:]   , 0 ), 1 )
    #rolled_rats = [ np.roll( rat, -i, -1 ) for i in range(7) ]
    #coefficienct = [ []
    #               , [-1,1]
    #               , [-3.0/2.0   , 2.0, -1.0/2.0 ]
    #               , [-11.0/6.0  , 3.0, -3.0/2.0 , 1.0/3.0 ]
    #               , [-25.0/12.0 , 4.0, -3.0     , 4.0/3.0 , -1.0/4.0 ]
    #               , [-137.0/60.0, 5.0, -5.0     , 10.0/3.0, -5.0/4.0 , 1.0/5.0 ]
    #               , [-49.0/20.0 , 6.0, -15.0/2.0, 20.0/3.0, -15.0/4.0, 6.0/5.0, -1.0/6.0 ] ]
    #ratcmp = -49.0/20.0*rolled_rats[0] + 6*rolled_rats[1] - 15.0/2.0*rolled_rats[2] + 20.0/3.0*rolled_rats[3] - 15.0/4.0*rolled_rats[4] + 6.0/5.0*rolled_rats[5] - 1.0/6.0*rolled_rats[6]
    #ratcmp = reduce( lambda a,b: a+b, [ c*d for c,d in zip(coefficienct[1],rolled_rats) ] )
    ratcmp = np.roll( rat, -1, -1 ) - rat
    correlator    = boot_object( rat   , time=np.arange(nt) )
    correlatorcmp = boot_object( ratcmp, time=0.5+np.arange(nt) )
    id         = cor_id(threept_set[0],keys=lam_id_keys)
    if 'fit' in args:
      correlator    = FitThreePointCorrelator( id_dict    , correlator   , save_path, args['fit'] )
      correlatorcmp = FitThreePointCorrelator( id_dict_cmp, correlatorcmp, save_path, args['fit'], flat=True )
    _ChromaPlot( [correlatorcmp], save=save_path+'png/threept_'+id+'_d'+str(delta)+'_cmp.png', xstride=12, xlims=ast.literal_eval(args.get('xlims','[0.0,24.0]')) )
    _ChromaPlot( [correlator]   , save=save_path+'png/threept_'+id+'_d'+str(delta)+'.png'    , xstride=12, xlims=ast.literal_eval(args.get('xlims','[0.0,24.0]')) )

def NonForwardThreePointRatio( nt, threept_cor, twopt_src_cor, twopt_sync_cor, time_reverse_avg=False ):
  if time_reverse_avg:
    threept_cor    = np.mean( threept_cor   , 1 )
    twopt_src_cor  = np.mean( twopt_src_cor , 1 )
    twopt_sync_cor = np.mean( twopt_sync_cor, 1 )
  else:
    threept_cor    = threept_cor[:,0,:,:]
    twopt_src_cor  = twopt_src_cor[:,0,:]
    twopt_sync_cor = twopt_sync_cor[:,0,:]

  rat = np.zeros( threept_cor.shape, dtype=threept_cor.dtype )
  for t in range(nt):
    for tau in range(nt):
      sqrt_term_num = twopt_src_cor[:,t-tau]  * twopt_sync_cor[:,tau] * twopt_sync_cor[:,t]  
      sqrt_term_den = twopt_sync_cor[:,t-tau] * twopt_src_cor[:,tau]  * twopt_src_cor[:,t] 
      sqrt_term     = sqrt_term_num / sqrt_term_den
      rat[:,t,tau]  = threept_cor[:,t,tau] / twopt_sync_cor[:,t] * np.sqrt( sqrt_term )

  correlator = boot_object( rat, time=np.arange(nt), tau=np.arange(nt) )
  return correlator

def TwoDimThreePoint( cors, args, save_path ):
  print 'Performing Two Dim fit'
  xlims = literal_eval(args.pop('xlims','[0.0,24.0]'))
  excluded_keys = data_keys + bookeeping_keys + [ 'mom', 'spin', 'q', 'match_vector' ]
  threept_sets = Partition( partial(EqualExcept, excluded_keys), cors )
  s_raw = np.asarray(ast.literal_eval(args['s']))
  s_raw = map(lambda (pp,p,a): (tuple(pp),tuple(p),a), s_raw)
  s_pp = map(Fst,s_raw)
  s = { pp:{} for pp in s_pp }
  for (pp,p,a) in s_raw:
    s[pp][p] = a

  keys = exclude( cors[0].keys(), excluded_keys )
  print keys
  for threept_set in threept_sets:
    print sub_dict( keys, threept_set[0] )

  cor_set = []
  for threept_set in threept_sets:
    nt = threept_set[0]['dim'][-1]
    data_type = args.get( 'data_type', 'cmpl' )
    coef_type = args['coef_type']
    pp = threept_set[0]['momentum']
    q  = threept_set[0]['q']
    p  = [ a - b for a,b in zip(pp,q) ]
    pp2 = sumsquares(pp)
    q2  = sumsquares(q)
    p2  = sumsquares(p)
    sunit = unit(s[tuple(pp)][tuple(p)])
    ff_kinematics = 'q{:02d}pp{:02d}p{:02d}'.format(q2,min(pp2,p2),max(pp2,p2)) + '_' + coef_type
    if Get('spin')(threept_set[0]) in [ 'up', 'down', 'avg' ]:
      threept_set_up, threept_set_down = SplitFilter( lambda x: Get('spin')(x) in [ 'up', 'avg' ], threept_set )
      twopt_sync_up    = np.mean( np.asarray([ cor['boot']      for cor in threept_set_up ]).real, 0 )
      twopt_src_up     = np.mean( np.asarray([ cor['boot2']     for cor in threept_set_up ]).real, 0 )
      threept_up       = np.mean( np.asarray([ cor['threeboot'] for cor in threept_set_up ])     , 0 )
      twopt_sync_down  = np.mean( np.asarray([ cor['boot']      for cor in threept_set_down ]).real, 0 )
      twopt_src_down   = np.mean( np.asarray([ cor['boot2']     for cor in threept_set_down ]).real, 0 )
      threept_down     = np.mean( np.asarray([ cor['threeboot'] for cor in threept_set_down ])     , 0 )
      twopt_src_unpol  = 0.5 * ( twopt_src_up  + twopt_src_down  )
      twopt_sync_unpol = 0.5 * ( twopt_sync_up + twopt_sync_down )
      threept_unpol    = 0.5 * ( threept_up    + threept_down    )
      threept_pol      = 0.5 * ( threept_up    - threept_down    )
      threept = sunit[0] * threept_unpol + sunit[-1] * theept_pol
    elif Get('spin')(threept_set[0]) in [ 'unpol', 'polx', 'poly', 'polz' ]:
      threept_set_unpol = filter( lambda x: Get('spin')(x) in [ 'unpol' ], threept_set )
      threept_set_polx  = filter( lambda x: Get('spin')(x) in [ 'polx'  ], threept_set )
      threept_set_poly  = filter( lambda x: Get('spin')(x) in [ 'poly'  ], threept_set )
      threept_set_polz  = filter( lambda x: Get('spin')(x) in [ 'polz'  ], threept_set )
      twopt_sync_unpol = np.mean( np.asarray([ cor['sync_boot']    for cor in threept_set_unpol ]), 0 )
      twopt_src_unpol  = np.mean( np.asarray([ cor['src_boot']     for cor in threept_set_unpol ]), 0 )
      threept_unpol    = np.mean( np.asarray([ cor['threept_boot'] for cor in threept_set_unpol ]), 0 )
      threept_polx     = np.mean( np.asarray([ cor['threept_boot'] for cor in threept_set_polx  ]), 0 )
      threept_poly     = np.mean( np.asarray([ cor['threept_boot'] for cor in threept_set_poly  ]), 0 )
      threept_polz     = np.mean( np.asarray([ cor['threept_boot'] for cor in threept_set_polz  ]), 0 )
      threept_pol_vec  = np.asarray( [threept_unpol,threept_polx,threept_poly,threept_polz] )
      threept = np.tensordot( sunit, threept_pol_vec, axes=(0,0) )

    correlator = NonForwardThreePointRatio( nt, threept, twopt_src_unpol, twopt_sync_unpol, time_reverse_avg=False )
    cor_set.append( correlator )

  
  combined_cors = boot_object( np.mean([ c['boot'] for c in cor_set ],0), time=np.arange(nt), tau=np.arange(nt) )
  #_ChromaPlot( [np.mean(cor_set,0)], plotter=MultiPlot, save=save_path+args.pop('save','png/twodim_threept_'+ff_kinematics+'_s'+id_dict['s_str']+'_'+data_type+'_'+id+'.png'), xstride=12, xlims=xlims, **plot_args )
  plot_str = '_'.join([ 'twodim_threept', ff_kinematics, data_type ])
  id = extend([ { 'ff' : ff_kinematics, 'op':'cloverExB_x_AB', 's':sunit, 'pp':pp, 'q':q, 'p':p }, sub_dict( ['dim','kappas','hadron'], threept_sets[0][0] ) ])
  if 'fit' in args:
    fit_obj = json.load(open(args['fit'],'r'))
    #print ff_kinematics
    #print sorted( fit_obj.keys() )
    print ff_kinematics
    if ff_kinematics in fit_obj:
      #print 'fitting: ', fit_obj[ff_kinematics]['trapfit']
      combined_cors, time_extent, tau_extent = FitThreePointTrapezoid( id, combined_cors, save_path, fit_obj[ff_kinematics]['trapfit'], data_type=data_type )
      args['time'] = args.get('time',time_extent)
      args['tau']  = args.get('tau' ,tau_extent)
  args.pop('s',None) #TO make the output string a little easier to read
  _ChromaPlot( [combined_cors], plotter=MultiPlot, save=save_path+'png/' + plot_str + '_combined.png', xstride=12, xlims=xlims, **args )
  _ChromaPlot( cor_set        , plotter=MultiPlot, save=save_path+'png/' + plot_str + '_seperate.png', xstride=12, xlims=xlims, **args )

def TwoDimThreePointR( cors, args, save_path ):
  print 'Performing Two Dim fit R'
  xlims = literal_eval(args.pop('xlims','[0.0,24.0]'))
  threept_sets = Partition( partial(EqualExcept, data_keys + bookeeping_keys + [ 'momentum', 'spin', 'q', 'match_vector' ]), cors )
  for threept_set in threept_sets:
    id_keys = exclude( lam_id_keys, ['type','op'] ) + ['label']
    id_dict = cor_id_dict( threept_set[0], fit='twodim_threeptr', keys=id_keys+['op'], extra_args=[('qsq',threept_set[0]['qsq'])] )
    id      = cor_id(      threept_set[0]                       , keys=id_keys )
    nt = id_dict['dim'][-1]
    correlator = boot_object( threept_set[0]['R'][:,0,:,:], time=np.arange(nt), tau=np.arange(nt) )
    plot_args = args.copy()
    if 'trapfit' in args:
      correlator, time_extent, tau_extent = FitThreePointTrapezoid( id_dict, correlator, save_path, args['trapfit'] )

      plot_args['time'] = str(time_extent)
      plot_args['tau']  = str(tau_extent)
    _ChromaPlot( [correlator], plotter=MultiPlot, save=save_path+'png/twodim_threept_'+id+'.png', xstride=12, xlims=xlims,  **plot_args )

twopt_ratios = {
  'mass'                 : MassRatio
, 'simple'               : SimpleRatio
, 'simple_odd'           : SimpleOddRatio
, 'simple_even'          : SimpleEvenRatio
, 'spin_odd'             : SpinDependentOddRatio
, 'spin_even'            : SpinDependentEvenRatio
, 'spin_partial_hilbert' : SpinDependentPartialHilbertRatio
, 'simple'               : SpinDependentSimpleRatio
, 'spin_odd_odd'         : SpinDependentOddOddRatio
, 'spin_independent_odd_odd_ratio' : SpinIndependentOddOddRatio
, 'spin_independent_phase_odd_odd_ratio' : SpinIndependentPhaseOddOddRatio
, 'polarised_phase_odd_odd_ratio' : PolarisedPhaseOddOddRatio
}  

threept_ratios = {
  'threepoint' : ThreePointMass
, 'twodimthreepoint' : TwoDimThreePoint
}

threeptR_ratios = {
  'twodimthreepointr' : TwoDimThreePointR
}
