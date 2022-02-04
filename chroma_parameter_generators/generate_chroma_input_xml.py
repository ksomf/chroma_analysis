import os
from math import *
from functools import partial
import stat

def create_executable_script( f, s ):
  open( f, 'w+' ).write( s )
  os.chmod( f, os.stat(f).st_mode | stat.S_IEXEC )

def Fst( xs ):
  return xs[0]

def Snd( xs ):
  return xs[1]

def Lst( xs ):
  return xs[-1]

def Negate( x ):
  if x != 0:
    return -x
  else:
    return x

def extend( objs ):
  return { k:v for obj in objs for k,v in obj.items() }

def cross( list, objlist ):
  return [ extend( [ obj2, obj1 ] ) for obj2 in list for obj1 in objlist ]

def collapse( obj ):
  return [ { k:v[i] for k,v in obj.items() } for i in range( len( obj[ obj.keys()[ 0 ] ] ) ) ]

def PartialCollapse( obj, donot_collapse=[] ):
  not_list = { k:v for k,v in obj.items() if not isinstance( v, list ) or k in donot_collapse }
  res = map( lambda s: extend([s,not_list]), collapse({ k:v for k,v in obj.items() if isinstance(v,list) and not k in donot_collapse }) )
  return res
def uncollapse( xs ):
  return { k:[ o[k] for o in xs ] for k in xs[0].keys() }
def listcross( accum, extra ):
  return [ ac + [ e ] for ac in accum for e in extra ]
def UncollapseCross( objss ):
  if len( objss ):
    return map( uncollapse, reduce( listcross, objss[1:], map( lambda s: [s], objss[0] ) ) )
  else:
    return []

def SetupLambda( obj ):
  items_to_append = []
  default_q = [[0,0,0]]*len(obj.get('lams',[1]))
  for i,q in enumerate( obj.get('q',default_q) ):
    if q != [0,0,0]:
      items_to_append.append( (i,map(Negate,q)) )

  for i,q in items_to_append:
    for k,v in obj.items():
      if k == 'q':
        obj[k].append(q)
      elif 'q_type' in obj and obj['q_type'][i] == 'sin' and k == 'lams' and obj[k][i] != 0.0:
        obj[k].append(-obj[k][i])
      else:
        obj[k].append(obj[k][i])

  op_dict = {
    'g1'    : 1
  , 'ig1'   : 1
  , 'g3'    : 4
  , 'ig3'    : 4
  , 'g4'    : 8
  , 'ig4'   : 8
  , 'ig1g5' : 14
  , 'g1g5'  : 14
  , 'ig3g5' : 11
  , 'g4g5' : 7
  , 'ig4g5' : 7
  }
  lams = [ (l,0.0) if not 'i' in o else (0.0,l) for l,o in zip( obj['lams'], obj['op'] ) ]
  obj['lam_real'] = map( Fst, lams )
  obj['lam_imag'] = map( Snd, lams )
  obj['op_num']   = [ op_dict[o] for o in obj['op'] ]
    
def HasNonzeroOps( obj, exclude ):
  op_matches = [ (i,v) for i,v in enumerate(obj['op']) if v in exclude ]
  res = all([ obj['lams'][i] != 0.0 for i,v in op_matches ])
  if res:
    print '[INFO] Excluded match:', op_matches, obj, exclude
  return res
    
def BuildLamObjects( objs, excluded_ops ):
  lists = map( partial(PartialCollapse,donot_collapse=['q']), objs )
  #res = filter( lambda s: not all( map( lambda t: t == 0, s['lams'] ) ), UncollapseCross( lists ) )
  res = filter( lambda s: sum(map( lambda t: t != 0, s['lams'] )) in [1,2], UncollapseCross( lists ) )
  res = filter( lambda s: all([ not HasNonzeroOps( s, excluded_op_group ) for excluded_op_group in excluded_ops ]), res )
  map( SetupLambda, res )
  return map( collapse, res )

def lennd( xs ):
  return sqrt( sum( [ x**2 for x in xs ] ) )
def dist_between_line_and_point( ls, le, pt ):
  if ls == le:
    le = [ 0.0 for l in ls ]
  vs = zip(ls,le,pt)
  sme = [ s - e for (s,e,p) in vs ]
  smp = [ s - p for (s,e,p) in vs ]
  pme = [ p - e for (s,e,p) in vs ]
  unit_line = [ s / lennd( sme ) for s in sme ]
  
  dot = lambda xs,ys: sum( [ x*y for (x,y) in zip(xs,ys) ] )
  coef = dot(smp,unit_line)

  projected_line = [ coef*l for l in unit_line ]
  
  normal_line = [ p - m for (p,m) in zip(smp,projected_line) ]
  #print 'Point:', pt, '(', ls, '--', le, ')', lennd( normal_line )
  
  return lennd(normal_line), copysign( lennd(projected_line), -coef ), lennd(sme), max(lennd(pme),lennd(smp))
def effective_dist_between_line_and_point( ls, le, pt ): 
   dist, lerp_dist, lerp_gap, metric = dist_between_line_and_point( ls, le, pt )
   #this metric determines the choice between initial vs linear initial guess roughly the best 
   return 1.1*dist + 0.4*metric 
#def negate_distance( xs, ys ):
#  return lennd( [ x+y for (x,y) in zip(xs,ys) ] )
def IsAvailable( available, needed ):
  return all( [ need in available for need in needed[1:-1] ] )

get_fh_param = lambda s: {'fh-params':s}
def AddShift( lams, res, objs, pick, picked, dependency ):
  #print 'PICKING:', pick
  picked.append(pick[0])
  if len( pick ) == 2:
    new_obj = get_fh_param( objs[pick[0]] )
    new_obj['initial_guess'] = { 'fh-params' : [] }
    dependency.append([])
  elif len( pick ) == 4:
    new_obj = get_fh_param( objs[pick[0]] )
    if pick[2] == pick[1]:
      new_obj['initial_guess'] = { 'fh-params' : [] }
    else:
      new_obj['initial_guess'] = get_fh_param( objs[pick[2]] )
    
    dist, lerp_dist, lerp_gap, dist = dist_between_line_and_point( lams[pick[1]], lams[pick[2]], lams[pick[0]] )
    new_obj['prior_perturbed_guess'] = get_fh_param( objs[pick[1]] )
    new_obj['prior_perturbed_guess'].update( {'factor': 1 + lerp_dist / lerp_gap} )
    dependency.append( [pick[1], pick[2]] )
    #print lams[pick[2]], '--', lams[pick[1]], '->', lams[pick[0]], '(', 1 + lerp_dist / lerp_gap, ')'
  res.append( new_obj )

def BuildPropagators( objs ):
  res = []
  if objs:
    lams = [ [ o['lams'] for o in obj ] for obj in objs ]
    #for a in enumerate(lams):
       #print a
    zero_shift_length = sorted( [ (i,lennd( ls )) for i,ls in enumerate(lams) ], key=Lst )
    #negate_distance   = [ [ negate_distance( ls1, ls2 ) for ls2 in lams ] for ls1 in lams ]
    lerp_distance     = [ (i,j,k,effective_dist_between_line_and_point( ls2, ls3, ls1 ))
                            for (k,ls3) in enumerate(lams)
                          for (j,ls2) in enumerate(lams)
                        for (i,ls1) in enumerate(lams) ]
    distances = sorted( zero_shift_length + lerp_distance, key=Lst )
    
    first_pick = zero_shift_length[ 0 ]
    picked_indicies = []
    dependency_indicies = []
    AddShift( lams, res, objs, first_pick, picked_indicies, dependency_indicies )
    while len( picked_indicies ) != len( lams ):
      available_options = filter( lambda s: IsAvailable( picked_indicies, s ), distances )
      available_options = filter( lambda s: Fst( s ) not in picked_indicies, available_options )
      AddShift( lams, res, objs, available_options[0], picked_indicies, dependency_indicies )
     
    dependencies = [ [a] + b for (a,b) in zip(picked_indicies,dependency_indicies) ]
    last_depend_indecies = [ max([k if i in v else 0 for (k,v) in enumerate(dependencies)]) for i in range( 1+max( map( max, dependencies ) ) ) ]
    erase_points = list( set( last_depend_indecies ) )
    erasers = [ (k,[ i for (i,v) in enumerate(last_depend_indecies) if v == k ]) for k in erase_points ]
    for k,erase in erasers:
      res[k].update( {'erase_objects': [ get_fh_param(objs[e]) for e in erase ]} )
      #print res[k]
  
  return res

def generate_invert_combine( render_ms, inverter_template, combiner_template, input_parameters, (jobber_name,jobber_template), (chroma_name, chroma_template), prefix='' ):
  q_type = input_parameters.get( 'q_type', 'cos' )
  excluded_op_combinations = input_parameters.get( 'excluded_op_combinations', [] )
  lam_q  = BuildLamObjects( input_parameters['shifts'], excluded_op_combinations )

  for shift in input_parameters['shifts']:
    shift['lams'] = map( Negate, shift['lams'] )
  lam_mq = BuildLamObjects( input_parameters['shifts'], excluded_op_combinations )
  for shift in input_parameters['shifts']:
    shift['lams'] = map( Negate, shift['lams'] )
  
  for shift in input_parameters['shifts']:
    if input_parameters['q'] == [0,0,0] and 'g5' in shift['op'] and 'q1-q2' in input_parameters['quarks']:
      shift['lams'] = list( set( map( Negate, shift['lams'] ) + shift['lams'] ) )
  lam_inv = BuildLamObjects( input_parameters['shifts'], excluded_op_combinations )

  input_parameters['propagators'] = BuildPropagators( lam_inv )
  print '[INFO] prop inversions: ', len(input_parameters['propagators'])
  if lam_q and len( lam_q[0] ) > 1:
    slice_lams = list( set( map( abs, map( lambda s: s[0]['lams'], lam_q ) ) ) ) 
    slices_q   = [ filter( lambda s: abs( s[0]['lams'] ) == sl, lam_q   ) for sl in slice_lams ]
    slices_mq  = [ filter( lambda s: abs( s[0]['lams'] ) == sl, lam_mq  ) for sl in slice_lams ]
    slices_inv = [ filter( lambda s: abs( s[0]['lams'] ) == sl, lam_inv ) for sl in slice_lams ]
    slices = len(slice_lams)
  else:
    slice_lams = [ 0.0     ]
    slices_q   = [ lam_q   ]
    slices_mq  = [ lam_mq  ]
    slices_inv = [ lam_inv ]
    slices = 1

  total_combinations = []
  readins = []
  combination_sets = []
  get_props = lambda prop1=[], prop2=[]: { 'prop1':{ 'fh-params':prop1 }, 'prop2':{ 'fh-params':prop2 } }
  for i, (sl, lq, lmq, linv) in enumerate( zip( slice_lams, slices_q, slices_mq, slices_inv ), 2 ):
    readins.append( [ {} ] + [ { 'fh-params' : [ lo for lo in los ] } for los in linv ] )
    combinations = []
    #print i, slices, sl
    if i == 2: #if it is first
      combinations = combinations + [ get_props() ]
    if 'q1' in input_parameters['quarks']:
      combinations = combinations + [ get_props(prop1=los) for los in lq ]
    if 'q2' in input_parameters['quarks']:
      combinations = combinations + [ get_props(prop2=los) for los in lq ]
    if 'q1+q2' in input_parameters['quarks']:
      combinations = combinations + [ get_props(prop1=los,prop2=los) for los in lq ]
    if 'q1-q2' in input_parameters['quarks']:
      combinations = combinations + [ get_props(prop1=los,prop2=losm) for los,losm in zip( lq, lmq ) ]
    total_combinations = total_combinations + combinations
    combination_sets.append( combinations )

  if input_parameters['unified']:
    input_parameters['unified_run'] = { 'combinations':total_combinations }
    open( prefix+'unified.xml.ms', 'w+' ).write( render_ms( inverter_template, input_parameters ) )
  else:
    input_parameters['combine_runs'] = [ { 'input_file_index' : i } for i,sl in enumerate( slice_lams, 1 ) ]
    input_parameters['slices'] = [ { 'slice' : sl } for sl in slice_lams ]
    open( prefix+'separate_1.xml.ms', 'w+' ).write( render_ms( inverter_template, input_parameters ) )
    for i,(sl,reads,combinations) in enumerate(zip(slice_lams,readins,combination_sets), 2):
      input_parameters['readins']      = reads
      input_parameters['combinations'] = combinations
      input_parameters['slice']        = sl
      open( prefix+'separate_' + str(i) + '.xml.ms'  , 'w+' ).write( render_ms( combiner_template, input_parameters ) )
  create_executable_script( jobber_name, render_ms(jobber_template,input_parameters) )
  create_executable_script( chroma_name, render_ms(chroma_template,input_parameters) )
  
def generate_wflow( render_ms, wflow_template, dir_template, shift_template, params ):
  wtimes = sorted( params['wilson_flows'] )
  wflows = [ { 'wtime'    : wtimes[0], 'wsteps'    : int(100 * wtimes[0])
             , 'incwtime' : wtimes[0], 'incwsteps' : int(100 * wtimes[0]) } ]
  for i in range(1,len(wtimes)):
    wflows.append( {
      'incwtime'  :       wtimes[i] - wtimes[i-1]
    , 'incwsteps' : int(( wtimes[i] - wtimes[i-1] ) * 100)
    , 'wtime' : wtimes[i]
    , 'wsteps' : int(wtimes[i] * 100)
    , 'start'  : wflows[i-1]
    } )
  params['wilson_flows'] = wflows
  open( 'input.xml.ms', 'w+' ).write( render_ms( wflow_template, params ) )
  open( 'directory.sh', 'w+' ).write( render_ms( dir_template  , params ) )
  open( 'shift.sh'    , 'w+' ).write( render_ms( shift_template, params ) )
  os.chmod( 'directory.sh', os.stat('directory.sh').st_mode | stat.S_IEXEC )
  os.chmod( 'shift.sh', os.stat('shift.sh').st_mode | stat.S_IEXEC )
