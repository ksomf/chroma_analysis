import json
import os.path
import numpy as np

from hl     import *
from hl_naming_convention import cor_path, chroma_str

def refresh_extensions( cor, file=None ):
  if file:
    cor['path'], cor['file']    = get_path_filename( file )
  cor['sync_cor_file']        = change_extension( cor['file'], 'dat'       )
  cor['src_cor_file']         = change_extension( cor['file'], 'dat2'      )
  cor['current_file']         = change_extension( cor['file'], 'three'     )
  cor['boot_sync_cor_file']   = change_extension( cor['file'], 'boot'      )
  cor['boot_src_cor_file']    = change_extension( cor['file'], 'boot2'     )
  cor['boot_threepoint_file'] = change_extension( cor['file'], 'threeboot' )

  cor['cor_file']      = cor['sync_cor_file']
  cor['boot_cor_file'] = cor['boot_sync_cor_file']

# take each q -q lambda pair and modify them into one
def cosine_redefine( prop ):
  if 'feynman_hellmann' not in prop:
    return
  feynman_hellmann_modifications = filter( lambda mod: mod['lambda'] != 0.0, prop.pop('feynman_hellmann') )
  res = []
  used_modifications = []
  for i,feynman_hellmann_modification in enumerate(feynman_hellmann_modifications):
    if i not in used_modifications:
      feynman_hellmann_modification_minus_q                  = dict( feynman_hellmann_modification
                                                                   , **{'op_q':map(Neg,feynman_hellmann_modification['op_q'])} )
      feynman_hellmann_modification_minus_q_and_minus_lambda = dict( feynman_hellmann_modification_minus_q
                                                                   , **{'lambda':-feynman_hellmann_modification_minus_q['lambda']})
      cosine_matches = [ j for j,v in enumerate(feynman_hellmann_modifications[i+1:],i+1) if v==feynman_hellmann_modification_minus_q ]
      sin_matches    = [ j for j,v in enumerate(feynman_hellmann_modifications[i+1:],i+1) if v==feynman_hellmann_modification_minus_q_and_minus_lambda ]
      if len(cosine_matches):
        used_modifications.append(cosine_matches[0])
        feynman_hellmann_modification['op'] = feynman_hellmann_modification['op']+'_cos'
      elif len(sin_matches):
        used_modifications.append(sin_matches[0])
        feynman_hellmann_modification['op'] = feynman_hellmann_modification['op']+'_sin'
      res.append(feynman_hellmann_modification)
  prop['feynman_hellmann'] = res
  #print prop

def u_d_classification( prop1, prop2 ):
  #assume that the operator ordering on the propagators is the same
  cosine_redefine( prop1 )
  cosine_redefine( prop2 )
  if 'feynman_hellmann' not in prop1 and 'feynman_hellmann' not in prop2:
    res = dict(prop1)
  elif 'feynman_hellmann' not in prop1:
    res = dict(prop2,**{'prop_modification_type':'quark2'})
  elif 'feynman_hellmann' not in prop2:
    res = dict(prop1,**{'prop_modification_type':'quark1'})
  else:
    fh_mod1 = prop1['feynman_hellmann']
    fh_mod2 =  prop2['feynman_hellmann']
    neg_fh_mod2 = [ dict(mod,**{'lambda':-mod['lambda']}) for mod in prop2['feynman_hellmann'] ]
    #print fh_mod1, fh_mod2, neg_fh_mod2
    if fh_mod1 == fh_mod2:
      res = dict(prop1,**{'prop_modification_type':'quark1_plus_quark2'})
    elif fh_mod1 == neg_fh_mod2:
      res = dict(prop1,**{'prop_modification_type':'quark1_minus_quark2'})
    else:
      print 'do not know what to do with pair: \n', prop1, '\n', prop2
      exit(1)

  return res

def load_correlator( file ):
  cor = json.load( open( file, 'r' ) )
  refresh_extensions( cor, file )
  if 'mom' not in cor and 'momentum' in cor:
    cor['mom'] = np.sum( np.asarray(cor['momentum'])**2 )
  if 'pp' not in cor and 'momentum' in cor:
    cor['pp'] = cor['momentum']
  if 'p' not in cor and cor.get('threepoint',False):
    cor['p'] = [ pp - q for (pp,q) in zip(cor['pp'],cor['q']) ]
  if 'propagators' in cor:
    props = cor.pop( 'propagators' )
    cor['extraction']['prop'] = props
    cor.update( u_d_classification( *props ) )
    for v in cor.get('feynman_hellmann',[]):
      if 'g' in v['op']:
        v['lambda'] = v['lambda'] / 2.0
  return cor

def LoadCorrelatorJSON( files ):
  cors = map( load_correlator, files )
  cors.sort( key=lambda s:len(s.get('op',[])), reverse=True )
  #cors = filter( lambda cor: max( map( abs, cor.get('lambdas',[0.0]) ) ) < 0.5, cors )
  #print [ s.get('type',None) for s in cors ]
  cors, zero_cors = SplitFilter( lambda s: 'feynman_hellmann' in s, cors )
  return cors, zero_cors

def LoadReweightJSONAndData( arg ):
  if '[' in arg:
    reweight_filenames = arg[1:-1].split(',')
  else:
    reweight_filenames = [ arg ]
  reweight_jsons     = map( LoadJSONFile, reweight_filenames )
  for file,jso in zip( reweight_filenames, reweight_jsons ):
    jso['shift'] = '.'.join( file.split('.')[:-1] ) + '_' +  jso['op']
  #loader = lambda rwfn: np.memmap( '.'.join( rwfn.split('.')[:-1] ) + '.dat', dtype=np.float64, mode='r' )
  loader = lambda rwfn: np.memmap( change_extension(rwfn,'dat'), dtype=np.float64, mode='r' )
  reweight_datas     = progmap( 'Mapping Reweight Files:', loader, reweight_filenames )
  reweight_datas     = [ rwght_data.reshape( (len(rwght_json['keys']),)+tuple(rwght_json['dim'][::-1]) ) for rwght_data, rwght_json in zip( reweight_datas, reweight_jsons ) ]
  print sng_report( 'Mapped Reweight Files:', reweight_datas[0].shape )
  return reweight_jsons, reweight_datas

def LoadCorrelatorData( args, cor, extract_keys=None ):
  if extract_keys is None:
    extract_keys = cor['keys']
  extract_keys = filter( Id, extract_keys )
  key_map = { k:i for i,k in enumerate( cor['keys'] ) }
  not_there_keys = [ k for k in extract_keys if k not in key_map ]
  if len(not_there_keys):
    print 'some keys were not found in', cor['path'], cor['file']
    print len(not_there_keys), not_there_keys
    exit(1)
  key_indicies = [ key_map[k] for k in extract_keys ]
  if key_indicies == range(len(extract_keys)):
    key_indicies = np.s_[0:len(extract_keys)]
  cor['used_keys'] = extract_keys
  if args.get('memmap',False):
    memory_loader       = np.memmap
    memory_loader_extra = { 'mode':'r' }
  else:
    memory_loader       = np.fromfile
    memory_loader_extra = {}
  unculled_shape = (len(cor['keys']),-1,cor['dim'][-1])
  cor['unculled_sync_correlator'] = memory_loader( cor_path(cor,'sync_cor_file'), dtype=np.complex128, **memory_loader_extra )
  cor['unculled_sync_correlator'] = cor['unculled_sync_correlator'].reshape( unculled_shape )
  cor['unculled_sync_correlator'] = cor['unculled_sync_correlator'][key_indicies,:,:]
  if os.path.isfile(cor_path(cor,'current_file')):
    cor['unculled_current'] = memory_loader( cor_path(cor,'current_file'), dtype=np.complex128 )
    cor['unculled_current'] = cor['unculled_current'].reshape( unculled_shape )
    cor['unculled_current'] = cor['unculled_current'][key_indicies,:,:]
    cor['unculled_source_correlator'] = memory_loader( cor_path(cor,'src_cor_file'), dtype=np.complex128 )
    cor['unculled_source_correlator'] = cor['unculled_source_correlator'].reshape( unculled_shape )
    cor['unculled_source_correlator'] = cor['unculled_source_correlator'][key_indicies,:,:]

#def LoadCorrelatorData( cor_json, keys=None, memmap=True, attach_key='cor' ):
#  nt = cor_json['dim'][-1]
#  if not keys:
#    keys = cor_json['keys']
#  cor_json['used_keys'] = keys
#  shape = (len(cor_json['keys']),-1,nt)
#  if memmap:
#    mapper = np.memmap
#    mapper_extra  = { 'mode':'r' }
#  else:
#    mapper = np.fromfile
#    mapper_extra  = {}
#  if 'raw' not in cor_json:
#    cor_json['key_map'] = { k:i for i,k in enumerate( cor_json['keys'] ) }
#    cor_json['raw']     = mapper( os.path.join(cor_json['path'],cor_json['sync_cor_file']), dtype=np.complex128, **mapper_extra )
#    cor_json['raw']     = cor_json['raw'].reshape( shape )
#  if 'rawinsertion' not in cor_json and os.path.isfile(os.path.join(cor_json['path'],cor_json['current_file'])):
#    cor_json['rawinsertion']     = np.fromfile( os.path.join(cor_json['path'],cor_json['current_file']), dtype=np.complex128 )
#    cor_json['rawinsertion']     = cor_json['rawinsertion'].reshape( shape )
#    cor_json['raw_source_twopt'] = mapper( os.path.join(cor_json['path'],cor_json['src_cor_file']), dtype=np.complex128, mode='r' )
#    cor_json['raw_source_twopt'] = cor_json['raw_source_twopt'].reshape( cor_json['raw'].shape )
#  indexes              = [ cor_json['key_map'][k] for k in keys ]
#  cor_json[attach_key] = cor_json['raw'][indexes,:,:]
#  if 'rawinsertion' in cor_json:
#    cor_json['insertion']    = cor_json['rawinsertion'][indexes,:,:]
#    cor_json['source_twopt'] = cor_json['raw_source_twopt'][indexes,:,:]
  
def LoadCorrelatorBootstraps( cor_json, attach_key='boot' ):
  #print cor_json['file']
  nt = cor_json['dim'][-1]
  if 'eff' in cor_json:
    return
  if cor_json['hadron'].count('g') == 2:
    time_directions = 1
  else:
    time_directions = 2
  cor_json['sync_boot'] = np.fromfile( cor_path(cor_json,'boot_cor_file'), dtype=np.complex128 )
  cor_json['sync_boot'] = cor_json['sync_boot'].reshape( (-1,time_directions,nt) )
  target_shape = cor_json['sync_boot'].shape
  if 'threepoint' in cor_json:
    cor_json['threept_boot'] = np.fromfile( cor_path(cor_json,'boot_threepoint_file'), dtype=np.complex128 )
    cor_json['threept_boot'] = cor_json['threept_boot'].reshape( target_shape + (nt,) )
    cor_json['src_boot']     = np.memmap( cor_path(cor_json,'boot_src_cor_file'), dtype=np.complex128 )
    cor_json['src_boot']     = cor_json['src_boot'].reshape( target_shape )

def _GetKeyAndSource( keysource ):
  tmp = keysource.split('.')
  return '.'.join(tmp[:-4]).replace('baryon_','').replace('meson_',''), np.asarray(map(int,tmp[-4:]))
def GetKeyAndSources( cor ):
  return zip( *map( _GetKeyAndSource, cor['keys'] ) )
 
data_keys = [ 'cor', 'key_map', 'used_keys', 'mom', 'type', 'lambdas', 'test', 'wgt', 'three_file' ]

def GetWeightPrefix( cor ):
  prefix = cor['weight_prefix']
  #print '[INFO] prefix is: ', prefix
  if 'lambdas' in cor:
    prefix = prefix + '_' + str(cor['lambdas'][0])
  if 'threepoint' in cor:
    prefix = prefix + '_threept'
  prefix = prefix + '_q' + threevec(cor.get('q',[0,0,0]))
  prefix = prefix + '_'
  return prefix

def GetThreepointFolderAndFile( cor, fit_key ):
  momfmter = lambda xs: ('{:+}'*len(xs)).format(*xs)
  extra_folder = fit_key
  #print cor['q'], cor['pp'], cor['p']
  file = '_'.join([cor['hadron'],cor['spin'],cor['op'],'q'+momfmter(cor['q']),'pp'+momfmter(cor['pp']),'p'+momfmter(cor['p'])] ) + '.json'
  return cor['path'] + '_' + extra_folder, file

def HasCurrentFile( cor, fit_key ):
  return os.path.isfile(change_extension( os.path.join(*GetThreepointFolderAndFile(cor,fit_key)), 'three' ))

def StoreCorrelatorJSON( cor, **kwargs ):
  if 'fit_key' in kwargs:
    cor['path'], cor['file'] = GetThreepointFolderAndFile( cor, kwargs['fit_key'] )
  refresh_extensions( cor )
  if not os.path.isdir( cor['path'] ):
    os.makedirs(cor['path'])
  after = []
  if 'current' in cor:
    current          = cor.pop('current')
    current_filename = cor['current_file']
    current.tofile( os.path.join(cor['path'],current_filename) )
    #print 'current: ', os.path.join(cor['path'],current_filename)
    source_twopt = cor.pop('unculled_src_correlator')
    source_twopt_file = cor['src_cor_file']
    source_twopt.tofile( os.path.join(cor['path'],source_twopt_file) )
    #print 'Src Twopt: ', os.path.join(cor['path'],source_twopt_file)
    after.append( ('current',current) )

  json_filename = cor.pop('file')
  data_filename = cor.pop('sync_cor_file')
  data          = cor.pop('unculled_sync_correlator')
  cor['propagators'] = cor['extraction'].pop('prop')

  [ cor.pop(k) for k in data_keys if k in cor ]
  json_filename = json_filename
  data_filename = data_filename 

  data.tofile( os.path.join(cor['path'],data_filename) )
  #print 'Sync Twopt: ', os.path.join(cor['path'],data_filename)
  json.dump(cor,open(os.path.join(cor['path'],json_filename),'w+'))
  #print 'JSON: ', os.path.join(cor['path'],json_filename)
  cor.update( dict( after ) )

#def StoreCombinationCorrelator( cor ):
#  refresh_extensions( cor, cor['file'] )
#  json_filename = cor.pop('file')
#  data_filename = cor.pop('data')
#  cor['extraction'].pop('prop')
#
#  [ cor.pop(k) for k in exclude(data_keys,['op','type']) if k in cor ]
#  json_filename = json_filename
#  data_filename = data_filename 
#
#  StoreBootData( cor )
#  [ cor.pop(s,'') for s in [ 'boot', 'threeboot', 'boot2', 'R', 'Runpol', 'Rpol' ] ]
#  #print cor
#  json.dump(cor,open(json_filename,'w+'))

def StoreFit( poly, save_path, save_append='', **kwargs ):
  res_json = kwargs
  res_json['eff'] = poly

  #print [ (k,chroma_str(k,kwargs[k])) for k in ['hadron','fit','prop_modification_type','data_type','pol_type','ff','mom','pp','feynman_hellmann'] if k in kwargs ]
  res_json['file_prefix'] = '_'.join( [ chroma_str(k,kwargs[k]) for k in ['hadron','fit','prop_modification_type','data_type','pol_type','ff','mom','pp','feynman_hellmann'] if k in kwargs ]  )
  print 'Saving Fit under: ', res_json['file_prefix']
  open( save_path+'fit/' + res_json['file_prefix'] + save_append + '.json', 'w+' ).write( json.dumps( res_json, indent=2, cls=NumpyJsonEncoder ) )
