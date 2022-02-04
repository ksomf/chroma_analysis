#!/bin/env python

import json
import os
import stat
import pystache as ms
from generate_chroma_input_xml import generate_chroma_input_xml, generate_wflow
import itertools
import operator

script_dir = os.path.dirname(os.path.abspath(__file__))
chroma_base_dir = os.path.dirname(script_dir)
src_base_dir = os.path.dirname(chroma_base_dir)
template_dir = '/'.join([script_dir,'templates'])

def get_factors( i ):
  res = [1]
  for j in range(2,i):
    if i % j == 0:
      res.append(j)
  res.append(i)
  return res
def load_ms( filepath ):
  file_contents = open( filepath, 'r' ).read().replace( '{{{', '[[[' ).replace( '}}}', ']]]' )
  return ms.parse( unicode( file_contents ) )

partial_templates = json.load( open( template_dir + 'partials.yml', 'r' ) )
ms_renderer = ms.Renderer(partials=partial_templates)
def render_ms( template, obj ):
  return ms_renderer.render( template, obj ).replace( '[[[', '{{' ).replace( ']]]', '}}' )

jobber_template   = load_ms( template_dir + 'jobber.sh.ms'          )
chroma_template   = load_ms( template_dir + 'chroma.sh.ms'          )
inverter_template = load_ms( template_dir + 'inverter.xml.ms.ms'    )
combiner_template = load_ms( template_dir + 'combiner.xml.ms.ms'    )
wflow_template    = load_ms( template_dir + 'wilson_flow.xml.ms.ms' )
dir_template      = load_ms( template_dir + 'directory.sh.ms'       )
shift_template    = load_ms( template_dir + 'traj.sh.ms'            )

hosts = {
  'raijin' : {
    'common'  : { 'raijin':True, 'jobfs':200, 'pbs':True, 'project':'<removed>', 'usencpus':True }
  , 'scripts' : [ { 'name':'kepler'        , 'queue':'gpu'      , 'gpn':8, 'nodes':1, 'cpn':24                      , 'walltime':'48:00:00', 'mpn':128, 'unified':False } 
                , { 'name':'pascal'        , 'queue':'gpupascal', 'gpn':4, 'nodes':1, 'cpn':24                      , 'walltime':'48:00:00', 'mpn':128, 'unified':False } 
                , { 'name':'sandybridge'   , 'queue':'normal'   ,          'nodes':4, 'cpn':16, 'geometry':'1 4 4 4', 'walltime':'48:00:00', 'mpn':32 , 'unified':True  }
                , { 'name':'sandybridgeexp', 'queue':'express'  ,          'nodes':4, 'cpn':16, 'geometry':'1 4 4 4', 'walltime':'24:00:00', 'mpn':32 , 'unified':True  }
                , { 'name':'broadwell'     , 'queue':'normalbw' ,          'nodes':5, 'cpn':28, 'geometry':'1 4 4 8', 'walltime':'48:00:00', 'mpn':128, 'unified':True  }
                , { 'name':'broadwellexp'  , 'queue':'expressbw',          'nodes':5, 'cpn':28, 'geometry':'1 4 4 8', 'walltime':'24:00:00', 'mpn':128, 'unified':True  }
                , { 'name':'knl'           , 'queue':'knl'      ,          'nodes':2, 'cpn':64, 'geometry':'1 4 4 8', 'walltime':'48:00:00', 'mpn':192, 'unified':True  } ]
  }
, 'hlogin' : {
  'common'  : { 'hlrn':True, 'pbs':True, 'cray':True, 'feature':'mpp', 'dds':True, 'SMALL_MESON':True }
  , 'scripts' : [ { 'name':'cray_account' , 'queue':'mppq', 'nodes':32, 'cpn':16, 'mpn':256, 'geometry':'4 4 4 8', 'unified':True, 'walltime':'12:00:00', 'account':'<removed>' }
                , { 'name':'cray_personal', 'queue':'mppq', 'nodes':32, 'cpn':16, 'mpn':256, 'geometry':'4 4 4 8', 'unified':True, 'walltime':'12:00:00', 'account':'<removed>' } ]
  }
, 'phoenix' : {
  'common'  : { 'phoenix':True, 'slurm':True, 'srcdir':src_base_dir }
  , 'scripts' : [ { 'name':'gpu', 'partition':'gpu', 'gpn':4, 'gres':'gpu:4,tmpfs:100GB', 'qos':'qxl', 'nodes':1, 'cpn':16, 'geometry':'1 2 2 4', 'mpc':3, 'walltime':'48:00:00', 'unified':False } ] 
  }
          
}

#hosts = {
#  'eddie'   : { 'name' : 'eddie' }
#, 'isaac'   : { 'name' : 'isaac', 'slurm' : { 'user' : '<removed>' }, 'partition' : 'batch', 'gpu_partition' : 'batch' }
#, 'hpn'     : { 'name' : 'isaacrc', 'slurm' : { 'user' : '<removed>' }, 'partition' : 'batch', 'gpu_partition' : 'gpu' }
#, 'phoenix' : { 'name' : 'phoenix', 'slurm' : { 'user' : '<removed>' }, 'partition' : 'batch', 'gpu_partition' : 'gpu' }
#, 'jrl'     : { 'name' : 'jureca', 'slurm' : { 'project' : '<removed>' }, 'partition' : 'batch', 'gpu_partition' : 'gpus' }
#, 'magnus'  : { 'name' : 'magnus', 'slurm' : { 'user' : '<removed>' } }
#, 'raijin'  : { 'name' : 'raijin', 'pbs' : { 'project' : '<removed>' }, 'partition' : 'normalbw', 'gpu_partition' : 'gpu' }
#}
hostname      = os.environ[ 'HOSTNAME' ]
if 'isaac' in hostname and 'phoenix' in hostname:
  server = 'hpn'
else:
  server = next( server for server in hosts.keys() if server in hostname )
server_params = map( lambda d: dict(hosts[server]['common'],**d), hosts[server]['scripts'] )

default_template_params = {
  "q":[0,0,0]
, "shifts" : []
, "quarks" : []
, "fermion" : {
    "action" : "SLRC"
    #"action" : "UNPRECONDITIONED_SLRC"
  , "csw"    : "2.65"
  }
, "cfg_type" : "ILDG"
, "cpu_precision" : "double"
, "gpu_precision" : "single"
}
default_job_params = {
  "resubmit" : False
, "sources" : 1
, "repeats" : 1
}

if 'wilson_flows' in params:
  generate_wflow( render_ms, wflow_template, dir_template, shift_template, params )
else:
  input_params = json.load( open( 'params.yml', 'r' ) )
  template_params = dict(default_template_params,**input_params)
  for submit_script_params in server_params:
    nodes = submit_script_params['nodes']
    if 'gpn' in submit_script_params:
      submit_script_params['gpus'] = nodes*submit_script_params['gpn']
    if not submit_script_params.get('cray',False):
      if 'mpn' in submit_script_params:
        submit_script_params['memory'] = nodes*submit_script_params['mpn']
    else:
      submit_script_params['memory'] = submit_script_params['mpn']
    submit_script_params['cpus']   = nodes*submit_script_params['cpn']

    dim = template_params['lattice']
    factors = map(get_factors,dim)
    possible_run_cpus = sorted(list(set(map(partial(reduce,operator.mul),itertools.product(*factors[2:])))))
    submit_script_params['run_cpus'] = filter( lambda x: x <= submit_script_params['cpus'], possible_run_cpus )[-1]
    print submit_script_params['name'], '({}->{})'.format(submit_script_params['cpus'],submit_script_params['run_cpus'])
    print submit_script_params

    params = dict( default_job_params, **dict( submit_script_params, **template_params ) )
    jobber_name = 'jobber_'+submit_script_params['name']
    chroma_name = 'chroma_'+submit_script_params['name']
    generate_invert_combine( render_ms, inverter_template, combiner_template, params, (jobber_name,jobber_template), (chroma_name,chroma_template), prefix='input_'+submit_script_params['name']+'_' )
