#!/bin/env python
import os
import copy
import pystache as ms
import json 
import stat
from hl import *

args = ArgvToArgsAndFiles({})

data = json.load( open( 'points_clover7_5.json', 'r+' ) )
seen_momenta = []
for group in data['operator_groups']:
  for op_group in group:
    for item in op_group['cors']:
      for i in item:
        item['pp'] = tuple(item['pp'])
        item['p']  = tuple(item['p'] )
        seen_momenta.append(item['pp'])
        seen_momenta.append(item['p'])
moms = list(set(seen_momenta))

list_fmter = lambda mom: '['+','.join(map(str,mom))+']'
mom_fmter = lambda mom: '{:+d}{:+d}{:+d}'.format(*mom)
qer = lambda pp, p: [ a-b for a,b in zip(pp,p) ]
mom2 = lambda pp: sum([a**2 for a in pp])
key_to_ms = lambda key: ''.join([ '{{#', key, '}}', '-', key, '={{.}} {{/', key, '}}' ])
def get_stuff(objs):
  for obj in objs:
    #print obj.keys()
    print obj, obj['pp']
  return map( lambda obj: ((obj.pop('pp'),obj.pop('p'),obj.pop('op'),obj.pop('s')),obj), objs )

#for q2_group in data['operator_groups']:
#  for coef_group in q2_group:
#    print coef_group
#    for (pp,p,op,s),rest in get_stuff(coef_group['cors']):
#      dict(rest.items() + {
#        'pp' : mom_fmter(pp)
#      , 'p'  : mom_fmter(p)
#      , 'p'     : mom_fmter(p)
#      , 'q'     : list_fmter( qer(pp,p) )
#      , 's'     : list_fmter( s )
#      , 'qf'    : mom_fmter(  qer(pp,p) )
#      , 'key'   : (mom2(pp),mom2(qer(pp,p)),mom2(p))
#      , 'as' : 'weights/as_' + data['actions'][op] + '.yml'
#      , 'op' : op
#      , 'longop' : '_'.join(data['actions'][op].split('_')[:-2])
#      , 'threepoint_prefix' : '_'.join( ['as'] + [data['actions'][op]] + data['actions'][op].split('_')[:-2] )
#      }.items())
templated = {
  'q2_groups' : [
    { 'same_q2' : [ { 's' : list_fmter([ list_fmter([ list_fmter(cor['pp']), list_fmter(cor['p']), list_fmter(cor['s']) ]) for cor in coef_group['cors'] ])
                    , 'coef_type' : coef_group['cors'][0]['coef_type']
                    , 'cor_set' : [ dict(rest.items() + {
        'pp' : mom_fmter(pp)
      , 'p'  : mom_fmter(p)
      , 'p'     : mom_fmter(p)
      , 'q'     : list_fmter( qer(pp,p) )
      #, 's'     : list_fmter( s )
      , 'qf'    : mom_fmter(  qer(pp,p) )
      , 'key'   : (mom2(pp),mom2(qer(pp,p)),mom2(p))
      , 'fit_key' : 'q{:02}pp{:02}p{:02}_'.format(mom2(qer(pp,p)),min(mom2(pp),mom2(p)),max(mom2(pp),mom2(p))) + rest['coef_type']
      , 'as' : 'weights/as_' + data['actions'][op] + '.yml'
      , 'op' : op
      , 'longop' : '_'.join(data['actions'][op].split('_')[:-2])
      , 'threepoint_prefix' : '_'.join( ['as'] + [data['actions'][op]] + data['actions'][op].split('_')[:-2] )
      }.items()) for (pp,p,op,s),rest in get_stuff(coef_group['cors']) ] }
    for coef_group in q2_group ] }
  for q2_group in data['operator_groups'] ]
, 'twopoint_syncs' : [ mom_fmter(pp) for pp in moms ]
, 'mom2s' : list(set([ mom2(pp) for pp in moms ]))
, 'kappa' : '0.132000'
, 'kappa2' : '0.132'
, 'dim'   : '24_24_24_48'
, 'dim2'  : '24x48'
, 'qgroups' : []
}
def filtercorsforfit( g ):
  res = [ copy.deepcopy( coef_group ) for coef_group in g if 'trapfit' in coef_group ]
  return res
if 'fit' in args:
  fit_obj = json.load(open(args['fit'],'r'))
  templated['groupsfitted'] = filter( lambda cor: cor.get('fit_key') in fit_obj, [ coef_set['cor_set'][0] for q2_group in templated['q2_groups'] for coef_set in q2_group['same_q2'] ] )

templated.update(args.items())
templated['twopt_folder']   = "twopt_cors"
templated['threept_folder'] = "clover7.5"
templated['all_threept_cors'] = [ cor for q2_group in templated['q2_groups'] for coef_group in q2_group['same_q2'] for cor in coef_group['cor_set'] ]
coef_types = list(set([ cor['coef_type'] for q2_group in templated['q2_groups'] for coef_group in q2_group['same_q2'] for cor in coef_group['cor_set'] ]))
templated['coef_sets'] = { coef_type : [ { 'same_q2' : filter( lambda coef_group: coef_group['cor_set'][0]['coef_type'] == coef_type, q2_group['same_q2'] ) }
  for q2_group in templated['q2_groups'] if coef_type in [cor['coef_type'] for coef_group in q2_group['same_q2'] for cor in coef_group['cor_set']] ] 
  for coef_type in coef_types }

partials = {
  'preamble'      : "#!/bin/bash -x\n\n set -e\n\n test -d {{threept_folder}}/png || mkdir {{threept_folder}}/png\n test -d {{threept_folder}}/fit || mkdir {{threept_folder}}/fit\n\n"
, 'four_kappa'    : "{{kappa}}_{{kappa}}_{{kappa}}_{{kappa}}"
, 'config_folder' : "{{dim}}_{{kappa}}_{{kappa}}"
, 'twopointer'    : "nucleon_{unpol,polx,poly,polz}_psync{{.}}.json"
, 'twopointerriders' : "{{> twopointfolder}}/nucleon_{unpol,polx,poly,polz}_psync{{.}}.{json,dat}"
, 'twopointfolder' : "{{> config_folder}}/_propshsh_propshsh"
, 'threepointer'  : "{{threept_folder}}/{{> config_folder}}/_propshsh_propshsh_{{fit_key}}/nucleon_{unpol,polx,poly,polz}_{{longop}}_q{{qf}}_pp{{pp}}_p{{p}}.json"
, 'threepointgen' : "hl_threepoint_gen {{threept_folder}}/{{> config_folder}}/_propshsh_propshsh/nucleon_{unpol,polx,poly,polz}_psync{{pp}}.json {{threept_folder}}/{{> config_folder}}/_propshsh_propshsh/nucleon_{unpol,polx,poly,polz}_psync{{p}}.json -action={{as}} -q={{q}} -save_dir={{threept_folder}} -fit_key={{fit_key}} -noupdate"
, 'efffit'        : "fit/effmass_nucleon_{{kappa2}}_{{kappa2}}_{{dim2}}_{{.}}_0.json"
, 'decompthreepoint' : "{{prefix}}_threept_{{label}}_nucleon_ff_{{dim}}_{{> four_kappa}}.json"
, 'decompthreepointfit' : "twodim_threeptr_nucleon_{{prefix}}_{{kappa2}}_{{kappa2}}_{{dim2}}_{{label}}_{{op}}.json"
, 'bootstrap_set' : 'hl_bootstrap {{#cor_set}} {{> threepointer}}{{/cor_set}} {{#twopoint_syncs}}{{threept_folder}}/{{> twopointfolder}}/{{> twopointer}} {{/twopoint_syncs}} -save_dir={{threept_folder}} -noupdate' 
, 'effmassargs' : ' '.join( map(key_to_ms,['ylims','trapfit','tau','time','data_type']) )
, 'job_preamble' : '#!/bin/bash -x\n#SBATCH -N 1\n#SBATCH -n 4\n#SBATCH --time=3:30:00\n#SBATCH --mem=72GB\n#SBATCH -p batch\n#SBATCH --output={{job_name}}_%A_%a.out\n#SBATCH -J {{job_name}}\n'
}
template = """{{> preamble}}
mkdir -p {{threept_folder}}/{{dim}}_{{kappa}}_{{kappa}}/_propshsh_propshsh
{{^justfit}}
{{^nogen}}
{{#twopoint_syncs}}
for file in {{twopt_folder}}/{{> twopointerriders}}; do test -L {{threept_folder}}/{{> config_folder}}/_propshsh_propshsh/${file##*/} || ln -s $(pwd)/$file {{threept_folder}}/{{> config_folder}}/_propshsh_propshsh/${file##*/}; done
{{/twopoint_syncs}}

{{#q2_groups}}
#Q2 GROUP 
{{#same_q2}}
##COEF GROUP
  {{#cor_set}}
{{> threepointgen}}\n
  {{/cor_set}}
{{/same_q2}}
{{/q2_groups}}
{{/nogen}}


{{#q2_groups}}
#Q2 GROUP 
{{#same_q2}}
##COEF GROUP
{{> bootstrap_set}}\n
{{/same_q2}}
{{/q2_groups}}

{{/justfit}}

hl_effmass -r=m {{#twopoint_syncs}} {{threept_folder}}/{{> twopointfolder}}/{{> twopointer}} {{/twopoint_syncs}} -xlims=[0.0,24.0] -ylims=[0.9,1.1] -fit=[9,21] -save_dir={{threept_folder}}

{{#q2_groups}}
#Q2 GROUP 
{{#same_q2}}
##COEF GROUP
hl_effmass {{#cor_set}} {{> threepointer}} {{/cor_set}} -r=2d3 {{#cors}}{{> effmassargs}}{{/cors}} -s={{s}} -save_dir={{threept_folder}} {{#fit}}-fit={{.}}{{/fit}} -coef_type={{coef_type}} {{#fast}}&{{/fast}}
{{/same_q2}}
{{/q2_groups}}
{{#fast}}wait{{/fast}}

hl_fit {{threept_folder}}/fit/effmass_nucleon_0.132_0.132_24x48*.json {{#groupsfitted}} {{threept_folder}}/fit/{{fit_key}}.json{{/groupsfitted}} -xlims=[0.0,0.3] -fs=28 -save_dir={{threept_folder}}

"""

gen_template = """{{> job_preamble}}
{{#twopoint_syncs}}
for file in {{twopt_folder}}/{{> twopointerriders}}; do test -L {{threept_folder}}/${file##*/} || ln -s $(pwd)/$file {{threept_folder}}/${file##*/}; done
{{/twopoint_syncs}}

#COEF GROUP
{{#cor_set}}
{{> threepointgen}}\n
{{/cor_set}}

{{> bootstrap_set}}
"""

job_submitter = """
#!/bin/bash -x
set -e
for script in seperate_jobs/*
do
  qsub $script
done
"""

ms_renderer = ms.Renderer(partials=partials)
open('do_clover7_5_threepoint_generation','w+').write(ms_renderer.render(template,templated))
os.chmod('do_clover7_5_threepoint_generation', os.stat('do_clover7_5_threepoint_generation').st_mode | stat.S_IEXEC)
if not os.path.exists('seperate_jobs'):
  os.mkdir('seperate_jobs')
#for i,q2_group in enumerate(templated['q2_groups']):
#  templated['same_q2'] = q2_group['same_q2']
for i,coef_group in enumerate([ t for q2_group in templated['q2_groups'] for t in q2_group['same_q2'] ]):
  templated['cor_set'] = coef_group['cor_set']
  templated['job_name'] = '3ptgen{:03}'.format(i)
  f = 'seperate_jobs/do_clover7_5_threepoint_generation_part'+'{:03}'.format(i)
  open(f,'w+').write(ms_renderer.render(gen_template,templated))
  os.chmod( f, os.stat(f).st_mode | stat.S_IEXEC )
open('do_all_jobs_again','w+').write(job_submitter)
os.chmod('do_all_jobs_again', os.stat('do_all_jobs_again').st_mode | stat.S_IEXEC)
for k,v in templated['coef_sets'].items():
  templated['q2_groups'] = v
  f = 'do_clover7_5_threepoint_generation_' + k
  open(f,'w+').write(ms_renderer.render(template,templated))
  os.chmod( f, os.stat(f).st_mode | stat.S_IEXEC )

