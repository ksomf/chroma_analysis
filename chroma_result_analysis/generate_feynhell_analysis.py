#!/bin/env python

from hl import *
import pystache as ms
import copy

args, params = ArgvToArgsAndFiles({}, 'yml')
if len(params):
  param_file = params[0]
else:
  param_file = 'params.yml'
params = LoadJSONFile( param_file )

if 'common_params' in params:
  params = extend( [params, LoadJSONFile( params['common_params'] )] )

momfmter = lambda xs: ('{:+}'*len(xs)).format(*xs)

q = params['q']
nq = map(Neg,q)

def neg_lam( obj ):
  res = copy.deepcopy(obj)
  for o in res['feynhell']:
    if 'lam' in o:
      l = o['lam']
      if l[0] == '-':
        o['lam'] = l[1:]
      else:
        o['lam'] = '-'+l
    if 'lam1' in o:
      l = o['lam1']
      if l[0] == '-':
        o['lam1'] = l[1:]
      else:
        o['lam1'] = '-'+l
    if 'lam2' in o:
      l = o['lam2']
      if l[0] == '-':
        o['lam2'] = l[1:]
      else:
        o['lam2'] = '-'+l
  return res
def do_prop( obj, type, sync ):
  if obj:
    for fhp in obj.get('feynhell',[]):
      fhp['qstr']  = momfmter(fhp['q' ])
      fhp['nqstr'] = momfmter(fhp['nq'])

  if obj == []:
    res = [[],[]]
  elif type == 'u':
    res = [obj,[]]
  elif type == 'd':
    res = [[],obj]
  elif type == 'upd':
    res = [obj,obj]
  elif type == 'umd':
    res = [obj,neg_lam(obj)]
  else:
    print 'unhandled type', type
    exit(1)
  return {'sync':sync,'props':res}
def make_set( obj, type ):
  simple_form = obj.get('simplified_ratio',False)
  ratio = obj['ratio']
  lams  = obj['lams']
  op    = obj['op']
  mom   = obj['mom']
  sync  = momfmter( mom )
  nsync = momfmter( map(Neg,mom) )
  momxfer_type = obj.get('momxfer_type','cos')
  obj['ff'] = 'pp{}'.format(sync)
  print sync, nsync
  res = []
  for lam1 in lams:
    nlam1 = '-'+lam1
    if ratio == 'se': #SIMPLE ODD
      if momxfer_type == 'cos':
        l1 = lam1
        l2 = lam1
        nl1 = nlam1
        nl2 = nlam1
      elif momxfer_type == 'sin':
        l1 = lam1
        l2 = nlam1
        nl1 = nlam1
        nl2 = lam1
      res.append(do_prop({'feynhell':[{'op':op,'lam1':l1,'lam2':l2,'q':q,'nq':nq}]},type,sync))
      if not simple_form:
        res.append(do_prop({'feynhell':[{'op':op,'lam1':nl1,'lam2':nl2,'q':q,'nq':nq}]},type,sync))
      res.append(do_prop([],type,sync))
      if obj.get('minus_mom',True) and not sync == nsync:
        res.append(do_prop({'feynhell':[{'op':op,'lam1':l1,'lam2':l2,'q':q,'nq':nq}]},type,nsync))
        if not simple_form:
          res.append(do_prop({'feynhell':[{'op':op,'lam1':nl1,'lam2':nl2,'q':q,'nq':nq}]},type,nsync))
        res.append(do_prop([],type,nsync))
    if ratio == 'sioo' or ratio == 'sipoo' or ratio == 'ppoo' : #SPIN INDEPEPEDANT ODD ODD
      for lam2 in lams:
        nlam2 = '-'+lam2
        if momxfer_type == 'cos':
          res.append(do_prop({'feynhell':[{'op':op[0],'lam': lam1,'q':q,'nq':nq},{'op':op[1],'lam': lam2,'q':q,'nq':nq}]},type,sync))
          res.append(do_prop({'feynhell':[{'op':op[0],'lam':nlam1,'q':q,'nq':nq},{'op':op[1],'lam':nlam2,'q':q,'nq':nq}]},type,sync))
          res.append(do_prop({'feynhell':[{'op':op[0],'lam': lam1,'q':q,'nq':nq},{'op':op[1],'lam':nlam2,'q':q,'nq':nq}]},type,sync))
          res.append(do_prop({'feynhell':[{'op':op[0],'lam':nlam1,'q':q,'nq':nq},{'op':op[1],'lam': lam2,'q':q,'nq':nq}]},type,sync))
          if not sync == nsync:
            res.append(do_prop({'feynhell':[{'op':op[0],'lam': lam1,'q':q,'nq':nq},{'op':op[1],'lam': lam2,'q':q,'nq':nq}]},type,nsync))
            res.append(do_prop({'feynhell':[{'op':op[0],'lam':nlam1,'q':q,'nq':nq},{'op':op[1],'lam':nlam2,'q':q,'nq':nq}]},type,nsync))
            res.append(do_prop({'feynhell':[{'op':op[0],'lam': lam1,'q':q,'nq':nq},{'op':op[1],'lam':nlam2,'q':q,'nq':nq}]},type,nsync))
            res.append(do_prop({'feynhell':[{'op':op[0],'lam':nlam1,'q':q,'nq':nq},{'op':op[1],'lam': lam2,'q':q,'nq':nq}]},type,nsync))
        elif momxfer_type == 'sin':                                                                                                       
          if 'g5' in op[0] or 'i' in op[0]:
            res.append(do_prop({'feynhell':[{'op':op[0],'lam1': lam2,'lam2':nlam2,'q':q,'nq':nq},{'op':op[1],'lam': lam1,'q':q,'nq':nq}]},type,sync))
            res.append(do_prop({'feynhell':[{'op':op[0],'lam1':nlam2,'lam2': lam2,'q':q,'nq':nq},{'op':op[1],'lam':nlam1,'q':q,'nq':nq}]},type,sync))
            res.append(do_prop({'feynhell':[{'op':op[0],'lam1':nlam2,'lam2': lam2,'q':q,'nq':nq},{'op':op[1],'lam': lam1,'q':q,'nq':nq}]},type,sync))
            res.append(do_prop({'feynhell':[{'op':op[0],'lam1': lam2,'lam2':nlam2,'q':q,'nq':nq},{'op':op[1],'lam':nlam1,'q':q,'nq':nq}]},type,sync))
            if not sync == nsync:
              res.append(do_prop({'feynhell':[{'op':op[0],'lam1': lam2,'lam2':nlam2,'q':q,'nq':nq},{'op':op[1],'lam': lam1,'q':q,'nq':nq}]},type,nsync))
              res.append(do_prop({'feynhell':[{'op':op[0],'lam1':nlam2,'lam2': lam2,'q':q,'nq':nq},{'op':op[1],'lam':nlam1,'q':q,'nq':nq}]},type,nsync))
              res.append(do_prop({'feynhell':[{'op':op[0],'lam1':nlam2,'lam2': lam2,'q':q,'nq':nq},{'op':op[1],'lam': lam1,'q':q,'nq':nq}]},type,nsync))
              res.append(do_prop({'feynhell':[{'op':op[0],'lam1': lam2,'lam2':nlam2,'q':q,'nq':nq},{'op':op[1],'lam':nlam1,'q':q,'nq':nq}]},type,nsync))
          else:
            res.append(do_prop({'feynhell':[{'op':op[0],'lam': lam1,'q':q,'nq':nq},{'op':op[1],'lam1': lam2,'lam2':nlam2,'q':q,'nq':nq}]},type,sync))
            res.append(do_prop({'feynhell':[{'op':op[0],'lam':nlam1,'q':q,'nq':nq},{'op':op[1],'lam1':nlam2,'lam2': lam2,'q':q,'nq':nq}]},type,sync))
            res.append(do_prop({'feynhell':[{'op':op[0],'lam': lam1,'q':q,'nq':nq},{'op':op[1],'lam1':nlam2,'lam2': lam2,'q':q,'nq':nq}]},type,sync))
            res.append(do_prop({'feynhell':[{'op':op[0],'lam':nlam1,'q':q,'nq':nq},{'op':op[1],'lam1': lam2,'lam2':nlam2,'q':q,'nq':nq}]},type,sync))
            if not sync == nsync:
              res.append(do_prop({'feynhell':[{'op':op[0],'lam': lam1,'q':q,'nq':nq},{'op':op[1],'lam1': lam2,'lam2':nlam2,'q':q,'nq':nq}]},type,nsync))
              res.append(do_prop({'feynhell':[{'op':op[0],'lam':nlam1,'q':q,'nq':nq},{'op':op[1],'lam1':nlam2,'lam2': lam2,'q':q,'nq':nq}]},type,nsync))
              res.append(do_prop({'feynhell':[{'op':op[0],'lam': lam1,'q':q,'nq':nq},{'op':op[1],'lam1':nlam2,'lam2': lam2,'q':q,'nq':nq}]},type,nsync))
              res.append(do_prop({'feynhell':[{'op':op[0],'lam':nlam1,'q':q,'nq':nq},{'op':op[1],'lam1': lam2,'lam2':nlam2,'q':q,'nq':nq}]},type,nsync))
  return res
extras = params['points']
#  { 'ratio':'se'   , 'lams':lams, 'op':'g3'            ,                        'mom':[0,0,0], 'xlims':'[0,18]', 'ylims':'[-5.0,-1.0]', 'fit':'[3,11]', 'unpol':True }
groups = []
for d in extras:
  if 'forced_fitting_range' in params and 'fit' in d:
    d['fit'] = params['forced_fitting_range']
  for t in d['contractions']:
    def _make_object(d):
      return extend([{ 'type':t
                     , 'meson' : d['hadron'].count('g') == 2
                     , 'ratio':d['ratio']
                     , 'set' : make_set(d,t)                 },d])
    groups.append(_make_object(d))
    #if d.get('minus_mom',True): #RUN DIAG SINGLE MOM RUNS TOO
    #  d1 = dict(d)
    #  d2 = dict(d)
    #  d1['minus_mom'] = False
    #  d2['minus_mom'] = False
    #  d1['save_append'] = 'single_mom'
    #  d2['save_append'] = 'single_mom'
    #  d2['mom'] = map(Neg,d['mom'])
    #  groups.append(_make_object(d1))
    #  groups.append(_make_object(d2))
data = {
  'kappa':params.get('kappa','0.120900')
, 'displace':'.5'
, 'dim' : params.get('dim',[32,32,32,64])
, 'xlims':'[0,24]'
, 'groups':groups
, 'freefield':params.get('freefield',False)
}
if 'common_keys' in params:
  data['common_keys'] = params['common_keys']
if 'remove_keys' in params:
  data['remove_keys'] = params['remove_keys']


data['groups2'] = [ {'sub_groups':s} for s in splitlist( data['groups'], 10 ) ]
data['hadron_sync_pairs'] = [{'hadron':h,'sync':s,'meson': h.count('g')==2 } for h,s in list(set( (d['hadron'],momfmter(d['mom'])) for d in extras )) ]
data.update(args)

key_to_ms = lambda key: ''.join([ '{{#', key, '}}', '-', key, '={{.}} {{/', key, '}}' ])
partials = {
  'preamble' : '#!/bin/bash -x\n\n set -e\n\n test -d png || mkdir png\n test -d fit || mkdir fit\n\n'
, 'pq' : '{{op}}_{{lam1}}{{lam}}_q{{qstr}}'
, 'mq' : '{{op}}_{{lam2}}{{lam}}_q{{nqstr}}'
, 'dim' : '{{#dim}}{{.}}_{{/dim}}'
, 'folder1' : '{{> dim}}{{kappa}}_{{kappa}}'
, 'folder2' : '{{#props}}_propshsh{{#feynhell}}_{{> pq}}{{/feynhell}}{{#feynhell}}_{{> mq}}{{/feynhell}}{{/props}}'
, 'folder2mass' : '_propshsh_propshsh'
, 'file'    : '{{hadron}}{{^meson}}_{{#pol}}{polx,poly,polz}{{/pol}}{{^pol}}unpol{{/pol}}{{/meson}}_psync{{sync}}.json'
, 'twopoint' : '{{> folder1}}/{{> folder2}}/{{> file}}'
, 'justsynckappa' : '{{> folder1}}/{{> folder2mass}}/{{> file}}'
, 'effmassargs' : ' '.join( map(key_to_ms,['simplified_ratio','ff','xlims','ylims','fit','trapfit','tau','time','pol_type','data_type','displace','save_append']) )
, 'bootstrapargs' : ' '.join( map(key_to_ms,['common_keys','remove_keys']) )
}
template = """{{> preamble}}

{{^justfit}}
{{^nobootstrap}}
{{^freefield}}
hl_bootstrap {{#groups}}{{#set}}{{> twopoint}} {{/set}}{{/groups}} {{#hadron_sync_pairs}}{{> justsynckappa}} {{/hadron_sync_pairs}} {{> bootstrapargs}} -noupdate 
{{/freefield}}
{{/nobootstrap}}

{{^freefield}}
hl_effmass {{#hadron_sync_pairs}}{{> justsynckappa}} {{/hadron_sync_pairs}} -r=m -xlims='[0,64]' -fit='[5,16]'
{{/freefield}}
{{#freefield}}
hl_freefield {{#hadron_sync_pairs}}{{> justsynckappa}} {{/hadron_sync_pairs}} -r=m -xlims='[0,64]' -fit='[5,16]'
{{/freefield}}

{{#groups2}}
{{#sub_groups}}
{{^freefield}}
hl_effmass -r={{ratio}} -t={{type}} {{#set}}{{> twopoint}} {{/set}} {{> effmassargs}} {{#fast}}&{{/fast}}
{{/freefield}}
{{#freefield}}
hl_freefield -r={{ratio}} -t={{type}} {{#set}}{{> twopoint}} {{/set}} {{> effmassargs}} {{#fast}}&{{/fast}}
{{/freefield}}
{{/sub_groups}}
{{#fast}}wait{{/fast}}
{{/groups2}}
{{/justfit}}


hl_fit fit/*.json
"""

ms_renderer = ms.Renderer(partials=partials)
create_executable_script( 'do_twopoint', ms_renderer.render(template,data) )
#open('do_twopoint','w+').write(ms_renderer.render(template,data))
