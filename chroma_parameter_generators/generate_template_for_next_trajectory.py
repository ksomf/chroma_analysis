#!/bin/env python
import os
import sys
import re
import pystache as ms
import random

def get_dim(file):
  fh = open(file,'r')
  for line in fh:
    if '<lx>' in line:
      nx = int(re.search('<lx>(.*)</lx>',line).group(1))
      ny = int(re.search('<ly>(.*)</ly>',line).group(1))
      nz = int(re.search('<lz>(.*)</lz>',line).group(1))
      nt = int(re.search('<lt>(.*)</lt>',line).group(1))
      return nx,ny,nz,nt
    elif 'latt_size' in line:
      dims = map( int, re.search('<latt_size>(.*)</latt_size>',line).group(1).split(' ') )
      return dims

params = {
  'TMP_DIR' : sys.argv[1]
, 'GEOMETRY' : sys.argv[4]
}
sources = int(sys.argv[2])
name    = sys.argv[3]

filelist   = open('filelist'  ,'r').read().splitlines()
donelist   = open('donelist'  ,'r').read().splitlines()
inproglist = open('inproglist','r').read().splitlines()
errorlist  = open('errorlist' ,'r').read().splitlines()

random.shuffle(filelist)
for file in filelist:
  file_name = file.split('/')[-1]
  if file_name[-1].isdigit():
    conf_name, conf_num = file_name.split('.lime')
    conf_gen = 1
  else:
    conf_name, conf_gen, conf_num, _ = file_name.split('.')
  conf_gen = int(conf_gen)
  conf_num = int(conf_num)
  ident = '{}.{}.{:05d}'.format(conf_name, conf_gen, conf_num)
  ongoing = filter( lambda x: x.startswith( ident ), donelist + inproglist + errorlist )
  ongoing_sources = map( lambda x: int(x.replace(ident,'').split('_')[1]), ongoing )
  if len(ongoing) < sources:
    source = next(x for x in range(sources) if x not in ongoing_sources)
    nx,ny,nz,nt = get_dim(file)
    
    random.seed(conf_gen+(38*conf_num)+(100*source+8))
    sx = random.randint(0,nx-1)
    sy = random.randint(0,ny-1)
    sz = random.randint(0,nz-1)
    st = random.randint(0,nt-1)
        
    tag = '_'.join( [ ident, str(source), str(sx), str(sy), str(sz), str(st) ] )
    open('inproglist','a').write(tag + '\n')
    params.update( { 'TRAJ':tag, 'CONF':file, 'NX':str(nx), 'NY':str(ny)
                   , 'NZ':str(nz), 'NT':str(nt), 'SX':str(sx)
                    , 'SY':str(sy), 'SZ':str(sz), 'ST':str(st) } )
    for template in filter( lambda x:x.endswith('.ms') and x.startswith('input_'+name), os.listdir('.') ):
      out_file = '.'.join( template.split('.')[:-2] ) + '_' + tag + '.' + template.split('.')[-2]
      open( params['TMP_DIR'] + '/' + out_file, 'w' ).write( ms.render( open(template,'r').read(), params ) )
    
    print tag
    exit( 0 )
    
exit(1)
