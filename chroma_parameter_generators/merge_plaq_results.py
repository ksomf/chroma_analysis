import sys
import json

output_name = sys.argv[1]
files       = sys.argv[2:]

jsos = [ json.load( open( f, 'r' ) ) for f in files ]
output = { k:v for k,v in jsos[0].items() if not isinstance( k, dict ) }
dict_key=None
for k,v in filter( lambda (a,b): isinstance( b, dict ), jsos[0].items() ):
  dict_key=k
  update = { k:{ k2:v2 for j in jsos for k2,v2 in j[k].items() } }
  output.update( update )

if 'keys' in jsos[0]:
  output['keys'] = [ v for j in jsos for v in j['keys'] ]
  b_out = open( '.'.join( output_name.split('.')[:-1] ) + '.dat', 'wb' )
  for f in files:
    b_file = '.'.join( f.split('.')[:-1] ) + '.dat'
    b_out.write( open( b_file, 'rb' ).read() )

mean = lambda xs: sum(xs)/len(xs)
for k,v in filter( lambda (a,b): 'mean' in a, jsos[0].items() ):
   #output.update( { k:mean([ x for j in jsos for x in [j[k]]*len(j[dict_key].keys()) ]) } )
   output.update( { k:mean([ x for j in jsos for x in [j[k]]*len(j['keys']) ]) } )

open( output_name, 'w+' ).write( json.dumps( output, indent=2 ) )
