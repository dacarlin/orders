from sys import stdin
import screed 
import os 
from core.db.amino_acid import aa, ecoli_codon, rc
import re 

for line in stdin:
  sseq, qseq, length, qstart, qend = line.split() 
  if int( length ) > 200:
    zipped = zip( sseq, qseq ) 
    diff = '+'.join( '%s%d%s' % ( j[0], i, j[1] ) for i, j in enumerate( zipped ) if j[0] != j[1] ) 

    for record in screed.open( os.getenv( 'SCAFFOLD' ) ):
      trim = record.sequence[ int( qstart ) - 1 : int( qend ) ].lower()  
      codon = [ trim[ i:i+3 ] for i in range( 0, 9900, 3 ) ] 

      for i in diff.split( '+' ):
        x = i[0].lower()
        y = int( i[1:-1] )
        z = i[-1].lower()
        
        if x == aa[ codon[ y ].lower() ]: #sanity check 
          codon[ y ] = ecoli_codon[ z ].upper()
        else:
          raise ValueError('mismatch')
  
    oligos = re.sub( r'([atcg]{15})[atcg]{0,}([atcg]{15})', r'\1 \2', ''.join( codon ) )
    print diff.lower(), oligos
