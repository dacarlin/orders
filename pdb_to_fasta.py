from __future__ import print_function 
from sys import argv 
#from core.db.amino_acid import THREE_to_one 

aa1_upper = list("ACDEFGHIKLMNPQRSTVWY")
aa1_lower = list("acdefghiklmnpqrstvwy")
aa3_upper = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
aa3_lower = "ala cys asp glu phe gly his ile lys leu met asn pro gln arg ser thr val trp tyr".split()

def THREE_to_one( one_letter ):
  return dict( zip( aa3_upper, aa1_lower ) )[one_letter]

print( '>%s' % argv[1] ) 
with open( argv[1] ) as fn:
  for line in fn:
    if ' CA ' in line and line.startswith( 'ATOM' ):
      print( THREE_to_one( line.split()[3] ), end='' )
