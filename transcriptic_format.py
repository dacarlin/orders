import sys 
import uuid

print 'mutant_label,sequencing_primer,oligo_label,sequence,scale,purification'

for line in sys.stdin:
  spl = line.strip( ).split( )
  mutant = spl[0]
  oligos = spl[1:]

  for index, oligo in enumerate( oligos ):
    print '%s,T7pro,%s,%s,25nm,standard' % ( mutant, uuid.uuid4(), oligo ) 
