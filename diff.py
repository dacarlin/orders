from sys import stdin
import screed 
import os 
import re 

def rc(seq):
  c = [{'a':'t', 't':'a', 'g':'c', 'c':'g', 'A':'T', 'C':'G', 'T':'A', 'G':'C', }[i] for i in str(seq)] 
  return ''.join(c)[::-1]

aa = { 'ata':'i', 'atc':'i', 'att':'i', 'atg':'m', 'aca':'t', 'acc':'t', 'acg':'t', 'act':'t', 'aac':'n', 'aat':'n', 'aaa':'k', 'aag':'k', 'agc':'s', 'agt':'s', 'aga':'r', 'agg':'r', 'cta':'l', 'ctc':'l', 'ctg':'l', 'ctt':'l', 'cca':'p', 'ccc':'p', 'ccg':'p', 'cct':'p', 'cac':'h', 'cat':'h', 'caa':'q', 'cag':'q', 'cga':'r', 'cgc':'r', 'cgg':'r', 'cgt':'r', 'gta':'v', 'gtc':'v', 'gtg':'v', 'gtt':'v', 'gca':'a', 'gcc':'a', 'gcg':'a', 'gct':'a', 'gac':'d', 'gat':'d', 'gaa':'e', 'gag':'e', 'gga':'g', 'ggc':'g', 'ggg':'g', 'ggt':'g', 'tca':'s', 'tcc':'s', 'tcg':'s', 'tct':'s', 'ttc':'f', 'ttt':'f', 'tta':'l', 'ttg':'l', 'tac':'y', 'tat':'y', 'taa':'_', 'tag':'_', 'tgc':'c', 'tgt':'c', 'tga':'_', 'tgg':'w'}
ecoli_codon = { 'g':'ggc', 'a':'gcg', 'v':'gtg', 'f':'ttt', 'e':'gaa', 'd':'gat', 'n':'aac', 'h':'cat', 'p':'ccg', 'q':'cag', 'w':'tgg', 'y':'tat', 'i':'att', 'm':'atg', 'c':'tgc', 'k':'aaa', 'l':'ctg', 'r':'cgt', 't':'acc', 's':'agc'}

for line in stdin:
  sseq, qseq, length, qstart, qend = line.split() 
  if int( length ) > 200:
    zipped = zip( sseq, qseq ) 
    diff = '+'.join( '%s%d%s' % ( j[0], i, j[1] ) for i, j in enumerate( zipped ) if j[0] != j[1] ) 

    for record in screed.open( os.getenv( 'SCAFFOLD' ) ):
      trim = record.sequence[ int( qstart ) - 1 : int( qend ) ].lower()  
      codon = [ trim[ i:i+3 ] for i in range( 0, len(trim), 3 ) ] 

      for i in diff.split( '+' ):
        x = i[0].lower()
        y = int( i[1:-1] )
        z = i[-1].lower()
        
        if x == aa[ codon[ y ].lower() ]: #sanity check 
          codon[ y ] = ecoli_codon[ z ].upper()
        else:
          raise ValueError('mismatch')
  
    oligos = re.sub( r'([atcg]{15})[atcg]{0,}([atcg]{15})', r'\1 \2', rc( ''.join( codon ) ) )
    print diff.lower(), oligos
