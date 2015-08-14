export SCAFFOLD='private/2azn.fasta'
for i in designs/*pdb; do 
  export alex_f="${i/pdb/fasta}"
  python pdb_to_fasta.py $i > $alex_f
  blastx -query ${SCAFFOLD} -subject $alex_f -outfmt "6 qseq sseq length qstart qend" | python diff.py  
done | python transcriptic_format.py 
