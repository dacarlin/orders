export SCAFFOLD='private/3orf.fasta'
for i in ~/Documents/boxi/etc/3orf_sandbox/*des*pdb; do 
  export alex_f="${i/pdb/fasta}"
  python pdb_to_fasta.py $i > $alex_f
  blastx -query ${SCAFFOLD} -subject $alex_f -outfmt "6 qseq sseq length qstart qend" | python diff.py  
done |  python transcriptic_format.py > kunkel_mutants.csv 

python kunkel_mutagenesis.py 
