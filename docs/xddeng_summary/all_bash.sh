
#First filtration
#Primary filteration "Presence - Absence"	Nucleotide - Nucleotide searches	≥ 0.7	≥ 100nt	≤ 1e-06	
mmseqs search $querydb $targetdb $resultsdb ./tmp/   --threads $THREADS --search-type 3 -s 1  --min-seq-id 0.70  --min-aln-len 100 -e 0.000001 --max-seqs 50  --max-accept 1 --split-memory-limit 374G
mmseqs convertalis $querydb $targetdb $resultsdb MMseqs2_C"$i"_vs_"$s".tsv --search-type 3 --format-output "query,target,evalue,pident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,mismatch"


#Diamond
#Primary filteration "Presence - Absence"	Translated nucleotide searches	≥ 0.5	≥ 33aa	≤ 1e-05
diamond makedb --in cated_sk100_reps.faa  --threads $THREADS  --db dmnd/"$i".dmnd
diamond blastx -b4.0  -p $THREADS -k 1 --query-cover $qcov  --evalue $evalue -q $input_seq -d $targetdb -o step_"$step"_vs_"$i".tsv -f 6 qseqid sseqid qstart qend sstart send evalue bitscore length pident mismatch slen #--subject-cover $scov

