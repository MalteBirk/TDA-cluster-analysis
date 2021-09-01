#conda activate TDA

#mkdir TDA
#cd TDA/

#get tropodithietic acid cluster from MiBIG
#https://mibig.secondarymetabolites.org/repository/BGC0000932/index.html#r1c1

mkdir referenceTDA
wget https://mibig.secondarymetabolites.org/repository/BGC0000932/BGC0000932.gbk -P referenceTDA/


#steal gbk splitter
#remember to credit https://github.com/davised/seqIO_extract
mkdir scripts
wget https://raw.githubusercontent.com/davised/seqIO_extract/master/seqIO_extract.py -P scripts
chmod 755 scripts/*

#extract fasta sequence of plasmid from genbank
scripts/seqIO_extract.py --fna referenceTDA/BGC0000932.gbk all > referenceTDA/BGC0000932.fna
faidx referenceTDA/BGC0000932.fna BGC0000932.1:103793-108389 > referenceTDA/tdaEA.fna 
 
#extract genes
echo -e "tdaA\ntdaB\ntdaC\ntdaD\ntdaE\ntdaF\n" > geneList.txt
scripts/seqIO_extract.py --outname gene --ffn --matchtype exact --searchname gene referenceTDA/BGC0000932.gbk geneList.txt > referenceTDA/BGC0000932.ffn  
rm geneList.txt
rm referenceTDA/BGC0000932.fna.fai

#make single fastas
#maybe superflous?
#cd referenceTDA
#cat BGC0000932.ffn  | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fasta")} print $0 > filename }'
#cd ..

#get maltes stuff
#must be root from maltes folder
#/home/s173696/Phaeobacter_WGS/Out_final_subset/
#find . -name '*ffn' -exec cp {} ../../../milst/TDA/Malte_backup/ \;
#find . -name "*fasta*" -exec cp {} ~/TDA/Malte_backup_genomes/ \;

#get all CDS
ncbi-genome-download -p 40 -F cds-fasta -l all  --flat-output -g Phaeobacter,Epibacterium,Pseudovibrio,Paracoccus,Stappia -o cds_fasta bacteria
gunzip cds_fasta/*

#get all genomes, only for naming purposes
ncbi-genome-download -p 40 -F fasta -l all  --flat-output -g Phaeobacter,Epibacterium,Pseudovibrio,Paracoccus,Stappia -o genomes bacteria
gunzip genomes/*

mkdir cds_fasta_rename/
cd cds_fasta
for i in *; do echo $i; suge=$(head -n1 ../genomes/$(echo $i | sed 's/cds_from_//g')); suge2=$(echo $suge | sed -r 's/>([a-zA-Z0-9._]* [a-zA-Z0-9]* [a-zA-Z0-9_\-.]* [a-zA-Z0-9_\-.]* [a-zA-Z0-9_\-.]*) .*/\1/g' | sed 's/ /_/g'); echo $suge2; bbrename.sh in=$i out=../cds_fasta_rename/$i prefix=$suge2 addprefix=t; done
cd ..	


#transfer to cds_fasta_rename
cp ../TDA_old/Malte_backup/ Malte_backup/ -r   
cp Malte_backup/* cds_fasta_rename/
rename 's/ffn/fna/' cds_fasta_rename/*  


#blast all genes
mkdir tdaBlast
ls cds_fasta_rename/ | parallel -j 20 'blastn  -perc_identity 50 -qcov_hsp_perc 50 -subject cds_fasta_rename/{} -query referenceTDA/BGC0000932.ffn -outfmt 6 -word_size 5 > tdaBlast/{.}.blast'

#blast against genomes
#ls genomes/   | parallel -j 20 'blastn  -perc_identity 50 -qcov_hsp_perc 50 -subject genomes/{} -query referenceTDA/tdaEA_full.fna -outfmt 6 -word_size 5 > tda_genome_blast/{.}.blast'                            

#have a count
find tdaBlast/ -mindepth 1 -type f,d   \( -empty -o -printf 'non-' \) -printf 'empty '   \( -type f -printf 'files' -o -printf 'directories' \)   -printf '\n' | sort | uniq -c

#   199 empty files
#	45 non-empty files 

#remove empties
find tdaBlast/ -size  0 -print -delete 

for i in tdaBlast/*; do if [[ $(wc -l <$i) -gt 4 ]]; then echo $i;  cat $i; fi; done 
for i in tdaBlast/*; do if [[ $(wc -l <$i) -lt 4 ]]; then echo $i;  cat $i; fi; done  



mkdir genes
mkdir aligns

for tda in tdaA tdaB tdaC tdaD tdaE; do
	echo $tda
	for blast in tdaBlast/*; do
		gene=$(grep $tda $blast | cut -f2 | cut -f2 -d'|')
		base=${blast##*/}
		cds="${base%.blast}.fna"
		faidx --regex $gene cds_fasta_rename/$cds
	done > genes/$tda.fasta

	mafft --quiet --adjustdirection genes/$tda.fasta > aligns/$tda.aln
done


for cds in cds_fasta_rename/*; do faidx --regex rpoB $cds;	done > genes/rpoB.fasta

mkdir prokka_out/
cd Malte_backup_genomes/
for i in *; do echo $i;j=$(basename $i .contigs.rename.fasta); echo $j; prokka $i -out ../prokka_out/$j -cpus 20 --quiet --locustag $j --prefix $j --usegenus --genus Phaeobacter --force; done

cd ../genomes/
for i in *; do echo $i;j=$(basename $i _genomic.fna); echo $j; prokka $i -out ../prokka_out/$j -cpus 20 --quiet --locustag $j --prefix $j --usegenus --genus Phaeobacter --force; done

roary -e --mafft -p 20 -cd 100 -i 95 -o roary-g100-i99 `find -name *gff`

for i in *; do  sed -i 's/ /_/g' $i/$i.ffn; faidx -f --regex "50S_ribosomal_protein_L11$" $i/$i.ffn ;done > ../genes/L11.ffn