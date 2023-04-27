#!/bash

echo "This scripts is to concatenate all abundance tsv to form count matrix table"
echo "*** Please be in the main directory which contains /samplename/abundance.tsv****"
#echo "read in all abundance tsv files"
abundance_tsv=`ls */abundance.tsv`


#echo "concatenate them"
#first create a concatenated file
paste -d' ' ${abundance_tsv} | awk 'NR>1' | awk '{for(i=5;i<=101;i+=5) printf "%s ",$i ;for(i=102;i<=119;i++) {printf "%s ",$i} ;print ""}' > .counts_only.tsv

#get the first abundance file 
first_abundance_tsv=`echo ${abundance_tsv} | cut --delimiter " " --fields 1`

#get transcript names
cat ${first_abundance_tsv} | awk 'NR>1' |  awk '{print $1}' > .transcript_id.tsv 

#paste together transcript id and counts table
paste -d " " .transcript_id.tsv .counts_only.tsv  > .transcript_counts_raw.tsv

#get the header name 
header="transcripts `ls -d *chr*/ | sed 's|/||g'` "

#add to final file
echo  ${header} > transcript_counts.csv

#add the transcript counts file to final oupit
cat .transcript_counts_raw.tsv  | tail -n+2>> transcript_counts.csv

#reformatting
sed -e 's/\s\+/,/g' transcript_counts.csv > .transcript_counts.new.tsv
sed '2,$s/.$//' .transcript_counts.new.tsv > transcript_counts.csv 

rm .transcript_counts_raw.tsv .counts_only.tsv .transcript_id.tsv .transcript_counts.new.tsv

