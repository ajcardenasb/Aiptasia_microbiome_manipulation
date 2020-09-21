

make.contigs(file=transplants.file, processors=18)
#Output File Names:
#transplants.trim.contigs.fasta
#transplants.trim.contigs.qual
#transplants.contigs.report
#transplants.scrap.contigs.fasta
#transplants.scrap.contigs.qual
#transplants.contigs.groups

##check your fasta file
summary.seqs(fasta=transplants.trim.contigs.fasta)
#Output File Names:
#transplants.trim.contigs.summary

##tabulate number of contigs per group
count.groups(group=transplants.contigs.groups)
#Output File Names:
#transplants.contigs.count.summary

##in case there are short sequencences
screen.seqs(fasta=transplants.trim.contigs.fasta, minlength=250, group=transplants.contigs.groups, processors=15, maxlength=550) 
#Output File Names:
#transplants.trim.contigs.good.fasta
#transplants.trim.contigs.bad.accnos
#transplants.contigs.good.groups


##counting your sequences per group after screen
count.groups(group=transplants.contigs.good.groups)
#Output File Names:
#transplants.contigs.good.count.summary

##remove redundant sequences 
unique.seqs(fasta=transplants.trim.contigs.good.fasta)
#Output File Names:
#transplants.trim.contigs.good.names
#transplants.trim.contigs.good.unique.fasta


##tabulate number of reads per unique sequence
count.seqs(name=transplants.trim.contigs.good.names, group=transplants.contigs.good.groups)
#Output File Names:
#transplants.trim.contigs.good.count_table


##cut off sequences that occur only once
split.abund(fasta=transplants.trim.contigs.good.unique.fasta, count=transplants.trim.contigs.good.count_table, cutoff=1)
#Output File Names:
#transplants.trim.contigs.good.rare.count_table
#transplants.trim.contigs.good.abund.count_table
#transplants.trim.contigs.good.unique.rare.fasta
#transplants.trim.contigs.good.unique.abund.fasta


## number of abundant seq per sample
count.groups(count=transplants.trim.contigs.good.abund.count_table)

#Output File Names:
#transplants.trim.contigs.good.abund.count.summary

## align sequences to silva reference 
align.seqs(fasta=transplants.trim.contigs.good.unique.abund.fasta, reference=~/database/16S_mothur/silva.bacteria.pcr.fasta, processors=16, flip=true)

#Output File Names:
#transplants.trim.contigs.good.unique.abund.align
#transplants.trim.contigs.good.unique.abund.align.report
#transplants.trim.contigs.good.unique.abund.flip.accnos

## get number of unique and aligned sequences.
summary.seqs(fasta=transplants.trim.contigs.good.unique.abund.align, count=transplants.trim.contigs.good.abund.count_table, processors=16)

#Output File Names:
#transplants.trim.contigs.good.unique.abund.summary

## cut the sequences that start before 25298 and end before 34113
screen.seqs(fasta=transplants.trim.contigs.good.unique.abund.align, count=transplants.trim.contigs.good.abund.count_table, start=25298, optimize=end, minlength=200, processors=16)
#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.align
#transplants.trim.contigs.good.unique.abund.bad.accnos
#transplants.trim.contigs.good.abund.good.count_table

## summary of sequences aligned to silva after start and end trimming
summary.seqs(fasta=transplants.trim.contigs.good.unique.abund.good.align, count=transplants.trim.contigs.good.abund.good.count_table)
#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.summary

## get number of sequences per group
count.groups(count=transplants.trim.contigs.good.abund.good.count_table)
#Output File Names:
#transplants.trim.contigs.good.abund.good.count.summary

## filter sequences --> get rid of uninformative column
filter.seqs(fasta=transplants.trim.contigs.good.unique.abund.good.align, vertical=T,trump=.)
transplants.filter
transplants.trim.contigs.good.unique.abund.good.filter.fasta

## get rid of redundancy due to filtering above
unique.seqs(fasta=transplants.trim.contigs.good.unique.abund.good.filter.fasta, count=transplants.trim.contigs.good.abund.good.count_table)
#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.count_table
#transplants.trim.contigs.good.unique.abund.good.filter.unique.fasta

## pre-clustering, allowing for up to 2 nt (1 per each 100 nt) difference between sequences
pre.cluster(fasta=transplants.trim.contigs.good.unique.abund.good.filter.unique.fasta, count=transplants.trim.contigs.good.unique.abund.good.filter.count_table, diffs=2)


#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.fasta
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.count_table


## get number of removed sequences per group
count.groups(count=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.count_table)
#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.count.summary


## identify chimeras from sequences, Vsearch (5 min). Needs either reference or counts
chimera.vsearch(fasta=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.fasta, count=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.count_table, processors= 16, dereplicate=t)

#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.pick.count_table
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.chimeras
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.accnos


## remove chimeras from sequences
remove.seqs(fasta=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.fasta, accnos=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.accnos)

#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.fasta

## get number of removed sequences per group
count.groups(count=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.pick.count_table)

#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.pick.count.summary

## classify sequences --> only want bacterial 16s, not 18s or from mitochondria etc.
classify.seqs(fasta=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.fasta, count=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference=~/database/16S_mothur/gg_13_5_99.fasta, taxonomy=~/database/16S_mothur/gg_13_5_99.pds.tax, cutoff=60, processors=16) 

#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pds.wang.taxonomy
transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pds.wang.tax.summary
#

## remove 18s, mitochondria, chloroplast etc. sequences
remove.lineage(fasta=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.fasta, count=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=chloroplast-Chloroplast-Mitochondria-mitochondria-unknown-Archaea-Eukaryota)

#Output File Names:
#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.fasta
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table


## get number of removed sequences per group
count.groups(count=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.pick.pick.count.summary


# cluster sequences into OTUs, taxlevel from kingdom (1) to species(7). need to have a higher cutoff here, can change it in next step. ### distance clustering was taking more than 24 h. ## For large databases like this one, mothur can split the process in two: cluster.split (took about 6 h) to generate the .file (clusters based on taxonomy?) and cluster.split using that file (clustering based on distances)

cluster.split(fasta=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.fasta, count=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, splitmethod=classify, taxlevel=4, cutoff=0.03, processors=16)

#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.dist
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.sensspec

# cutoff of 0.03 for OTUs
make.shared(list=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)

#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared



# classify OTUs
classify.otu(list=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=0.03)

#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.tax.summary

# get representative OTU
get.oturep(fasta=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.fasta, count=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, list=transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, label=0.03, method=abundance)

#Output File Names:
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.count_table
#transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.fasta



cp transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared transplants.final.OTU_table

cp transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy  transplants.final.taxonomy

grep -v ">" transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.fasta  | sed 's/-//g' > temp_transplants.final.fasta

awk '{print ">"$1}' transplants.trim.contigs.good.unique.abund.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy  | tail -n +2 > headers_fasta

paste headers_fasta temp_transplants.final.fasta |  tr "\t" "\n" > transplants.final.fasta

scp cardena@10.74.186.160:/home/cardena/Rubens_Miseq/Only50Samples/transplants.final* .

scp cardena@10.74.186.160:/home/cardena/Rubens_Miseq/Only50Samples/transplants.*.count.summary .

 
remove.seqs(accnos=contaminants, fasta=transplants.final.fasta)




# classify OTUs to latest SILVA release
classify.seqs(fasta=transplants.final.fasta, template=~/databases/Ribosomal/SILVA_138/silva.nr_v138.align, taxonomy=~/databases/Ribosomal/SILVA_138/silva.nr_v138.tax)

