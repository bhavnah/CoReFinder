# CoReFinder

The pipeline uses SOAP mappings of paired-end (PE) reads to identify collapsed and repeated/non-collapsed regions in genome assemblies based on the -r parameter (i.e. how  to  report repeat hits). The values used for -r are as follows:

-r 0: map uniquely to the genome

-r 1: map reads to more that one location in the genome, but only one hit is reported at random

-r 2: map reads to all locations in the genome

The collapsed and repeated/non-collapsed regions are defined as follows:

Repeated/non-collapsed region:

If we map reads to the reference with –r0, there will be regions where reads do not map. This could happen due to of a number of reasons: the region was not sequenced, and therefore not present in the reads, the region in the genome where we want to map the reads contain N's or tiny repetitive regions or more likely that they represent large duplicated regions in the assembly. Because the reads can map to more than one place, using -r0 will discard them. If we map the same reads to the reference assembly and identify regions empty with –r 0 but covered with –r1 then these are likely to be repeated in the assembly. Also, if we map reads with -r2, a duplicated, non-collapsed region will roughly have half of the reads mapped at -r1 compared to -r2, a triplicated region should have roughly a third, etc...

Collapsed region:

Coverage higher than the average coverage with -r 0, -r 1 and -r 2 will be repetitive/collapsed regions of the assembly. This is because these reads have nowhere else to map in the genome, so they will pile up at the collapsed portion of the assembly.

# Prerequisites

1. PE reads: A coverage of at least 30x is required for more accurate estimation of candidate regions. The reads should also be assessed for quality; if after quality control the read depth or length are considerably reduced then the data should not be used.
2. Reference genome: The pipeline has been tested using the same reads that were used to build the corresponding reference genome. Since presence/absence variation occurs in individuals of the same species, using reads from one individual to look for collapsed or repeated regions in the reference genome of another individual will not yield accuarate results and should be avoided.
3. SOAPaligner/soap2
4. SAMtools (lastest version)
5. bedtools (at least version v2.25.0, so that genomecov -d reports positions with 0 reads aligning)
6. R (any version)

# Methodology

## 1. SOAP alignments of reads

The steps below are described for -r 0 mappings, but the same steps must be repeated for -r 1 and -r 2 such that for each dataset/library there are 3 sets of mappings (-r 0,1,2).

1. Format reference genome sequence:

```
./2bwt-builder genome.fasta
```

2. Align reads as follows:

```
soap -p 4 -r 0 -v 2 -m 0 -x 1000 -a R1.fastq -b R2.fastq -o r0.paired.soap -2 r0.single.soap
```

The above will run soap2 with 4 CPUs, an insert size from 0 to 1000 using the 2 input files R1.fastq and R2.fastq. We will only use the paired soap output for the next steps.

3. Prepare soap output for bedtools:

```
soap2sam.pl r0.paired.soap > r0.paired.sam
samtools faidx genome.fasta
samtools view -bt genome.fasta.fai r0.paired.sam > r0.paired.bam
samtools sort r0.paired.bam r0.paired.sorted
samtools index r0.paired.sorted.bam
```

At this stage you should have one paired.sorted.bam file for each of -r 0, -r 1, and -r 2.

4. Calculate per base coverage using bedtools:

```
bedtools genomecov -d -ibam r0.paired.sorted.bam > r0.cov
bedtools genomecov -d -ibam r1.paired.sorted.bam > r1.cov
bedtools genomecov -d -ibam r2.paired.sorted.bam > r2.cov
bedtools unionbedg -i r0.cov r1.cov r2.cov -header -names r0 r1 r2 > all.cov
```

In case unionbedg does not work, use bash paste:

```
paste r0.cov r1.cov r2.cov | cut -f1-3,6,9 > all.cov
```

At this stage the output will look like this:

```
chr1	1	0	0	0
chr1	2	0	0	0
chr1	3	0	0	0
chr1	4	0	0	0
chr1	5	0	0	0
chr1	6	0	0	0
chr1	7	0	0	0
chr1	8	1	1	1
chr1	9	1	1	1
```

Since you have mapped reads to a reference genome which consists of chromosomes/scaffolds, you should split all.cov by chromosome/scaffold as tthis will reduce the run time. For e.g. extract the results for chr1:

```
grep -w chr1 all.cov > chr1.cov
```
 
Once you have extracted all chromosomes/scaffolds of interest, add a header each resulting file so that it looks like this:

```
chrom	pos	r0	r1	r2
chr1	1	0	0	0
chr1	2	0	0	0
chr1	3	0	0	0
chr1	4	0	0	0
chr1	5	0	0	0
chr1	6	0	0	0
chr1	7	0	0	0
chr1	8	1	1	1
chr1	9	1	1	1
```

5. Identify collapsed and repetitive/non-collapsed regions:

Run the script on each chromosome.cov file:

```
Rscript compareCoverage.R chr1.cov chr1
```

The resulting output (tab-delimited) will look like this:

```
region.number   start   end     region.length   r2      r1      r0      copy number     region.type
1       18380   18466   87      19      7       0       2.714   rnc
2       18541   18598   58      17      8       0       2.125   rnc
3       31260   31875   616     23      7       0       3.286   rnc
4       31943   32085   143     7       3       0       2.333   rnc
5       32093   32616   524     20      5       0       4       rnc
6       32669   34310   1642    21      7       0       3       rnc
7       37098   37162   65      26      26      26      1       coll
8       39069   39195   127     28      28      28      1       coll
9       53162   53278   117     31      31      31      1       coll
10      72267   72366   100     29      29      28      1       coll
```

# Visualisation (optional)

Plot the median coverage (using user specified bin sizes) for -r 0, -r 1, -r 2 in regions on interest across the genome:

```
python plotMeanCoverage.py chr1.cov 1000 1 100000
```

The above will plot the mean coverage on chr1 in regions 1 bp to 100 Kbp using a bin size of 1000.
