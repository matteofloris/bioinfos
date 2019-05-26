#!/usr/bin/env python
#title           :MomentEstimator.py
#description     :This python script will predict the age of a recent allele following Slatkin & Rannala, Annu. Rev. Genomics Hum. Genet. 2000.
#author          :Matteo Floris, Univ. of Sassari, Italy, matteo.floris@gmail.com
#date            :20190529
#version         :1.0
#usage           :python MomentEstimator.py
#notes           :
#python_version  :2.7  
#==============================================================================

# first of all, you should download the impute map from http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/
# and then run plink to interpolate the genetic positions for you vcf file: 
# plink --vcf YOUR_VCF.vcf --maf 0.001 --from-bp 18209262 --to-bp 19209262 --cm-map genetic_map_chr19_combined_b37.txt 19 --make-just-bim --out SS3514 --chr 19
# 
# prepara la mappa genetica

gmap = {}
c = 0
f = open ("SS3514.bim", "r")
o = open("SS3514.2.bim", "w")
precedente_prima = 0
precedente_terza = 0
o.write("position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n")
for l in f:
    cols = l.split("\t")
    terza = float(cols[2].strip()) - precedente_terza
    prima = (float(cols[3].strip()) - precedente_prima)/1000000
    try: valore = str(terza/prima)
    except: valore = "0"
    s = "\t".join([  cols[3].strip(), valore, cols[2].strip() ]) + "\n"
    o.write(s)
    c += 1
    if c > 0:
        precedente_prima = float(cols[3].strip())
        precedente_terza = float(cols[2].strip())
f.close()
o.close()

import sys, os, math
from scipy.stats.mstats import gmean

### here load the new genetic map
gmap = {}
f = open("SS3514.2.bim", "r")
for l in f:
    if "position" not in l:
        posi, COMBINED_rate = int(l.split()[0].strip()), float(l.split()[1].strip())
        gmap[posi] = COMBINED_rate
f.close()

### position of your mutation
target_mut = 18709262 

cromosomes_with_mutation = []
snps = {}
alleles = {}

vcf = open("YOUR_VCF.vcf", "r")
for l in vcf:
    if "#" not in l[0]:
        phased = l.split()[9:] ### colonne degli individui nel VCF
        posi = int( l.split()[1] )
        if posi > 18209262 and posi < 19209262:
            if posi == target_mut:
                for sample in phased:
                    ### la mutazione deve essere biallellica
                    ### per ogni sample, memorizza il genotipo
                    genotype = sample.split(":")[0].strip()
                    if genotype in ["0|0", "1|0", "0|1", "1|1"]:
                        cromosomes_with_mutation.append(genotype)
            else:
                if posi not in snps:
                    snps[posi] = []
                    alleles[posi] = [l.split()[3].strip(), l.split()[4].strip()]
                    for sample in phased:
                        genotype = sample.split(":")[0].strip()
                        if genotype in ["0|0", "1|0", "0|1", "1|1"]:
                            snps[posi].append(genotype)
                        else: snps[posi].append("-|-")

N_CHROMOSOMES = len(phased) * 2

# x is the frequency of chromosomes carrying both the mutation and the SNP allele 
# that has a higher frequency of the mutation on its background 

# y is the frequency of chromosomes carrying the mutation and the other SNP allele

x = 0
y = 0

values = []
values_uncorr = []
counts = 0
#### t si calcola per ogni snp:
for posi in snps:
    ref, alt = 0.0, 0.0
    i = 0
    while i < len(cromosomes_with_mutation):
        if cromosomes_with_mutation[i] == "1|0" and snps[posi][i] == "1|0": alt += 1.0
        if cromosomes_with_mutation[i] == "0|1" and snps[posi][i] == "0|1": alt += 1.0
        if cromosomes_with_mutation[i] == "0|1" and snps[posi][i] == "1|0": ref += 1.0
        if cromosomes_with_mutation[i] == "1|0" and snps[posi][i] == "0|1": ref += 1.0
        if cromosomes_with_mutation[i] == "1|0" and snps[posi][i] == "1|1": alt += 1.0
        if cromosomes_with_mutation[i] == "0|1" and snps[posi][i] == "1|1": alt += 1.0
        if cromosomes_with_mutation[i] == "1|0" and snps[posi][i] == "0|0": ref += 1.0
        if cromosomes_with_mutation[i] == "0|1" and snps[posi][i] == "0|0": ref += 1.0

        i += 1

    ref = ref / N_CHROMOSOMES
    alt = alt / N_CHROMOSOMES

    r = 0.487
    x = ref
    y = alt
    if alt > ref: 
        x = alt
        y = ref
    if posi in gmap:
        theta = abs(gmap[posi] - gmap[target_mut])
        #
        if x + y > 0 and posi != target_mut and x != y and x > 0 and y > 0 and theta < 1:
            t = ( 1 / math.log(1 - theta) ) * math.log( (x - y) / (1 - y) )
            counts += 1
            correction = (1/r) * math.log( ( theta*(math.e**r) )/( math.e**r - 1 ) )
            t_new = t - correction
            values.append(t_new)
            values_uncorr.append(t)
            print posi, x, y, theta, t, t_new

result = gmean(values) ## media geometrica
result_uncorr = gmean(values_uncorr) ## media geometrica

print result_uncorr, result

