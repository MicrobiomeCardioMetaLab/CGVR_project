#### main codes for GMR project 
# 1.Viral genome identification
# 2.Compare with other database
# 3.Pan-genome generation
# 4.SNV calling
# 5.Genomic dissimilarity

##############################################
####### 1. Viral genome identification #######
##############################################

# We utilized a single sample viral genome identification pipeline, code see https://github.com/RasmussenLab/phamb

##############################################
######## 2. Compare with other db #########
##############################################
# 1. create a blast+ database:
makeblastdb -in merged.fa -dbtype nucl -out merged_db

# 2. perform all-vs-all blastn of sequences
blastn -query merged.fa \
       -db merge_db \
       -outfmt '6 std qlen slen' \
       -max_target_seqs 10000 \
       -out merged_blastn.tsv

# 3. calculate pairwise ANI by combining local alignments between sequence pairs
python /public/home/bioinfo_wang/00_software/checkv/scripts/anicalc.py \
       -i merged_blastn.tsv \
       -o merged_blastn_ani.tsv

# 4. perform CD-HIT-like clustering using the MIUVIG recommended-parameters (95% ANI + 85% AF)
python /public/home/bioinfo_wang/00_software/checkv/scripts/aniclust.py \
       --fna merged.fa \
       --ani merged_blastn_ani.tsv \
       --out merged_blastn_cluster.tsv \
       --min_ani 95 --min_tcov 85 --min_qcov 0

##############################################
######## 3. Pangenome generation #############
##############################################

roary -e --mafft -i 90 -cd 90 -s gff/*.gff -f roary/ -p 2

##############################################
############ 4. SNV Calling ##################
##############################################

snippy --reference ${ref} --ctgs ${query} --outdir ${outdir} --prefix ${sample} --cpus 1 --mincov 5

##############################################
####### 5. Genomic dissimilarity #############
##############################################

dfs = []
for file in os.listdir('snippy/'):
    vcf_file = pd.read_csv(os.path.join(path, 'snippy', file, file + '.csv'))
    snp_file = vcf_file[vcf_file['TYPE'] == 'snp']     
    snp_file['snv_site'] = snp_file['CHROM'].astype(str) + '__' + snp_file['POS'].astype(str) + '__' + snp_file['REF']
    df = snp_file[['snv_site', 'ALT']]
    df.rename(columns={'ALT': file}, inplace = True)
    dfs.append(df)
result_df = reduce(lambda left, right: pd.merge(left, right, on='snv_site', how='outer'), dfs)

filter_df = result_df.set_index('snv_site')
mask = filter_df.count(axis=1) / filter_df.shape[1] >= 0.2
dat = filter_df[mask]
print(dat)
a = dat.apply(pd.value_counts, axis=1)
b = a.isnull().sum(axis=1).tolist()

# the case of 2 kind of base
e = [i for i, x in enumerate(b) if x == 2]
snv2 = dat.iloc[e, :]

dis_matrix = []
for a in list(snv2.columns):
    print(a)
    for b in list(snv.columns):
        dis_matrix.append(calcu_dis(a, b))
df = pd.concat(dis_matrix, axis=0)
dat = df.pivot_table(index='sample1', columns='sample2', values='dis')
dat.to_csv(path + 'dis.xls', sep='\t')
