import os

# shortcuts
pj = os.path.join

# input directories
data = config.get('data', 'data')
genomes_dir = pj(data, 'genomes')
phenotypes_dir = pj(data, 'phenotypes')
templates_dir = config.get('templates', 'templates')

# input files
k12_genome = pj(data, 'genome.faa')
phenotypes = pj(phenotypes_dir, 'phenotypes.tsv')
input_file = pj(data, 'inputs.tsv')
strains = [x.split('.')[0] for x in os.listdir(genomes_dir)
           if x.endswith('.fasta')]
genomes = [pj(genomes_dir, x + '.fasta') for x in strains]
rna_samples_file = pj(data, 'rna_samples.tsv')
rna_samples = {x.rstrip().split('\t')[0] for x in open(rna_samples_file)
               if x.rstrip().split('\t')[0] != 'strain'}
rna_reads = [(x.rstrip().split('\t')[0],
	      x.rstrip().split('\t')[1],
              x.rstrip().split('\t')[2])
             for x in open(rna_samples_file)
             if x.rstrip().split('\t')[0] != 'strain']
html_template = pj(templates_dir, 'html.tpl')
# urls for online data
ecoref_phenotypes = 'https://evocellnet.github.io/ecoref/data/phenotypic_data.tsv'
ecoref_strains = 'https://evocellnet.github.io/ecoref/data/strains.tsv'
# reports
report1_template = pj(templates_dir, 'plots.ipynb')
report2_template = pj(templates_dir, 'odds_ratio.ipynb')
report3_template = pj(templates_dir, 'simulations.ipynb')
report4_template = pj(templates_dir, 'virulence_genes.ipynb')
report5_template = pj(templates_dir, 'chemical.ipynb')

# configurable stuff
# edit at will or change these settings with --config
uniref50 = config.get('uniref50', 'db/uniref50')
power_ogs = 'pks2,group_2650'
simulated_power_ogs = 'group_7955,fabG'
trimmomatic_dir = config.get('trimmomatic_dir',
                             'software/trimmomatic-0.36-5/share/trimmomatic')

# output directories
out = config.get('out', 'out')
annotations_dir = pj(out, 'annotations')
parsnp_tree_dir = pj(out, 'parsnp')
gubbins_tree_dir = pj(out, 'gubbins')
roary_dir = pj(out, 'roary')
associations_dir = pj(out, 'associations')
kmer_counts_dir = pj(associations_dir, 'kmer_counts')
kmer_mappings_dir = pj(associations_dir, 'kmer_mappings')
refseq_dir = pj(out, 'refseq')
rna_dir = pj(out, 'rna')
plots_dir = pj(out, 'plots')
notebooks_dir = config.get('notebooks', 'notebooks')
# simulations
simulated_parsnp_tree_dir = pj(refseq_dir, 'parsnp')
simulated_gubbins_tree_dir = pj(refseq_dir, 'gubbins')
simulated_annotations_dir = pj(refseq_dir, 'annotations')
simulated_roary_dir = pj(refseq_dir, 'roary')

# inputs but generated manualy during the analysis
virulence = pj(out, 'virulence_genes.tsv')
# output files
annotations = [pj(annotations_dir, x, x + '.gff')
               for x in strains]
parsnp_tree = pj(parsnp_tree_dir, 'parsnp.tree')
polished_parsnp_tree = pj(parsnp_tree_dir, 'tree.nwk')
parsnp_xmfa = pj(parsnp_tree_dir, 'parsnp.xmfa')
parsnp_alignment = pj(parsnp_tree_dir, 'parsnp.fasta')
baps = pj(out, 'baps.tsv')
gubbins_tree = pj(gubbins_tree_dir, 'gubbins.final_tree.tre')
polished_gubbins_tree = pj(gubbins_tree_dir, 'tree.nwk')
gubbins_prefix = pj(gubbins_tree_dir, 'gubbins')
kmers = pj(out, 'kmers.gz')
sketches_base = pj(out, 'sketches')
sketches = sketches_base + '.msh'
mash_distances = pj(out, 'mash.tsv')
parsnp_similarities = pj(out, 'parsnp.tsv')
gubbins_similarities = pj(out, 'gubbins.tsv')
roary = pj(roary_dir, 'gene_presence_absence.Rtab')
roarycsv = pj(roary_dir, 'gene_presence_absence.csv')
sampled_pangenome = pj(roary_dir, 'sampled_pangenome.faa')
transcripts = [pj(annotations_dir, x, x + '.transcripts')
               for x in rna_samples]
indexes = [pj(annotations_dir, x, x + '.index')
           for x in rna_samples]
rna_counts = [pj(rna_dir, x[0], x[1], 'abundance.tsv')
              for x in rna_reads]
de_genes = pj(rna_dir, 'overall.csv')
fold_changes = pj(rna_dir, 'fold_changes.tsv')
refseq = pj(refseq_dir, 'refseq.tsv')
# associations
# kmers
associations_cont_lmm_kmer = pj(associations_dir, 'associations_cont_lmm_kmer.tsv')
patterns_cont_lmm_kmer = pj(associations_dir, 'patterns_cont_lmm_kmer.txt')
lineage_cont_lmm_kmer = pj(associations_dir, 'lineage_cont_lmm_kmer.tsv')
h2_cont_lmm_kmer = pj(associations_dir, 'h2_cont_lmm_kmer.txt')
# pangenome
associations_cont_lmm_rtab = pj(associations_dir, 'associations_cont_lmm_rtab.tsv')
patterns_cont_lmm_rtab = pj(associations_dir, 'patterns_cont_lmm_rtab.txt')
lineage_cont_lmm_rtab = pj(associations_dir, 'lineage_cont_lmm_rtab.tsv')
h2_cont_lmm_rtab = pj(associations_dir, 'h2_cont_lmm_rtab.txt')
# associations - downstream
references = pj(associations_dir, 'references.txt')
# kmers
filtered_cont_lmm_kmer = pj(associations_dir, 'filtered_cont_lmm_kmer.tsv')
qq_cont_lmm_kmer = pj(associations_dir, 'qq_cont_lmm_kmer.png')
annotated_cont_lmm_kmer = pj(associations_dir, 'annotated_cont_lmm_kmer.tsv')
kmer_counts_lmm = [pj(kmer_counts_dir, x + '.txt')
                   for x in strains]
kmer_mappings_lmm = [pj(kmer_mappings_dir, x + '.tsv')
                     for x in strains]
kmer_count_lmm = pj(associations_dir, 'kmer_counts_lmm.txt')
summary_cont_lmm_kmer = pj(associations_dir, 'summary_cont_lmm_kmer.tsv')
summary_lineage_cont_lmm_kmer = pj(associations_dir, 'summary_lineage_cont_lmm_kmer.tsv')
kmer_gene_count_lmm = pj(associations_dir, 'kmer_gene_counts_lmm.txt')
kmer_lineage_gene_count_lmm = pj(associations_dir, 'kmer_lineage_gene_counts_lmm.txt')
# pangenome
filtered_cont_lmm_rtab = pj(associations_dir, 'filtered_cont_lmm_rtab.tsv')
qq_cont_lmm_rtab = pj(associations_dir, 'qq_cont_lmm_rtab.png')
rtab_gene_count_lmm = pj(associations_dir, 'rtab_gene_counts_lmm.txt')
# associations - downstream - genes
gene_distances = pj(associations_dir, 'gene_distances.tsv')
associated_ogs = pj(associations_dir, 'associated_ogs.txt')
sampled_ogs = pj(associations_dir, 'associated_ogs.faa')
uniref = pj(associations_dir, 'associated_ogs.faa.uniref50.tsv')
unirefnames = pj(associations_dir, 'associated_ogs.faa.uniref50.names.tsv')
# restricted strains analysis
restricted_covariates = pj(associations_dir, 'restricted_covariates.tsv')
associations_restricted = pj(associations_dir, 'associations_restricted.tsv')
patterns_restricted = pj(associations_dir, 'patterns_restricted.txt')
lineage_restricted = pj(associations_dir, 'lineage_restricted.tsv')
h2_restricted = pj(associations_dir, 'h2_restricted.txt')
filtered_restricted = pj(associations_dir, 'filtered_restricted.tsv')
qq_restricted = pj(associations_dir, 'qq_restricted.png')
annotated_restricted = pj(associations_dir, 'annotated_restricted.tsv')
summary_restricted = pj(associations_dir, 'summary_restricted.tsv')
# offline annotation
eggnog = pj(associations_dir, 'associated_ogs.faa.emapper.annotations')
unified_annotations = pj(associations_dir, 'associated_ogs.final.tsv')
# power analysis
odds_ratio = pj(associations_dir, 'odds_ratio.tsv')
power_analysis = pj(associations_dir, 'power_analysis.tsv')
simulated_parsnp_tree = pj(simulated_parsnp_tree_dir, 'parsnp.tree')
simulated_parsnp_xmfa = pj(simulated_parsnp_tree_dir, 'parsnp.xmfa')
simulated_parsnp_alignment = pj(simulated_parsnp_tree_dir, 'parsnp.fasta')
simulated_gubbins_tree = pj(simulated_gubbins_tree_dir, 'gubbins.final_tree.tre')
simulated_polished_gubbins_tree = pj(simulated_gubbins_tree_dir, 'tree.nwk')
simulated_gubbins_prefix = pj(simulated_gubbins_tree_dir, 'gubbins')
simulated_gubbins_similarities = pj(refseq_dir, 'gubbins.tsv')
simulated_roary = pj(simulated_roary_dir, 'gene_presence_absence.Rtab')
simulated_power_analysis = pj(refseq_dir, 'power_analysis.tsv')
# plots
viz_tree = pj(plots_dir, '4_tree.pdf')
# reports
# associtions results
report1_nb = pj(notebooks_dir, 'plots.ipynb')
report1 = pj(notebooks_dir, 'plots.html')
# odds ratio
report2_nb = pj(notebooks_dir, 'odds_ratio.ipynb')
report2 = pj(notebooks_dir, 'odds_ratio.html')
# power analysis
report3_nb = pj(notebooks_dir, 'simulations.ipynb')
report3 = pj(notebooks_dir, 'simulations.html')
# virulence genes
report4_nb = pj(notebooks_dir, 'virulence_genes.ipynb')
report4 = pj(notebooks_dir, 'virulence_genes.html')
# comparison with ecoref phenotypes
report5_nb = pj(notebooks_dir, 'chemical.ipynb')
report5 = pj(notebooks_dir, 'chemical.html')
# all reports
reports = [report1, report2,
           report3, report4,
           report5]

rule annotate:
  input: annotations

rule:
  input:
    ref=k12_genome,
    genome=pj(genomes_dir, '{strain}.fasta')
  output: pj(annotations_dir, '{strain}', '{strain}.gff')
  params: pj(annotations_dir, '{strain}')
  threads: 5
  shell:
    'prokka --outdir {params} --force --prefix {wildcards.strain} --addgenes --locustag {wildcards.strain} --mincontiglen 200 --genus Escherichia -species coli --strain {wildcards.strain} --proteins {input.ref} --cpus {threads} {input.genome}'

rule:
  output: parsnp_tree
  params:
    genomes_dir=genomes_dir,
    outdir=parsnp_tree_dir,
    focus=pj(genomes_dir, 'IAI01.fasta')
  threads: 40
  shell:
    'parsnp -d {params.genomes_dir} -r {params.focus} -p {threads} -o {params.outdir} -v -c'

rule make_tree:
  input: parsnp_tree
  output: polished_parsnp_tree
  shell:
    'src/fix_tree_labels {input} {output}'

rule:
  input: parsnp_tree
  output: parsnp_alignment
  params: parsnp_xmfa
  shell:
    'harvesttools -x {params} -M {output}'

rule:
  input: parsnp_alignment
  output: gubbins_tree
  params: gubbins_prefix
  threads: 40
  shell:
    'run_gubbins.py --verbose --threads {threads} {input} --prefix {params}'

rule structure:
  input: parsnp_alignment
  output: baps
  threads: 20
  shell:
    'src/baps.R {input} {output} --cores {threads}'

rule make_gubbins_tree:
  input: gubbins_tree
  output: polished_gubbins_tree
  shell:
    'src/fix_tree_labels {input} {output}'

rule do_kmers:
  input: input_file
  output: kmers
  shell:
    'fsm-lite -l {input} -t tmp.txt -m 9 -M 100 -s 1 -v | gzip > {output}'

rule:
  output: sketches
  threads: 5
  params:
    gdir=genomes_dir,
    base=sketches_base
  shell:
    'mash sketch -p {threads} -s 10000 -o {params.base} {params.gdir}/*.fasta'

rule mash:
  input: sketches
  output: mash_distances
  threads: 5
  shell:
    'mash dist -p {threads} {input} {input} | square_mash > {output}'

rule similarity:
  input:
    parsnp=polished_parsnp_tree,
    gubbins=polished_gubbins_tree
  output:
    pout=parsnp_similarities,
    gout=gubbins_similarities
  shell:
    '''
    python src/phylogeny_distance.py --calc-C {input.parsnp} > {output.pout}
    python src/phylogeny_distance.py --calc-C {input.gubbins} > {output.gout}
    '''

rule:
  input:
    roary=roary
  params:
    pangenome=roarycsv,
    annotations=annotations_dir
  output:
    sampled_pangenome
  shell:
    'src/sample_pangenome {params.pangenome} {params.annotations} --focus-strain IAI39 --focus-strain IAI01 > {output}'

rule:
  input:
    pj(annotations_dir, '{strain}', '{strain}.gff'),
    roary
  output: pj(annotations_dir, '{strain}', '{strain}.transcripts')
  params:
    anndir=annotations_dir,
    roarycsv=roarycsv
  shell:
    'src/pangenome2transcripts {params.roarycsv} {params.anndir} {wildcards.strain} --reference IAI39 > {output}'

rule:
  input: annotations
  output: roary
  params: roary_dir
  threads: 40
  shell:
    'rm -rf {params} && roary -p {threads} -f {params} -s -v -g 100000 {input}'

rule pangenome:
  input:
    roary, sampled_pangenome, transcripts

rule:
  input:
    pj(annotations_dir, '{strain}', '{strain}.transcripts')
  output:
    pj(annotations_dir, '{strain}', '{strain}.index')
  shell:
    'kallisto index -i {output} {input}'

rule:
  input:
    index=pj(annotations_dir, '{strain}', '{strain}.index'),
    rf=rna_samples_file,
  output:
    pj(rna_dir, '{strain}', '{replica}', 'abundance.tsv')
  params:
    odir1=pj(rna_dir, '{strain}'),
    odir=pj(rna_dir, '{strain}', '{replica}'),
    tdir=trimmomatic_dir,
    average=130,
    sd=70
  shell:
    'mkdir -p {params.odir1} && trimmomatic SE $(awk \'{{if ($1 == "{wildcards.strain}" && $2 == "{wildcards.replica}") print $3}}\' {input.rf}) /dev/stdout ILLUMINACLIP:{params.tdir}/adapters/TruSeq3-SE.fa:2:30:10 | kallisto quant -b 100 -i {input.index} -o {params.odir} --bias --single -l {params.average} -s {params.sd} /dev/stdin'

rule:
  input:
    rna_counts,
    rf=rna_samples_file
  output:
    de_genes
  params:
    rna_dir
  threads: 40
  shell:
    'Rscript src/deseq.R {input.rf} {params} --cores {threads} --pvalue 0.01 --foldchange 0.0 --reference IAI55'

rule:
  input:
    de_genes
  output:
    fold_changes
  params:
    rna_dir
  shell:
    'src/merge_fold_changes $(find {params} -type f -name \'*.csv\' ! -wholename \'{input}\') > {output}'

rule transcriptomics:
  input:
    fold_changes

rule:
  input:
    phenotype=phenotypes,
    kmers=kmers,
    dist=mash_distances,
    sim=gubbins_similarities
  output:
    associations=associations_cont_lmm_kmer,
    patterns=patterns_cont_lmm_kmer,
    lineage=lineage_cont_lmm_kmer,
    h2=h2_cont_lmm_kmer
  threads: 40
  params:
    dimensions=3
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column killed --kmers {input.kmers} --max-dimensions {params.dimensions} --lineage --lineage-file {output.lineage} --cpu {threads} --output-patterns {output.patterns} --distance {input.dist} --lmm --similarity {input.sim} 2>&1 > {output.associations} | grep \'h^2\' > {output.h2}'

rule associate_kmers:
  input:
    associations_cont_lmm_kmer

rule:
  input:
    phenotype=phenotypes,
    rtab=roary,
    dist=mash_distances,
    sim=gubbins_similarities
  output:
    associations=associations_cont_lmm_rtab,
    patterns=patterns_cont_lmm_rtab,
    lineage=lineage_cont_lmm_rtab,
    h2=h2_cont_lmm_rtab
  threads: 5
  params:
    dimensions=3
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column killed --pres {input.rtab} --max-dimensions {params.dimensions} --lineage --lineage-file {output.lineage} --cpu {threads} --output-patterns {output.patterns} --distance {input.dist} --lmm --similarity {input.sim} 2>&1 > {output.associations} | grep \'h^2\' > {output.h2}'

rule associate_pangenome:
  input:
    associations_cont_lmm_rtab

rule:
  input:
    aclk=associations_cont_lmm_kmer,
    aclr=associations_cont_lmm_rtab,
    pk=patterns_cont_lmm_kmer,
    pr=patterns_cont_lmm_rtab
  output:
    fclk=filtered_cont_lmm_kmer,
    fclr=filtered_cont_lmm_rtab
  shell:
    '''
    cat <(head -1 {input.aclk}) <(awk -v pval=$(python src/count_patterns.py {input.pk} | tail -n 1 | awk '{{print $2}}') '$4<pval {{print $0}}' {input.aclk}) > {output.fclk}
    cat <(head -1 {input.aclr}) <(awk -v pval=$(python src/count_patterns.py {input.pr} | tail -n 1 | awk '{{print $2}}') '$4<pval {{print $0}}' {input.aclr}) > {output.fclr}
    ''' 

rule:
  input:
    aclk=associations_cont_lmm_kmer,
    aclr=associations_cont_lmm_rtab
  output:
    qclk=qq_cont_lmm_kmer,
    qclr=qq_cont_lmm_rtab
  shell:
    '''
    python src/qq_plot.py {input.aclk} --output {output.qclk}
    python src/qq_plot.py {input.aclr} --output {output.qclr}
    '''

rule:
  output: references
  params:
    genomes_dir,
    annotations_dir 
  shell:
    'src/prepare_ref_file {params} --ref IAI01 --ref IAI39 > {output}'

rule:
  input:
    fclk=filtered_cont_lmm_kmer,
    ref=references,
    roary=roary
  output:
    dclk=annotated_cont_lmm_kmer
  params:
    gd=genomes_dir,
    pangenome=roarycsv
  shell:
    '''
    python src/annotate_hits.py {input.fclk} {input.ref} {output.dclk} --tmp-prefix /tmp/ --roary {params.pangenome}
    rm -f {params.gd}/*.pac {params.gd}/*.sa {params.gd}/*.amb {params.gd}/*.ann {params.gd}/*.bwt
    '''

rule:
  input:
    aclk=annotated_cont_lmm_kmer,
    lclk=lineage_cont_lmm_kmer
  output:
    sclk=summary_cont_lmm_kmer,
    slclk=summary_lineage_cont_lmm_kmer
  params:
    minsize=20
  shell:
    '''
    python src/summarise_annotations.py {input.aclk} --min-size {params} > {output.sclk}
    python src/summarise_annotations.py {input.aclk} --lineage $(head -n 2 {input.lclk} | tail -n 1 | awk '{{print $1}}') --min-size {params} > {output.slclk}
    '''

rule:
  input:
    fclk=filtered_cont_lmm_kmer,
    genome=pj(genomes_dir, '{strain}.fasta')
  output:
    pj(kmer_counts_dir, '{strain}.txt')
  shell:
    'src/map_back {input.fclk} {input.genome} --bwa-algorithm fastmap | wc -l > {output}'

rule:
  input:
    kmer_counts_lmm
  output:
    kmer_count_lmm
  params:
    kmer_counts_dir 
  shell:
    'for i in $(ls {params}); do echo $(basename $i .txt) $(cat {params}/$i); done > {output}'

rule:
  input:
    summary=summary_cont_lmm_kmer,
    pangenome=roary
  output:
    kmer_gene_count_lmm
  params:
    roarycsv
  shell:
    'src/summary2genes {input.summary} {params} > {output}'

rule:
  input:
    summary=summary_lineage_cont_lmm_kmer,
    pangenome=roary
  output:
    kmer_lineage_gene_count_lmm
  params:
    roarycsv
  shell:
    'src/summary2genes {input.summary} {params} > {output}'

rule:
  input:
    summary=filtered_cont_lmm_rtab,
    pangenome=roary
  output:
    rtab_gene_count_lmm
  params:
    roarycsv
  shell:
    'src/summary2genes {input.summary} {params} > {output}'

rule:
  input:
    rtab=filtered_cont_lmm_rtab,
    kmer=summary_cont_lmm_kmer
  output:
    associated_ogs
  shell:
    ''
    '(tail -n+2 {input.rtab} | awk \'{{print $1}}\' && tail -n+2 {input.kmer} | awk \'{{print $1}}\') | sort | uniq > {output}'

rule:
  input:
    ogs=associated_ogs,
    roary=roary
  params:
    pangenome=roarycsv,
    annotations=annotations_dir
  output:
    sampled_ogs
  shell:
    'src/sample_pangenome {params.pangenome} {params.annotations} --focus-strain IAI39 --focus-strain IAI01 --groups {input.ogs} > {output}'
    

rule downstream:
  input:
    qq_cont_lmm_kmer,
    qq_cont_lmm_rtab,
    kmer_count_lmm,
    kmer_gene_count_lmm,
    kmer_lineage_gene_count_lmm,
    rtab_gene_count_lmm,
    sampled_ogs

rule:
  input:
    roary,
    ogs=associated_ogs
  params:
    pangenome=roarycsv,
    annotations=annotations_dir
  output:
    gene_distances
  shell:
    'python3 src/ogs_distances {params.pangenome} {params.annotations} --groups {input.ogs} > {output}'

rule:
  input:
    sampled_ogs
  params:
    uniref50
  output:
    uniref
  threads: 40
  shell:
    'blastp -query {input} -num_threads {threads} -db {params} -outfmt 6 > {output}'

rule:
  input:
    uniref
  output:
    unirefnames
  shell:
    'src/uniref2genes {input} > {output}'

rule:
  input:
    sampled_ogs
  output:
    eggnog
  shell:
    'echo "OFFLINE ANNOTATION: please submit {input} to http://eggnogdb.embl.de/#/app/emapper and download the output to {output}"'

rule:
  input:
    unirefnames,
    sampled_ogs,
    eggnog
  output:
    unified_annotations
  shell:
    'src/unify_annotations {input} > {output}'

rule:
  input:
    fclk=filtered_cont_lmm_kmer,
    genome=pj(genomes_dir, '{strain}.fasta'),
    gff=pj(annotations_dir, '{strain}', '{strain}.gff'),
    roary=roary
  params:
    roarycsv
  output:
    pj(kmer_mappings_dir, '{strain}.tsv')
  shell:
    'mkdir -p tmp_map_back_{wildcards.strain} && src/map_back {input.fclk} {input.genome} --bwa-algorithm fastmap --print-details --tmp-prefix tmp_map_back_{wildcards.strain} --gff {input.gff} --roary {params} > {output} && rm -rf tmp_map_back_{wildcards.strain}'

rule annotate_hits:
  input:
    gene_distances,
    unified_annotations,
    kmer_mappings_lmm

rule:
  input:
    p=phenotypes,
    s=summary_cont_lmm_kmer,
    r=roary
  params:
    r=roarycsv,
    t=2
  output:
    restricted_covariates
  shell:
    'src/restrict_strains {params.r} {input.s} {input.p} --threshold {params.t} > {output}'

rule:
  input:
    phenotype=phenotypes,
    kmers=kmers,
    dist=mash_distances,
    sim=gubbins_similarities,
    cov=restricted_covariates
  output:
    associations=associations_restricted,
    patterns=patterns_restricted,
    lineage=lineage_restricted,
    h2=h2_restricted
  threads: 40
  params:
    dimensions=3
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column killed --kmers {input.kmers} --max-dimensions {params.dimensions} --lineage --lineage-file {output.lineage} --cpu {threads} --output-patterns {output.patterns} --distance {input.dist} --lmm --similarity {input.sim} --covariates {input.cov} --use-covariates 2 2>&1 > {output.associations} | grep \'h^2\' > {output.h2}'

rule:
  input:
    ar=associations_restricted,
    pr=patterns_restricted
  output:
    filtered_restricted
  shell:
    '''cat <(head -1 {input.ar}) <(awk -v pval=$(python src/count_patterns.py {input.pr} | tail -n 1 | awk '{{print $2}}') '$4<pval {{print $0}}' {input.ar}) > {output}'''

rule:
  input:
    associations_restricted
  output:
    qq_restricted
  shell:
    'python src/qq_plot.py {input} --output {output}'

rule:
  input:
    fr=filtered_restricted,
    ref=references,
    roary=roary
  output:
    annotated_restricted
  params:
    gd=genomes_dir,
    pangenome=roarycsv
  shell:
    '''
    python src/annotate_hits.py {input.fr} {input.ref} {output} --tmp-prefix /tmp/ --roary {params.pangenome}
    rm -f {params.gd}/*.pac {params.gd}/*.sa {params.gd}/*.amb {params.gd}/*.ann {params.gd}/*.bwt
    '''

rule:
  input:
    annotated_restricted
  output:
    summary_restricted
  shell:
    'python src/summarise_annotations.py {input} > {output}'

rule restricted:
  input:
    filtered_restricted,
    qq_restricted,
    summary_restricted

rule:
  input:
    associations_cont_lmm_rtab,
    phenotypes,
    roary
  output:
    odds_ratio
  shell:
    'src/get_odds_ratio {input} > {output}'

rule:
  input:
    k12_genome
  output:
    refseq
  params:
    r=refseq_dir,
    a=simulated_annotations_dir
  threads: 40
  shell:
    '''
    ncbi-genome-download bacteria --format fasta --output-folder {params.r} --parallel {threads} --retries 10 --genus "Escherichia coli" -vvv --metadata-table {output} --human-readable --assembly-level complete
    mkdir -p {params.r}/unpacked
    src/unpack_genomes {params.r}/human_readable {params.r}/unpacked > unpack.sh && bash unpack.sh && rm -rf {params.r}/human_readable {params.r}/refseq
    for genome in $(ls {params.r}/unpacked); do prokka --outdir {params.a}/$genome --force --prefix $genome --addgenes --locustag $genome --mincontiglen 200 --genus Escherichia -species coli --strain $genome --proteins {input} --cpus {threads} {params.r}/unpacked/$genome; done
    '''

rule:
  input: refseq
  output: simulated_roary
  params:
    r=simulated_roary_dir,
    a=simulated_annotations_dir
  threads: 40
  shell:
    'rm -rf {params.r} && roary -p {threads} -f {params.r} -s -v -g 100000 {params.a}/*/*.gff'

rule:
  input: refseq
  output: simulated_parsnp_tree
  params:
    genomes_dir=pj(refseq_dir, 'unpacked'),
    outdir=simulated_parsnp_tree_dir,
    focus=pj(refseq_dir, 'unpacked', 'genome001')
  threads: 40
  shell:
    'parsnp -d {params.genomes_dir} -r {params.focus} -p {threads} -o {params.outdir} -v -c'

rule:
  input: simulated_parsnp_tree
  output: simulated_parsnp_alignment
  params: simulated_parsnp_xmfa
  shell:
    'harvesttools -x {params} -M {output}'

rule:
  input: simulated_parsnp_alignment
  output: simulated_gubbins_tree
  params: simulated_gubbins_prefix
  threads: 40
  shell:
    'run_gubbins.py --verbose --threads {threads} {input} --prefix {params}'

rule:
  input: simulated_gubbins_tree
  output: simulated_polished_gubbins_tree
  shell:
    'src/fix_tree_labels {input} {output}'

rule:
  input:
    simulated_polished_gubbins_tree
  output:
    simulated_gubbins_similarities
  shell:
    '''
    python src/phylogeny_distance.py --calc-C {input} > {output}
    '''

rule:
  input:
    roary,
    gubbins_similarities
  output:
    power_analysis
  params:
    power_ogs
  shell:
    'src/power_simulation {input} {params} > {output}'

rule:
  input:
    simulated_roary,
    simulated_gubbins_similarities
  output:
    simulated_power_analysis
  params:
    simulated_power_ogs
  shell:
    'src/power_simulation {input} {params} > {output}'

rule simulations:
  input:
    odds_ratio,
    power_analysis,
    simulated_power_analysis

rule:
  input:
    phenotypes,
    kmer_count_lmm,
    kmer_gene_count_lmm,
    rtab_gene_count_lmm
  output:
    viz_tree
  shell:
    'bash src/plot_trees.sh'

rule:
  input:
    rt=report1_template,
    ht=html_template,
    gd=gene_distances,
    sk=summary_cont_lmm_kmer,
    ua=unified_annotations
  output:
    report1
  params:
    report1_nb
  shell:
    'python3 src/run_notebook.py {input.rt} {params} -k dists=../{input.gd} -k kmer_hits=../{input.sk} -k names=../{input.ua} && jupyter nbconvert --to html --template {input.ht} {params} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=6000'

rule:
  input:
    rt=report2_template,
    ht=html_template,
    odds=odds_ratio,
    f=filtered_cont_lmm_rtab
  output:
    report2
  params:
    report2_nb
  shell:
    'python3 src/run_notebook.py {input.rt} {params} -k odds_ratio=../{input.odds} -k filtered=../{input.f} && jupyter nbconvert --to html --template {input.ht} {params} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=600'

rule:
  input:
    rt=report3_template,
    ht=html_template,
    s1=power_analysis,
    s2=simulated_power_analysis
  output:
    report3
  params:
    report3_nb
  shell:
    'python3 src/run_notebook.py {input.rt} {params} -k simulations1=../{input.s1} -k simulations2=../{input.s2} && jupyter nbconvert --to html --template {input.ht} {params} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=600'

rule:
  input:
    rt=report4_template,
    ht=html_template,
    o=odds_ratio,
    v=virulence,
    f=filtered_cont_lmm_rtab,
    n=unified_annotations,
    p=phenotypes,
    t=polished_gubbins_tree,
    r=roary,
    spangenome=sampled_pangenome,
    mapping=kmer_mappings_dir,
    fold=fold_changes
  output:
    report4
  params:
    report4_nb
  shell:
    'python3 src/run_notebook.py {input.rt} {params} -k odds_ratio=../{input.o} -k virulence=../{input.v} -k filtered=../{input.f} -k tnames=../{input.n} -k phenotypes=../{input.p} -k tree=../{input.t} -k rtab=../{input.r} -k spangenome=../{input.spangenome} -k mapping=../{input.mapping} -k fold_changes=../{input.fold} && jupyter nbconvert --to html --template {input.ht} {params} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=600'
 
rule:
  input:
    rt=report5_template,
    ht=html_template,
    p2=phenotypes,
    g=genomes_dir,
    f=summary_cont_lmm_kmer,
    r=roary
  output:
    report5
  params:
    r=report5_nb,
    s=ecoref_strains,
    p1=ecoref_phenotypes,
  threads: 20
  shell:
    'python3 src/run_notebook.py {input.rt} {params.r} -k cores={threads} -k strains="{params.s}" -k filtered=../{input.f} -k phenotypes="{params.p1}" -k pathogenicity=../{input.p2} -k gdir=../{input.g} -k rtab=../{input.r} && jupyter nbconvert --to html --template {input.ht} {params.r} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=600'
 
rule plots:
  input:
    viz_tree,
    reports
    
