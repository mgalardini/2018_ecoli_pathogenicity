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
report_template = pj(templates_dir, 'plots.ipynb')
html_template = pj(templates_dir, 'html.tpl')

# binaries/databases
# edit at will or change these settings with --config
uniref50 = config.get('uniref50', 'db/uniref50')

# output directories
out = config.get('out', 'out')
annotations_dir = pj(out, 'annotations')
parsnp_tree_dir = pj(out, 'parsnp')
gubbins_tree_dir = pj(out, 'gubbins')
roary_dir = pj(out, 'roary')
associations_dir = pj(out, 'associations')
kmer_counts_dir = pj(associations_dir, 'kmer_counts')
kmer_mappings_dir = pj(associations_dir, 'kmer_mappings')
plots_dir = pj(out, 'plots')
notebooks_dir = config.get('notebooks', 'notebooks')

# output files
annotations = [pj(annotations_dir, x, x + '.gff')
               for x in strains]
parsnp_tree = pj(parsnp_tree_dir, 'parsnp.tree')
polished_parsnp_tree = pj(parsnp_tree_dir, 'tree.nwk')
parsnp_xmfa = pj(parsnp_tree_dir, 'parsnp.xmfa')
parsnp_alignment = pj(parsnp_tree_dir, 'parsnp.fasta')
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
# plots
viz_tree = pj(plots_dir, '4_tree.pdf')
report_nb = pj(notebooks_dir, 'plots.ipynb')
report = pj(notebooks_dir, 'plots.html')

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
  threads: 20
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
  threads: 20
  shell:
    'run_gubbins.py --verbose --threads {threads} {input} --prefix {params}'

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

rule pangenome:
  input: annotations
  output: roary
  params: roary_dir
  threads: 20
  shell:
    'rm -rf {params} && roary -p {threads} -f {params} -s -v -g 100000 {input}'

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
    rm {params.gd}/*.pac {params.gd}/*.sa {params.gd}/*.amb {params.gd}/*.ann {params.gd}/*.bwt
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
    genome=pj(genomes_dir, '{strain}.fasta')
  output:
    pj(kmer_mappings_dir, '{strain}.tsv')
  shell:
    'src/map_back {input.fclk} {input.genome} --bwa-algorithm fastmap --print-details > {output}'

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
    rm {params.gd}/*.pac {params.gd}/*.sa {params.gd}/*.amb {params.gd}/*.ann {params.gd}/*.bwt
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
    rt=report_template,
    ht=html_template,
    gd=gene_distances,
    sk=summary_cont_lmm_kmer,
    ua=unified_annotations
  output:
    report
  params:
    report_nb
  shell:
    'python src/run_notebook.py {input.rt} {params} -k dists=../{input.gd} -k kmer_hits=../{input.sk} -k names=../{input.ua} && jupyter nbconvert --to html --template {input.ht} {params} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=600'
 
rule plots:
  input:
    viz_tree,
    report
    
