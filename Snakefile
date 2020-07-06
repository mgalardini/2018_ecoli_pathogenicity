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
k12_gbk = pj(data, 'genome.gbk')
k12_references = pj(data, 'references.txt')
k12_uniprot = pj(data, 'locus_uniprot.tsv')
phenotypes = pj(phenotypes_dir, 'phenotypes.tsv')
phenotypes_ecoli = pj(phenotypes_dir, 'phenotypes_ecoli.tsv')
og_names = pj(data, 'hpi.tsv')
og_names_other = pj(data, 'others.tsv')
phylogroups = pj(data, 'phylogroups.tsv')
phylogroups_ecoli = pj(data, 'phylogroups_ecoli.tsv')
input_file = pj(data, 'inputs.tsv')
input_file_ecoli = pj(data, 'inputs_ecoli.tsv')
strains = [x.split('.')[0] for x in os.listdir(genomes_dir)
           if x.endswith('.fasta')]
strains_ecoli = [x.split()[0] for x in open(input_file_ecoli)
                 if x.split()[0] != 'strain']
genomes = [pj(genomes_dir, x + '.fasta') for x in strains]
survival1 = pj(data, 'liste_souris_NILS46.csv')
survival2 = pj(data, 'liste_souris_NILS9.csv')
ybactin = pj(data, 'yersiniabactin.tsv')
virulence_uniprot = pj(data, 'virulence_genes.faa')
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
report6_template = pj(templates_dir, 'gene_view.ipynb')
report7_template = pj(templates_dir, 'survival.ipynb')
report8_template = pj(templates_dir, 'yersiniabactin.ipynb')
report9_template = pj(templates_dir, 'cog_enrich.ipynb')

# configurable stuff
# edit at will or change these settings with --config
uniref50 = config.get('uniref50', 'db/uniref50')
min_kmer_size = 20

# output directories
out = config.get('out', 'out')
annotations_dir = pj(out, 'annotations')
snps_dir = pj(out, 'snippy')
parsnp_tree_dir = pj(out, 'parsnp')
gubbins_tree_dir = pj(out, 'gubbins')
unitigs_dir = pj(out, 'unitigs')
unitigs_ecoli_dir = pj(out, 'unitigs_ecoli')
roary_dir = pj(out, 'roary')
associations_dir = pj(out, 'associations')
kmer_counts_dir = pj(associations_dir, 'kmer_counts')
kmer_mappings_dir = pj(associations_dir, 'kmer_mappings')
kmer_mappings_ecoli_dir = pj(associations_dir, 'kmer_mappings_ecoli')
refseq_dir = pj(out, 'refseq')
plots_dir = pj(out, 'plots')
staramr_dir = pj(out, 'staramr')
maps_dir = pj(plots_dir, 'maps')
notebooks_dir = config.get('notebooks', 'notebooks')
# simulations
simulated_parsnp_tree_dir = pj(refseq_dir, 'parsnp')
simulated_gubbins_tree_dir = pj(refseq_dir, 'gubbins')
simulated_annotations_dir = pj(refseq_dir, 'annotations')
simulated_roary_dir = pj(refseq_dir, 'roary')

# output files
sift = pj(out, 'sift.tab.gz')
annotations = [pj(annotations_dir, x, x + '.gff')
               for x in strains]
snps = [pj(snps_dir, x, 'snps.vcf')
        for x in strains]
annotated_snps = [pj(snps_dir, x, 'annotated_snps.tsv')
                  for x in strains]
staramrs = [pj(staramr_dir, x, 'summary.tsv')
            for x in strains]
staramr = pj(out, 'staramr.tsv')
parsnp_tree = pj(parsnp_tree_dir, 'parsnp.tree')
polished_parsnp_tree = pj(parsnp_tree_dir, 'tree.nwk')
parsnp_xmfa = pj(parsnp_tree_dir, 'parsnp.xmfa')
parsnp_alignment = pj(parsnp_tree_dir, 'parsnp.fasta')
baps = pj(out, 'baps.csv')
baps_clusters = pj(out, 'baps.tsv')
gubbins_tree = pj(gubbins_tree_dir, 'gubbins.final_tree.tre')
polished_gubbins_tree = pj(gubbins_tree_dir, 'tree.nwk')
gubbins_prefix = pj(gubbins_tree_dir, 'gubbins')
kmers = pj(unitigs_dir, 'unitigs.txt')
kmers_ecoli = pj(unitigs_ecoli_dir, 'unitigs.txt')
sketches_base = pj(out, 'sketches')
sketches = sketches_base + '.msh'
mash_distances = pj(out, 'mash.tsv')
mash_distances_ecoli = pj(out, 'mash_ecoli.tsv')
parsnp_similarities = pj(out, 'parsnp.tsv')
gubbins_similarities = pj(out, 'gubbins.tsv')
parsnp_similarities_ecoli = pj(out, 'parsnp_ecoli.tsv')
gubbins_similarities_ecoli = pj(out, 'gubbins_ecoli.tsv')
roary = pj(roary_dir, 'gene_presence_absence.Rtab')
roarycsv = pj(roary_dir, 'gene_presence_absence.csv')
sampled_pangenome = pj(roary_dir, 'sampled_pangenome.faa')
virulence = pj(out, 'virulence_genes.tsv')
refseq = pj(refseq_dir, 'refseq.tsv')
# associations
# kmers
associations_cont_lmm_kmer = pj(associations_dir, 'associations_cont_lmm_kmer.tsv')
patterns_cont_lmm_kmer = pj(associations_dir, 'patterns_cont_lmm_kmer.txt')
lineage_cont_lmm_kmer = pj(associations_dir, 'lineage_cont_lmm_kmer.tsv')
h2_cont_lmm_kmer = pj(associations_dir, 'h2_cont_lmm_kmer.txt')
# kmers ecoli
associations_cont_lmm_kmer_ecoli = pj(associations_dir, 'associations_cont_lmm_kmer_ecoli.tsv')
patterns_cont_lmm_kmer_ecoli = pj(associations_dir, 'patterns_cont_lmm_kmer_ecoli.txt')
lineage_cont_lmm_kmer_ecoli = pj(associations_dir, 'lineage_cont_lmm_kmer_ecoli.tsv')
h2_cont_lmm_kmer_ecoli = pj(associations_dir, 'h2_cont_lmm_kmer_ecoli.txt')
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
binary_kmer_annotations = pj(associations_dir, 'annotated_binary_cont_lmm_kmer.tsv')
kmer_counts_lmm = [pj(kmer_counts_dir, x + '.txt')
                   for x in strains]
kmer_mappings_lmm = [pj(kmer_mappings_dir, x + '.tsv')
                     for x in strains]
kmer_count_lmm = pj(associations_dir, 'kmer_counts_lmm.txt')
summary_cont_lmm_kmer = pj(associations_dir, 'summary_cont_lmm_kmer.tsv')
summary_lineage_cont_lmm_kmer = pj(associations_dir, 'summary_lineage_cont_lmm_kmer.tsv')
kmer_gene_count_lmm = pj(associations_dir, 'kmer_gene_counts_lmm.txt')
kmer_lineage_gene_count_lmm = pj(associations_dir, 'kmer_lineage_gene_counts_lmm.txt')
# kmers ecoli
filtered_cont_lmm_kmer_ecoli = pj(associations_dir, 'filtered_cont_lmm_kmer_ecoli.tsv')
qq_cont_lmm_kmer_ecoli = pj(associations_dir, 'qq_cont_lmm_kmer_ecoli.png')
annotated_cont_lmm_kmer_ecoli = pj(associations_dir, 'annotated_cont_lmm_kmer_ecoli.tsv')
kmer_mappings_lmm_ecoli = [pj(kmer_mappings_ecoli_dir, x + '.tsv')
                           for x in strains_ecoli]
summary_cont_lmm_kmer_ecoli = pj(associations_dir, 'summary_cont_lmm_kmer_ecoli.tsv')
summary_lineage_cont_lmm_kmer_ecoli = pj(associations_dir, 'summary_lineage_cont_lmm_kmer_ecoli.tsv')
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
# associations - downstream - snps - k12
sampled_ogs_k12 = pj(associations_dir, 'associated_ogs.k12.tsv')
annotated_kmer_k12 = pj(associations_dir, 'annotated_kmer_k12.txt')
associated_kmer_k12 = pj(associations_dir, 'annotated_kmer_k12.bed')
collated_sift = pj(associations_dir, 'nonsyn.k12.tsv')
associations_cont_lmm_nonsyn = pj(associations_dir, 'associations.nonsyn.k12.tsv')
# associations - downstream - genes - ecoli
associated_ogs_ecoli = pj(associations_dir, 'associated_ogs_ecoli.txt')
sampled_ogs_ecoli = pj(associations_dir, 'associated_ogs_ecoli.faa')
uniref_ecoli = pj(associations_dir, 'associated_ogs.faa.uniref50.ecoli.tsv')
unirefnames_ecoli = pj(associations_dir, 'associated_ogs.faa.uniref50.names.ecoli.tsv')
unified_annotations_ecoli = pj(associations_dir, 'associated_ogs_ecoli.final.tsv')
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
eggnog = pj(associations_dir, 'associated.eggnogg.tsv')
eggnog_ref = pj(associations_dir, 'IAI39.tsv')
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
# gene cassette across all strains
report6_nb = pj(notebooks_dir, 'gene_view.ipynb')
report6 = pj(notebooks_dir, 'gene_view.html')
# survival curve
report7_nb = pj(notebooks_dir, 'survival.ipynb')
report7 = pj(notebooks_dir, 'survival.html')
# yersiniabactin
report8_nb = pj(notebooks_dir, 'yersiniabactin.ipynb')
report8 = pj(notebooks_dir, 'yersiniabactin.html')
# enrichment
report9_nb = pj(notebooks_dir, 'cog_enrich.ipynb')
report9 = pj(notebooks_dir, 'cog_enrich.html')
# all reports (minus the ones relying on offline files
reports = [report1, report2,
           report3, report4,
           report5, report6,
           report7, report8]
offline_reports = [report9]

rule download_sift:
  output: sift
  shell:
    'wget -O {output} http://ftp.ebi.ac.uk/pub/databases/mutfunc/mutfunc_v2/ecoli/sift.tab.gz'

rule staramr:
  input: staramrs
  output: staramr
  params: staramr_dir
  shell:
    '''
    head -n 1 $(find {params} -type f -name '*.tsv' | head -n 1) > {output}
    sed -s 1d {input} >> {output}
    '''

rule:
  input: pj(genomes_dir, '{strain}.fasta')
  output: pj(staramr_dir, '{strain}', 'summary.tsv')
  params: pj(staramr_dir, '{strain}')
  shell:
    'staramr search -n 1 -o {params} {input}'

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

rule snps:
  input: snps

rule:
  input:
    ref=k12_gbk,
    genome=pj(genomes_dir, '{strain}.fasta')
  output: pj(snps_dir, '{strain}', 'snps.vcf')
  params:
    d=pj(snps_dir, '{strain}'),
    t='{strain}'
  threads: 1
  shell:
    'mkdir -p /tmp/{params.t} && snippy --force --outdir {params.d} --ref {input.ref} --ctgs {input.genome} --cpus {threads} --ram 8 --tmpdir /tmp/{params.t}'

rule annotate_snps:
  input: annotated_snps
  output: collated_sift 
  params: snps_dir
  shell:
    'python3 src/collate_sift.py {params} > {output}'

rule:
  input:
    annotated_kmer_k12,
    sampled_ogs_k12
  output: associated_kmer_k12
  shell:
    'python3 src/annotated2bed.py {input} > {output}'

rule:
  input:
    s=sift,
    b=associated_kmer_k12,
    v=pj(snps_dir, '{strain}', 'snps.vcf'),
    u=k12_uniprot
  output:pj(snps_dir, '{strain}', 'annotated_snps.tsv') 
  shell:
    'bedtools intersect -a {input.v} -b {input.b} | python3 src/vcf2tsv - {input.s} --ids {input.u} > {output}'

rule associate_nonsyn:
  input:
    phenotype=phenotypes,
    rtab=collated_sift,
    sim=gubbins_similarities,
    baps=phylogroups
  output:
    associations=associations_cont_lmm_nonsyn
  threads: 1
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column killed --pres {input.rtab} --uncompressed --cpu {threads} --lmm --similarity {input.sim} --lineage-clusters {input.baps} > {output}'


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

rule:
  input: parsnp_alignment
  output: baps
  threads: 20
  shell:
    'src/baps.R {input} {output} --cores {threads}'

rule structure:
  input: baps
  output: baps_clusters
  shell:
    '''sed -e 's/,/\t/g' -e 's/.fasta//g' -e 's/.ref//g' {input} | tail -n+2 | awk '{{print $1"\t"$2}}' > {output}'''

rule make_gubbins_tree:
  input: gubbins_tree
  output: polished_gubbins_tree
  shell:
    'src/fix_tree_labels {input} {output}'

rule do_kmers:
  input: input_file
  output: kmers
  params: unitigs_dir
  threads: 40
  shell:
    'rm -rf {params} && unitig-counter -strains {input} -output {params} -nb-cores {threads}'

rule do_kmers_ecoli:
  input: input_file_ecoli
  output: kmers_ecoli
  params: unitigs_ecoli_dir
  threads: 40
  shell:
    'rm -rf {params} && unitig-counter -strains {input} -output {params} -nb-cores {threads}'

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

rule mash_ecoli:
  input:
    m=mash_distances,
    i=input_file_ecoli
  output: mash_distances_ecoli
  shell:
    'python3 src/reduce_squared_matrix.py {input.i} {input.m} > {output}'

rule similarity:
  input:
    parsnp=polished_parsnp_tree,
    gubbins=polished_gubbins_tree
  output:
    pout=parsnp_similarities,
    gout=gubbins_similarities
  shell:
    '''
    python3 src/phylogeny_distance.py --calc-C {input.parsnp} > {output.pout}
    python3 src/phylogeny_distance.py --calc-C {input.gubbins} > {output.gout}
    '''

rule similarity_ecoli:
  input:
    parsnp=parsnp_similarities,
    gubbins=gubbins_similarities,
    strains=input_file_ecoli
  output:
    pout=parsnp_similarities_ecoli,
    gout=gubbins_similarities_ecoli
  shell:
    '''
    python3 src/reduce_squared_matrix.py {input.strains} {input.parsnp} > {output.pout}
    python3 src/reduce_squared_matrix.py {input.strains} {input.gubbins} > {output.gout}
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
    q=virulence_uniprot,
    s=sampled_pangenome
  output:
    virulence
  shell:
    '''echo -e "gene\\tog" > {output};blastp -evalue 1E-4 -query {input.q} -subject {input.s} -outfmt "6 qaccver saccver pident qcovs evalue" | awk '{{if ($4 > 80 && $3 > 50) print $0}}' | awk '{{print $1"\\t"$2}}' >> {output}'''

rule:
  input: annotations
  output: roary
  params: roary_dir
  threads: 40
  shell:
    'rm -rf {params} && roary -p {threads} -f {params} -s -v -g 100000 {input}'

rule pangenome:
  input:
    roary

rule:
  input:
    phenotype=phenotypes,
    kmers=kmers,
    dist=mash_distances,
    sim=gubbins_similarities,
    baps=phylogroups
  output:
    associations=associations_cont_lmm_kmer,
    patterns=patterns_cont_lmm_kmer,
    lineage=lineage_cont_lmm_kmer,
    h2=h2_cont_lmm_kmer
  threads: 40
  params:
    dimensions=3
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column killed --kmers {input.kmers} --uncompressed --max-dimensions {params.dimensions} --lineage --lineage-file {output.lineage} --cpu {threads} --output-patterns {output.patterns} --distance {input.dist} --lmm --similarity {input.sim} --lineage-clusters {input.baps} 2>&1 > {output.associations} | grep \'h^2\' > {output.h2}'

rule associate_kmers:
  input:
    associations_cont_lmm_kmer

rule:
  input:
    phenotype=phenotypes_ecoli,
    kmers=kmers_ecoli,
    dist=mash_distances_ecoli,
    sim=gubbins_similarities_ecoli,
    baps=phylogroups_ecoli
  output:
    associations=associations_cont_lmm_kmer_ecoli,
    patterns=patterns_cont_lmm_kmer_ecoli,
    lineage=lineage_cont_lmm_kmer_ecoli,
    h2=h2_cont_lmm_kmer_ecoli
  threads: 40
  params:
    dimensions=3
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column killed --kmers {input.kmers} --uncompressed --max-dimensions {params.dimensions} --lineage --lineage-file {output.lineage} --cpu {threads} --output-patterns {output.patterns} --distance {input.dist} --lmm --similarity {input.sim} --lineage-clusters {input.baps} 2>&1 > {output.associations} | grep \'h^2\' > {output.h2}'

rule associate_kmers_ecoli:
  input:
    associations_cont_lmm_kmer_ecoli

rule:
  input:
    phenotype=phenotypes,
    rtab=roary,
    dist=mash_distances,
    sim=gubbins_similarities,
    baps=phylogroups
  output:
    associations=associations_cont_lmm_rtab,
    patterns=patterns_cont_lmm_rtab,
    lineage=lineage_cont_lmm_rtab,
    h2=h2_cont_lmm_rtab
  threads: 5
  params:
    dimensions=3
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column killed --pres {input.rtab} --max-dimensions {params.dimensions} --lineage --lineage-file {output.lineage} --cpu {threads} --output-patterns {output.patterns} --distance {input.dist} --lmm --similarity {input.sim} --lineage-clusters {input.baps} 2>&1 > {output.associations} | grep \'h^2\' > {output.h2}'

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
    aclk=associations_cont_lmm_kmer_ecoli,
    pk=patterns_cont_lmm_kmer_ecoli,
  output:
    fclk=filtered_cont_lmm_kmer_ecoli,
  shell:
    '''
    cat <(head -1 {input.aclk}) <(awk -v pval=$(python src/count_patterns.py {input.pk} | tail -n 1 | awk '{{print $2}}') '$4<pval {{print $0}}' {input.aclk}) > {output.fclk}
    ''' 

rule:
  input:
    aclk=associations_cont_lmm_kmer,
  output:
    qclk=qq_cont_lmm_kmer_ecoli
  shell:
    '''
    python src/qq_plot.py {input.aclk} --output {output.qclk}
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
    fclk=filtered_cont_lmm_kmer_ecoli,
    ref=references,
    roary=roary
  output:
    dclk=annotated_cont_lmm_kmer_ecoli
  params:
    gd=genomes_dir,
    pangenome=roarycsv
  shell:
    '''
    python src/annotate_hits.py {input.fclk} {input.ref} {output.dclk} --tmp-prefix /tmp/ --roary {params.pangenome}
    rm -f {params.gd}/*.pac {params.gd}/*.sa {params.gd}/*.amb {params.gd}/*.ann {params.gd}/*.bwt
    '''

rule annotated:
  input:
    aclk=annotated_cont_lmm_kmer,
    lclk=lineage_cont_lmm_kmer,
    r=roary
  output:
    sclk=summary_cont_lmm_kmer,
    slclk=summary_lineage_cont_lmm_kmer
  params:
    minsize=min_kmer_size
  shell:
    '''
    python src/summarise_annotations.py {input.aclk} {input.r} --min-size {params} > {output.sclk}
    python src/summarise_annotations.py {input.aclk} {input.r} --lineage $(head -n 2 {input.lclk} | tail -n 1 | awk '{{print $1}}') --min-size {params} > {output.slclk}
    '''

rule annotated_ecoli:
  input:
    aclk=annotated_cont_lmm_kmer_ecoli,
    lclk=lineage_cont_lmm_kmer_ecoli,
    r=roary
  output:
    sclk=summary_cont_lmm_kmer_ecoli,
    slclk=summary_lineage_cont_lmm_kmer_ecoli
  params:
    minsize=min_kmer_size
  shell:
    '''
    python src/summarise_annotations.py {input.aclk} {input.r} --min-size {params} > {output.sclk}
    python src/summarise_annotations.py {input.aclk} {input.r} --lineage $(head -n 2 {input.lclk} | tail -n 1 | awk '{{print $1}}') --min-size {params} > {output.slclk}
    '''

rule:
  input:
    annotated_cont_lmm_kmer,
    kmers,
    mash_distances
  output:
    binary_kmer_annotations
  params:
    minsize=min_kmer_size
  shell:
    'python src/annotated2binary.py {input} --min-size {params} > {output}'

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

rule:
  input:
    kmer_mappings_lmm
  params:
    kdir=kmer_mappings_dir,
    minsize=min_kmer_size
  output:
    associated_ogs
  shell:
    '''awk '{{if (length($2) >= {params.minsize} && $8 != "") print $8}}' {params.kdir}/*.tsv | sort | uniq -c | sort -n | awk '{{print $2"\\t"$1}}' > {output}'''

rule:
  input:
    k12=k12_genome,
    ogs=sampled_ogs
  output:
    sampled_ogs_k12
  shell:
    '''blastp -query {input.ogs} -subject {input.k12} -evalue 1E-4 -outfmt "6 qseqid sseqid pident qcovs" | awk '{{if ($3 > 70 && $4 > 70) print $1"\\t"$2}}' > {output}'''

rule:
  input:
    kmer=filtered_cont_lmm_kmer,
    ref=k12_references
  output:
     annotated_kmer_k12
  shell:
     'python3 src/annotate_hits.py {input} {output} --id locus_tag' 

rule unitigs_snps:
  input:
    sampled_ogs_k12,
    annotated_kmer_k12

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

rule:
  input:
    fclk=filtered_cont_lmm_kmer_ecoli,
    genome=pj(genomes_dir, '{strain}.fasta'),
    gff=pj(annotations_dir, '{strain}', '{strain}.gff'),
    roary=roary
  params:
    roarycsv
  output:
    pj(kmer_mappings_ecoli_dir, '{strain}.tsv')
  shell:
    'mkdir -p tmp_map_back_{wildcards.strain} && src/map_back {input.fclk} {input.genome} --bwa-algorithm fastmap --print-details --tmp-prefix tmp_map_back_{wildcards.strain} --gff {input.gff} --roary {params} > {output} && rm -rf tmp_map_back_{wildcards.strain}'

rule:
  input:
    kmer_mappings_lmm_ecoli
  params:
    kdir=kmer_mappings_ecoli_dir,
    minsize=min_kmer_size
  output:
    associated_ogs_ecoli
  shell:
    '''awk '{{if (length($2) >= {params.minsize} && $8 != "") print $8}}' {params.kdir}/*.tsv | sort | uniq -c | sort -n | awk '{{print $2"\\t"$1}}' > {output}'''

rule:
  input:
    ogs=associated_ogs_ecoli,
    roary=roary
  params:
    pangenome=roarycsv,
    annotations=annotations_dir
  output:
    sampled_ogs_ecoli
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
    kmer_mappings_lmm,
    sampled_ogs,
    binary_kmer_annotations

rule downstream_ecoli:
  input:
    summary_cont_lmm_kmer_ecoli,
    qq_cont_lmm_kmer_ecoli,
    kmer_mappings_lmm_ecoli,
    sampled_ogs_ecoli

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
    f1=unirefnames,
    f2=sampled_ogs,
    f3=og_names
  output:
    unified_annotations
  shell:
    'src/unify_annotations {input.f1} {input.f2} --names {input.f3} > {output}'

rule annotate_hits:
  input:
    gene_distances,
    unified_annotations

rule:
  input:
    sampled_ogs_ecoli
  params:
    uniref50
  output:
    uniref_ecoli
  threads: 40
  shell:
    'blastp -query {input} -num_threads {threads} -db {params} -outfmt 6 > {output}'

rule:
  input:
    uniref_ecoli
  output:
    unirefnames_ecoli
  shell:
    'src/uniref2genes {input} > {output}'

rule:
  input:
    f1=unirefnames_ecoli,
    f2=sampled_ogs_ecoli,
    f3=og_names
  output:
    unified_annotations_ecoli
  shell:
    'src/unify_annotations {input.f1} {input.f2} --names {input.f3} > {output}'

rule annotate_hits_ecoli:
  input:
    unified_annotations_ecoli

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
  shell:
    'src/power_simulation {input} > {output}'

rule:
  input:
    simulated_roary,
    simulated_gubbins_similarities
  output:
    simulated_power_analysis
  shell:
    'src/power_simulation {input} > {output}'

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
    rtab_gene_count_lmm, 
    phylogroups,
    roary,
    virulence,
    og_names,
    og_names_other
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
    ske=summary_cont_lmm_kmer_ecoli,
    ua=unified_annotations,
    hp=og_names,
    hp1=og_names_other,
    roary=roary
  output:
    report1
  params:
    report=report1_nb,
    roarycsv=roarycsv
  shell:
    'python3 src/run_notebook.py {input.rt} {params.report} -k dists=../{input.gd} -k kmer_hits=../{input.sk} -k kmer_hits_ecoli=../{input.ske} -k names=../{input.ua} -k hpi=../{input.hp} -k others=../{input.hp1} -k rtab=../{params.roarycsv} && jupyter nbconvert --to html --template {input.ht} {params.report} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=6000'

rule:
  input:
    rt=report2_template,
    ht=html_template,
    odds=odds_ratio,
    f=associated_ogs,
    n=unified_annotations,
    s=summary_cont_lmm_kmer,
    hp=og_names,
    hp1=og_names_other
  output:
    report2
  params:
    report2_nb
  shell:
    'python3 src/run_notebook.py {input.rt} {params} -k odds_ratio=../{input.odds} -k filtered=../{input.f} -k names=../{input.n} -k kmer_hits=../{input.s} -k hpi=../{input.hp} -k others=../{input.hp1} && jupyter nbconvert --to html --template {input.ht} {params} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=600'

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
  output:
    report4
  params:
    report4_nb
  shell:
    'python3 src/run_notebook.py {input.rt} {params} -k odds_ratio=../{input.o} -k virulence=../{input.v} -k filtered=../{input.f} -k tnames=../{input.n} -k phenotypes=../{input.p} -k tree=../{input.t} -k rtab=../{input.r} -k spangenome=../{input.spangenome} -k mapping=../{input.mapping} && jupyter nbconvert --to html --template {input.ht} {params} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=600'
 
rule:
  input:
    rt=report5_template,
    ht=html_template,
    p2=phenotypes,
    g=genomes_dir,
    f=summary_cont_lmm_kmer,
    r=roary,
    hp=og_names,
    hp1=og_names_other,
    ba=binary_kmer_annotations,
    samr=staramr
  output:
    report5
  params:
    r=report5_nb,
    s=ecoref_strains,
    p1=ecoref_phenotypes,
  threads: 40
  shell:
    'python3 src/run_notebook.py {input.rt} {params.r} -k cores={threads} -k strains="{params.s}" -k filtered=../{input.f} -k phenotypes="{params.p1}" -k pathogenicity=../{input.p2} -k gdir=../{input.g} -k rtab=../{input.r} -k hpi=../{input.hp} -k others=../{input.hp1} -k binary=../{input.ba} -k staramr=../{input.samr} && jupyter nbconvert --to html --template {input.ht} {params.r} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=6000'

rule:
  input:
    rt=report6_template,
    ht=html_template,
    f=summary_cont_lmm_kmer,
    n=unified_annotations,
    r=roary,
    a=annotations_dir,
    h=og_names,
    o=maps_dir
  output:
    report6
  params:
    r=report6_nb,
    roary=roarycsv,
    s='IAI39'
  shell:
    'python3 src/run_notebook.py {input.rt} {params.r} -k kmer_hits=../{input.f} -k names=../{input.n} -k pangenome=../{params.roary} -k ref_strain={params.s} -k ref_annotation_dir=../{input.a} -k hpi=../{input.h} -k outdir=../{input.o} && jupyter nbconvert --to html --template {input.ht} {params.r} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=600'

rule:
  input:
    rt=report7_template,
    ht=html_template,
    s1=survival1,
    s2=survival2
  output:
    report7
  params:
    report7_nb,
  shell:
    'python3 src/run_notebook.py {input.rt} {params} -k nils46=../{input.s1} -k nils9=../{input.s2} && jupyter nbconvert --to html --template {input.ht} {params} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=600'

rule:
  input:
    rt=report8_template,
    ht=html_template,
    y=ybactin
  output:
    report8
  params:
    report8_nb,
  shell:
    'python3 src/run_notebook.py {input.rt} {params} -k data=../{input.y} && jupyter nbconvert --to html --template {input.ht} {params} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=600'

rule:
  input:
    rt=report9_template,
    ht=html_template,
    ref=eggnog_ref,
    ann=eggnog,
    roary=roary,
    names=unified_annotations,
    hpi=og_names,
    others=og_names_other
  output:
    report9
  params:
    r=report9_nb,
    roary=roarycsv
  shell:
    'python3 src/run_notebook.py {input.rt} {params.r} -k ref=../{input.ref} -k associated=../{input.ann} -k roary=../{params.roary} -k names=../{input.names} -k hpi=../{input.hpi} -k others=../{input.others} && jupyter nbconvert --to html --template {input.ht} {params.r} --ExecutePreprocessor.enabled=True --ExecutePreprocessor.timeout=600'
 
rule all:
  input:
    viz_tree,
    reports

rule all_offline:
  input:
    viz_tree,
    reports,
    offline_reports
