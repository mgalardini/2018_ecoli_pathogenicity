import os

# shortcuts
pj = os.path.join

# input directories
data = config.get('data', 'data')
genomes_dir = pj(data, 'genomes')
phenotypes_dir = pj(data, 'phenotypes')

# input files
k12_genome = pj(data, 'genome.faa')
phenotypes = pj(phenotypes_dir, 'phenotypes.tsv')
input_file = pj(data, 'inputs.tsv')
strains = [x.split('.')[0] for x in os.listdir(genomes_dir)
           if x.endswith('.fasta')]
genomes = [pj(genomes_dir, x + '.fasta') for x in strains]

# output directories
out = config.get('out', 'out')
annotations_dir = pj(out, 'annotations')
parsnp_tree_dir = pj(out, 'parsnp')
gubbins_tree_dir = pj(out, 'gubbins')
roary_dir = pj(out, 'roary')
associations_dir = pj(out, 'associations')

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
# associations
# kmers
associations_kmer = pj(associations_dir, 'associations_kmer.tsv')
patterns_kmer = pj(associations_dir, 'patterns_kmer.txt')
lineage_kmer = pj(associations_dir, 'lineage_kmer.tsv')
associations_cont_kmer = pj(associations_dir, 'associations_cont_kmer.tsv')
patterns_cont_kmer = pj(associations_dir, 'patterns_cont_kmer.txt')
lineage_cont_kmer = pj(associations_dir, 'lineage_cont_kmer.tsv')
associations_lmm_kmer = pj(associations_dir, 'associations_lmm_kmer.tsv')
patterns_lmm_kmer = pj(associations_dir, 'patterns_lmm_kmer.txt')
lineage_lmm_kmer = pj(associations_dir, 'lineage_lmm_kmer.tsv')
h2_lmm_kmer = pj(associations_dir, 'h2_lmm_kmer.txt')
associations_cont_lmm_kmer = pj(associations_dir, 'associations_cont_lmm_kmer.tsv')
patterns_cont_lmm_kmer = pj(associations_dir, 'patterns_cont_lmm_kmer.txt')
lineage_cont_lmm_kmer = pj(associations_dir, 'lineage_cont_lmm_kmer.tsv')
h2_cont_lmm_kmer = pj(associations_dir, 'h2_cont_lmm_kmer.txt')
# pangenome
associations_rtab = pj(associations_dir, 'associations_rtab.tsv')
patterns_rtab = pj(associations_dir, 'patterns_rtab.txt')
lineage_rtab = pj(associations_dir, 'lineage_rtab.tsv')
associations_cont_rtab = pj(associations_dir, 'associations_cont_rtab.tsv')
patterns_cont_rtab = pj(associations_dir, 'patterns_cont_rtab.txt')
lineage_cont_rtab = pj(associations_dir, 'lineage_cont_rtab.tsv')
associations_lmm_rtab = pj(associations_dir, 'associations_lmm_rtab.tsv')
patterns_lmm_rtab = pj(associations_dir, 'patterns_lmm_rtab.txt')
lineage_lmm_rtab = pj(associations_dir, 'lineage_lmm_rtab.tsv')
h2_lmm_rtab = pj(associations_dir, 'h2_lmm_rtab.txt')
associations_cont_lmm_rtab = pj(associations_dir, 'associations_cont_lmm_rtab.tsv')
patterns_cont_lmm_rtab = pj(associations_dir, 'patterns_cont_lmm_rtab.txt')
lineage_cont_lmm_rtab = pj(associations_dir, 'lineage_cont_lmm_rtab.tsv')
h2_cont_lmm_rtab = pj(associations_dir, 'h2_cont_lmm_rtab.txt')

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
  input: genomes_dir
  output: parsnp_tree
  params:
    outdir=parsnp_tree_dir,
    focus=pj(genomes_dir, 'IAI01.fasta')
  threads: 20
  shell:
    'parsnp -d {input} -r {params.focus} -p {threads} -o {params.outdir} -v -c'

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
  input: genomes_dir
  output: sketches
  threads: 5
  params: sketches_base
  shell:
    'mash sketch -p {threads} -s 10000 -o {params} {input}/*.fasta'

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
    dist=mash_distances
  output:
    associations=associations_kmer,
    patterns=patterns_kmer,
    lineage=lineage_kmer
  threads: 40
  params:
    dimensions=3
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column phenotype --kmers {input.kmers} --max-dimensions {params.dimensions} --lineage --lineage-file {output.lineage} --cpu {threads} --output-patterns {output.patterns} --distance {input.dist} > {output.associations}'

rule:
  input:
    phenotype=phenotypes,
    kmers=kmers,
    dist=mash_distances
  output:
    associations=associations_cont_kmer,
    patterns=patterns_cont_kmer,
    lineage=lineage_cont_kmer
  threads: 40
  params:
    dimensions=3
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column killed --kmers {input.kmers} --max-dimensions {params.dimensions} --lineage --lineage-file {output.lineage} --cpu {threads} --output-patterns {output.patterns} --distance {input.dist} > {output.associations}'

rule:
  input:
    phenotype=phenotypes,
    kmers=kmers,
    dist=mash_distances,
    sim=gubbins_similarities
  output:
    associations=associations_lmm_kmer,
    patterns=patterns_lmm_kmer,
    lineage=lineage_lmm_kmer,
    h2=h2_lmm_kmer
  threads: 40
  params:
    dimensions=3
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column phenotype --kmers {input.kmers} --max-dimensions {params.dimensions} --lineage --lineage-file {output.lineage} --cpu {threads} --output-patterns {output.patterns} --distance {input.dist} --lmm --similarity {input.sim} 2>&1 > {output.associations} | grep \'h^2\' > {output.h2}'

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
    associations_kmer,
    associations_cont_kmer,
    associations_lmm_kmer,
    associations_cont_lmm_kmer

rule:
  input:
    phenotype=phenotypes,
    rtab=roary,
    dist=mash_distances
  output:
    associations=associations_rtab,
    patterns=patterns_rtab,
    lineage=lineage_rtab
  threads: 5
  params:
    dimensions=3
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column phenotype --pres {input.rtab} --max-dimensions {params.dimensions} --lineage --lineage-file {output.lineage} --cpu {threads} --output-patterns {output.patterns} --distance {input.dist} > {output.associations}'

rule:
  input:
    phenotype=phenotypes,
    rtab=roary,
    dist=mash_distances
  output:
    associations=associations_cont_rtab,
    patterns=patterns_cont_rtab,
    lineage=lineage_cont_rtab
  threads: 5
  params:
    dimensions=3
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column killed --pres {input.rtab} --max-dimensions {params.dimensions} --lineage --lineage-file {output.lineage} --cpu {threads} --output-patterns {output.patterns} --distance {input.dist} > {output.associations}'

rule:
  input:
    phenotype=phenotypes,
    rtab=roary,
    dist=mash_distances,
    sim=gubbins_similarities
  output:
    associations=associations_lmm_rtab,
    patterns=patterns_lmm_rtab,
    lineage=lineage_lmm_rtab,
    h2=h2_lmm_rtab
  threads: 5
  params:
    dimensions=3
  shell:
    'pyseer --phenotypes {input.phenotype} --phenotype-column phenotype --pres {input.rtab} --max-dimensions {params.dimensions} --lineage --lineage-file {output.lineage} --cpu {threads} --output-patterns {output.patterns} --distance {input.dist} --lmm --similarity {input.sim} 2>&1 > {output.associations} | grep \'h^2\' > {output.h2}'

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
    associations_rtab,
    associations_cont_rtab,
    associations_lmm_rtab,
    associations_cont_lmm_rtab
