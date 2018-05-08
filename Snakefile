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
roary_dir = pj(out, 'roary')

# output files
annotations = [pj(annotations_dir, x, x + '.gff')
               for x in strains]
parsnp_tree = pj(parsnp_tree_dir, 'parsnp.tree')
polished_parsnp_tree = pj(parsnp_tree_dir, 'tree.nwk')
parsnp_xmfa = pj(parsnp_tree_dir, 'parsnp.xmfa')
parsnp_alignment = pj(parsnp_tree_dir, 'parsnp.fasta')
kmers = pj(out, 'kmers.gz')
sketches_base = pj(out, 'sketches')
sketches = sketches_base + '.msh'
mash_distances = pj(out, 'mash.tsv')
roary = pj(roary_dir, 'gene_presence_absence.Rtab')

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

rule make_gubbins_tree:
  input: parsnp_tree
  output: parsnp_alignment
  params: parsnp_xmfa
  shell:
    'harvesttools -x {params} -M {output}'

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

rule pangenome:
  input: annotations
  output: roary
  params: roary_dir
  threads: 20
  shell:
    'rm -rf {params} && roary -p {threads} -f {params} -s -v -g 100000 {input}'
