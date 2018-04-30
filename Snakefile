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

# output files
annotations = [pj(annotations_dir, x, x + '.gff')
               for x in strains]
parsnp_tree = pj(parsnp_tree_dir, 'parsnp.tree')
kmers = pj(out, 'kmers.gz')

rule annotate:
  input: annotations

rule:
  input:
    ref=k12_genome,
    genome=pj(genomes_dir, '{strain}.fasta')
  output: pj(annotations_dir, '{strain}', '{strain}.gff')
  params: pj(annotations_dir, '{strain}')
  shell:
    'prokka --outdir {params} --force --prefix {wildcards.strain} --addgenes --locustag {wildcards.strain} --mincontiglen 200 --genus Escherichia -species coli --strain {wildcards.strain} --proteins {input.ref} --cpus 1 {input.genome}'

rule make_tree:
  input: genomes_dir
  output: parsnp_tree
  params:
    outdir=parsnp_tree_dir,
    focus=pj(genomes_dir, 'IAI01.fasta')
  threads: 20
  shell:
    'parsnp -d {input} -r {params.focus} -p {threads} -o {params.outdir} -v -c'

rule do_kmers:
  input: input_file
  output: kmers
  shell:
    'fsm-lite -l {input} -t tmp.txt -m 9 -M 100 -s 1 -v | gzip > {output}'
