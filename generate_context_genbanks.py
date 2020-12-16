import os, sys
import pandas as pd
from Bio.Alphabet import generic_dna, generic_protein
from Bio import SeqIO, SeqFeature
from Bio.SeqFeature import FeatureLocation
import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(
    description='Given a directory of protein fasta files containing the ORFs of interest and the nucleotide sequences from which they originate, generates optionally annotated genbank files for use with Clinker.')
requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument('-nucdir', metavar='[NUCLEOTIDE FASTA DIRECTORY]',
                    help='Directory containing nucleotide sequences.', required=True)
requiredNamed.add_argument('-protdir', metavar='[PRODIGAL PROTEIN FASTA DIRECTORY]',
                    help='Directory of prodigal output for your fastas in amino acid format.', default=None, required=True)
parser.add_argument('-pfamscandir', metavar='[Pfamscan annotation directory]',
                    help='Directory of pfamscan annotations (filename must follow the convention [genome id].faa.out)', default=None)
parser.add_argument('-kofamscandir', metavar='[Kofamscan annotation directory]',
                    help='Directory of kofamscan annotations (filename must follow the convention [genome id].faa.out.parsed.good)')
parser.add_argument('-goososfile', metavar='[GOOSOS output file]', help='path to all_hits_evalues_df.tsv from GOOSOS.py',
                    default=None)
parser.add_argument('-custom_annotationfile', metavar='Custom annotation two-column tsv file',
                    help='Two column no header file containing labeled ORF IDs and their corresponding custom annotation.', default=None)

requiredNamed.add_argument('-outdir', metavar='[OUTPUT DIRECTORY]', help='Output directory')

args = parser.parse_args()

if not os.path.exists(args.outdir):
	os.system('mkdir ' + args.outdir)
outdir = args.outdir

protein_dir = os.path.abspath(args.protdir)
contigs_dir = os.path.abspath(args.nucdir)

goososfile = args.goososfile

#Parse custom annotations file
if args.custom_annotationfile != None:
    annotationfile = os.path.abspath(args.custom_annotationfile)
    annotation_df = pd.read_csv(annotationfile, sep='\t', names=['orf_id', 'annotation'])
    annotation_dict = dict(zip(annotation_df.orf_id.tolist(), annotation_df.annotation.tolist()))
    test_proteinfile = os.listdir(protein_dir)[0]
    for rec in SeqIO.parse(os.path.join(protein_dir, test_proteinfile), 'fasta'):
        genome_id = test_proteinfile.split('.faa')[0]
        if '|' in rec.id and '|' not in annotation_df.iloc[0].orf_id:
            unlabeled_annotation = True
        else:
            unlabeled_annotation = False
        break
else:
    annotation_dict = None

#Just need pfamscan annotations for hydrogenases here
pfam_dir = args.pfamscandir

kofam_dir = args.kofamscandir

kofam_header = ['orf_id', 'knum', 'threshold', 'bitscore', 'evalue', 'description']

genbank_dir = os.path.abspath(args.outdir)

def parse_pfam_outfile(test_outfile):
    header = ['orf_id', 'aln_start', 'aln_end', 'envelope_start', 'envelope_end', 'pfam_acc', 'pfam_name',
         'hmm_type', 'hmm_start', 'hmm_end', 'hmm_length', 'bitscore', 'evalue', 'significance', 'clan']
    with open(test_outfile, 'r') as infile:
        lines = [x.rstrip() for x in infile.readlines()]

    split_lines = [x.split() for x in list(filter(lambda x: not x.startswith('#'), lines))]
    test_pfam_df = pd.DataFrame(split_lines[1:], columns=header)
    test_pfam_df['scaffold_id'] = test_pfam_df['orf_id'].apply(lambda x: '_'.join(x.split('_')[:-1]).split('|')[1])
    test_pfam_df['orfnum'] = test_pfam_df['orf_id'].apply(lambda x: int(x.split('_')[-1]))
    return test_pfam_df

def grab_scaffold(header):
    if '|' in header:
        return '_'.join(header.split('|')[1].split('_')[:-1])
    else:
        return '_'.join(header.split('_')[:-1])

def construct_scaffold_genbank(protein_recs, protein_file, scaffold_id, outdir=outdir):

    genome_id = protein_file.split('.faa')[0]

    #Get the scaffold rec
    scaffold_filter = lambda x: x.id == scaffold_id
    contigs = list(filter(scaffold_filter, SeqIO.parse(os.path.join(contigs_dir, protein_file.replace('.faa', '.fna')), 'fasta')))

    genbank_rec = contigs[0]

    #Scaffold-only pfamscan annotations
    if pfam_dir != None:
        pfam_file = list(filter(lambda x: x == protein_file + '.out', os.listdir(pfam_dir)))[0]
        pfam_df = parse_pfam_outfile(os.path.join(pfam_dir, pfam_file))
    else:
        pfam_df = None

    if kofam_dir != None:
        try:
            kofam_file = list(filter(lambda x: x == protein_file + '.out.parsed.good', os.listdir(kofam_dir)))[0]
        except:
            kofam_scafdf = None
            

        if kofam_scafdf != None:
            kofam_df = pd.read_csv(os.path.join(kofam_dir, kofam_file), sep='\t', names=kofam_header)

            kofam_df['scaffold_id'] = kofam_df.orf_id.apply(lambda x: '_'.join(x.split('_')[:-1]))
            kofam_df['scaffold_id'] = kofam_df.scaffold_id.apply(lambda x: x.split('|')[1] if '|' in x else x)
            kofam_scafdf = kofam_df[kofam_df.scaffold_id == scaffold_id]
    else:
        kofam_scafdf = None


    #Start the nucleotide sequence of the genbank file at the start position of the first CDS
    total_start = int(protein_recs[0].description.split(' # ')[1])
    total_end = int(protein_recs[-1].description.split(' # ')[2])
    genbank_rec.seq = genbank_rec.seq[total_start:total_end+1]
    genbank_rec.seq.alphabet = generic_dna
    count = 0
    if goososfile != None:
        goosos_df = pd.read_csv(goososfile, sep='\t')
    else:
        goosos_df = None

    for protein_rec in protein_recs:
        #Prepare location info and construct SeqFeature object
        start = int(protein_rec.description.split(' # ')[1]) - total_start
        startpos = SeqFeature.ExactPosition(start)
        end = int(protein_rec.description.split(' # ')[2]) - total_start
        endpos = int(SeqFeature.ExactPosition(end))
        strand = int(protein_rec.description.split(' # ')[3])
        rec_location = FeatureLocation(startpos, endpos)
        rec_feature = SeqFeature.SeqFeature(rec_location, type="CDS", strand=strand)

        #Add ORF name without genome ID
        if '|' in protein_rec.id:
            rec_feature.qualifiers['protein_id'] = protein_rec.id.split('|')[1]
        else:
            rec_feature.qualifiers['protein_id'] = protein_rec.id
        rec_feature.qualifiers['translation'] = protein_rec.seq

        #Get pfam info
        if pfam_dir != None:
            red_pfam_df = pfam_df[pfam_df.orf_id == protein_rec.id]
            domains = '+'.join(red_pfam_df.pfam_name.tolist())
            rec_feature.qualifiers['name'] = domains

        #Get kofam info
        if kofam_dir != None:
            red_kofam_df = kofam_scafdf[kofam_scafdf.orf_id == protein_rec.id]

            kofam_annotations = '+'.join(red_kofam_df.description.tolist())

            if kofam_annotations != '':
                rec_feature.qualifiers['locus_tag'] = kofam_annotations

        if goosos_df != None:
            red_goosos_df = goosos_df[goosos_df.orf_id == protein_rec.id]
            goosos_annotations = '+'.join(goosos_df.family_hmm.tolist())
            rec_feature.qualifiers['label'] = goosos_annotations

        if annotation_dict != None:

            if protein_rec.id in annotation_dict:
                rec_feature.qualifiers['id'] = annotation_dict[protein_rec.id]
            elif '|' in protein_rec.id and unlabeled_annotation:
                if protein_rec.id.split('|') in annotation_dict:
                    rec_feature.qualifiers['id'] = annotation_dict[protein_rec.id.split('|')]
            elif '|' not in protein_rec.id and not unlabeled_annotation:
                if genome_id + '|' + protein_rec.id in annotation_dict:
                    rec_feature.qualifiers['id'] = annotation_dict[genome_id + '|' + protein_rec]
        genbank_rec.features.append(rec_feature)


    SeqIO.write(genbank_rec, os.path.join(outdir, genome_id + '_' + scaffold_id + '.gbk'), 'genbank')
    return




def main():
    for protein_file in tqdm(os.listdir(protein_dir)):
        #Get the protein recs
        protein_recs = list(SeqIO.parse(os.path.join(protein_dir, protein_file), 'fasta'))
        rec_ids = [rec.id for rec in protein_recs]
        scaffolds = list(map(grab_scaffold, rec_ids))
        if len(set(scaffolds)) > 1:
            #There are multiple scaffolds in your dataset
            scaffolds = list(set(scaffolds))
            for scaffold in scaffolds:
                scaffold_proteins = list(filter(lambda x: grab_scaffold(x.id) == scaffold))
                construct_scaffold_genbank(scaffold_proteins, protein_file, scaffold)
        else:
            construct_scaffold_genbank(protein_recs, protein_file, scaffolds[0])

if __name__ == "__main__":
    main()
