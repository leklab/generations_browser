import time
import pymongo
import gzip
#import argparse


CHROMOSOMES = ['chr%s' % x for x in range(1, 23)]
CHROMOSOMES.extend(['chrX', 'chrY', 'chrM'])
CHROMOSOME_TO_CODE = { item: i+1 for i, item in enumerate(CHROMOSOMES) }


def get_single_location(chrom, pos):
	"""
	Gets a single location from chromosome and position
	chr must be actual chromosme code (chrY) and pos must be integer
	Borrowed from xbrowse
	"""
	return CHROMOSOME_TO_CODE[chrom] * int(1e9) + pos


def get_xpos(chrom, pos):
	"""
	Borrowed from xbrowse
	"""
	if not chrom.startswith('chr'):
		chrom = 'chr{}'.format(chrom)
	return get_single_location(chrom, int(pos))

def connect_db(host = 'localhost', port = 27017, db_name = 'exac'):
    """
    Connects to the specific database.
    """
    client = pymongo.MongoClient(host=host, port=port)
    return client[db_name]


def get_genes_from_gencode_gtf(gtf_file):
	"""
	Parse gencode GTF file;
	Returns iter of gene dicts
	"""
	for line in gtf_file:
		        
		if line.startswith('#'):
			continue
		
		fields = line.strip('\n').split('\t')

		if fields[2] != 'gene':
			continue

		chrom = fields[0][3:]
		start = int(fields[3]) + 1  # bed files are 0-indexed
		stop = int(fields[4]) + 1
		info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
		info = {k: v.strip('"') for k, v in info.items()}
		gene_id = info['gene_id'].split('.')[0]

		gene = {
			'gene_id': gene_id,
			'gene_name': info['gene_name'],
			'gene_name_upper': info['gene_name'].upper(),
			'chrom': chrom,
			'start': start,
			'stop': stop,
			'strand': fields[6],
			'xstart': get_xpos(chrom, start),
			'xstop': get_xpos(chrom, stop),
		}

		yield gene


def get_transcripts_from_gencode_gtf(gtf_file):

	"""
	Parse gencode GTF file;
	Returns iter of transcript dicts
	"""
	for line in gtf_file:
		if line.startswith('#'):
			continue
	
		fields = line.strip('\n').split('\t')

		if fields[2] != 'transcript':
			continue

		chrom = fields[0][3:]
		start = int(fields[3]) + 1  # bed files are 0-indexed
		stop = int(fields[4]) + 1
		info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
		info = {k: v.strip('"') for k, v in info.items()}
		gene_id = info['gene_id'].split('.')[0]
		transcript_id = info['transcript_id'].split('.')[0]

		gene = {
			'transcript_id': transcript_id,
			'gene_id': gene_id,
			'chrom': chrom,
			'start': start,
			'stop': stop,
			'strand': fields[6],
			'xstart': get_xpos(chrom, start),
			'xstop': get_xpos(chrom, stop),
		}

		yield gene


def get_exons_from_gencode_gtf(gtf_file):
	"""
	Parse gencode GTF file;
	Returns iter of transcript dicts
	"""
	for line in gtf_file:
		if line.startswith('#'):
			continue
		fields = line.strip('\n').split('\t')

		if fields[2] not in ['exon', 'CDS', 'UTR']:
			continue

		chrom = fields[0][3:]
		feature_type = fields[2]
		start = int(fields[3]) + 1  # bed files are 0-indexed
		stop = int(fields[4]) + 1
		info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
		info = {k: v.strip('"') for k, v in info.items()}
		gene_id = info['gene_id'].split('.')[0]
		transcript_id = info['transcript_id'].split('.')[0]

		exon = {
			'feature_type': feature_type,
			'transcript_id': transcript_id,
			'gene_id': gene_id,
			'chrom': chrom,
			'start': start,
			'stop': stop,
			'strand': fields[6],
			'xstart': get_xpos(chrom, start),
			'xstop': get_xpos(chrom, stop),
		}
		yield exon


def get_canonical_transcripts(canonical_transcript_file):
    for line in canonical_transcript_file:
        gene, transcript = line.strip().split()
        yield gene, transcript

def get_omim_associations(omim_file):
    for line in omim_file:
        fields = line.strip().split('\t')
        if len(fields) == 4:
            yield fields
        else:
            yield None


def get_dbnsfp_info(dbnsfp_file):
    """
    Parse dbNSFP_gene file;
    Returns iter of transcript dicts
    """
    header = dbnsfp_file.readline().split('\t')
    fields = dict(zip(header, range(len(header))))
    for line in dbnsfp_file:
        line = line.split('\t')
        other_names = line[fields["Gene_old_names"]].split(';') if line[fields["Gene_old_names"]] != '.' else []
        if line[fields["Gene_other_names"]] != '.':
            other_names.extend(line[fields["Gene_other_names"]].split(';'))
        gene_info = {
            'gene_name': line[fields["Gene_name"]],
            'ensembl_gene': line[fields["Ensembl_gene"]],
            'gene_full_name': line[fields["Gene_full_name"]],
            'gene_other_names': other_names
        }
        yield gene_info



def load_gene_models(gencode_gtf = 'gencode.gtf.gz', canonical_transcript = 'canonical_transcripts.txt.gz', dbnsfp= 'dbNSFP2.6_gene.gz', omim = 'omim_info.txt.gz'):

	db = connect_db()

	db.genes.drop()
	db.transcripts.drop()
	db.exons.drop()
	print('Dropped db.genes, db.transcripts, and db.exons.')

	start_time = time.time()

	canonical_transcripts = {}
	
	with gzip.open(canonical_transcript,'rt') as canonical_transcript_file:
		for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
			canonical_transcripts[gene] = transcript

	omim_annotations = {}
	with gzip.open(omim,'rt') as omim_file:
		for fields in get_omim_associations(omim_file):
			if fields is None:
				continue
			gene, transcript, accession, description = fields
			omim_annotations[gene] = (accession, description)

	dbnsfp_info = {}
	with gzip.open(dbnsfp, 'rt') as dbnsfp_file:
		for dbnsfp_gene in get_dbnsfp_info(dbnsfp_file):
			other_names = [other_name.upper() for other_name in dbnsfp_gene['gene_other_names']]
			dbnsfp_info[dbnsfp_gene['ensembl_gene']] = (dbnsfp_gene['gene_full_name'], other_names)

	print('Done loading metadata. Took %s seconds' % int(time.time() - start_time))

	# grab genes from GTF
	start_time = time.time()
	with gzip.open(gencode_gtf, 'rt') as gtf_file:
		for gene in get_genes_from_gencode_gtf(gtf_file):
			gene_id = gene['gene_id']
			if gene_id in canonical_transcripts:
				gene['canonical_transcript'] = canonical_transcripts[gene_id]
			if gene_id in omim_annotations:
				gene['omim_accession'] = omim_annotations[gene_id][0]
				gene['omim_description'] = omim_annotations[gene_id][1]
			if gene_id in dbnsfp_info:
				gene['full_gene_name'] = dbnsfp_info[gene_id][0]
				gene['other_names'] = dbnsfp_info[gene_id][1]
			db.genes.insert(gene, w=0)

	print('Done loading genes. Took %s seconds' % int(time.time() - start_time))

	start_time = time.time()
	db.genes.ensure_index('gene_id')
	db.genes.ensure_index('gene_name_upper')
	db.genes.ensure_index('gene_name')
	db.genes.ensure_index('other_names')
	db.genes.ensure_index('xstart')
	db.genes.ensure_index('xstop')
	print('Done indexing gene table. Took %s seconds' % int(time.time() - start_time))

	# and now transcripts
	start_time = time.time()
	with gzip.open(gencode_gtf, 'rt') as gtf_file:
		for transcript in get_transcripts_from_gencode_gtf(gtf_file):
			db.transcripts.insert(transcript)
	#db.transcripts.insert((transcript for transcript in get_transcripts_from_gencode_gtf(gtf_file)), w=0)
	print('Done loading transcripts. Took %s seconds' % int(time.time() - start_time))

	start_time = time.time()
	db.transcripts.ensure_index('transcript_id')
	db.transcripts.ensure_index('gene_id')
	print('Done indexing transcript table. Took %s seconds' % int(time.time() - start_time))

	# Building up gene definitions
	start_time = time.time()
	with gzip.open(gencode_gtf, 'rt') as gtf_file:
		#db.exons.insert((exon for exon in get_exons_from_gencode_gtf(gtf_file)), w=0)
		db.exons.insert((exon for exon in get_exons_from_gencode_gtf(gtf_file)))

	print('Done loading exons. Took %s seconds' % int(time.time() - start_time))

	start_time = time.time()
	db.exons.ensure_index('exon_id')
	db.exons.ensure_index('transcript_id')
	db.exons.ensure_index('gene_id')
	print('Done indexing exon table. Took %s seconds' % int(time.time() - start_time))

	return []



if __name__ == '__main__':
	load_gene_models()













