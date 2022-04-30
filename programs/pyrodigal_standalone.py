'''
This script manages the behavior of pyrodigal for protein prediction.
'''
import pyrodigal as pd
import gzip
import argparse
	
#Iterator for agnostic reader
class agnostic_reader_iterator:
	def __init__(self, reader):
		self.handle_ = reader.handle
		self.is_gz_ = reader.is_gz
		
	def __next__(self):
		if self.is_gz_:
			line = self.handle_.readline().decode()
		else:
			line = self.handle_.readline()
		
		#Ezpz EOF check
		if line:
			return line
		else:
			raise StopIteration

#File reader that doesn't care if you give it a gzipped file or not.
class agnostic_reader:
	def __init__(self, file):
		self.path = file
		
		with open(file, 'rb') as test_gz:
			#Gzip magic number
			is_gz = (test_gz.read(2) == b'\x1f\x8b')
		
		self.is_gz = is_gz
		
		if is_gz:
			self.handle = gzip.open(self.path)
		else:
			self.handle = open(self.path)
			
	def __iter__(self):
		return agnostic_reader_iterator(self)
		
	def close(self):
		self.handle.close()
	
class pyrodigal_manager:
	def __init__(self, file = None, aa_out = None, nt_out = None, is_meta = False, full_headers = True, trans_table = 11,
				num_bp_fmt = True, verbose = True, do_compress = False, compare_against = None):
		#Input NT sequences
		self.file = file
		
		#List of seqs read from input file.
		self.sequences = None
		#Concatenation of up to first 32 million bp in self.sequences - prodigal caps at this point.
		self.training_seq = None
		
		#Predicted genes go here
		self.predicted_genes = None
		#Record the translation table used.
		self.trans_table = trans_table
		
		#This is the pyrodigal manager - this does the gene predicting.
		self.manager = pd.OrfFinder(meta=is_meta)
		self.is_meta = is_meta
		
		#Full prodigal header information includes more than just a protein number.
		#If full_headers is true, protein deflines will match prodigal; else, just protein ID.
		self.full_headers = full_headers
		
		#Prodigal prints info to console. I enhanced the info and made printing default, but also allow them to be totally turned off.
		self.verbose = verbose
		
		#Prodigal formats outputs with 70 bases per line max
		self.num_bp_fmt = num_bp_fmt
		
		#File names for outputs
		self.aa_out = aa_out
		self.nt_out = nt_out
		
		#List of proteins in excess of 100K base pairs (HMMER's limit) and their lengths. This is also fastAAI specific.
		self.excluded_seqs = {}
		
		#Gzip outputs if asked.
		self.compress = do_compress
		
		#Normally, we don't need to keep an input sequence after it's had proteins predicted for it - however
		#For FastAAI and MiGA's purposes, comparisons of two translation tables is necessary.
		#Rather than re-importing sequences and reconstructing the training sequences, 
		#keep them for faster repredict with less I/O
		self.compare_to = compare_against
		if self.compare_to is not None:
			self.keep_seqs = True
			self.keep_after_train = True
		else:
			self.keep_seqs = False
			self.keep_after_train = False
	
	#Imports a fasta as binary.
	def import_sequences(self):
		if self.sequences is None:
			self.sequences = {}
			
		#check for zipped and import as needed.
		with open(self.file, 'rb') as test_gz:
			#Gzip magic number
			is_gz = (test_gz.read(2) == b'\x1f\x8b')
		
		if is_gz:
			fh = gzip.open(self.file)
		else:
			fh = open(self.file, "rb")
		
		imp = fh.readlines()
		
		#imp = []
		#for line in fh:
		#	imp.append(line.decode().strip())
		
		fh.close()
		
		cur_seq = None
		for s in imp:
			s = s.decode().strip()
			#> is 62 in ascii. This is asking if the first character is '>'
			if s.startswith(">"):
				#Skip first cycle, then do for each after
				if cur_seq is not None:
					self.sequences[cur_seq] = ''.join(self.sequences[cur_seq])
					self.sequences[cur_seq] = self.sequences[cur_seq].encode()
					#print(cur_seq, len(self.sequences[cur_seq]))
				cur_seq = s[1:]
				cur_seq = cur_seq.split()[0]
				cur_seq = cur_seq.encode('utf-8')
				self.sequences[cur_seq] = []
			else:
				#Remove the newline character.
				#bases = s[:-1]
				self.sequences[cur_seq].append(s)
		
		#Final set
		self.sequences[cur_seq] = ''.join(self.sequences[cur_seq])
		self.sequences[cur_seq] = self.sequences[cur_seq].encode()
		
		#Now we have the data, go to training.
		if not self.is_meta:
			self.train_manager()
		
	#Collect up to the first 32 million bases for use in training seq.
	def train_manager(self):
		running_sum = 0
		seqs_added = 0
		if self.training_seq is None:
			self.training_seq = []
			for seq in self.sequences:
				running_sum += len(self.sequences[seq])
				if seqs_added > 0:
					#Prodigal interleaving logic - add this breaker between sequences, starting at sequence 2
					self.training_seq.append(b'TTAATTAATTAA')
					running_sum += 12
					
				seqs_added += 1
					
				#Handle excessive size
				if running_sum >= 32000000:					
					print("Warning:  Sequence is long (max 32000000 for training).")
					print("Training on the first 32000000 bases.")
				
					to_remove = running_sum - 32000000
					
					#Remove excess characters
					cut_seq = self.sequences[seq][:-to_remove]
					#Add the partial seq
					self.training_seq.append(cut_seq)
					
					#Stop the loop and move to training
					break
				
				#add in a full sequence
				self.training_seq.append(self.sequences[seq])

			if seqs_added > 1:
				self.training_seq.append(b'TTAATTAATTAA')
				
			self.training_seq = b''.join(self.training_seq)
		
		if len(self.training_seq) < 20000:
			if self.verbose:
				print("Can't train on 20 thousand or fewer characters. Switching to meta mode.")
			self.manager = pd.OrfFinder(meta=True)
			self.is_meta = True
		else:
			if self.verbose:
				print("")
				#G is 71, C is 67; we're counting G + C and dividing by the total.
				gc = round(((self.training_seq.count(67) + self.training_seq.count(71))/ len(self.training_seq)) * 100, 2)
				print(len(self.training_seq), "bp seq created,", gc, "pct GC")
				
			#Train
			self.manager.train(self.training_seq, translation_table = self.trans_table)
		
		if not self.keep_after_train:
			#Clean up
			self.training_seq = None
		
	def predict_genes(self):
		if self.is_meta:
			print("Finding genes in metagenomic mode")
		else:
			print("Finding genes with translation table", self.trans_table)
			print("")
			
		self.predicted_genes = {}
		for seq in self.sequences:
			
			if self.verbose:
				print("Finding genes in sequence", seq.decode(), "("+str(len(self.sequences[seq]))+ " bp)... ", end = '')
				
			self.predicted_genes[seq] = self.manager.find_genes(self.sequences[seq])
				
			#If we're comparing multiple tables, then we want to keep these for re-prediction.
			if not self.keep_seqs:
				#Clean up
				self.sequences[seq] = None
			
			if self.verbose:
				print("done!")
		
		#Run alt comparisons in gene predict.
		if self.compare_to is not None:
			while len(self.compare_to) > 0:
				try:
					next_table = int(self.compare_to.pop(0))
					self.compare_alternative_table(next_table)
				except:
					print("Alternative table comparison failed! Skipping.")
	
	#Predict genes with an alternative table, compare results, and keep the winner.	
	def compare_alternative_table(self, table):
			if table == self.trans_table:
				print("You're trying to compare table", table, "with itself.")
			else:
				if self.verbose:
					print("Comparing translation table", self.trans_table, "against table", table)

				old_table = self.trans_table
				old_genes = self.predicted_genes
				old_size = 0
				for seq in self.predicted_genes:
					for gene in self.predicted_genes[seq]:
						old_size += (gene.end - gene.begin)
				
				self.trans_table = table
				self.train_manager()
				self.predict_genes()
				
				new_size = 0
				for seq in self.predicted_genes:
					for gene in self.predicted_genes[seq]:
						new_size += (gene.end - gene.begin)

				if (old_size / new_size) > 1.1:
					if self.verbose:
						print("Translation table", self.trans_table, "performed better than table", old_table, "and will be used instead.")
				else:
					if self.verbose:
						print("Translation table", self.trans_table, "did not perform significantly better than table", old_table, "and will not be used.")
					self.trans_table = old_table
					self.predicted_genes = old_genes
				
				#cleanup
				old_table = None
				old_genes = None
				old_size = None
				new_size = None
			
	#Break lines into size base pairs per line. Prodigal's default for bp is 70, aa is 60.
	def num_bp_line_format(self, string, size = 70):
		#ceiling funciton without the math module
		ceiling = int(round((len(string)/size)+0.5, 0))
		formatted = '\n'.join([string[(i*size):(i+1)*size] for i in range(0, ceiling)])
		return formatted
			
	#Writeouts
	def write_nt(self):
		if self.nt_out is not None:
			if self.verbose:
				print("Writing nucleotide sequences... ")
			if self.compress == '1' or self.compress == '2':
				out_writer = gzip.open(self.nt_out+".gz", "wb")
				
				content = b''
				
				for seq in self.predicted_genes:
					seqname = b">"+ seq + b"_"
					#Gene counter
					count = 1
					for gene in self.predicted_genes[seq]:
						#Full header lines
						if self.full_headers:
							content += b' # '.join([seqname + str(count).encode(), str(gene.begin).encode(), str(gene.end).encode(), str(gene.strand).encode(), gene._gene_data.encode()])
						else:
							#Reduced headers if we don't care.
							content += seqname + str(count).encode()
							
						content += b'\n'
							
						if self.num_bp_fmt:
							#60 bp cap per line
							content += self.num_bp_line_format(gene.sequence(), size = 70).encode()
						else:
							#One-line sequence.
							content += gene.sequence().encode()
							
						content += b'\n'
						count += 1
				
				out_writer.write(content)
				out_writer.close()
			
			if self.compress == '0' or self.compress == '2':
				out_writer = open(self.nt_out, "w")
			
				for seq in self.predicted_genes:
					#Only do this decode once.
					seqname = ">"+ seq.decode() +"_"
					#Gene counter
					count = 1
					
					for gene in self.predicted_genes[seq]:
						#Full header lines
						if self.full_headers:
							#Standard prodigal header
							print(seqname + str(count), gene.begin, gene.end, gene.strand, gene._gene_data, sep = " # ", file = out_writer)
						else:
							#Reduced headers if we don't care.
							print(seqname + str(count), file = out_writer)
							
						if self.num_bp_fmt:
							#60 bp cap per line
							print(self.num_bp_line_format(gene.sequence(), size = 70), file = out_writer)
						else:
							#One-line sequence.
							print(gene.sequence(), file = out_writer)
							
						count += 1
							
				out_writer.close()
		

		
	def write_aa(self):
		if self.aa_out is not None:
			if self.verbose:
				print("Writing amino acid sequences...")
			if self.compress == '1' or self.compress == '2':
				out_writer = gzip.open(self.aa_out+".gz", "wb")
				
				content = b''
				
				for seq in self.predicted_genes:
					seqname = b">"+ seq + b"_"
					#Gene counter
					count = 1
					for gene in self.predicted_genes[seq]:
						#Full header lines
						if self.full_headers:
							content += b' # '.join([seqname + str(count).encode(), str(gene.begin).encode(), str(gene.end).encode(), str(gene.strand).encode(), gene._gene_data.encode()])
						else:
							#Reduced headers if we don't care.
							content += seqname + str(count).encode()
							
						content += b'\n'	
						
						if self.num_bp_fmt:
							#60 bp cap per line
							content += self.num_bp_line_format(gene.translate(), size = 60).encode()
						else:
							#One-line sequence.
							content += gene.sequence().translate().encode()
							
						content += b'\n'
						count += 1
				
				out_writer.write(content)
				out_writer.close()
			
			if self.compress == '0' or self.compress == '2':
				out_writer = open(self.aa_out, "w")
				
				for seq in self.predicted_genes:
					#Only do this decode once.
					seqname = ">"+ seq.decode() +"_"
					#Gene counter
					count = 1
					
					for gene in self.predicted_genes[seq]:
						#Full header lines
						if self.full_headers:
						#Standard prodigal header
							print(seqname + str(count), gene.begin, gene.end, gene.strand, gene._gene_data, sep = " # ", file = out_writer)
						else:
							#Reduced headers if we don't care.
							print(seqname + str(count), file = out_writer)
							
						if self.num_bp_fmt:
							#60 bp cap per line
							print(self.num_bp_line_format(gene.translate(), size = 60), file = out_writer)
						else:
							#One-line sequence of amino acid translations.
							print(gene.translate(), file = out_writer)
						
						count += 1
							
				out_writer.close()

	
def options():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description = '''
	Replacer for Prodigl gene prediction using Pyrodigal.
	
	Currently only allows output of nucleotide and amino acid FASTA sequences, but has a few advantages:
		30-50% faster than prodigal.
		Takes gzipped inputs and (optionally) compresses outputs.
		
	If compressing outputs, ".gz" will be added to the name(s) given in --aa_out and --nt_out.''')
	parser.add_argument('-i', '--input',  dest = 'input', default = None, help = 'Input nucleotide file in FASTA format. Required.')
	
	parser.add_argument('-a', '--aa_out',  dest = 'aa_out', default = None, help = 'Output proteins AA translations to this file.')
	parser.add_argument('-d', '--nt_out',  dest = 'nt_out', default = None, help = 'Output proteins NT sequences to this file.')
	
	parser.add_argument('-g', '--trans_table', dest = 'table', default = 11, help = 'Select translation table to use (default 11)')
	
	parser.add_argument('--meta_mode', dest = 'meta', action = 'store_true', help = 'Run in metagenomic mode rather than training on your inputs. Off by default.')
	
	parser.add_argument('--compress',  dest = 'do_compress', default = 0, help = 'Gzip protein outputs. 0 = uncompressed, 1 = compressed, 2 = write both compressed and uncompressed.')
	parser.add_argument('--fixed_width',  dest = 'format_out', action='store_false', help = 'Print AA or NT sequences 1/line, instead of a fixed number of 70 BP or 60 AA per line.')

	parser.add_argument('--quiet',  dest = 'verbose', action='store_false', help = 'Print updates to console.')
	parser.add_argument('--short_header', dest = 'trunc', action='store_false', help = 'Print minimal sequence identifiers instead of normal Prodigal info. Only use if you only care about the protein sequences themselves. Off by default.')

	parser.add_argument('--compare_alts', dest = 'compare_against', default = None, help = "Comma-separated list of alternative translations tables to compare the table given by --trans_table against. Only the winner will be written.")
	
	opts = parser.parse_known_args()[0]
	return opts, parser
	
def main():
	opts, parser = options()
	do_continue = True
	input, nt_out, aa_out = opts.input, opts.nt_out, opts.aa_out
	if input is None:
		print("I need an input file.")
		parser.print_help()
		do_continue = False
		
	if nt_out is None and aa_out is None and do_continue:
		print("I need at least one output type.")
		parser.print_help()
		do_continue = False
		
	if do_continue:
		trans_table = int(opts.table)
		comp = opts.do_compress
		if comp not in ['0', '1', '2']:
			print("Compression option:", comp, "not recognized. Must be '0', '1', or '2' Defaulting to uncompressed output.")
			comp = '0'
		fmt_lines = opts.format_out
		verbose = opts.verbose
		short_head = opts.trunc
		meta = opts.meta
		alts = opts.compare_against
		if alts is not None:
			alts = alts.split(",")

		manager = pyrodigal_manager(file = input, aa_out = aa_out, nt_out = nt_out, is_meta = meta, full_headers = short_head, num_bp_fmt = fmt_lines, verbose = verbose, do_compress = comp, trans_table = trans_table, compare_against = alts)
		manager.import_sequences()
		manager.predict_genes()
		manager.write_nt()
		manager.write_aa()
	
if __name__=="__main__":
	main()
	
