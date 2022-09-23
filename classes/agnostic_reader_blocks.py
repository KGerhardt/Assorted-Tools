import gzip
from io import BufferedReader
has_xopen = False

try:
	from xopen import xopen
	has_xopen = True
except:
	pass

#Iterator for agnostic reader
class agnostic_reader_iterator:
	def __init__(self, reader):
		self.handle_ = reader.handle
		self.is_gz_ = reader.is_gz
		
		self.remnant = ""
		self.nextlines = []
		self.nextindex = 0
		
	def next_block(self, size=65536):
		b = self.handle_.read(size)
		if b:
			b = self.remnant + b
			self.remnant = ""
		
			if not b.endswith("\n"):
				backtrack = b.rindex("\n")+1
				self.remnant = b[backtrack:]
				b = b[:backtrack]
			
			b = b.splitlines(True)
			
		return b
				
		
	def __next__(self):
		if self.is_gz_:
			if has_xopen:
				line = self.handle_.readline()
			else:
				line = self.handle_.readline().decode()
				
		else:
			line = self.next_block()
		
		#return or stop iteration	
		if line:
			return line
		else:
			raise StopIteration
			
class agnostic_reader:
	def __init__(self, file):
		self.path = file
		
		with open(file, 'rb') as test_gz:
			#Gzip magic number
			is_gz = (test_gz.read(2) == b'\x1f\x8b')
		
		self.is_gz = is_gz
		
		if self.is_gz:
			if has_xopen:
				self.handle = xopen(self.path, threads = 0)
			else:
				ref = gzip.open(self.path)
				self.handle = BufferedReader(ref)
		else:
			self.handle = open(self.path)
			
	def __iter__(self):
		return agnostic_reader_iterator(self)
		
	def close(self):
		if self.is_gz:
			if has_xopen:
				self.handle.close()
			else:
				self.handle.close()
				ref.close()
		else:
			self.handle.close()
	
def agnostic_open(file):
	ar = agnostic_reader(file)
	if ar.is_gz:
		for line in ar:
			yield line
	else:
		for block in ar:
			for line in block:
				yield line
	
	ar.close()

#Example usage
'''	
import sys
file = sys.argv[1]

for line in agnostic_open(file):
	print(line)
'''
