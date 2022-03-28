import datetime
from math import ceil
import sys

class progress_tracker:
	def __init__(self, total, step_size = 2, message = None):
		self.current_count = 0
		self.max_count = total
		#Book keeping.
		self.start_time = None
		self.end_time = None
		#Show progrexx every [step] percent
		self.step = step_size
		self.justify_size = ceil(100/self.step)
		self.last_percent = 0
		self.message = message
		
		self.start()

	def curtime(self):
		time_format = "%d/%m/%Y %H:%M:%S"
		timer = datetime.datetime.now()
		time = timer.strftime(time_format)
		return time
		
	def start(self):
		if self.message is not None:
			print(self.message)
		
		print("")
		try:
			percentage = (self.current_count/self.max_count)*100
			sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/self.step)).ljust(self.justify_size)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(self.current_count) + " of " + str(self.max_count) + ' ) at ' + self.curtime() + "\n")
			sys.stdout.flush()
		except:
			#It's not really a big deal if the progress bar cannot be printed.
			pass
	
	def update(self):
		self.current_count += 1
		percentage = (self.current_count/self.max_count)*100
		try:
			if percentage // self.step > self.last_percent:
				sys.stdout.write('\033[A')
				sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*int(percentage/self.step)).ljust(self.justify_size)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(self.current_count) + " of " + str(self.max_count) + ' ) at ' + self.curtime() + "\n")
				sys.stdout.flush()
				self.last_percent = percentage // self.step
			#Bar is always full at the end.
			if count == self.max_count:
				sys.stdout.write('\033[A')
				sys.stdout.write("Completion".rjust(3)+ ' |'+('#'*self.justify_size).ljust(self.justify_size)+'| ' + ('%.2f'%percentage).rjust(7)+'% ( ' + str(self.current_count) + " of " + str(self.max_count) + ' ) at ' + self.curtime() + "\n")
				sys.stdout.flush()
				#Add space at end.
				print("")
		except:
			#It's not really a big deal if the progress bar cannot be printed.
			pass
					
					

		
			