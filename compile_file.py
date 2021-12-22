#coding = utf-8
#!/usr/bin/python3
import os

class compile(object):
	def __init__(self):
		self.input_gate = 'output'
		self.gates = self.get_gate()
		for gate in self.gates:
			paths = self.get_source(gate)
			self.add_in(paths,gate)

	def get_gate(self):
		gates = []
		for i in os.listdir(self.input_gate):
			gates.append(i)
		return gates

	def get_source(self,gate):
		paths = []
		for root, dirs, files in os.walk(gate):
			for name in files:
				path = os.path.join(root,name)
				paths.append(path)
		return paths

	def add_in(self,paths,gate):
		for path in paths:
			with open(path,'r') as f:
				content = f.read()
			with open(gate + '.fasta','a',encoding='utf8') as g:
				g.write(content)

if __name__ == '__main__':
	compile()