#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

# return string of specified number of bytes in hex
def random_vector(length):
	return str(os.urandom(length).encode('hex'))	

# read input as an integer, provided it is of type
# - integer, long
# - hexstring starting with or without '0x' 
def parse_input(x):
	if isinstance(x, str):
		y = x.replace('0x', '')
		try:
			z = int(y,16)
			return True, z
		except Exception as e:
			return False, -1
	try:
		y = int(x)
		return True, y
	except:
		return False, -2

# return the minimal size of binary vectors needed to describe all integers
# in the basis
def vector_length(basis):
	length = 0
    	for v in basis:
    		length = max(length, v.bit_length() )
    	return length

# transorm integer tot a binary array
def num2vec(num, length):
	v = []
	t = num
	for i in xrange(length):
		v.append(t%2)
		t = (t>>1)
	v.reverse()
	return v

# transorm a binary array tot an integer
def vec2num(vec):
	n = 0
	l = len(vec)
	for i in xrange(l):
		n *= 2
		n += int( vec[i] )
	return n

# write input x as a linear combination
# of the basis vectors
def linear_combination(x, basis):
	res, val = parse_input(x)
    	if res:
    		x_int = val
    	else:
    		return False
    	M = []
	size_basis = len(basis)
	dimension = x_int.bit_length()
	for v in basis:
		dimension = max(v.bit_length(), dimension) 
    	for v in basis:
		vec = num2vec(v, dimension)
		M.append(vec)
	A = matrix(GF(2), size_basis, dimension, M).transpose()
	r = A.rank()
	vec = num2vec(x_int, dimension)
	M.append(vec)
	A = matrix(GF(2), size_basis+1, dimension, M).transpose()
	print r, A.rank()
	if A.rank() > r: #in this case the target was not a lin combo of the basis, so we fail
		return []
	else:
		B = A.echelon_form()
		combo_index = B.transpose()[size_basis][:size_basis]
		combo_list = []
		res = 0
		for i in xrange(size_basis):
			if combo_index[i] == 1:
				combo_list.append(basis[i])
		return combo_list
	

def random_affine_map(input_dimension, output_dimension):
	basis = []
	A = affine_gf2()
	while len(basis) < min(input_dimension, output_dimension):
		L = []
		for v in basis:
			vec = num2vec(v, output_dimension)
			L.append(vec)
		r = int(random_vector(8*output_dimension),16)%(2**output_dimension)
		vec = num2vec(r, output_dimension)
		L.append(vec)
		M = matrix(GF(2), len(L), output_dimension, L)
		if M.rank() > len(basis):
			A.L.map[1<<len(basis)] = r
			basis.append(r)
	A.constant = int(random_vector(8*output_dimension),16)%(2**output_dimension)
	A.Linv = A.L.construct_inverse()
	return A

class linear_gf2():
	def __init__(self):
		self.observed = {}
		self.map = {}
	
	def get_size_observed_inputs(self):
		return vector_length(self.observed.keys())

	def get_size_observed_outputs(self):
		return vector_length(self.observed.values())

	def normalise(self):
		normal = {}
		size = vector_length(self.map.keys())
		if len(self.map) < size:
			print("error: not enough vectors in basis")
			return
		for i in xrange(size):
			x = 1<<i
			normal[x] = self.evaluate(x)
		self.map.clear()
		self.map = dict(normal)

	def evaluate(self, x):
		if len(self.map) == 0:
			print("error: map not initialized")
			return None
		f,x = parse_input(x)
		if not f:
			print("error: input badly formatted")
			return None
		if x == 0:
			return 0
		combo = linear_combination(x, self.map.keys())
		if combo == []:
			print("error: input not in support of map")
			return None
		res = 0
		for v in combo:
			res ^= self.map[v]
		return res

	def construct_inverse(self):
		Linv = linear_gf2()
		Linv.map = {y:x for x,y in self.map.items()}
		return Linv
	
	def fixed_points(self):
		self.normalise()
		M = []
		size_basis = len(self.map)
		size = vector_length(self.map.values())
		if size != size_basis:
			print("error: input and outputvectors are not of same length")
			return None
		for i in xrange(size):
			v = num2vec( self.map[1<<(size-i-1)], size)
			M.append(v)	
		A = matrix(GF(2), size_basis, size, M).transpose()
		AA = A + matrix.identity(GF(2), size)
		K = []
		for v in AA.right_kernel().basis():
			K.append(vec2num(v))
		E = []
		n = len(K)
		for i in xrange(2**n): 
			t = i
			res = 0
			for j in xrange(n):
				if t&1 == 1:
					res ^= K[j]
				t = t>>1
			E.append(res)
		return E				

class affine_gf2():
	def __init__(self):
		self.observed = {}
		self.L = linear_gf2()
		self.Linv = linear_gf2()
		self.constant = None

	def get_size_observed_inputs(self):
		return vector_length(self.observed.keys())

	def get_size_observed_outputs(self):
		return vector_length(self.observed.values())

	def construct_inverse_lin_map(self):
		self.Linv = self.L.construct_inverse()

	def evaluate(self, x):
		res = self.L.evaluate(x)
		if res is not None:
			res ^= self.constant
		return res
		
	def preimage(self, x):
		f,res = parse_input(x)
		if f:
			return self.Linv.evaluate(res ^ self.constant)
		return None

   	def status(self):
		print("number of vectors in basis = %d" % len(self.L.map))
		print("input vectors consist of %d bits" % self.get_size_observed_inputs() )
		print("output vectors consist of %d bits" % self.get_size_observed_outputs() )
	
	# add 
	def observe(self, x, y):
		res, val = parse_input(x)
		if res:
			x_int = val
		else:
			return False
		res, val = parse_input(y)
		if res:
			y_int = val
		else:
			return False
		self.observed[x_int] = y_int
		return True
	
	# find the fixes points, i.e. x such that A(x) = x
	def fixed_points(self):
		# use the orthonormal representation
		self.transform()
		self.L.normalise()
		M = []
		size_basis = len(self.L.map)
		size = vector_length(self.L.map.values())
		# check for input and output domains to match
		if size != size_basis:
			print("error: input and outputvectors are not of same length")
			return None
		# start constructing the augmented matrix
		for i in xrange(size):
			v = num2vec(self.L.map[1<<(size-i-1)], size)
			v.append(0)
			M.append(v)
		# construct the last column of the augmented matrix
		v = num2vec(self.constant, size)
		v.append(1)
		M.append(v)
		# convert the arrau tot a matrix over GF(2)
		A = matrix(GF(2), size_basis+1, size+1, M).transpose()
		# construct AA = A + I, fixed points of A can be derived from the kernel of AA
		AA = A + matrix.identity(GF(2), size+1)
		K = []
		# get a basis of the kernel of AA
		for v in AA.right_kernel().basis():
			K.append(vec2num(v))
		E = []
		n = len(K)
		# get all the vectors in the kernel of AA and check them, 
		# since not all of them correspond to fixed points of A 
		for i in xrange(2**n): 
			t = i
			res = 0
			for j in xrange(n):
				if t&1 == 1:
					res ^= K[j]
				t = t>>1
			# we are only interested in vectors in the kernel with its last coordinate
			# equal to '1'
			if res&1 == 1:
				E.append( (res-1)/2 )
		return E				

	def add_to_basis(self, x,y, size):
		M = []
		size_basis = len(self.L.map)
		for vec in self.L.map.keys():
			v = num2vec(vec, size)
			M.append(v)	
		v = num2vec(x, size)
		M.append(v)
		A = matrix(GF(2), size_basis+1, size, M).transpose()
		if A.rank() > size_basis: 
			#in this case the target was not a lin combo of the basis
			self.L.map[x] = y
			return True
		else:
			# we will use the extra relation to determine the constant or check for consistency
			combo = linear_combination(x, self.L.map.keys())
			res = y
			for v in combo:
				res ^= self.L.map[v]
			if res != 0:
				if self.constant == None:
					print("[+] recovered constant")
					self.constant = res
					return True	
				if self.constant != res:
					print("[-] inconsistent equation found")
					return False						
		
	def transform(self):
		num_observed = len(self.observed)
		observed_input_vector_size = vector_length(self.observed.keys())
		observed_output_vector_size = vector_length(self.observed.values())
		for x in self.observed.keys():
			self.add_to_basis(x, self.observed[x], observed_input_vector_size)
		for x in self.L.map.keys():
			y = self.L.map[x]
			self.L.map[x] = y^self.constant
		self.construct_inverse_lin_map()
		return





