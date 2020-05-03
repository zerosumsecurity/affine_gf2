#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

def random_vector(length):
	return str(os.urandom(length).encode('hex'))	

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

def vector_length(basis):
	length = 0
    	for v in basis:
    		length = max(length, v.bit_length() )
    	return length

def num2vec(num, length):
	v = []
	t = num
	for i in xrange(length):
		v.append(t%2)
		t = (t>>1)
	v.reverse()
	return v

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


D = [['1c72dfaf68f4b9bb', '003e4014add801a2'], ['963ead99d67b1370', '9ac5a128895204de'], ['ea39fccec802b2ff', '05951f3a4bff981d'], ['ae30c81c7b277ca9', '32ec7fc6cb575680'], ['fd74e161b8aa8dcf', '60bf126cc51e2529'], ['df4c576d46c52886', '6f549e89d232fa17'], ['b44c33c234b364f1', 'ca79c4179138c593'], ['9f65ccc3c9884d3e', '35153d8ffc129afb'], ['4032b8283df5f441', 'cc62df84a2a60385'], ['34f89684a9f71a8c', '9050b0eb2e3c1c03'], ['bc6dc24c33bb215e', '92a236f4dfba59e6'], ['3b87afc5cd2fbf5a', '18e06c3891fa3d43'], ['d8c74f3340d71dad', '33612f3b6fa223ab'], ['b5e565e34289b5fb', '5bf27c80e0d34fcd'], ['6b37c776d164899c', '4866465a93ea54e5'], ['361ed7d4449b8aca', '2d50bad52741668b'], ['9e20b5df7f28fb76', '19cb57857aaebfac'], ['3f01b0050fc13576', '11989734d6ccf08a'], ['5cb5f4ec01812b15', '3a74bea1e41ff69f'], ['8bb77e7c25d9a9cc', '819e6f2acbe13fcd'], ['99c8a621a9f85b90', '163e30c5d7d9f231'], ['301f91d08923d0e1', 'ad29964a7db48912'], ['a57a64e26c1e1f8f', '7d908e75301deaff'], ['12eb185fb7763175', '47b5b6f651976dd3'], ['fca00dd8286bc9fe', '699fa0b10ba531fa'], ['e9ebdefea208e6e6', '5ecc9dfd3ecaca04'], ['8166910d721de1ba', 'c8b20f8dfe995550'], ['9210bf7121ae8af8', 'df992c6a72541cd2'], ['536998f7f8495fcd', '7c0914287272254a'], ['e2ac7c978f9179e7', '5a7a5b749cfa1a2f'], ['6567e5dede78eb31', '6e7c9c3711954116'], ['54a69a3c0cb3e5c7', '72703319dfe0e862'], ['3ca21bfbc677daad', '9012458496c64503'], ['06e137939fb4834e', '0694b29d06b52f68'], ['549ef9139dc1b0bd', 'bf88a8bfe1f56981'], ['0897f434169f8133', 'b15367c4aa44a49e'], ['d996b5ceb705bb9b', '0f9496fb23858fd0'], ['dbcf8714b4898c09', '1788a0b89f38bf82'], ['2bc8fa5cb1ad52eb', 'b767c2d1a74d6dfb'], ['b3e71db2a5f9d4d9', '295a06a13c818f6c'], ['fd30b46d5f2786f3', '9f51fdc257398e51'], ['1c474f342c013ec5', '2967e6aed645bd1b'], ['0b019107c26aa083', '979944808c2eb4fa'], ['b7d74085989b2ea5', '92a4d8c37f046457'], ['3b6dacbcb51b4be8', '4889c0e9199a6abe'], ['b202c3677d25d1b8', '734ab506c90ac5bc'], ['08f9e37d39b59d01', '6be0ae27b6ce3133'], ['9ca9613f484d0ba6', '745eb6b0c66f3268'], ['43f366bed467aedd', '4f8b93cd973aa51e'], ['7e0098e44b5485d7', 'fe235e39ce5b7c10'], ['8bbe2e14b76a64eb', 'd00845043aa569b0'], ['8ffe5eb43e492681', 'a3fd00011a384c0c'], ['9c970cd823358011', '70907c28fed405b7'], ['d5876f62efc802ed', '29dc4b2995132757'], ['40f4c099c8c723b2', '69a381efabc33240'], ['e5208e70b556b25b', 'c7f77cf16772b507'], ['f1209ccddcf64b9c', 'ee5e1ce144f3addf'], ['34b2d8fcfadc39fe', 'd3d3455d380bdb6d'], ['365d496ab29d68f6', '6525465767d67f76'], ['5d4675d3f86a8f13', 'ffd14dcbe4239927'], ['b3671c810f9dd0e3', '93ab79e537a315fc'], ['f51755a1400cf3b5', '63cb13b229b232da'], ['1a414a01acbecda7', '254adfc7d15196fb'], ['ed4fac9fb727b0d7', '3113760cf49e7d0d'], ['5c532ef16e941dc3', '276298f09e6fb861'], ['34a1034ba9b7acfa', '6b56934aeb840000'], ['587e166c47827c31', '53b81364452a4adb'], ['bb75a92cfe72fb36', 'ac3903e37f5e0209'], ['02d7819e28fe05b7', 'c932c5b7a5433b9b'], ['db138d35d55b645c', '513b8c4110ee8319'], ['8ee84fc93ddbf811', 'e75ec4722f7e2ebe'], ['4e5e16d943e85516', '107fa4fdfc474a5c'], ['fbdb545b9a868613', '835968d7daa270dc'], ['fc89aa7d5aa5f144', '54b357a9f71c1cf1'], ['8c5da5d210692371', 'bd63235cb3463306'], ['c6699a04230335a8', '1d59e0a0a134066a'], ['f84470f50b84602c', '8fb1af6501789874'], ['33b77d25c5b81f1f', '484cfcc45a77c56c'], ['efc5c144e5aad0ad', 'eea5eb85b04b2c24'], ['ad8ab5b37577c096', '341aded55cd9810b']]




A = affine_gf2()
for t in D:
	A.observe(t[0],t[1])
A.transform()



