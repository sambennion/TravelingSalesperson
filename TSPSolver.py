#!/usr/bin/python3

import queue
from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))




import time
import numpy as np
from TSPClasses import *
import copy
import heapq
import itertools



class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution,
		time spent to find solution, number of permutations tried during search, the
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''

	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this
		algorithm</returns>
	'''
	#Time Complexity is O(n) where n is the len(cities)
	#Space Complexity is O(n) where n is the size of cities
	def short_edge(self, currentCity, cities):
		shortest = float("Inf")
		index = 0
		for i in range(len(cities)):
			cost = currentCity.costTo(cities[i])
			if cost < shortest:
				shortest = cost
				index = i
		return cities.pop(index)


	def greedy( self,time_allowance=60.0 ):
		# Set cities to shallow copy of getCities(). 
		# Copy because we'll be popping them off in my algorithm 
		# so we don't affect the original cities object.
		cities = copy.deepcopy(self._scenario.getCities())
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		route = []
		while not foundTour and time.time()-start_time < time_allowance:
			
			currentCity = cities.pop(count)
			route = [currentCity]
			while not foundTour and time.time()-start_time < time_allowance and len(cities) > 0:
				currentCity = self.short_edge(currentCity, cities)
				route.append(currentCity)
			bssf = TSPSolution(route)
			# print(route)
			# bssf in the greedy sense will just be the first one that 
			# find a tour.
			if(bssf.cost < math.inf):
				foundTour = True
			else:
				cities = copy.deepcopy(self._scenario.getCities())
			count += 1
		end_time = time.time()
		results = {}
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results



	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints:
		max queue size, total number of states created, and number of pruned states.</returns>
	'''

	def reduce_matrix(self, matrix):
		n = len(matrix)
		# print(matrix)
		min_of_rows = np.amin(matrix, axis=1)
		# print(min_of_rows)
		#reduce rows first
		for i in range(n):
			if(min_of_rows[i] < np.inf):
				matrix[i:] -= min_of_rows[i]
			else:
				min_of_rows[i] = 0 #This is for summing at the end.
		matrix[matrix<0] = 0
		# reduce each column
		min_of_columns = np.amin(matrix, axis=0)
		# print(np.shape(min_of_columns))
		# print(matrix)
		# print(min_of_columns)
		for i in range(n):
			if(min_of_columns[0, i] < np.inf):
				matrix[:,i] = matrix[:,i] - min_of_columns[0, i]
			else:
				min_of_columns[0, i] = 0
		matrix[matrix<0] = 0
		cost = np.sum(min_of_rows) + np.sum(min_of_columns)
		# print(cost)
		return matrix, cost


	def branchAndBound( self, time_allowance=60.0 ):
		cities = copy.deepcopy(self._scenario.getCities())
		n = len(cities)
		#Matrix uses O(n^2 space)
		matrix = np.matrix(np.ones((n,n)) * np.inf)
		#O(n^2) time to initialize this array with distances.
		for i in range(n):
			for j in range(n):
				matrix[i, j] = cities[i].costTo(cities[j])
		#Running a greedy first thing off the bat will create a baseline for upperbound.
		#It tells us a few things that will save us time, it gives us a starting node that has
		#a solution, and it will allow us to prune out routes that aren't fruitful early on.
		
		results = self.greedy(time_allowance=60.0)
		state_count = n
		pruned = 0
		lower_bound = results['cost']
		starting_node = results['count'] - 1
		bssf = results['soln']
		#Runs in n^2 time
		start_time = time.time()
		reduced_matrix, initial_cost = self.reduce_matrix(copy.deepcopy(matrix))
		# currentCity = cities[starting_node]

		# reduced_matrix[:,0] += np.inf
		# print(reduced_matrix)
		pq = []
		#Somehow need to keep track of parent cost, 
		# added matrix cost, and the reducucing cost.
		
		#Initial loop to add nodes and costs to pq. Runs O(n) time
		for i in range(len(cities)):
			if start_time - time_allowance <= 0:
					break
			edge_cost = reduced_matrix[starting_node, i]
			if edge_cost == np.inf:
				continue
			new_reduced = copy.deepcopy(reduced_matrix)
			new_reduced[starting_node] += np.inf
			new_reduced[:,i] += np.inf
			new_reduced[i, starting_node] = np.inf
			new_reduced, reduction_cost = self.reduce_matrix(copy.deepcopy(new_reduced))
			# print(new_reduced)
			node_cost = edge_cost + initial_cost + reduction_cost
			
			if node_cost >= lower_bound:
				
				# print("Node: " + str(node_cost) + " Upper: " + str(lower_bound))
				pruned += 1
			else:
				# print(i)
				heapq.heappush(pq, (node_cost, i, copy.copy(new_reduced), [starting_node, i]))
				# print(len(pq))
		# print(pq[2])
		while len(pq) > 0 and start_time - time_allowance > 0:
			initial_cost, index, reduced_matrix, path = heapq.heappop(pq)
			if initial_cost >= lower_bound:
				pruned += 1
				continue
			for i in range(len(cities)):
				
				if start_time - time_allowance <= 0:
					break
				#We don't want to cycle back to the starting node
				if i == starting_node:
					continue
				edge_cost = reduced_matrix[index, i]
				if edge_cost == np.inf:
					continue
				#same as above
				new_reduced = copy.deepcopy(reduced_matrix)
				new_reduced[index] += np.inf
				new_reduced[:,i] += np.inf
				new_reduced[i, index] = np.inf
				new_reduced, reduction_cost = self.reduce_matrix(copy.copy(new_reduced))
				#Same
				node_cost = edge_cost + initial_cost + reduction_cost
				if node_cost >= lower_bound:
					pruned += 1
				else:
					path.append(i)
					if len(path) < n:
						heapq.heappush(pq, (node_cost, i, copy.copy(new_reduced), path))
					else:
						#create new lower bound bssf
						lower_bound = node_cost
						# print(path)
						# print(matrix)
						# route = [cities[j] for j in path]
						bssf = TSPSolution([cities[j] for j in path])
						# print(bssf.route)


		end_time = time.time()
		results['cost'] = lower_bound
		results['time'] = end_time - start_time
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = pruned

			# else:
				# new_reduced[currentCity] = 

		
		return results



	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution,
		time spent to find best solution, total number of solutions found during search, the
		best solution found.  You may use the other three field however you like.
		algorithm</returns>
	'''

	def fancy( self,time_allowance=60.0 ):
		pass
