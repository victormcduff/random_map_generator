import math
import numpy as np
from random import *
from scipy.spatial import Voronoi, voronoi_plot_2d


class DrawPoligons: #class to draw polygons 

	def __init__(self, seed_points, x_size, y_size, norm_steps):
		self.seed_points = seed_points
		self.x_size = x_size
		self.y_size = y_size
		self.norm_steps = norm_steps


	def starter_points(self):

		self.points = np.zeros((self.seed_points,2))

		for i in range(self.seed_points):
			x_coord = randint(0,self.x_size)
			y_coord = randint(0,self.y_size)

			self.points[i,0] = x_coord
			self.points[i,1] = y_coord


	def region_extractor(self):

		x = []
		y = []

		vor = Voronoi(self.points)
		regions = vor.regions
		vertices = vor.vertices

		regions.remove([])

		for i in range(np.size(regions)):

			line_region = np.asarray(regions[i])
	
			x_line = []
			y_line = []

			if np.all(line_region >= 0):

				for k in range(np.size(line_region)):
					a = line_region[k]
					x_line.append(vertices[a,0])
					y_line.append(vertices[a,1])

				x.append(x_line)
				y.append(y_line)

		return x, y


	def find_centroids(self):

		centroids = np.zeros((len(self.regionNodsX),2))

		for i in range(len(self.regionNodsX)):
			centre_x = sum(self.regionNodsX[i])/len(self.regionNodsX[i]) 
			centre_y = sum(self.regionNodsY[i])/len(self.regionNodsY[i])

			centroids[i,0] = centre_x
			centroids[i,1] = centre_y

		return centroids


	def normalise(self):

		self.starter_points()

		for i in range(self.norm_steps):
			self.regionNodsX, self.regionNodsY = self.region_extractor()
			self.points = self.find_centroids()

		self.regionNodsX, self.regionNodsY = self.region_extractor()
		self.points = self.find_centroids()

		return self.points, self.regionNodsX, self.regionNodsY

