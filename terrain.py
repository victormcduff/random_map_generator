##Random map generation
##Use opensimplex noise with voronoi polygons 

##Arkadiy Ukolov (2019) 

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from random import *
from scipy.spatial import Voronoi, voronoi_plot_2d
from opensimplex import OpenSimplex
from poligons_drawing import DrawPoligons

from scipy.spatial import Delaunay
from scipy.spatial import KDTree
from scipy.optimize import fsolve


def noise_gen_simplex(nx, ny):
	return gen.noise2d(nx, ny) / 2 + 0.5


def dictionaryGenerator(points, x, y):
	dictionary = {}

	for i in range(len(x)):
		for j in range(len(x[i])):

			dictionary[i,j] = {
					'Neibouring poligons': [],
					'Neibouring nods': []
				}

			adjacentNods, nodsNumber = findAdjacentNods(i, j, x, y)

			for k in range(nodsNumber):
				dictionary[i,j]['Neibouring poligons'].append(adjacentNods[k][0])
				dictionary[i,j]['Neibouring nods'].append(adjacentNods[k][1])

	return dictionary


def generateRegionNoise(points, x, y):
	global x_size, y_size

	noiseRegionNods = []
	noiseRegionCentres = np.zeros(np.size(points[:,0]))

	for i in range(np.size(x)):

		line_noise = []

		for k in range(np.size(x[i])):

			nx = x[i][k]/x_size
			ny = y[i][k]/y_size

			line_noise.append(0.4 * noise_gen_simplex(nx, ny) + 0.9 * noise_gen_simplex(2 * nx, 2 * ny) + 0.7 * noise_gen_simplex(4 * nx, 4 * ny))

		noiseRegionNods.append(line_noise) #the noise is arranged with the regions


	for i in range(np.size(noiseRegionNods)):
		noiseRegionCentres[i] = np.sum(noiseRegionNods[i])/np.size(noiseRegionNods[i])

	noiseRegionCentres = (noiseRegionCentres - np.min(noiseRegionCentres)) / (np.max(noiseRegionCentres) - np.min(noiseRegionCentres))

	for i in range(np.size(noiseRegionNods)):
			noiseRegionNods[i] = (noiseRegionNods[i] - np.min(noiseRegionNods[i])) / (np.max(noiseRegionNods[i]) - np.min(noiseRegionNods[i]))

	return noiseRegionNods, noiseRegionCentres



def land_distribution(points, noiseRegionCentres, x, y, levels): #points and regions assumed to be ARRANGED

	poligonType = []
	nodType = []

	for i in range(np.size(points[:,0])):

		line = []

		if noiseRegionCentres[i] >= levels[3]:
			plt.fill(x[i], y[i], '#fcfcfc') #snow peaks

			poligonType.append('snow')

			for k in range(np.size(x[i])):
				line.append('snow')

			nodType.append(line)


		elif noiseRegionCentres[i] >= levels[2]:
			plt.fill(x[i], y[i], '#7f7f7f') #mountains

			poligonType.append('mountain')

			for k in range(np.size(x[i])):
				line.append('mountain')

			nodType.append(line)

		elif noiseRegionCentres[i] >= levels[1]:
			plt.fill(x[i], y[i], '#2f8200') #land

			poligonType.append('land')

			for k in range(np.size(x[i])):
				line.append('land')

			nodType.append(line)

		elif noiseRegionCentres[i] >= levels[0]:
			plt.fill(x[i], y[i], '#004cff') #sea

			poligonType.append('water')

			for k in range(np.size(x[i])):
				line.append('water')

			nodType.append(line)

		else:
			plt.fill(x[i], y[i], '#0036b7') #ocean

			poligonType.append('water')

			for k in range(np.size(x[i])):
				line.append('water')

			nodType.append(line)

	return poligonType, nodType



def findAdjacentNods(k, n, x, y):
	adjacentNods = [] #here will lie the arranged coordinates in form: [row, column]
	nodsNumber = 0

	for a in range(len(x)):
		for b in range(len(x[a])):
			if (x[a][b] == x[k][n]) and (y[a][b] == y[k][n]):

				if (len(x[a]) > (b+1)):
					adjacentNods.append([a, b+1])

				else:
					adjacentNods.append([a, 0])

				if ((b-1) >= 0):
					adjacentNods.append([a, b-1])
				else:
					adjacentNods.append([a, len(x[a])-1])

				nodsNumber += 2

	return adjacentNods, nodsNumber

def findAdjacentCentres(points, x, y):

	adjacentCentres = []

	for i in range(np.size(points[:,0])):

		adjacentCentresLine = []
		centresNumber = 0
		centreID = 0

		for k in range(len(x[i])):
			adjacentNods, nodsNumber = findAdjacentNods(i, k, x, y)

			for n in range(nodsNumber):
				adjacentCentresLine.append(adjacentNods[n][0])

		adjacentCentres.append(adjacentCentresLine)

	return adjacentCentres #it IS consistent with the, well, other stuff. ID of a line is, well, ID, y'know


def coastGeneration(points, noiseRegionCentres, x, y, poligonType, nodType):
	global poligonStructure

	for i in range(np.size(points[:,0])):
		if (poligonType[i] == 'land'):
			
			water = False
			land = False

			for j in range(len(x[i])):
				for k in range(len(poligonStructure[i,j]['Neibouring poligons'])):
					ID = poligonStructure[i,j]['Neibouring poligons'][k]

					if (poligonType[ID] == 'water'):
						water = True

					if (poligonType[ID] != 'water'):
						land = True

			if water and land:
				plt.fill(x[i], y[i], '#fff8af')
				poligonType[i] = 'coast'

				for k in range(np.size(x[i])):
					nodType[i][k] = 'coast'


def river_generator(noiseRegionCentres, x, y, poligonType, nodType, riversNumber):
	global poligonStructure

	numberOfRivers = 0

	while numberOfRivers < riversNumber:

		k = randint(0,len(x)-1)
		n = randint(0,len(x[k])-1)

		trigger = 0
		flow = 1.5

		while (poligonType[k] != 'water') and (trigger != 1):
			
			neibouringPoligons = poligonStructure[k,n]['Neibouring poligons']
			neibouringNods = poligonStructure[k,n]['Neibouring nods']

			pr = 0
			for p in range(len(neibouringPoligons)):

				minimum = noiseRegionCentres[neibouringPoligons[0]]

				if noiseRegionCentres[neibouringPoligons[p]] <= minimum:
					minimum = noiseRegionCentres[neibouringPoligons[p]]
					pr = p

			rowToGo = neibouringPoligons[pr]
			columnToGo = neibouringNods[pr]

			if (nodType[rowToGo][columnToGo] == 'river_%d' % (numberOfRivers)): #so if you face YOURSELF

				poligonType[rowToGo] = 'water'
				plt.fill(x[rowToGo], y[rowToGo], '#2867fc', edgecolor = '#004cff', linewidth = flow)
				trigger = 1

				numberOfRivers += 1
			
			elif (poligonType[rowToGo] != 'water'):

				plt.plot([x[k][n], x[rowToGo][columnToGo]], [y[k][n], y[rowToGo][columnToGo]], '#004cff', linewidth = flow)
				nodType[rowToGo][columnToGo] = 'river_%d' % (numberOfRivers) #give all rivers its own number

				k = rowToGo
				n = columnToGo

				flow += 0.15 #to grow as you flow

			else:
				trigger = 1

				numberOfRivers += 1
 

def cityGenerator(x, y, poligonType, citiesNumber):
	global poligonStructure

	previous_city_number = -1

	numberOfCities = 0

	citiesID = [] #would store city's row in the poligon matrix and its unique COLOUR

	while numberOfCities != citiesNumber: #need to make some dependance on the terrain type...

		k = randint(0,len(x)-1)

		if (poligonType[k] != 'water'):
			plt.fill(x[k], y[k], 'k')
			poligonType[k] = 'city'

			citiesID.append([k, randint(0,255)/255, randint(0,255)/255, randint(0,255)/255])

			if (previous_city_number != -1):
				roadGenerator(k, previous_city_number, x, y, citiesID) #numer of the row

			previous_city_number = k
			numberOfCities += 1

	return citiesID


def roadGenerator(currentN, destinationN, x, y, citiesID):
	global poligonStructure, points, poligonType

	destinationX = points[destinationN,0]
	destinationY = points[destinationN,1]

	while currentN != destinationN:

		minimumX = abs(points[currentN,0] - destinationX)
		minimumY = abs(points[currentN,1] - destinationY)
		minimumLength = np.sqrt(minimumX ** 2 + minimumY ** 2)
		minimumN = currentN

		for j in range(len(x[currentN])):
			for k in range(len(poligonStructure[currentN,j]['Neibouring poligons'])):
				ID = poligonStructure[currentN,j]['Neibouring poligons'][k]

				lengthX = abs(points[ID,0] - destinationX)
				lengthY = abs(points[ID,1] - destinationY)

				length = np.sqrt(lengthX ** 2 + lengthY ** 2)

				if length < minimumLength:
					minimumLength = length
					minimumN = ID
		
		if (poligonType[minimumN] == 'water'):
			plt.plot([points[currentN,0], points[minimumN,0]], [points[currentN,1], points[minimumN,1]], 'k--', linewidth = 1.5)

			if (poligonType[currentN] != 'water'):
				plt.fill(x[currentN],y[currentN],'k')
				poligonType[currentN] = 'city'

				citiesID.append([currentN, randint(0,255)/255, randint(0,255)/255, randint(0,255)/255])

		elif (poligonType[currentN] == 'water') and (poligonType[minimumN] != 'water'):
			plt.fill(x[minimumN],y[minimumN],'k') 
			poligonType[minimumN] = 'city'

			citiesID.append([minimumN, randint(0,255)/255, randint(0,255)/255, randint(0,255)/255])

			plt.plot([points[currentN,0], points[minimumN,0]], [points[currentN,1], points[minimumN,1]], 'k--', linewidth = 1.5)

		else:
			plt.plot([points[currentN,0], points[minimumN,0]], [points[currentN,1], points[minimumN,1]], 'k-',  linewidth = 2)

		currentN = minimumN


def countries(x, y, citiesID):
	global poligonStructure, poligonType

	countryID = np.zeros(len(x)) #the array of ALL the poligons, with the ID of the countries
	countryID = countryID - 1 #everything is -1, that means the land is empty
	occupiedID = []
	
	occupation = 0

	for i in range(len(citiesID)):
		for j in range(len(x[citiesID[i][0]])):
			for k in range(len(poligonStructure[citiesID[i][0],j]['Neibouring poligons'])):

				ID = poligonStructure[citiesID[i][0],j]['Neibouring poligons'][k]

				if (countryID[ID] == -1) and (poligonType[ID] != 'city'):
					
					countryID[ID] = citiesID[i][0]

					colour = citiesID[i][1:4]

					plt.fill(x[ID], y[ID], color = colour, alpha = 0.5)

					occupation += 1
					occupiedID.append([ID,citiesID[i][1],citiesID[i][2],citiesID[i][3]])

	for i in range(len(occupiedID)):
		for j in range(len(x[occupiedID[i][0]])):
			for k in range(len(poligonStructure[occupiedID[i][0],j]['Neibouring poligons'])):

				ID = poligonStructure[occupiedID[i][0],j]['Neibouring poligons'][k]

				if (ID not in occupiedID) and (poligonType[ID] != 'city'):

					colour = occupiedID[i][1:4]

					plt.fill(x[ID], y[ID], color = colour, alpha = 0.5)
					occupiedID.append([ID,occupiedID[i][1],occupiedID[i][2],occupiedID[i][3]])



def show_poligons(points, xAxis, yAxis, x, y):
	#plt.plot(points[:, 0], points[:, 1], 'ro')
	
	x = np.asarray(x)
	y = np.asarray(y)

	"""
	for i in range(np.size(x)): #just draws poligons. Replacable by colouring and vice versa
		plt.plot(x[i], y[i], 'k-') #this replaces rivers, so DON'T

		fix_end_linesX = []
		fix_end_linesX.append(x[i][0])
		fix_end_linesX.append(x[i][-1])

		fix_end_linesY = []
		fix_end_linesY.append(y[i][0])
		fix_end_linesY.append(y[i][-1])

		plt.plot(fix_end_linesX, fix_end_linesY, 'k-') #replaces rivers as well
	"""

	plt.axis((10,xAxis,10,yAxis))
	plt.show()

##NEED TO STORE NEIBOURS FOR THE CENTRES OF POLIGONS, NOT THE NODS
##NEED TO CREATE BIOMES! (temperature, etc)
##MAKE SMALL ISLANDS IN OCEANS?
##COUNTRIES
##NAMES FOR EVERYTHING (at least for the cities, or I need to create LOTS of new entities, basically rewriting the entire code)
#need to create different layers to show

seed_points = 1000

x_size = 500
y_size = 500

xAxis = 490 #axis to show
yAxis = 490

levels = [0.3, 0.5, 0.7, 0.9] #limits: ocean - 0.3 - sea - 0.4 - land - 0.7 - mountains - 0.9 - snow peaks
citiesNumber = seed_points/100 #numer of cities in the world
riversNumber = seed_points/100 #numer of rivers in the world

norm_steps = 5

print('Enter seed:')
seed = int(input())
gen = OpenSimplex(seed)

print('Normalising poligons...')
poligons = DrawPoligons(seed_points, x_size, y_size, norm_steps)
points, regionNodsX, regionNodsY = poligons.normalise()

print('Generating dictionary...')
poligonStructure = dictionaryGenerator(points, regionNodsX, regionNodsY)

print('Evaluating noise...')
noiseRegionNods, noiseRegionCentres = generateRegionNoise(points, regionNodsX, regionNodsY) #generate some noise for terrain elevation

print('Distributing land...')
poligonType, nodType = land_distribution(points, noiseRegionCentres, regionNodsX, regionNodsY, levels) #distribute land/water according to the evaluation

print('Finding coast...')
coastGeneration(points, noiseRegionCentres, regionNodsX, regionNodsY, poligonType, nodType)

print('Drawing rivers...')
river_generator(noiseRegionCentres, regionNodsX, regionNodsY, poligonType, nodType, riversNumber)

print('Building cities...')
citiesID = cityGenerator(regionNodsX, regionNodsY, poligonType, citiesNumber)

#print('Making countries...')
#countries(regionNodsX, regionNodsY, citiesID)

print('Processing image...')
show_poligons(points, xAxis, yAxis, regionNodsX, regionNodsY)