import numpy as np
from graph_tool.all import *
from collections import Counter
import os
import sys
import time

def Vertex_Layer(hierarchy):
	if hierarchy == 0:
		return [1, 1]
	else:
		prevList = Vertex_Layer(hierarchy - 1)
		l = 2*len(prevList) - 1 
		newList = [0] * l
		for x in range(l):
			if(x%2 == 0):
				newList[x] = prevList[int(x/2)]
			else:
				newList[x] = 2**hierarchy
		return newList

def Rudin_Shapiro(hierarchy):
	letter = ['A', 'B', 'C', 'D']
	inputSequence = np.random.choice(letter)
	for _ in range(int(hierarchy)):
		outputSequence = ""
		for char in inputSequence:
			if char == 'A':
				outputSequence += 'A' + 'B'
			elif char == 'B':
				outputSequence += 'A' + 'C'
			elif char == 'C':
				outputSequence += 'D' + 'B'
			elif char == 'D':
				outputSequence += 'D' + 'C'
		inputSequence = outputSequence
	return outputSequence

def Edge_Generator(g, initialVertex, finalVertex, edgeType, edgeProbability):
	
	pA3 = 0.8376738751139919
	pA1 = 0.07757137896493078

	pB3 = 0.599877104432572
	pB1 = 0.1746402319986869

	pC3 = 0.32922402701244635
	pC1 = 0.24455644234185844

	pD3 = 0.20089067267992552
	pD1 = 0.24731761894342835
	
	if(edgeType=="A"):
		p1 = pA1
		p3 = pA3
	elif(edgeType=="B"):
		p1 = pB1
		p3 = pB3
	elif(edgeType=="C"):
		p1 = pC1
		p3 = pC3
	else:
		p1 = pD1
		p3 = pD3
		
	if edgeProbability < p3:
		g.add_edge(g.vertex_index[initialVertex], g.vertex_index[finalVertex])
		g.add_edge(g.vertex_index[finalVertex], g.vertex_index[initialVertex])
			
	elif p3 <= edgeProbability < p3 + p1:
		g.add_edge(g.vertex_index[initialVertex], g.vertex_index[finalVertex])
			
	elif p3 + p1 <= edgeProbability < p3 + 2 * p1:
		g.add_edge(g.vertex_index[finalVertex], g.vertex_index[initialVertex])
		
	return

def Lattice_Generator(typeLayer, vertexLayer):
	
	np.random.seed()
	g = Graph()
	g.add_vertex(sum(vertexLayer))
	edgeProbability = np.random.uniform(0,1,int(sum(vertexLayer)/2+ sum(vertexLayer)-2))
	edgeIndex = 0
	vertexPerLayer = 0
	
	for layerIndex in range(len(vertexLayer)-1):
		
		if(vertexLayer[layerIndex] < vertexLayer[layerIndex + 1]):
			
			numberConnection = int(vertexLayer[layerIndex+1]/ vertexLayer[layerIndex])
			
			for vertexLayerIndex in range(vertexLayer[layerIndex]):
				
				for connectionIndex in range(numberConnection):
					
					initialVertex = vertexPerLayer + vertexLayerIndex
					finalVertex = vertexPerLayer + vertexLayer[layerIndex] + vertexLayerIndex*numberConnection + connectionIndex
					Edge_Generator(g, initialVertex, finalVertex, typeLayer[layerIndex], edgeProbability[edgeIndex])
					edgeIndex += 1
					
			vertexPerLayer += vertexLayer[layerIndex]
			
		else:
			
			numberConnection = int(vertexLayer[layerIndex]/ vertexLayer[layerIndex + 1])
			
			for vertexLayerIndex in range(int(vertexLayer[layerIndex]/numberConnection)):
				
				for connectionIndex in range(numberConnection):
					
					initialVertex = vertexPerLayer + connectionIndex + vertexLayerIndex*numberConnection
					finalVertex = vertexPerLayer + vertexLayer[layerIndex] + vertexLayerIndex
					Edge_Generator(g, initialVertex, finalVertex, typeLayer[layerIndex], edgeProbability[edgeIndex])
					edgeIndex += 1
					
			vertexPerLayer += vertexLayer[layerIndex]
			
	return g

def Percolation_Check(g, vertexLayer):
	
	compDown = graph_tool.topology.label_out_component(g, 0, label=None)
	compUp = graph_tool.topology.label_out_component(g, sum(vertexLayer)-1, label=None)
	
	aDown = compDown.a[0]
	bDown = compDown.a[-1]
	
	aUp = compUp.a[0]
	bUp = compUp.a[-1]
	if(aDown == bDown & aUp == bUp):
		return int(3)
	elif(aDown == bDown):
		return int(1)
	elif(aUp == bUp):
		return int(2)
	return int(0)

def Find_Index(arr, element):
	indices = np.where(arr == element)[0]
	if indices.size > 0:
		return indices[0]
	else:
		return -1  # Return -1 if the element is not found in the array
	
def GSCC_GOUT(g, vertexLayer):
	comp, hist, is_attractor = graph_tool.topology.label_components(g, attractors=True)
	gscc = max(hist)
	fractionGSCC = gscc/sum(vertexLayer)
	
	clusterIndex = Find_Index(hist, gscc)
	vertexIndex = Find_Index(comp.a, clusterIndex)
	
	compOUT = graph_tool.topology.label_out_component(g, vertexIndex, label=None)
	gout = np.sum(compOUT.a == 1)
	fractionGOUT = gout/sum(vertexLayer)
	
	return fractionGSCC, fractionGOUT

def Result_Generator(typeLayer, vertexLayer):
	g = Lattice_Generator(typeLayer, vertexLayer)
	percolation = Percolation_Check(g, vertexLayer)
	fractionGSCC, fractionGOUT = GSCC_GOUT(g, vertexLayer)
	return [percolation, fractionGSCC, fractionGOUT]

def Check_File(index, numberFiles, hierarchy):
	
	dirPath = "data/"
	folderName = f"N{hierarchy}/"
	fullDirPath = dirPath + folderName
	
	if not os.path.exists(fullDirPath):
		os.mkdir(fullDirPath)
	
	fileName = f"N{hierarchy}_results_{numberFiles + 1 + index}.npy"
		
	if os.path.exists(fullDirPath + fileName):
		#print("File exists!")
		return None
		
	else:
		emptyArray = np.array([])
		np.save(fullDirPath + fileName, emptyArray)
		return fileName

def Save_Data(hierarchy, fileName, results):
	dirPath = "data/"
	folderName = f"N{hierarchy}/"
	fullDirPath = dirPath + folderName
	results = np.array(results)
	np.save(fullDirPath + fileName, results)
	return
def Simulation(n, hierarchy, vertexLayer):
	
	results = []
	
	for _ in range(n):
		typeLayer = Rudin_Shapiro(hierarchy)
		result = Result_Generator(typeLayer, vertexLayer)
		results.append(result)
		
	return results

def Save_Data(hierarchy, fileName, results):
	dirPath = "data/"
	folderName = f"N{hierarchy}/"
	fullDirPath = dirPath + folderName
	results = np.array(results)
	np.save(fullDirPath + fileName, results)
	return

def main():
	
	n = 100
	generations = [7]
	numberFiles = 0
	numberRuns = int(sys.argv[1])
	
	for hierarchy in generations:	
		vertexLayer = Vertex_Layer(hierarchy)
	
		for runIndex in range(numberRuns):
			fileName = Check_File(runIndex, numberFiles, hierarchy)
			
			if(fileName == None):
				continue
		
			results = Simulation(n, hierarchy, vertexLayer)
			Save_Data(hierarchy, fileName, results)

	return

if __name__ == "__main__":
    main()
