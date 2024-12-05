import numpy as np
import matplotlib.pyplot as plt

class TrussStructure:

    def __init__(self):
        pass

    def DefineGeometry(self, baseShape, radius, height, rings):
        self.columns = []
        ringPitch = height / (rings+1)

        if baseShape == 'square':
            self.baseNodes = [[radius, radius, 0], [-radius, radius, 0], [-radius, -radius, 0], [radius, -radius, 0]]
        else:
            basePitch = (2 * np.pi) / baseShape
            self.baseNodes = []
            for i in range(baseShape):
                baseNode = [radius * np.cos(basePitch*i), radius * np.sin(basePitch*i), 0]
                self.baseNodes.append(baseNode)
            
        for baseNode in self.baseNodes:
            column = []
            for i in range(rings+2):
                column.append([baseNode[0], baseNode[1], ringPitch*i])
            self.columns.append(column)

        self.elements = []

        # #vertical connections between nodes in the same column
        for column in self.columns:
            for i in range(len(column)-1):
                element = [column[i], column[i+1]]
                self.elements.append(element)

        #horizonal connections between nodes on the same vertical level
        for i in range(len(self.columns)):
            for j in range(len(self.columns[0])):
                element = [self.columns[i][j], self.columns[(i+1) % (len(self.columns))][j]]
                self.elements.append(element)
        
        # #diagonal connections
        for i in range(len(self.columns)):
            for j in range(len(self.columns[0])-1):
                element = [self.columns[i][j], self.columns[(i+1) % len(self.columns)][j+1]]
                self.elements.append(element)


truss = TrussStructure()
truss.DefineGeometry(20, 1, 6, 2)

fig = plt.figure()
ax = plt.axes(projection ='3d')

for element in truss.elements:
    x = [element[0][0], element[1][0]]
    y = [element[0][1], element[1][1]]
    z = [element[0][2], element[1][2]]
    ax.plot3D(x, y, z, 'black')

plt.show()