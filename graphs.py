import numpy as np

class graph:

    def __init__(self):
        self.edgeSet = {

            0 : []
        }
        ##
        ##  Instantiating this as [0] is mathematically inaccurate but it is the only way
        ##  that Python syntax will allow us to do this. This will not be a problem in practice 
        ##  because we will never use a graph instantiated directly from this without further 
        ##  modification in an actual algorithm. 
        ##
        self.incidenceMatrix = np.array([0])

        self.adjacencyMatrix = np.array([0])

        self.vertexNames = []

        self.vertexCount = 1

        self.edgeCount = 0

        self.edgesWeighted = False

        self.vertexNames = []


    ##  
    ##  Function Description: This function allows the user to see the data members and methods of the 
    ##  graph class. 
    ##  
    def classDescription(self):
        print("The graph class is a class used to define data members and functionality of graphs: ")
        print("Included in this class are the following data members: \n")

        print("edgeSet: a dictionary in which the keys are numbers, representing vertices, while the values are ")
        print("arrays of numbers, representing the vertices connected to the key vertex.\n")
        
        print("vertexNames: An array of strings which contains the names associated with each vertex ")
        print("of the graph. \n")

        print("edgesWeighted: a boolean variable which is equal to True if the edges of the graph have weights and ")
        print("which is equal to False if they do not. \n")

        print("vertexCount: a positive integer which is equal to the number of vertices in the graph. \n")

        print("edgeCount: a positive integer which is equal to the number of edges in the graph. \n")

        print("incidenceMatrix: a two-dimensional numpy array where each collumn represents an edge. The ")
        print("coefficient incidenceMatrix[i][j] is 1 if the edge represented leaves from vertex i and ")
        print("is -1 if the edge represented arrives at vertex j. Otherwise it is 0. There are edgeCount ")
        print("columns and vertexCount rows in this matrix.\n")

        print("adjacencyMatrix: a two-dimensional numpy which is square and of size vertexCount. The coefficient ")
        print("adjacencyMatrix[i][j] is zero if no edge connects from vertex i to vertex j. If the graph has ")
        print("weighted edges and there is an edge from vertex i to vertex j then adjacencyMatrix[i][j] will equal ")
        print("the weight of the edge. If the graph is unweighted and there is an edge connected from vertex i to ")
        print("vertex j then the coefficient adjacencyMatrix[i][j] will equal 1.\n")

        print("Included in this class are the following methods: \n")

        print("addEdge(self, source, sink, weight): This function adds an edge of the specified weight which connects ")
        print("from the specified (by index) source vertex to the specified (by index) sink vertex. \n")



    def fromUserInput(self):

        print("State the number of vertices in the graph: ")

        vertexCountSet = False

        while vertexCountSet == False:

            self.vertexCount = int(input())

            if self.vertexCount <= 0:
                print("ERROR: Only non-trivial defined graphs may be constructed: the number of vertices must be greater thanor equal to zero: ") 
            else:
                vertexCountSet = True

        self.edgeCount = 0

        for j in range(self.vertexCount):

            print(f"What is the name of vertex {j}")

            self.vertexNames.append(input())

            edges = []

            ##
            ##  There is no input validation here: the user's responsibility is to 
            ##  do their own input correctly. 
            ##
            for k in range(self.vertexCount):
                print(f"Does vertex {j} connect to vertex {k}?")

                print(">Y")
                print(">N")

                selection = input()

                if selection == 'Y' or selection == 'y':
                    edges.append(k)

            self.edgeSet.update({j : edges})

            self.edgeCount += len(edges)

            edges = []

        self.incidenceMatrix = np.zeros((self.vertexCount, self.edgeCount))

        currentColumn = 0

        for j in range(self.vertexCount):

            edges = self.edgeSet[j]

            for k in range(len(edges)):

                self.incidenceMatrix[j][currentColumn] = 1

                self.incidenceMatrix[edges[k]][currentColumn] = -1

                currentColumn += 1 

        print("Are the edges of the graph weighted? ")

        print(">Y")
        print(">N") 

        selection = input()

        self.adjacencyMatrix = np.zeros((self.vertexCount, self.vertexCount)) 
    
        if (selection == 'n' or selection == 'N'):

            self.edgesWeighted = False

            row = -1
            column = -1

            for j in range(self.edgeCount):
            
                i = 0

                while row == -1 or column == -1:
                
                    if self.incidenceMatrix[i][j] == 1:
                        row = i

                    elif self.incidenceMatrix[i][j] == -1:
                        column = i
                
                    i += 1 
            
                self.adjacencyMatrix[row][column] = 1

                row = column = -1

                i = 0

        else:

            decoyMatrix = np.zeros((self.vertexCount, self.vertexCount))

            self.edgesWeighted = True

            row = -1
            column = -1

            for j in range(self.edgeCount):
            
                i = 0

                while row == -1 or column == -1:
                
                    if self.incidenceMatrix[i][j] == 1:
                        row = i

                    elif self.incidenceMatrix[i][j] == -1:
                        column = i
                
                    i += 1 
            
                decoyMatrix[row][column] = 1

                row = column = -1

                i = 0

            for j in range(self.vertexCount):

                for k in range(self.vertexCount):

                    if decoyMatrix[j][k] == 1:

                        print(f"Select the weight on the edge from vertex {j} to vertex {k} ")

                        self.adjacencyMatrix[j][k] = float(input())
    ##
    ##  Function Description: 
    ##
    def addEdge(self, source, sink, weight):

        self.edgeCount += 1

        self.adjacencyMatrix[source][sink] = weight

        edges = self.edgeSet[source]

        edges.append(sink)

        edges.sort()

        newIncidenceMatrix = np.zeros((self.vertexCount, self.edgeCount))

        currentColumn = 0

        for j in range(self.vertexCount):

            edges = self.edgeSet[j]

            for k in range(len(edges)):

                newIncidenceMatrix[j][currentColumn] = 1

                newIncidenceMatrix[edges[k]][currentColumn] = -1

                currentColumn += 1 

        self.incidenceMatrix = newIncidenceMatrix

            

        
