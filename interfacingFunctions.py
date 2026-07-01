##  
##  This file is inteded to contain the files which will be used to manage the execution of 
##  algorithms that cannot always be assumed to be constant in their inputs, outputs, or 
##  steps. Examples include user input functions, data management functions, or compact 
##  algorithms. 
##
import numpy as np
import matrixFunctions
import graphs
import pulp

class linearProgram:

    ##
    ##  Function Description: This is the default constructor of the linearProgram class. 
    ##  this is used only to instantiate memory for a linearProgram object and should not,
    ##  under any circumstances, be used without later calling one of the assignment functions
    ##  of the class. 
    ##
    def __init__(self):
        self.rowCount = 2

        self.columnCount = 2

        self.constraintMatrix = np.array([[1, 0], [0, 1]])

        self.constraintVector = np.array([0, 0])

        self.costFunction = np.array([0, 0])

        self.sucessfullAllocation = False

    ##
    ##  linearProgram Class User Input Constructor: 
    ##
    ##  Method Description: This is the user input constructor of the linearProgram class. 
    ##  This will be used to allow for the construction of a linearProgram object using 
    ##  only user input
    ##
    ##  Precondition: The method is called. It will be assumed that the user only argues
    ##  integral type arguments for the requested inputs: a better language than Python  
    ##  would have some way to implement safety features for the event in which they 
    ##  don't, but we are not writting in a better language than Python. 
    ##
    ##  Postcondition: An object of the linearProgram class is constructed to the users
    ##  specification. 
    ##
    def createProgramFromUserInput(self):
        
        rowCountSpecified = False

        while rowCountSpecified == False:

            print("State the number of rows of the matrix: ")    

            self.rowCount = int(input()) 

            if self.rowCount <= 0:

                print("ERROR: The number of rows of the matrix must be a positive integer: ")

            else:
                rowCountSpecified = True

        columnCountSpecified = False

        while columnCountSpecified == False:

            print("State the number of columns of the matrix: ")

            self.columnCount = int(input())

            if self.columnCount < 0:

                print("ERROR: The number of columns of the matrix must be a positive integer: ")

            else:
                columnCountSpecified = True

        print("Define the matrix of constraints for your standard form linear program: the left-hand-side of the equation Ax = b: ")

        matrix = createMatrixFromUserInput(self.rowCount, self.columnCount)

        self.constraintMatrix = matrix

        print("Define the vector of constraints for your standard form linear program: the right-hand-side of the equation Ax = b: ")

        self.constraintVector = createVector(len(self.constraintMatrix))

        print("Define the cost function for your standard form linear program: ")

        self.costFunction = createVector(len(self.constraintMatrix[0]))

        self.sucessfullAllocation = True

    ##
    ##  Function descripion: This function prompts the user for information on the LP and creates the standard
    ##  form version of said LP given this. 
    ##
    def createFromUserInputGeneralForm(self):
        print("State the dimension of the vector space as an integer: ")

        tempColumnCount = int(input())

        tempCostFunction = np.zeros(tempColumnCount, dtype=float)

        for j in range(tempColumnCount):
            print(f"Please state the cost of variable {j}")

            tempCostFunction[j] = float(eval(input()))

        print("State the number of constraints: ")

        self.rowCount = int(input())

        self.constraintVector = np.zeros(self.rowCount, dtype=float)

        tempConstraintMatrix = np.zeros((self.rowCount, tempColumnCount), dtype=float)

        geqRows = []

        eqRows = []

        leqRows = []

        for j in range(self.rowCount):

            print(f"What form is the constraint {j} of?:")
            print("1: >=")
            print("2: = ")
            print("3: =<")

            selection = int(input())

            match selection:
                case 1:
                    geqRows.append(j)
                case 2:
                    eqRows.append(j)
                case 3:
                    leqRows.append(j)

            print(f"State the coefficients on the left hand side of the equation a_{j}x (>=, =, =<) b_{j}")

            for k in range(tempColumnCount):
                print(f"Coefficient {k} is : ")
                tempConstraintMatrix[j][k] = float(eval(input()))

            print(f"State the right hand side of the same equation: ")
            self.constraintVector[j] = float(eval(input())) 

        self.columnCount = tempColumnCount + len(geqRows) + len(leqRows)

        self.constraintMatrix = np.zeros((self.rowCount, self.columnCount), dtype=float)

        for j in range(len(tempConstraintMatrix)):
            for k in range(len(tempConstraintMatrix[0])):
                self.constraintMatrix[j][k] = tempConstraintMatrix[j][k]   
            ##
            ##  This if-else feels like it flows into the opposite cases that it should flow to, 
            ##  but it works, so I'm going to keep it like this. 
            ##
            if j in geqRows:
                self.constraintMatrix[j][j + tempColumnCount] = -1
            elif j in leqRows:
                self.constraintMatrix[j][j + tempColumnCount] = 1
        
        self.costFunction = np.zeros(self.columnCount, dtype=float)

        for j in range(len(tempCostFunction)):
            self.costFunction[j] = tempCostFunction[j]

        self.sucessfullAllocation = True
        


    def createProgramFromGraph(self, G : graphs.graph):
        ##
        ##  If the graph used to create this program is weighted then the program created will be a 
        ##  minimum cost flow problem. 
        ##
        if G.edgesWeighted:

            print("Is the flow accross the arcs of the graph weighted from below? ")

            print(">Y")
            print(">N")

            coconut = input()
            ##
            ##  This part will create the linear program for if there are both an upper and lower bound to the 
            ##  flow accross an edge. If there is only an upper bound to flow accross an edge, go to the else 
            ##  half of this flow control. 
            ##
            if coconut == 'y' or coconut == 'Y':

                ##
                ##  Prompt for and scan in the weight on each edge into a constraint vector for the lower bound.
                ##  This will probably need to be a doubly nested for loop using the edges in G.edgeSet. We will
                ##  need a separate counter variable for the current coefficient of the vector we are working with. 
                ##

                ##
                ##  You will need three loops of this form: one for each of the lowerBoundVector, one for the upperBoundVector,
                ##  and one for the demandFunction. 
                ##

                demandFunction = np.zeros(G.vertexCount)

                for j in range(G.vertexCount):

                    print(f"State the demand at vertex {j}, \"{G.vertexNames[j]}\"")

                    demandFunction[j] = float(eval(input()))

                upperBound = np.zeros(G.edgeCount)

                lowerBound = np.zeros(G.edgeCount)

                index = 0

                for j in range(G.vertexCount):
                    
                    edges = G.edgeSet[j]

                    for k in range(len(edges)):

                        print(f"State the upper bound for flow cross the edge ({j}, {edges[k]}), connected vertex {j}, {G.vertexNames[j]} to vertex {edges[k]}, {G.vertexNames[edges[k]]}")

                        upperBound[index] = float(eval(input()))

                        print(f"State the lower bound for flow cross the edge ({j}, {edges[k]}), connected vertex {j}, {G.vertexNames[j]} to vertex {edges[k]}, {G.vertexNames[edges[k]]}")
                          
                        lowerBound[index] = float(eval(input()))

                        index += 1

                index = 0

                self.constraintVector = np.zeros(G.vertexCount + 2*G.edgeCount)

                for j in range(G.vertexCount): 

                    self.constraintVector[j] = demandFunction[j]

                for j in range(G.edgeCount):

                    self.constraintVector[j + G.vertexCount] = upperBound[j]

                    self.constraintVector[j + G.vertexCount + G.edgeCount] = lowerBound[j]

                self.constraintMatrix = np.zeros((G.vertexCount + 2*G.edgeCount, 3*G.edgeCount))

                for j in range(len(G.incidenceMatrix)):
                    for k in range(len(G.incidenceMatrix[0])):

                        self.constraintMatrix[j][k] = G.incidenceMatrix[j][k]

                for j in range(G.edgeCount):
                    ##
                    ##  upper bound
                    ##
                    self.constraintMatrix[j + G.vertexCount][j] = 1

                    self.constraintMatrix[j + G.vertexCount][j + G.edgeCount] = 1
                    ##
                    ##  lower bound
                    ##
                    self.constraintMatrix[j + G.vertexCount + G.edgeCount][j] = 1 

                    self.constraintMatrix[j + G.vertexCount + G.edgeCount][j + 2*G.edgeCount] = -1

                self.costFunction = np.zeros(len(self.constraintMatrix[0]))

                index = 0

                for j in range(G.vertexCount):
                    
                    edges = G.edgeSet[j]

                    for k in range(len(edges)):

                        self.costFunction[index] = G.adjacencyMatrix[j][edges[k]]

                        index += 1

                self.sucessfullAllocation = True


            ##
            ##  This is for if the weights accross the edges are only bounded from above. This will mostly
            ##  be the same as the if half of this if statement. 
            ##
            else:

                ##
                ##  Prompt for and scan in the weight on each edge into a constraint vector for the lower bound.
                ##  This will probably need to be a doubly nested for loop using the edges in G.edgeSet. We will
                ##  need a separate counter variable for the current coefficient of the vector we are working with. 
                ##

                ##
                ##  You will need three loops of this form: one for each of the lowerBoundVector, one for the upperBoundVector,
                ##  and one for the demandFunction. 
                ##

                demandFunction = np.zeros(G.vertexCount)

                for j in range(G.vertexCount):

                    print(f"State the demand at vertex {j}, \"{G.vertexNames[j]}\"")

                    demandFunction[j] = float(eval(input()))

                upperBound = np.zeros(G.edgeCount)

                index = 0

                for j in range(G.vertexCount):
                    
                    edges = G.edgeSet[j]

                    for k in range(len(edges)):

                        print(f"State the upper bound for flow cross the edge ({j}, {edges[k]}), connected vertex {j}, {G.vertexNames[j]} to vertex {edges[k]}, {G.vertexNames[edges[k]]}")

                        upperBound[index] = float(eval(input()))

                        index += 1

                index = 0

                self.constraintVector = np.zeros(G.vertexCount + G.edgeCount)

                for j in range(G.vertexCount): 

                    self.constraintVector[j] = demandFunction[j]

                for j in range(G.edgeCount):

                    self.constraintVector[j + G.vertexCount] = upperBound[j]



                self.constraintMatrix = np.zeros((G.vertexCount + G.edgeCount, 2*G.edgeCount))

                for j in range(len(G.incidenceMatrix)):
                    for k in range(len(G.incidenceMatrix[0])):

                        self.constraintMatrix[j][k] = G.incidenceMatrix[j][k]

                for j in range(G.edgeCount):
                    ##
                    ##  upper bound
                    ##
                    self.constraintMatrix[j + G.vertexCount][j] = 1

                    self.constraintMatrix[j + G.vertexCount][j + G.edgeCount] = 1

            index = 0

            self.costFunction = np.zeros(len(self.constraintMatrix[0]))

            for j in range(G.vertexCount):
                    
                edges = G.edgeSet[j]

                for k in range(len(edges)):

                    self.costFunction[index] = G.adjacencyMatrix[j][edges[k]]

                    index += 1

            self.sucessfullAllocation = True

        ##  
        ##  If the graph used to create this LP is not weighted then the program created will be a 
        ##  maximal flow problem, rather than a minimal cost flow problem. 
        ##  
        else:

            index = 0

            capacity = np.zeros(G.edgeCount)

            for j in range(G.vertexCount):

                edges = G.edgeSet[j]

                for k in range(len(edges)):

                    print(f"State the constraint on the flow accross the edge ({j}, {edges[k]}), between the edges connecting vertex {G.vertexNames[j]} to vertex {G.vertexNames[k]} (as a nonnegative real number):")

                    capacity[index] = float(eval(input()))

                    index += 1

            print(G.edgeSet)

            print("Which of these vertices is to be the sink node?")

            sink = int(input())

            print("Which of these vertices is to be the source")

            source = int(input())

            G.addEdge(sink, source, 1)

            ##  
            ##  You should now have enough to create the linear program. 
            ##  
            self.constraintMatrix = np.zeros((G.vertexCount + G.edgeCount, 2* G.edgeCount))

            for j in range(len(G.incidenceMatrix)):
                for k in range(len(G.incidenceMatrix[0])):

                    self.constraintMatrix[j][k] = G.incidenceMatrix[j][k]

            for j in range(G.edgeCount - 1):
                self.constraintMatrix[j + G.vertexCount][j] = self.constraintMatrix[j + G.vertexCount][j + G.edgeCount] = 1

            self.constraintVector = np.zeros(len(self.constraintMatrix))
            ##
            ##  We only have capacities defined for edges prior to the connecting between the sink and source, which
            ##  has already happened. If we do not subtract one here we will get an off by one error. 
            ##
            for j in range(G.edgeCount - 1):

                self.constraintVector[j + G.vertexCount] = capacity[j]

            self.costFunction = np.zeros(len(self.constraintMatrix[0]))

            self.costFunction[G.edgeCount -1] = -1

            self.sucessfullAllocation = True     

    ##  
    ##  Function Description: This function will create the LP for the currency exchange problem (or for any problem analagous).
    ##  
    def createCurrencyExchange(self, G :graphs.graph):

        index = 0

        capacity = np.zeros(G.edgeCount)

        for j in range(G.vertexCount):

            edges = G.edgeSet[j]

            for k in range(len(edges)):

                print(f"State the constraint on the flow accross the edge ({j}, {edges[k]}), between the edges connecting vertex {G.vertexNames[j]} to vertex {G.vertexNames[edges[k]]} (as a nonnegative real number):")

                capacity[index] = float(eval(input()))

                index += 1

        print(G.edgeSet)

        print("Which of these vertices is to be the sink node?")

        sink = int(input())

        print("Which of these vertices is to be the source")

        source = int(input())

        G.addEdge(sink, source, 1)

        ##  
        ##  You should now have enough to create the linear program. 
        ##  

        for j in range(G.vertexCount):
            for k in range(G.vertexCount):
                if not G.adjacencyMatrix[j][k] == 0:

                    for l in range(G.edgeCount):
                        if G.incidenceMatrix[j][l] > 0 and G.incidenceMatrix[k][l] < 0:
                            
                            G.incidenceMatrix[j][l] = G.adjacencyMatrix[j][k]

                            G.incidenceMatrix[k][l] = -1*G.adjacencyMatrix[j][k] 

        self.constraintMatrix = np.zeros((G.vertexCount + G.edgeCount, 2* G.edgeCount))

        for j in range(len(G.incidenceMatrix)):
            for k in range(len(G.incidenceMatrix[0])):

                self.constraintMatrix[j][k] = G.incidenceMatrix[j][k]

        for j in range(G.edgeCount - 1):
            self.constraintMatrix[j + G.vertexCount][j] = self.constraintMatrix[j + G.vertexCount][j + G.edgeCount] = 1


        self.constraintVector = np.zeros(len(self.constraintMatrix))
        ##
        ##  We only have capacities defined for edges prior to the connecting between the sink and source, which
        ##  has already happened. If we do not subtract one here we will get an off by one error. 
        ##
        for j in range(G.edgeCount - 1):

            self.constraintVector[j + G.vertexCount] = capacity[j]

        self.costFunction = np.zeros(len(self.constraintMatrix[0]))

        self.costFunction[G.edgeCount -1] = -1

        self.sucessfullAllocation = True     


##
##  Function description: This function uses pulp to solve a currency exchange problem
##
def solveCurrencyExchange(program : linearProgram):

    ##
    ##  Problem definition
    ##
    prob = pulp.LpProblem("Currency_Exchange", pulp.LpMinimize)

    ##
    ##  define decission variables:
    ##
    x = [pulp.LpVariable(f"x{i}", lowBound = 0, cat = "pulp.LpContinuous") for i in range(program.rowCount)]


    ##
    ##  Defining objective function
    ##
    prob += pulp.lpSum(program.costFunction[j] * x[j] for j in range(program.rowCount)), "Total Cost"

    ##
    ##  Constraints Ax = b
    ##
    for j in range(program.columnCount):
        prob += pulp.lpSum(program.constraintMatrix[j][k] * x[k] for k in range(len(x))) == program.constraintVector[j], f"Resource_{j}"

    ##
    ##  Solve:
    ##
    prob.solve()

    ##
    ##  Results:
    ##
    print(f"Status: {pulp.LpStatus[prob.status]}")
    print(f"Optimal x:", [x[j].varValue for j in range(program.columnCount)])
    print(f"Optimal Cost: ", pulp.value(prob.objective))





def createMatrixFromUserInput(): 
    ##  
    ##  This can, and should, have more type safety related exceptions. That will come 
    ##  later. Given that this function should ever need to be called I hope that this   
    ##  will not end up being a problem in practice. 
    ##

    columnCountFound = False
    rowCountFound = False
    ##  
    ##  This loop control structure, repeated twice, should be sufficient. 
    ##  
    while rowCountFound == False: 
        print("How many rows are there in this matrix? Enter your answer as a positive integer")
        rowCount = int(input())

        if rowCount <= 0:
            print("Only matrices with more than zero rows are defined: please enter a valid value.")
        else:
            rowCountFound = True

    while columnCountFound == False: 
        print("How many columns are there in this matrix? Enter your answer as a positive integer")
        columnCount = int(input())

        if columnCount <= 0:
            print("Only matrices with more than zero columns are defined: please enter a valid value.")
        else:
            columnCountFound = True

    matrix = np.zeros((rowCount, columnCount))

    for i in range(rowCount):
        for j in range(columnCount):
            print(f"Write the coefficient at row {i} and at column {j} of this matrix")
            matrix[i][j] = float(eval(input()))

    return matrix


def createMatrixFromUserInput(rowCount, columnCount):

    if rowCount <= 0:
        print("Only matrices with more than zero rows are defined: please enter a valid value.")

        return False
    
    elif columnCount <= 0:
        print("Only matrices with more than zero column are defined: please enter a valid value.")        

        return False
    
    else:
        matrix = np.zeros((rowCount, columnCount))

        for i in range(rowCount):
            for j in range(columnCount):
                print(f"Write the coefficient at row {i} and at column {j} of this matrix")
                matrix[i][j] = float(eval(input()))

    return matrix

##
##  Function Description: This function takes the length of a vector as an input from the 
##  user, prompts the user for the coefficients of the vector as arguments, then returns 
##  the vector. 
##
def createVector(length):

    if length <= 0:
        print("ERROR: Only vectors with positive integer dimension are defined.")

        return False

    returnVector = np.zeros(length)

    for j in range(length):
        print(f"Write coefficient {j} of the vector")
        returnVector[j] = float(eval(input()))

    return returnVector


