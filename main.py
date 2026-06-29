import numpy as np
import matrixFunctions
import interfacingFunctions
import graphs

program = interfacingFunctions.linearProgram()

program.createFromUserInputGeneralForm()

print("Constraints:")
print(program.constraintMatrix)
print("Constraint Vector:")
print(program.constraintVector)
print("Cost Function: ")
print(program.costFunction)

print("Test of simplex method phase one program: ")
print(matrixFunctions.simplexPhaseOne(program))


'''constraintMatrix = np.array([[1, -1, 1, 0, 0], [1, 1, 0, 1, 0], [0, 1, 0, 0, 1]], dtype=float)
cost = np.array([-1, 0, 0, 0, 0])
BFS = np.array([0, 0, 3, 8, 4])
basis = [2, 3, 4]

print(matrixFunctions.revisedSimplexWithBlandRuleDebug(constraintMatrix, cost, basis, BFS))'''

'''
constraintMatrix = np.array([[1, 2, 2, 1, 0, 0], [2, 1, 2, 0, 1, 0], [2, 2, 1, 0, 0, 1]])

costFunction = np.array([-10, -12, -12, 0, 0, 0])

BFS = np.array([0, 0, 0, 20, 20, 20])

basis = [3, 4, 5]

print(matrixFunctions.simplexTableauWithBlandsRuleDebug(basis, BFS, constraintMatrix, costFunction))
'''

'''augmentedBasisMatrix = np.array([[1, 2, 3, -4], [-2, 3, 1, 2], [4, -3, -2, 2]], dtype=float)

reducedAugmentedBasisMatrix = matrixFunctions.reduceAugmentedBasisMatrix(augmentedBasisMatrix, 2)

print("Augmented Basis Matrix:")
print(augmentedBasisMatrix)
print("Reduced Augmented Basis Matrix")
print(reducedAugmentedBasisMatrix)'''


'''tableau = matrixFunctions.createTableau(basis, BFS, constraintMatrix, costFunction)

print(tableau)

reducedTableau = matrixFunctions.reduceTableau(tableau, 3, 2)

print(reducedTableau)'''

""" G = graphs.graph()

G.fromUserInput()

print(G.edgeSet)

print(G.edgeCount)

print(G.incidenceMatrix)

print(G.adjacencyMatrix)





##  
##  We will now be testing the interfacing function default constructor for the linearProgram type. 
##  

optimimationProblem = interfacingFunctions.linearProgram()

optimimationProblem.createProgramFromGraph(G)

print(optimimationProblem.constraintMatrix)

print(optimimationProblem.constraintVector)

print(optimimationProblem.costFunction)
 """
##  
##  
##  These were test variables for the standard simplex method function. 
##  
##  

##constraintMatrix = np.array([[1, -1, 1, 0], [1, 1, 0, 1]])

##costFunction = np.array([-2, -1, 0, 0])

##basis = [2, 3]

##constraintVector = np.array([2, 6])

##initialBFS = np.array([0, 0, 2, 6])

##solution = matrixFunctions.standardSimplexMethod(initialBFS, basis, constraintMatrix, constraintVector, costFunction)

##print(solution)