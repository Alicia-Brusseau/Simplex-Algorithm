import numpy as np
import interfacingFunctions

##  
##  Tested:
##  

##  
##  Function Description: This function takes the constraint equations of a 
##  program as well as a basis and returns the basis matrix corresponding to
##  a BFS and returns the corresponding basis matrix.
##  
##
def createBasisMatrix(basis, constraints):
    if len(basis) >= len(constraints[0]):
        return False
    
    ##
    ##  Remember that a solution is only a basic feasible solution if it 
    ##  has a number of linearly independent constraints equal to the 
    ##  dimension of the vector space we are working in. This requires 
    ##  that our basis matrix must be invertible and, as such, square. 
    ##
    ##  This is one of the many failsafes built into this code to 
    ##  ensure that the solutions used in this algorithm remain as 
    ##  basic feasible solutions. 
    ##  
    if len(basis) != len(constraints):
        return False

    for j in range(len(basis)):
        if basis[j] >= len(constraints[0]):
            return False
        
    basisMatrix = np.zeros((len(constraints), len(constraints)))

    for j in range(len(constraints)):
        for k in range(len(basis)):
            basisMatrix[j][k] = constraints[j][basis[k]]

    return basisMatrix

##  
##  Function Description: This function takes, as its arguments, a basis matrix 
##  (basisMatrix), a constraint matrix which the basisMatrix is from (constraints),  
##  and a column index in the range of consraints. The function returns to the 
##  scope which called it the inverse image of the column vector of the constraint
##  matrix under the transformation of basisMatrix. 
##
def findBasisRestrictedDirectionVector(basisMatrix, constraints, columnIndex):
    if columnIndex >= len(constraints[0]):
        return False
    
    elif len(basisMatrix) != len(basisMatrix[0]):
        return False
    ##  
    ##  Remember, basis matrices must be invertible. It is not sufficient
    ##  for them to be square.   
    ##
    elif np.linalg.det(basisMatrix) == 0:
        return False
    
    else:

        columnVector = np.zeros(len(constraints))

        for j in range(len(constraints)):
            columnVector[j] = constraints[j][columnIndex]

    return np.matmul(np.linalg.inv(basisMatrix), columnVector)

##
##  Function Description: This function takes a basis matrix, a matrix of constraints,
##  a cost function, and a basis (also corresponding to the basis matrix) and returns 
##  the reduced costs in the direction of each basis. 
##
##  Precondition: All arguments must be numpy arrays, except for basis, which is a 
##  list of indices. The formal parameter constraints must be the constraint matrix    
##  of a linear program, while constFunction must be the cost function of the same 
##  linear program and basisMatrix must be a basisMatrix derived from constraints 
##  with the basis indices in basis: all dimensions must match accordingly.
##
##  Postcondition: A numpy array of equal length to constFunction is returned to 
##  the scope which called this function, containing the reduced cost in each  
##  respective index direction.  
##
def computeReducedCosts(basisMatrix, constraints, costFunction, basis): 

    basisCost = np.zeros(len(basis))

    for j in range(len(basis)):
        basisCost[j] = costFunction[basis[j]]

    reducedCosts = np.zeros(len(constraints[0]))

    columnVector = np.zeros(len(constraints))

    for j in range(len(constraints[0])):
        for k in range(len(constraints)):
            columnVector[k] = constraints[k][j]

        reducedCosts[j] = costFunction[j] - np.dot(np.matmul(np.linalg.inv(basisMatrix), columnVector), basisCost)

    return reducedCosts

##
##  Function Description: This function takes as its argument a basis restricted
##  solution vector to a linear program along with a basis restricted direction 
##  vector and returns the minimal ratio between the coefficients of the former 
##  with the latter (at respective indices) if the direction vector is not 
##  all negative, otherwise it returns False. 
##
def findMinimalRatio(basisRestrictedSolution, basisRestrictedDirectionVector):

    if len(basisRestrictedDirectionVector) != len(basisRestrictedSolution):
        return False

    else:
        
        minimalRatio = np.inf
        minimalIndex = 0

        for j in range(len(basisRestrictedSolution)):
            if basisRestrictedDirectionVector[j] > 0:
                if basisRestrictedSolution[j] / basisRestrictedDirectionVector[j] < minimalRatio:
                    minimalIndex = j
                    minimalRatio = basisRestrictedSolution[j] / basisRestrictedDirectionVector[j]

        if minimalRatio == np.inf:
            return False
        else: 
            return (minimalIndex, minimalRatio)


##  
##  Function Description: This function performs the standard simplex method using Blands rule for 
##  cycle avoidance. 
##
def standardSimplexMethod(initialBFS, initialBasis, constraintMatrix, constraintVector, costFunction):

    ##  
    ##  All of this code is intended to be the first computation of the reduced costs
    ##  along edges of the polyhedron. The variables introduced by doing this will be 
    ##  used later if the initial bfs used were not optimal. 
    ##
    basis = initialBasis

    bfs = initialBFS

    basisMatrix = createBasisMatrix(basis, constraintMatrix)

    basisSolution = np.matmul(np.linalg.inv(basisMatrix), constraintVector)

    reducedCosts = computeReducedCosts(basisMatrix, constraintMatrix, costFunction, basis)

    minimalReducedCostIndex = 0

    minimalReducedCost = reducedCosts[minimalReducedCostIndex]

    for j in range(len(reducedCosts)):
        if reducedCosts[j] < minimalReducedCost:
            minimalReducedCost = reducedCosts[j]
            ##  
            ##  Recall: we are to have j enter the basis. 
            ##  
            minimalreducedCostIndex = j

    ##  
    ##  We will now account for two cases: The first being the case where the first 
    ##  basic feasible solution was optimal and the second for if the first basic  
    ##  feasible solution was not optimal. 
    ##
    if minimalReducedCost >= 0: 
        return (True, initialBFS)
    else:
        ##
        ##  This will be our loop condition. 
        ##
        solutionOptimal = False

        while solutionOptimal != True:
        
        
            basisSolution = np.matmul(np.linalg.inv(basisMatrix), constraintVector)


            unitDirectionVector = findBasisRestrictedDirectionVector(basisMatrix, constraintMatrix, minimalReducedCostIndex)

            if not (findMinimalRatio(basisSolution, unitDirectionVector) == False):

                (minimalRatioIndex, minimalRatio) = findMinimalRatio(basisSolution, unitDirectionVector)

            else:
                ##  
                ##  Unbounded solution:    
                ##
                return (False, bfs, basis, unitDirectionVector)



            ##  
            ##  Create new BFS here.  
            ##

            directionVector = np.zeros(len(initialBFS))

            for j in range(len(basis)):
                directionVector[basis[j]] = -1* unitDirectionVector[j]

            directionVector[minimalReducedCostIndex] = 1

            for j in range(len(directionVector)):
                directionVector[j] *= minimalRatio


            for j in range(len(bfs)):
                bfs[j] += directionVector[j]


            basis.remove(basis[minimalRatioIndex])

            basis.append(minimalReducedCostIndex)

            basis.sort()

            ##  
            ##  Check termination condition here.
            ##  

            basisMatrix = createBasisMatrix(basis, constraintMatrix)

            reducedCosts = computeReducedCosts(basisMatrix, constraintMatrix, costFunction, basis)

            solutionOptimal = True

            minimalReducedCost = reducedCosts[0]
            minimalRatioIndex = 0

            for j in range(len(reducedCosts)):
                if reducedCosts[j] < minimalReducedCost:
                    minimalReducedCost = reducedCosts[j]

                    minimalReducedCostIndex = j

            if minimalReducedCost < 0:
                solutionOptimal = False

        return (solutionOptimal, bfs)
   
    
##  
##  
##  Debugging version of the standard simplex method algorithm. This does the same thing as the standard simplex
##  method algorithm except for that it prints at each step along the way. 
##  
def standardSimplexMethodPrint(initialBFS, initialBasis, constraintMatrix, constraintVector, costFunction):

    ##  
    ##  All of this code is intended to be the first computation of the reduced costs
    ##  along edges of the polyhedron. The variables introduced by doing this will be 
    ##  used later if the initial bfs used were not optimal. 
    ##
    basis = initialBasis

    print("Our initial basis is:")
    print(basis)

    bfs = initialBFS

    print("Our initial BFS is: ")
    print(bfs)

    basisMatrix = createBasisMatrix(basis, constraintMatrix)

    print("The first basis matrix of ")
    print(constraintMatrix)
    print("corresponding to our basis is")
    print(basisMatrix)

    basisSolution = np.matmul(np.linalg.inv(basisMatrix), constraintVector)

    print("Our solution with respect to this basis is")
    print(basisSolution)

    reducedCosts = computeReducedCosts(basisMatrix, constraintMatrix, costFunction, basis)

    print("The first reduced costs in each index direction are")
    print(reducedCosts)

    minimalReducedCostIndex = 0

    minimalReducedCost = reducedCosts[minimalReducedCostIndex]

    for j in range(len(reducedCosts)):
        if reducedCosts[j] < minimalReducedCost:
            minimalReducedCost = reducedCosts[j]
            ##  
            ##  Recall: we are to have j enter the basis. 
            ##  
            minimalreducedCostIndex = j

    ##  
    ##  We will now account for two cases: The first being the case where the first 
    ##  basic feasible solution was optimal and the second for if the first basic  
    ##  feasible solution was not optimal. 
    ##
    if minimalReducedCost >= 0: 
        return (True, initialBFS)
    else:
        ##
        ##  This will be our loop condition. 
        ##
        solutionOptimal = False

        while solutionOptimal != True:
        
        
            basisSolution = np.matmul(np.linalg.inv(basisMatrix), constraintVector)

            print("The solution to the lp restricted to the basis is")
            print(basisSolution)

            unitDirectionVector = findBasisRestrictedDirectionVector(basisMatrix, constraintMatrix, minimalReducedCostIndex)

            print("The direction vector in the direction of the greatest reduced cost is")
            print(unitDirectionVector)


            if not (findMinimalRatio(basisSolution, unitDirectionVector) == False):

                (minimalRatioIndex, minimalRatio) = findMinimalRatio(basisSolution, unitDirectionVector)

                print("The minimal ratio of these coefficients at positive direction values is")
                print(minimalRatio)

            else:
                ##  
                ##  Unbounded solution:  
                ##     
                
                return (False, bfs, basis, unitDirectionVector)



            ##  
            ##  Create new BFS here. I'm worried I might cause a deep/shalow copy
            ##  error here by mistake.   
            ##

            directionVector = np.zeros(len(initialBFS))

            for j in range(len(basis)):
                directionVector[basis[j]] = -1* unitDirectionVector[j]

            directionVector[minimalReducedCostIndex] = 1

            for j in range(len(directionVector)):
                directionVector[j] *= minimalRatio


            print("The unrestricted direction vector is")
            print(directionVector)


            for j in range(len(bfs)):
                bfs[j] += directionVector[j]

            print("The new BFS is")
            print(bfs)


            basis.remove(basis[minimalRatioIndex])

            basis.append(minimalReducedCostIndex)

            basis.sort()

            print("The new basis is")
            print(basis)

            ##  
            ##  Check termination condition here.
            ##  

            basisMatrix = createBasisMatrix(basis, constraintMatrix)

            print("Our new basis matrix is")
            print(basisMatrix)

            reducedCosts = computeReducedCosts(basisMatrix, constraintMatrix, costFunction, basis)

            print("Our new reduced costs are")
            print(reducedCosts)

            solutionOptimal = True

            minimalReducedCost = reducedCosts[0]
            minimalRatioIndex = 0

            for j in range(len(reducedCosts)):
                if reducedCosts[j] < minimalReducedCost:
                    minimalReducedCost = reducedCosts[j]

                    minimalReducedCostIndex = j

            if minimalReducedCost < 0:
                solutionOptimal = False

        return (solutionOptimal, bfs)

##
##  Function Description: This function performs the simplex tableau version of the simplex
##  algorithm.
##
def simplexTableauWithBlandsRule(basis, basicFeasibleSolution, constraintMatrix, costFunction):
    
    BFS = basicFeasibleSolution.copy()

    dimension = len(BFS)

    nonBasicVariables = []

    for j in range(dimension):
        if not j in basis:

            nonBasicVariables.append(j)

    tableau = createTableau(basis, BFS, constraintMatrix, costFunction)

    reducedCost = 0

    ratios = [] 

    leaveIndex = -1

    enterIndex = -1

    while True:
        ##
        ##  Check if the solution is optimal 
        ##
        reducedCost = np.inf

        for j in range(len(basicFeasibleSolution)):
            if tableau[1][j + 2] < 0 and reducedCost >= 0:
                reducedCost = tableau[1][j + 2]  

                enterIndex = j
        ##
        ##  What is returned if this is the optimal solution. 
        ##
        if reducedCost >= 0:
            return (BFS, np.matmul(costFunction, BFS), True)
        
        else:

            ##
            ##  Create the direction vector
            ##
            basisReducedDirectionVector = np.zeros(len(basis))

            for j in range(len(basisReducedDirectionVector)):

                basisReducedDirectionVector[j] = tableau[2 + j][enterIndex+2]

            ##  
            ##  Check if the solution is unbounded
            ##  

            isUnbounded = True

            for j in range(len(basisReducedDirectionVector)):
                if basisReducedDirectionVector[j] > 0:
                    isUnbounded = False
            
            ##
            ##  What is returned if the solution is unbounded. 
            ##
            if isUnbounded == True:

                directionVector = np.zeros(len(BFS))

                for j in range(len(basis)):

                    directionVector[basis[j]] = basisReducedDirectionVector[j]

                return (BFS, directionVector, False)

            ##
            ##  Find the index to leave the basis
            ##
            minimumRatio = np.inf

            for j in range(len(basis)):

                ratios.append(BFS[basis[j]]/basisReducedDirectionVector[j])
                if BFS[basis[j]]/basisReducedDirectionVector[j] < minimumRatio and basisReducedDirectionVector[j] > 0:
                    minimumRatio = BFS[basis[j]]/basisReducedDirectionVector[j]

                    leaveIndex = j
            ##
            ##  Update the tableau
            ##
            tableau = reduceTableau(tableau, leaveIndex + 2, enterIndex + 2)
            
            ##
            ##  Update the basis
            ##
            for j in range(len(basis)):

                if tableau[j + 2][0] == basis[leaveIndex]:
                    tableau[j + 2][0] = enterIndex

            for j in range(len(basis)):
                if j == leaveIndex:
                    basis[j] = enterIndex
            for j in range(len(nonBasicVariables)):
                if nonBasicVariables[j] == enterIndex:
                    nonBasicVariables[j] = enterIndex

            ##
            ##  Update the basis
            ##
            BFS = np.zeros(len(BFS))

            for j in range(len(basis)):
                BFS[int(tableau[j + 2][0])] = tableau[j + 2][1]


##
##  Function Description: This function creates the tableau used in the simplex tableau version
##  of the simplex algorithm. 
##
##  Precondition: the variable basicFeasibleSolution is not reduced to the basis. basis is a list,
##  not a numpy array. 
##  Postcondition: 
##
def createTableau(basis, basicFeasibleSolution, constraintMatrix, costFunction):

    ##
    ##  This will be the basicFeasibleSolution reduced to its basis. 
    ##
    basisFeasibleSolution = np.zeros(len(basis))

    for j in range(len(basis)):
        basisFeasibleSolution[j] = basicFeasibleSolution[basis[j]]

    tableau = np.zeros((len(basis) + 2, len(basicFeasibleSolution) + 2))

    basisMatrix = np.zeros((len(constraintMatrix), len(basis)))

    sortedBasis = sorted(basis)

    for j in range(len(constraintMatrix)):
        for k in range(len(sortedBasis)):

            basisMatrix[j][k] = constraintMatrix[j][basis[k]]

    for j in range(len(basis)):

        tableau[j + 2][0] = basis[j]

        tableau[j + 2][1] = basisFeasibleSolution[j]

    basisCostFunction = np.zeros(len(basis))

    for j in range(len(basis)):
        basisCostFunction[j] = costFunction[basis[j]]

    rowVector = np.zeros(len(constraintMatrix))

    for j in range(len(constraintMatrix[0])):

        tableau[0][j + 2] = j

        for k in range(len(rowVector)):

            rowVector[k] = constraintMatrix[k][j]

        tableau[1][j + 2] = costFunction[j] - np.matmul(basisCostFunction, np.matmul(np.linalg.inv(basisMatrix), rowVector))

        directionVector = np.matmul(np.linalg.inv(basisMatrix), rowVector)

        for k in range(len(directionVector)):

            tableau[k + 2][j + 2] = directionVector[k]

    return tableau

def reduceTableau(tableau, row, column):

    reducedTableau = tableau.copy()

    constant = tableau[row][column]

    for j in range(1, len(tableau)):

        if not j == row:

            rowMultiple = -1*(tableau[j][column] / constant)

            for k in range(1, len(tableau[0])):
                reducedTableau[j][k] += rowMultiple * tableau[row][k]

    for j in range(1, len(tableau[0])):

        reducedTableau[row][j] *= (1/constant)

    return reducedTableau

##
##  Function Description: This function performs the simplex tableau version of the simplex
##  algorithm.
##
def simplexTableauWithBlandsRuleDebug(basis, basicFeasibleSolution, constraintMatrix, costFunction):
    
    BFS = basicFeasibleSolution.copy()

    print("The First BFS is")

    print(BFS)

    dimension = len(BFS)

    nonBasicVariables = []

    for j in range(dimension):
        if not j in basis:

            nonBasicVariables.append(j)

    tableau = createTableau(basis, BFS, constraintMatrix, costFunction)

    print("corresponding tableau:")

    print(tableau)

    reducedCost = 0

    ratios = [] 

    leaveIndex = -1

    enterIndex = -1

    while True:
        ##
        ##  Check if the solution is optimal 
        ##
        reducedCost = np.inf

        for j in range(len(basicFeasibleSolution)):
            if tableau[1][j + 2] < 0 and reducedCost >= 0:
                reducedCost = tableau[1][j + 2]  

                enterIndex = j
        ##
        ##  What is returned if this is the optimal solution. 
        ##
        if reducedCost >= 0:
            return (BFS, np.matmul(costFunction, BFS), True)
        
        else:

            print("Index entering basis:")
            print(enterIndex)

            ##
            ##  Create the direction vector
            ##
            basisReducedDirectionVector = np.zeros(len(basis))

            for j in range(len(basisReducedDirectionVector)):

                basisReducedDirectionVector[j] = tableau[2 + j][enterIndex+2]

            ##  
            ##  Check if the solution is unbounded
            ##  

            isUnbounded = True

            for j in range(len(basisReducedDirectionVector)):
                if basisReducedDirectionVector[j] > 0:
                    isUnbounded = False

            print("Direction Vector: ")
            
            
            ##
            ##  What is returned if the solution is unbounded. 
            ##
            if isUnbounded == True:

                directionVector = np.zeros(len(BFS))

                for j in range(len(basis)):

                    directionVector[basis[j]] = basisReducedDirectionVector[j]

                return (BFS, directionVector, False)
            
            print(basisReducedDirectionVector)

            ##
            ##  Find the index to leave the basis
            ##
            minimumRatio = np.inf

            for j in range(len(basis)):

                ratios.append(BFS[basis[j]]/basisReducedDirectionVector[j])
                if BFS[basis[j]]/basisReducedDirectionVector[j] < minimumRatio and basisReducedDirectionVector[j] > 0:
                    minimumRatio = BFS[basis[j]]/basisReducedDirectionVector[j]

                    leaveIndex = j
            ##
            ##  Update the tableau
            ##
            tableau = reduceTableau(tableau, leaveIndex + 2, enterIndex + 2)
            
            ##
            ##  Update the basis
            ##
            for j in range(len(basis)):

                if tableau[j + 2][0] == basis[leaveIndex]:
                    tableau[j + 2][0] = enterIndex

            for j in range(len(basis)):
                if j == leaveIndex:
                    basis[j] = enterIndex
            for j in range(len(nonBasicVariables)):
                if nonBasicVariables[j] == enterIndex:
                    nonBasicVariables[j] = enterIndex

            print("New Tableau:")

            print(tableau)

            ##
            ##  Update the basis
            ##
            BFS = np.zeros(len(BFS))

            for j in range(len(basis)):
                BFS[int(tableau[j + 2][0])] = tableau[j + 2][1]

            print("New Basis:")

            print(basis)




##
##  Function description: helper function to the revised simplex method which reduces the augmented
##  basis matrix around its pivot.
##
##  Precondition: The augmentedBasisMatrix must be a float. There is no way to ensure that the copy is    
##  a float from within the function scope because Python is a bad language. 
##
def reduceAugmentedBasisMatrix(augmentedBasisMatrix, exitRow):

    newAugmentedBasisMatrix = augmentedBasisMatrix.copy()

    for j in range(len(newAugmentedBasisMatrix)):

        if not j == exitRow:
            rowMultiple = float(newAugmentedBasisMatrix[j][-1] / newAugmentedBasisMatrix[exitRow][-1])

            for k in range(len(newAugmentedBasisMatrix[0])):
                newAugmentedBasisMatrix[j][k] -= (rowMultiple * newAugmentedBasisMatrix[exitRow][k] )

    constant = newAugmentedBasisMatrix[exitRow][-1]

    for k in range(len(augmentedBasisMatrix[0])):
        newAugmentedBasisMatrix[exitRow][k] *= float((1/constant))

    return newAugmentedBasisMatrix

##
##  Function Description: This function is intended to perform phase one of the simplex method: 
##
def simplexPhaseOne(program : interfacingFunctions.linearProgram):

    tempMatrix = np.zeros((program.rowCount, program.columnCount + program.rowCount), dtype = float)

    basicFeasibleSolution = np.zeros(program.columnCount + program.rowCount, dtype=float)

    costFunction = np.zeros(program.columnCount + program.rowCount, dtype=float)

    for j in range(program.rowCount):
        for k in range(program.columnCount):
            tempMatrix[j][k] = program.constraintMatrix[j][k]

        if program.constraintVector[j] >= 0:
            tempMatrix[j][j + program.columnCount] = 1

            basicFeasibleSolution[j + program.columnCount] = program.constraintVector[j]

        else:
            tempMatrix[j][j + program.columnCount] = -1

            basicFeasibleSolution[j + program.columnCount] = -1*program.constraintVector[j]

            for k in range(len(program.constraintMatrix[0])):
                tempMatrix[j][k] *= -1

    basis = []

    for j in range(program.columnCount, len(costFunction)):
        costFunction[j] = 1

        basis.append(j)

    ##
    ##  you still need the basis. 
    ##
    (basicFeasibleSolution, costOrDireciton, status, basis) = simplexTableauWithBlandsRuleBasis(basis, basicFeasibleSolution, tempMatrix, costFunction)

    return (basicFeasibleSolution, costOrDireciton, status, basis)

##
##  Function Description: This function performs the simplex tableau version of the simplex
##  algorithm.
##
##
##  Note: For whatever reason this code does not work on some graphs. I have only found it to not work on graphs.
##  specifically, it does not work on the graph from figure 2. 
##
def simplexTableauWithBlandsRuleBasis(basis, basicFeasibleSolution, constraintMatrix, costFunction):
    
    BFS = basicFeasibleSolution.copy()

    dimension = len(BFS)

    nonBasicVariables = []

    for j in range(dimension):
        if not j in basis:

            nonBasicVariables.append(j)

    tableau = createTableau(basis, BFS, constraintMatrix, costFunction)

    reducedCost = 0

    ratios = [] 

    leaveIndex = -1

    enterIndex = -1

    while True:
        ##
        ##  Check if the solution is optimal 
        ##
        reducedCost = np.inf

        for j in range(len(basicFeasibleSolution)):
            if tableau[1][j + 2] < 0 and reducedCost >= 0:
                reducedCost = tableau[1][j + 2]  

                enterIndex = j
        ##
        ##  What is returned if this is the optimal solution. 
        ##
        if reducedCost >= 0:
            return (BFS, np.matmul(costFunction, BFS), True, basis)
        
        else:

            ##
            ##  Create the direction vector
            ##
            basisReducedDirectionVector = np.zeros(len(basis))

            for j in range(len(basisReducedDirectionVector)):

                basisReducedDirectionVector[j] = tableau[2 + j][enterIndex+2]

            ##  
            ##  Check if the solution is unbounded
            ##  

            isUnbounded = True

            for j in range(len(basisReducedDirectionVector)):
                if basisReducedDirectionVector[j] > 0:
                    isUnbounded = False
            
            ##
            ##  What is returned if the solution is unbounded. 
            ##
            if isUnbounded == True:

                directionVector = np.zeros(len(BFS))

                for j in range(len(basis)):

                    directionVector[basis[j]] = basisReducedDirectionVector[j]

                return (BFS, directionVector, False, basis)

            ##
            ##  Find the index to leave the basis
            ##
            minimumRatio = np.inf

            for j in range(len(basis)):

                ratios.append(BFS[basis[j]]/basisReducedDirectionVector[j])
                if BFS[basis[j]]/basisReducedDirectionVector[j] < minimumRatio and basisReducedDirectionVector[j] > 0:
                    minimumRatio = BFS[basis[j]]/basisReducedDirectionVector[j]

                    leaveIndex = j
            ##
            ##  Update the tableau
            ##
            tableau = reduceTableau(tableau, leaveIndex + 2, enterIndex + 2)
            
            ##
            ##  Update the basis
            ##
            for j in range(len(basis)):

                if tableau[j + 2][0] == basis[leaveIndex]:
                    tableau[j + 2][0] = enterIndex

            for j in range(len(basis)):
                if j == leaveIndex:
                    basis[j] = enterIndex
            for j in range(len(nonBasicVariables)):
                if nonBasicVariables[j] == enterIndex:
                    nonBasicVariables[j] = enterIndex

            ##
            ##  Update the basis
            ##
            BFS = np.zeros(len(BFS))

            for j in range(len(basis)):
                BFS[int(tableau[j + 2][0])] = tableau[j + 2][1]

##
##  Function Description: This function is intended to perform phase one of the simplex method: 
##
def simplexPhaseOneDebug(program : interfacingFunctions.linearProgram):

    tempMatrix = np.zeros((program.rowCount, program.columnCount + program.rowCount), dtype = float)

    basicFeasibleSolution = np.zeros(program.columnCount + program.rowCount, dtype=float)

    costFunction = np.zeros(program.columnCount + program.rowCount, dtype=float)

    for j in range(program.rowCount):
        for k in range(program.columnCount):
            tempMatrix[j][k] = program.constraintMatrix[j][k]

        if program.constraintVector[j] >= 0:
            tempMatrix[j][j + program.columnCount] = 1

            basicFeasibleSolution[j + program.columnCount] = program.constraintVector[j]

        else:
            tempMatrix[j][j + program.columnCount] = -1

            basicFeasibleSolution[j + program.columnCount] = -1*program.constraintVector[j]

            for k in range(len(program.constraintMatrix[0])):
                tempMatrix[j][k] *= -1

    basis = []

    for j in range(program.columnCount, len(costFunction)):
        costFunction[j] = 1

        basis.append(j)

    ##
    ##  you still need the basis. 
    ##
    (basicFeasibleSolution, costOrDireciton, status, basis) = simplexTableauWithBlandsRuleBasisDebug(basis, basicFeasibleSolution, tempMatrix, costFunction)

    return (basicFeasibleSolution, costOrDireciton, status, basis)

##
##  Function Description: This function performs the simplex tableau version of the simplex
##  algorithm.
##
def simplexTableauWithBlandsRuleBasisDebug(basis, basicFeasibleSolution, constraintMatrix, costFunction):
    
    BFS = basicFeasibleSolution.copy()

    print("The First BFS is")

    print(BFS)

    dimension = len(BFS)

    nonBasicVariables = []

    for j in range(dimension):
        if not j in basis:

            nonBasicVariables.append(j)

    tableau = createTableau(basis, BFS, constraintMatrix, costFunction)

    print("corresponding tableau:")

    print(tableau)

    reducedCost = 0

    ratios = [] 

    leaveIndex = -1

    enterIndex = -1

    while True:
        ##
        ##  Check if the solution is optimal 
        ##
        reducedCost = np.inf

        for j in range(len(basicFeasibleSolution)):
            if tableau[1][j + 2] < 0 and reducedCost >= 0:
                reducedCost = tableau[1][j + 2]  

                enterIndex = j
        ##
        ##  What is returned if this is the optimal solution. 
        ##
        if reducedCost >= 0:
            return (BFS, np.matmul(costFunction, BFS), True, basis)
        
        else:

            print("Index entering basis:")
            print(enterIndex)

            ##
            ##  Create the direction vector
            ##
            basisReducedDirectionVector = np.zeros(len(basis))

            for j in range(len(basisReducedDirectionVector)):

                basisReducedDirectionVector[j] = tableau[2 + j][enterIndex+2]

            ##  
            ##  Check if the solution is unbounded
            ##  

            isUnbounded = True

            for j in range(len(basisReducedDirectionVector)):
                if basisReducedDirectionVector[j] > 0:
                    isUnbounded = False

            print("Direction Vector: ")
            
            
            ##
            ##  What is returned if the solution is unbounded. 
            ##
            if isUnbounded == True:

                directionVector = np.zeros(len(BFS))

                for j in range(len(basis)):

                    directionVector[basis[j]] = basisReducedDirectionVector[j]

                return (BFS, directionVector, False, basis)
            
            print(basisReducedDirectionVector)

            ##
            ##  Find the index to leave the basis
            ##
            minimumRatio = np.inf

            for j in range(len(basis)):

                ratios.append(BFS[basis[j]]/basisReducedDirectionVector[j])
                if BFS[basis[j]]/basisReducedDirectionVector[j] < minimumRatio and basisReducedDirectionVector[j] > 0:
                    minimumRatio = BFS[basis[j]]/basisReducedDirectionVector[j]

                    leaveIndex = j
            ##
            ##  Update the tableau
            ##
            tableau = reduceTableau(tableau, leaveIndex + 2, enterIndex + 2)
            
            ##
            ##  Update the basis
            ##
            for j in range(len(basis)):

                if tableau[j + 2][0] == basis[leaveIndex]:
                    tableau[j + 2][0] = enterIndex

            for j in range(len(basis)):
                if j == leaveIndex:
                    basis[j] = enterIndex
            for j in range(len(nonBasicVariables)):
                if nonBasicVariables[j] == enterIndex:
                    nonBasicVariables[j] = enterIndex

            print("New Tableau:")

            print(tableau)

            ##
            ##  Update the basis
            ##
            BFS = np.zeros(len(BFS))

            for j in range(len(basis)):
                BFS[int(tableau[j + 2][0])] = tableau[j + 2][1]

            print("New Basis:")

            print(basis)

##  
##  Function Description: This function performs the revised simplex algorithm. 
## 
##
##
##  This does not work because at the last pivot the interpreter decides to reinterpret some variables as 
##  integers rather than floats. I don't know what is causing this to happen and every method I have tried
##  to explicity cast a type has failed. I will, in the very near future, fix this by rewritting this 
##  entire program in a language that doesn't decide types of its own volition.  
##
def revisedSimplexWithBlandsRule(constraintMatrix, costFunction, basis, basicFeasibleSolution):
    basisMatrix = np.zeros((len(constraintMatrix), len(basis)))

    for j in range(len(constraintMatrix)):
        for k in range(len(basis)):

            basisMatrix[j][k] = constraintMatrix[j][basis[k]]

    basisFeasibleSolution = np.zeros(len(basis))

    basisCost = np.zeros(len(basis))

    for j in range(len(basis)):

        basisFeasibleSolution[j] = basicFeasibleSolution[basis[j]]

        basisCost[j] = costFunction[basis[j]]

    inverseBasisMatrix = np.linalg.inv(basisMatrix)

    priceVector = np.matmul(basisCost, inverseBasisMatrix)

    reducedCosts= np.zeros(len(costFunction))

    columnVector = np.zeros(len(constraintMatrix))

    for j in range(len(costFunction)):

        for k in range(len(constraintMatrix)):

            columnVector[k] = constraintMatrix[k][j]

        reducedCosts[j] = costFunction[j] - np.matmul(priceVector, columnVector)

    solutionOptimal = True

    enterIndex = -1

    for j in range(len(reducedCosts)):

        if reducedCosts[j] < 0:

            enterIndex = j

            solutionOptimal = False

    if solutionOptimal:

        return (basicFeasibleSolution, np.matmul(costFunction, basicFeasibleSolution), solutionOptimal)
    ##
    ##  All of the rest of the algorithm will be done in here. 
    ##
    else:

        newBasicFeasibleSolution = basicFeasibleSolution.copy()

        while not solutionOptimal:

            for j in range(len(constraintMatrix)):

                columnVector[j] = constraintMatrix[j][enterIndex]

            basisDirectionVector = np.matmul(inverseBasisMatrix, columnVector)

            for j in range(len(basis)):
                basisDirectionVector[j] *= -1

            solutionUnbounded = True

            for j in range(len(basisDirectionVector)):

                if basisDirectionVector[j] < 0:

                    solutionUnbounded = False

            if solutionUnbounded:
                ##
                ##  This is where I will return unbounded solutions. 
                ##
                ##  Must: 
                ##      > Find basis unrestricted direction vector
                ##      > Return (BFS, basis unrestricted direction vector, false)

                directionVector = np.zeros(len(newBasicFeasibleSolution))

                for j in range(len(basisDirectionVector)):

                    directionVector[basis[j]] = basisDirectionVector[j]

                return (newBasicFeasibleSolution, directionVector, False)

                

            else:
                ##  
                ##  now to pivot
                ##  
                exitIndex = -1

                minimalRatio = np.inf

                for j in range(len(basis)):

                    if basisDirectionVector[j] < 0:

                        if -1*(newBasicFeasibleSolution[basis[j]] / basisDirectionVector[j]) < minimalRatio:

                            minimalRatio = -1*(newBasicFeasibleSolution[basis[j]] / basisDirectionVector[j])

                            exitIndex = j

                directionVector = np.zeros(len(newBasicFeasibleSolution))

                for j in range(len(basisDirectionVector)):

                    directionVector[basis[j]] = basisDirectionVector[j]

                directionVector[exitIndex] = 1

                for j in range(len(newBasicFeasibleSolution)):

                    newBasicFeasibleSolution[j] += minimalRatio * directionVector[j]

                ##
                ##  Create new basis matrix. 
                ##
                augmentBasisMatrix = np.zeros((len(constraintMatrix), len(basis) + 1), dtype = float)

                for j in range(len(inverseBasisMatrix)):

                    for k in range(len(inverseBasisMatrix[0])):

                        augmentBasisMatrix[j][k] = inverseBasisMatrix[j][k]

                    augmentBasisMatrix[j][-1] = -1*basisDirectionVector[j]

                augmentBasisMatrix = reduceAugmentedBasisMatrix(augmentBasisMatrix, exitIndex)                    

                for j in range(len(inverseBasisMatrix)):

                    for k in range(len(inverseBasisMatrix[0])):

                        inverseBasisMatrix[j][k] = augmentBasisMatrix[j][k]

                ##
                ##  It is now time to update the basis. 
                ##
                basis[exitIndex] = enterIndex

                basis.sort()

                for j in range(len(basisCost)):

                    basisCost[j] = costFunction[basis[j]]


                priceVector = np.matmul(basisCost, inverseBasisMatrix)

                for j in range(len(reducedCosts)):

                    for k in range(len(constraintMatrix)):

                        columnVector[k] = constraintMatrix[k][j]

                    reducedCosts[j] = costFunction[j] - np.matmul(priceVector, columnVector)

                ##
                ##  Check termination condition: 
                ##
                solutionOptimal = True

                for j in range(len(reducedCosts)):

                    if reducedCosts[j] < 0:

                        solutionOptimal = False

                        enterIndex = j

                if solutionOptimal:

                    return (newBasicFeasibleSolution, np.matmul(costFunction, newBasicFeasibleSolution), solutionOptimal)

##  
##  Function Description: This function performs the revised simplex algorithm. 
##  
def revisedSimplexWithBlandsRuleDebug(constraintMatrix, costFunction, basis, basicFeasibleSolution):

    basisMatrix = np.zeros((len(constraintMatrix), len(basis)))

    for j in range(len(constraintMatrix)):
        for k in range(len(basis)):

            basisMatrix[j][k] = constraintMatrix[j][basis[k]]

    print("Basis matrix")
    print(basisMatrix)

    basisFeasibleSolution = np.zeros(len(basis))

    basisCost = np.zeros(len(basis))

    for j in range(len(basis)):

        basisFeasibleSolution[j] = basicFeasibleSolution[basis[j]]

        basisCost[j] = costFunction[basis[j]]

    print("Basis restricted basic feasible solution")
    print(basicFeasibleSolution)

    inverseBasisMatrix = np.linalg.inv(basisMatrix)

    print("Inverse of basis matrix")
    print(inverseBasisMatrix)

    priceVector = np.matmul(basisCost, inverseBasisMatrix)

    print("Price vector")
    print(priceVector)

    reducedCosts= np.zeros(len(costFunction))

    columnVector = np.zeros(len(constraintMatrix))

    for j in range(len(costFunction)):

        for k in range(len(constraintMatrix)):

            columnVector[k] = constraintMatrix[k][j]

        reducedCosts[j] = costFunction[j] - np.matmul(priceVector, columnVector)

    print("Reduced costs")
    print(reducedCosts)

    solutionOptimal = True

    enterIndex = -1

    for j in range(len(reducedCosts)):

        if reducedCosts[j] < 0:

            enterIndex = j

            solutionOptimal = False

    if solutionOptimal:

        return (basicFeasibleSolution, np.matmul(costFunction, basicFeasibleSolution), solutionOptimal)
    ##
    ##  All of the rest of the algorithm will be done in here. 
    ##
    else:

        newBasicFeasibleSolution = basicFeasibleSolution.copy()

        while not solutionOptimal:

            for j in range(len(constraintMatrix)):

                columnVector[j] = constraintMatrix[j][enterIndex]

            basisDirectionVector = np.matmul(inverseBasisMatrix, columnVector)

            print("Basis restricted direction vector")
            print(basisDirectionVector)

            for j in range(len(basis)):
                basisDirectionVector[j] *= -1

            solutionUnbounded = True

            for j in range(len(basisDirectionVector)):

                if basisDirectionVector[j] < 0:

                    solutionUnbounded = False

            if solutionUnbounded:
                ##
                ##  This is where I will return unbounded solutions. 
                ##
                ##  Must: 
                ##      > Find basis unrestricted direction vector
                ##      > Return (BFS, basis unrestricted direction vector, false)

                directionVector = np.zeros(len(newBasicFeasibleSolution))

                for j in range(len(basisDirectionVector)):

                    directionVector[basis[j]] = basisDirectionVector[j]

                return (newBasicFeasibleSolution, directionVector, False)

                

            else:
                ##  
                ##  now to pivot
                ##  
                exitIndex = -1

                minimalRatio = np.inf

                for j in range(len(basis)):

                    if basisDirectionVector[j] < 0:

                        if -1*(newBasicFeasibleSolution[basis[j]] / basisDirectionVector[j]) < minimalRatio:

                            minimalRatio = -1*(newBasicFeasibleSolution[basis[j]] / basisDirectionVector[j])

                            exitIndex = j

                print("Minimum ratio")
                print(minimalRatio)

                directionVector = np.zeros(len(newBasicFeasibleSolution))

                for j in range(len(basisDirectionVector)):

                    directionVector[basis[j]] = basisDirectionVector[j]

                directionVector[exitIndex] = 1

                print("basis unrestricted direction vector")
                print(directionVector)

                for j in range(len(newBasicFeasibleSolution)):

                    newBasicFeasibleSolution[j] += minimalRatio * directionVector[j]
                
                print("New basic feasible solution")
                print(newBasicFeasibleSolution)

                ##
                ##  Time to create new basis matrix. 
                ##
                augmentBasisMatrix = np.zeros((len(constraintMatrix), len(basis) + 1), dtype = float)

                for j in range(len(inverseBasisMatrix)):

                    for k in range(len(inverseBasisMatrix[0])):

                        augmentBasisMatrix[j][k] = inverseBasisMatrix[j][k]

                    augmentBasisMatrix[j][-1] = -1*basisDirectionVector[j]

                print("Augmented inverse basis matrix")
                print(augmentBasisMatrix)

                augmentBasisMatrix = reduceAugmentedBasisMatrix(augmentBasisMatrix, exitIndex)                    

                print("Reduced augmented basis matrix")
                print(augmentBasisMatrix)

                for j in range(len(inverseBasisMatrix)):

                    for k in range(len(inverseBasisMatrix[0])):

                        inverseBasisMatrix[j][k] = augmentBasisMatrix[j][k]

                print("new inverse basis matrix")
                print(inverseBasisMatrix)

                ##
                ##  It is now time to update the basis. 
                ##
                basis[exitIndex] = enterIndex

                basis.sort()

                print("new basis")
                print(basis)

                for j in range(len(basisCost)):

                    basisCost[j] = costFunction[basis[j]]

                print("new basis restricted cost")
                print(basisCost)

                priceVector = np.matmul(basisCost, inverseBasisMatrix)

                print("new price vector")
                print(priceVector)

                for j in range(len(reducedCosts)):

                    for k in range(len(constraintMatrix)):

                        columnVector[k] = constraintMatrix[k][j]

                    reducedCosts[j] = costFunction[j] - np.matmul(priceVector, columnVector)

                print("new reduced costs")
                print(reducedCosts)

                ##
                ##  Check termination condition: 
                ##
                solutionOptimal = True

                for j in range(len(reducedCosts)):

                    if reducedCosts[j] < 0:

                        solutionOptimal = False

                        enterIndex = j

                if solutionOptimal:

                    return (newBasicFeasibleSolution, np.matmul(costFunction, newBasicFeasibleSolution), solutionOptimal)

                
