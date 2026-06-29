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
##  Preconditions: 
##
##  Postconditions: 
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
##  Untested:
##

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
##  Untested
##

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
##  I think I will start from scratch here. I want the indices of the basis to be 
##  kept in order. I also want the indices of the basis and the indices not in the  
##  basis to be tracked the entire way through the execution of this algorithm. 
##

##  
##  This function takes a linear program in standard form and performs the standard
##  simplex method on it. It returns the optimal solution to the linear program if   
##  it is finite, along with True, and returns the final direction vector if the 
##  solution is unbounded, along with False.
##
##  Precondtion: A standard form linear program, along with an initial basic feasible
##  solution, are used to call the function. The initial BFS, constraintMatrix, and 
##  constraintVector are all numpy arrays. No error checking is done within this  
##  function to ensure that the linear program is in standard form or that the feasible
##  region is non-empty.
##
##  Postcondition: The optimal BFS is returned along with True if it exists. The final
##  direction vector along with False is returned otherwise. 
##


##  
##  
##  This will be the standard version of the standard simplex method. I will also 
##  create a version of the standard simplex method that prints at each step in the 
##  agorithm. 
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
                ##  Unbounded solution: You must also return the direction vector 
                ##  and where you are directioning from, but that can come later.   
                ##
                ##  
                ##  It is now later: fuck, what do I do?
                ##  
                ##  Actually, this is probably fine as is. I'll just need to   
                ##  test it. 
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
##  Debugging version of the above code. This does the same thing as the standard simplex
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
                ##  Unbounded solution: You must also return the direction vector 
                ##  and where you are directioning from, but that can come later.   
                ##
                ##  
                ##  It is now later: fuck, what do I do?
                ##  
                ##  Actually, this is probably fine as is. I'll just need to   
                ##  test it. 
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
##  Precondition: 
##
##  Postcondition: 
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
##  Precondition: 
##
##  Postcondition: 
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

    basis = []

    for j in range(program.columnCount, len(costFunction)):
        costFunction[j] = 1

        basis.append(j)

    ##
    ##  you still need the basis. 
    ##
    (basicFeasibleSolution, costOrDireciton, status) = simplexTableauWithBlandsRule(basis, basicFeasibleSolution, tempMatrix, costFunction)

    return (basicFeasibleSolution, costOrDireciton, status)