import numpy as np
import interfacingFunctions
import matrixFunctions
import graphs
import pulp

def testForMatrixInstantiationInPulp():
    ##
    ##  These will be the arrays used to create the LP. 
    ##
    A = [[1, -1, 1, 0, 0], [1, 1, 0, 1, 0], [0, 1, 0, 0, 1]]

    ##
    ##  Objective array
    ##
    C = [-1, 0, 0, 0, 0]

    ##
    ##  Problem definition
    ##
    prob = pulp.LpProblem("Prodution_Planning", pulp.LpMinimize)

    ##
    ##  define decission variables:
    ##
    x = [pulp.LpVariable(f"x{i}", lowBound = 0, cat = "pulp.LpContinuous") for i in range(len(A[0]))]


    ##
    ##  Defining objective function
    ##
    prob += pulp.lpSum(C[j] * x[j] for j in range(len(A[0]))), "Total Cost"

    b = [3, 8, 4]

    ##
    ##  Constraints Ax = b
    ##
    for j in range(len(A)):
        prob += pulp.lpSum(A[j][k] * x[k] for k in range(len(x))) == b[j], f"Resource_{j}"

    ##
    ##  Solve:
    ##
    prob.solve()

    ##
    ##  Results:
    ##
    print(f"Status: {pulp.LpStatus[prob.status]}")
    print(f"Optimal x:", [x[j].varValue for j in range(len(A[0]))])
    print(f"Optimal Cost: ", pulp.value(prob.objective))

def testForCurrencyExchangeProblem():

    G = graphs.graph()

    G.fromUserInput()

    print(G.edgeSet)

    print(G.edgeCount)

    print(G.incidenceMatrix)

    print(G.adjacencyMatrix)

    optimimationProblem = interfacingFunctions.linearProgram()

    optimimationProblem.createCurrencyExchange(G)

    print(optimimationProblem.constraintMatrix)

    print(optimimationProblem.constraintVector)

    print(optimimationProblem.costFunction)

    interfacingFunctions.solveCurrencyExchange(optimimationProblem)

def pulpTestI():

    ourProblem = pulp.LpProblem() 

    ##  
    ##  Problem decision variables:
    ##  
    x1 = pulp.LpVariable("Variable1", lowBound=0, cat=pulp.LpContinuous)

    x2 = pulp.LpVariable("Variable2", lowBound=0, cat = pulp.LpContinuous)

    x3 = pulp.LpVariable("Variable3", lowBound=0, cat=pulp.LpContinuous)

    x4 = pulp.LpVariable("Variable4", lowBound=0, cat=pulp.LpContinuous)

    x5 = pulp.LpVariable("Variable5", lowBound=0, cat=pulp.LpContinuous)

    x6 = pulp.LpVariable("Variable6", lowBound=0, cat=pulp.LpContinuous)

    x7 = pulp.LpVariable("Variable7", lowBound=0, cat=pulp.LpContinuous)

    x8 = pulp.LpVariable("Variable8", lowBound=0, cat=pulp.LpContinuous)

    ##
    ##  Objective:
    ##
    ourProblem += 1 * x1 + 2*x2 + 3*x3 + 4*x4 + 5*x5 + 6*x6 + 7*x7 + 8*x8

    ##  
    ##  Row One Constraint:
    ##  
    ourProblem += x1 + 4*x2 + 5*x3 + 6*x4 == 8

    ##
    ##  Row Two Constraint:
    ##
    ourProblem += 6*x1 + 3*x2 + 5*x3 + x4 == 7

    ##  
    ##  Row Three Constraint:
    ##  
    ourProblem += x1 + x3 + x5 == 6

    ##
    ##  Row Four Constraint:
    ##
    ourProblem += x2 + x4 + x6 == 5

    ##
    ##  Row Five Constraint: 
    ##
    ourProblem += x1 + x2 + x3 + x4 - x7 == 4

    ##
    ##  Row Six Constraint: 
    ##
    ourProblem += x1 + x3 -x8 == 3

    ## write LP
    ourProblem.writeLP("sampleLinearProgram.lp")
    ourProblem.solve()

    print(pulp.LpStatus[ourProblem.status])

    ##
    ##  get variable values
    ##
    varvals = {}

    for v in ourProblem.variables():

        varvals[v.name] = v.varValue
        print(v.name, "=", v.varValue)

def simplexMethodPhaseOneTest():

    program = interfacingFunctions.linearProgram()

    program.createFromUserInputGeneralForm()

    print("Constraints:")
    print(program.constraintMatrix)
    print("Constraint Vector:")
    print(program.constraintVector)
    print("Cost Function: ")
    print(program.costFunction)

    print("Test of simplex method phase one program: ")

    (BFS, costOrDirection, status, basis) = matrixFunctions.simplexPhaseOne(program)

    print(f"({BFS}, {costOrDirection}, {status}, {basis})")

def revisedSimplexTest():
    constraintMatrix = np.array([[1, -1, 1, 0, 0], [1, 1, 0, 1, 0], [0, 1, 0, 0, 1]], dtype=float)
    cost = np.array([-1, 0, 0, 0, 0])
    BFS = np.array([0, 0, 3, 8, 4])
    basis = [2, 3, 4]

    print(matrixFunctions.revisedSimplexWithBlandsRuleDebug(constraintMatrix, cost, basis, BFS))

def simplexTableauTest():
    constraintMatrix = np.array([[1, 2, 2, 1, 0, 0], [2, 1, 2, 0, 1, 0], [2, 2, 1, 0, 0, 1]])

    costFunction = np.array([-10, -12, -12, 0, 0, 0])

    BFS = np.array([0, 0, 0, 20, 20, 20])

    basis = [3, 4, 5]

    print(matrixFunctions.simplexTableauWithBlandsRuleDebug(basis, BFS, constraintMatrix, costFunction))

def revisedSimplexRowReductionTest():
    
    augmentedBasisMatrix = np.array([[1, 2, 3, -4], [-2, 3, 1, 2], [4, -3, -2, 2]], dtype=float)

    reducedAugmentedBasisMatrix = matrixFunctions.reduceAugmentedBasisMatrix(augmentedBasisMatrix, 2)

    print("Augmented Basis Matrix:")
    print(augmentedBasisMatrix)
    print("Reduced Augmented Basis Matrix")
    print(reducedAugmentedBasisMatrix)

def graphUserInputInstantiationTest():
    G = graphs.graph()

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

def standardSimplexMethodTest():
    
    constraintMatrix = np.array([[1, -1, 1, 0], [1, 1, 0, 1]])

    costFunction = np.array([-2, -1, 0, 0])

    basis = [2, 3]

    constraintVector = np.array([2, 6])

    initialBFS = np.array([0, 0, 2, 6])

    solution = matrixFunctions.standardSimplexMethod(initialBFS, basis, constraintMatrix, constraintVector, costFunction)

    print(solution)