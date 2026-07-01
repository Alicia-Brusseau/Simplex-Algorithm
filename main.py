import numpy as np
import matrixFunctions
import interfacingFunctions
import graphs
import pulp
import tests
import time
import matplotlib.pyplot as plt

'''
##  
##  The following are the three numpy arrays used to create the linear programs represented by the given graphs: 
##  

A1 = np.array([[1, 1, 1, 0, 0, 0, 0, 0, 0, 0], [-1, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, -1, 0, -1, 1, 0, 0, 0, 0, 0], [0, 0, -1, 0, -1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0, 0 ], [0, 0, 0, 1, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 1]], dtype=float)

b1 = np.array([8, 0, -2, -6, 0, 0, 0, 0, 0], dtype=float)

C1 = np.array([15, 4, 2, 8, 31, 0, 0, 0, 0, 0], dtype=float)


prob1 = pulp.LpProblem("graph_A", pulp.LpMinimize)

x1 = [pulp.LpVariable(f"x{i}", lowBound = 0, cat = "pulp.LpContinuous") for i in range(len(A1[0]))]

prob1 += pulp.lpSum(C1[j] * x1[j] for j in range(len(A1[0]))), "Total Cost"

for j in range(len(A1)):
    prob1 += pulp.lpSum(A1[j][k] * x1[k] for k in range(len(x1))) == b1[j], f"Constraint_{j}"


prob1.solve()

  
print(f"Status: {pulp.LpStatus[prob1.status]}")
print(f"Optimal x:", [x1[j].varValue for j in range(len(A1[0]))])
print(f"Optimal Cost: ", pulp.value(prob1.objective))

A2 = np.array([[1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [-1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, -1, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]], dtype=float)

b2 = np.array([8, 0, -2, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=float)

C2 = np.array([15, 4, 2, 8, 31, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], dtype=float)

prob2 = pulp.LpProblem("graph_B", pulp.LpMinimize)

x2 = [pulp.LpVariable(f"x{i}", lowBound = 0, cat = "pulp.LpContinuous") for i in range(len(A2[0]))]

prob2 += pulp.lpSum(C2[j] * x2[j] for j in range(len(A2[0]))), "Total Cost"

for j in range(len(A2)):
    prob2 += pulp.lpSum(A2[j][k] * x2[k] for k in range(len(x2))) == b2[j], f"Constraint_{j}"


prob2.solve()

  
print(f"Status: {pulp.LpStatus[prob2.status]}")
print(f"Optimal x:", [x2[j].varValue for j in range(len(A2[0]))])
print(f"Optimal Cost: ", pulp.value(prob2.objective))

A3 = np.array([[1, 0, 0, 1, 0, 0, 0, 0], [-1, 1, 0, 0, 0, 0, 0, 0], [0, -1, 1, 0, 0, 0, 0, 0], [0, 0, -1, -1, 0, 0, 0, 0], [1, 0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 0, 1, 0, 0], [0, 0, 1, 0, 0, 0, 1, 0], [0, 0, 0, 1, 0, 0, 0, 1]], dtype=float)

b3 = np.array([16, 0, 0, -16, 3, 3, 3, 9], dtype=float)

C3 = np.array([5, 5, 4, 16, 0, 0, 0, 0], dtype=float)

prob3 = pulp.LpProblem("graph_B", pulp.LpMinimize)

x3 = [pulp.LpVariable(f"x{i}", lowBound = 0, cat = "pulp.LpContinuous") for i in range(len(A3[0]))]

prob3 += pulp.lpSum(C3[j] * x3[j] for j in range(len(A3[0]))), "Total Cost"

for j in range(len(A3)):
    prob3 += pulp.lpSum(A3[j][k] * x3[k] for k in range(len(x3))) == b3[j], f"Constraint_{j}"


prob3.solve()

  
print(f"Status: {pulp.LpStatus[prob3.status]}")
print(f"Optimal x:", [x3[j].varValue for j in range(len(A3[0]))])
print(f"Optimal Cost: ", pulp.value(prob3.objective))


'''


##  
##  This was intended to be the prompts which take in the user's input in order to 
##  to create the specified graphs then which performs the required analysis with them. When this was
##  attempted the functions designed to perform the simplex method in tableau form failed, which did not occur at 
##  any point during testing. 
##  
'''
graph3 = graphs.graph()

graph3.fromUserInput()

program3 = interfacingFunctions.linearProgram()

program3.createProgramFromGraph(graph3)

start = time.time()

(basicFeasibleSolution, costOrDireciton, status, basis) = matrixFunctions.simplexPhaseOne(program3)

if status == True:

    matrixFunctions.simplexTableauWithBlandsRule(basis, basicFeasibleSolution, program3.constraintMatrix, program3.costFunction)

end = time.time()

print("Runtime: ", end - start)
'''

##
##  Used to produce the graphs shown in the report. 
##
'''
decisionVariablesCount = np.array([8, 10, 15])

constraintCount = np.array([8, 9, 14])

runtime = np.array([0.02, 0, 0])

plt.scatter(constraintCount, runtime)
plt.title("Scatter Plot: runtime vs number of constraints")
plt.xlabel("Number of constraints")
plt.ylabel("Runtime ")
plt.show()
'''