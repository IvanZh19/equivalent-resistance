import numpy as np
import pprint

class Network:
    """
    Represents a network of passive resistors. Nodes are numbered starting at 0.
    """

    def __init__(self, num_nodes):
        """
        To init a network. None represents no edge, nodes are indexed by ints,
        Cell values are resistances.
        """
        self.matrix = [[None for i in range(num_nodes)] for j in range(num_nodes)]
        self.num_nodes = num_nodes

    def add_resistor(self, node1, node2, resistance):
        self.matrix[node1][node2] = resistance
        self.matrix[node2][node1] = resistance

    def del_resistor(self, node1, node2):
        self.matrix[node1][node2] = None
        self.matrix[node2][node1] = None

    def get_resistor(self, node1, node2):
        return self.matrix[node1][node2]

    def get_num_nodes(self):
        return self.num_nodes

def ohms_law_equations(network):
    """
    Returns a list of equations v-ir = 0 for each edge in the resistive network.
    Input: network is a Network instance.
    Output: a list of equations, which are lists of elements that sum to 0
        elements in each equation list can be numerical values, or (variable, coefficient).
    """
    n = network.get_num_nodes()
    equations = []
    for row_ix in range(n):
        for col_ix in range(row_ix+1):
            resistance = network.get_resistor(row_ix, col_ix)
            if resistance != None:
                equation = []
                equation.append(((col_ix, row_ix), -resistance))
                equation.append((col_ix, 1))
                equation.append((row_ix, -1))
                equations.append(equation)
    return equations

def kcl_equations(network):
    """
    Returns a list of all KCL equations from the resistive network.
    Input: network is a Network instance.
    Output: a list of equations, which are lists of elements that sum to 0
        elements in each equation list can be numerical values, or (variable, coefficient).
        current variables are tuples of vertices (smaller vtx, larger vtx)
        note that these are ordered in increasing node number.
    """
    n = network.get_num_nodes()
    equations = []
    for row_ix in range(n):
        equation = []
        for col_ix in range(n):
            if network.get_resistor(row_ix, col_ix) != None:
                # for sign change since we want currents to be uniquely (smaller_num, larger_num)
                if row_ix <= col_ix:
                    equation.append(((row_ix, col_ix), 1))
                else:
                    equation.append(((col_ix, row_ix), -1))
        equations.append(equation)
    return equations

def equivalent_resistance(network, pos_node, neg_node, voltage):
    """
    Performs a node analysis on the network to determine the equivalent resistance of the network.
    Input: network is a Network instance.
    Output: the equivalent resistance of the network.
    """
    # collect variables
    vars = set()
    kcl_eqns = kcl_equations(network)
    ohms_law_eqns = ohms_law_equations(network)
    equations = kcl_eqns + ohms_law_eqns
    for eqn in equations:
        vars = vars.union(set([clause[0] for clause in eqn]))
    vars_vector = list(vars)
    n = len(vars_vector)
    var_columns = {vars_vector[i] : i for i in range(n)}

    # now create and fill in the matrix constraints
    b = [0 for i in range(n)]
    A = [[0 for i in range(n)] for j in range(n)]
    for i, eqn in enumerate(equations):
        for clause in eqn:
            var_column = var_columns[clause[0]]
            A[i][var_column] = clause[1]

    # replace the positive and negative node KCL eqns with voltage constraints
    # note this relies on equations[i] = kcl_eqns[i] being KCL at node i
    neg_ix = var_columns[neg_node]
    ground_constraint = [0 for i in range(n)]
    ground_constraint[neg_ix] = 1
    A[neg_node] = ground_constraint

    pos_ix = var_columns[pos_node]
    positive_constraint = [0 for i in range(n)]
    positive_constraint[pos_ix] = 1
    A[pos_node] = positive_constraint
    b[pos_node] = voltage

    pprint.pprint(A)
    print(b)

    x = np.linalg.solve(A, b)
    is_exact = np.allclose(A @ x, b)
    print(is_exact)

    # use total output current to determine equivalent resistance
    # for var in vars_vector:
    #     if type(var) == tuple and pos_node in var:




    return x, vars_vector









if __name__=="__main__":
    network = Network(3)
    network.add_resistor(0, 1, 10)
    network.add_resistor(1, 2, 20)
    network.add_resistor(2, 0, 30)

    kcl = kcl_equations(network)
    print(kcl)

    ohms = ohms_law_equations(network)
    print(ohms)

    print(equivalent_resistance(network, 0, 1, 10))
