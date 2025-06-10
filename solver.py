import numpy as np
import pprint

class Network:
    """
    Represents a network of passive resistors.
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
    Returns a list of KCL equations from the resistive network.
    Input: network is a Network instance.
    Output: a list of equations, which are lists of elements that sum to 0
        elements in each equation list can be numerical values, or (variable, coefficient).
        current variables are tuples of vertices (smaller vtx, larger vtx)
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

def node_analysis(equations, pos_node, neg_node, voltage):
    """
    Returns an assignment of variables such that the input equations are satisfied.
    Input: equations is a list of equations, which are lists of elements that sum to 0
        elements in each equation list can be numerical, or (variable, coefficient)
        pos_node and neg_node are the nodes to be treated as the high and low ends
        of the input voltage.
    Output: a valid assignment to all variables to satisfy all equations.
    """
    # collect variables
    vars = set()
    for eqn in equations:
        vars = vars.union(set([clause[0] for clause in eqn]))
    vars_vector = list(vars)
    n = len(vars_vector)
    ix_lookup = {vars_vector[i] : i for i in range(n)}

    # now create and fill in the matrix constraints
    # we must also enforce both pos_node = voltage and neg_node = 0
    b = [0 for i in range(n)]
    b.append(0)
    b.append(voltage)

    print(len(b), 'b len')
    print(b)

    A = [[0 for i in range(n)] for j in range(len(equations)+2)]

    print(len(A), 'vert size')
    print(len(A[0]), 'hor size')

    for i, eqn in enumerate(equations):
        for clause in eqn:
            var_column = ix_lookup[clause[0]]
            A[i][var_column] = clause[1]
    pos_ix = ix_lookup[pos_node]
    neg_ix = ix_lookup[neg_node]
    A[-2][neg_ix] = 1
    A[-1][pos_ix] = 1

    pprint.pprint(A)


    x = np.linalg.solve(A, b)

    # use total output current to determine equivalent resistance









if __name__=="__main__":
    network = Network(3)
    network.add_resistor(0, 1, 10)
    network.add_resistor(1, 2, 20)
    network.add_resistor(2, 0, 30)

    kcl = kcl_equations(network)
    print(kcl)

    ohms = ohms_law_equations(network)
    print(ohms)

    equations = kcl + ohms
    print(node_analysis(equations, 0, 1, 10))
