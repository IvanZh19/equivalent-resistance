def ohms_law_equations(network):
    """
    Returns a list of equations v-ir = 0 for each edge in the resistive network.
    Input: network is a square matrix representing a resistive network graph.
        rows/cols 0 and 1 are for the source and sink.
        ints/floats are taken as resistances, and None is taken as no edge.
        note that network is by definition symmetric.
    Output: a list of equations, which are lists of elements that sum to 0
        elements in each equation list can be numerical values, or (variable, coefficient).
    """
    n = len(network)
    equations = []
    for row_ix in range(n):
        for col_ix in range(row_ix+1):
            if network[row_ix][col_ix] != None:
                equation = []
                equation.append(((col_ix, row_ix), -network[row_ix][col_ix]))
                equation.append((col_ix, 1))
                equation.append((row_ix, -1))
                equations.append(equation)
    return equations

def kcl_equations(network):
    """
    Returns a list of KCL equations from the resistive network.
    Input: network is a square matrix representing a resistive network graph.
        rows/cols 0 and 1 are for the source and sink.
        ints/floats are taken as resistances, and None is taken as no edge.
        note that network is by definition symmetric.
    Output: a list of equations, which are lists of elements that sum to 0
        elements in each equation list can be numerical values, or (variable, coefficient).
        current variables are tuples of vertices (smaller vtx, larger vtx)
    """
    n = len(network)
    equations = []
    for row_ix in range(n):
        equation = []
        for col_ix in range(n):
            if network[row_ix][col_ix] != None:
                if row_ix <= col_ix:
                    equation.append(((row_ix, col_ix), 1))
                else:
                    equation.append(((col_ix, row_ix), -1))
        equations.append(equation)
    return equations








if __name__=="__main__":
    network = [[None, 2, 3], [4, None, 6], [7, 8, None]]
    print(kcl_equations(network))

    print(ohms_law_equations(network))
