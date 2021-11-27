import numpy as np

def is_numeric(x):
    "Checks if the input is a number."
    return isinstance(x, (int, float))

def is_negative(x):
    "Checks if the input is negative."
    return x < 0


class LinearCombination:
    "A linear combination of monomial terms."

    def __init__(self, *coefs):
        self.coefs = list(coefs)

    def append(self, coef):
        "Adds a new monomial with the given coefficient."
        assert is_numeric(coef)
        self.coefs.append(coef)

    @property
    def vector_form(self):
        return np.array(self.coefs, dtype=np.float64)

    def __len__(self):
        return len(self.coefs)

    def __repr__(self):
        buffer = ''
        for index, coef in enumerate(self.coefs):
            if coef == 0:
                continue

            if index > 0:
                if coef >= 0:
                    buffer += ' + '
                else:
                    buffer += ' - '
                    coef = -coef

            buffer += f'{coef}·'

            variable = f'x_{index + 1}'
            buffer += variable

        return buffer


class LinearConstraint:
    def __init__(self, *parts):
        assert len(parts) >= 3, 'Need at least one coefficient, an operator and the free term'

        coefs = parts[:-2]
        assert all(map(is_numeric, coefs))
        operator = parts[-2]
        assert operator in ('<=', '>=', '='), f'Invalid operator: `{operator}`'
        free_term = parts[-1]

        self.combination = LinearCombination(*coefs)
        self.operator = operator
        self.free_term = free_term

    @property
    def is_equality(self):
        return self.operator == '='

    @property
    def is_less_than(self):
        return self.operator == '<='

    @property
    def is_greater_than(self):
        return self.operator == '>='

    def make_equality(self):
        "Turns this inequality constraint into an equality constraint by adding a slack variable."
        assert not self.is_equality

        if self.is_less_than:
            self.add_variable(1)
        else:
            self.add_variable(-1)

        self.operator = '='

    def add_variable(self, coef):
        self.combination.append(coef)

    def add_free_variable(self):
        "Adds a new variable which doesn't affect this contraint."
        self.combination.append(0)

    @property
    def vector_form(self):
        return np.hstack([self.combination.vector_form, self.free_term])

    def __len__(self):
        return len(self.combination)

    def __repr__(self):
        return f'{self.combination} {self.operator} {self.free_term}'


class LinearProgram:
    def __init__(self, constraints, objective, maximize=True):
        self.constraints = constraints
        self.objective = objective
        self.maximize = maximize
        self.num_free_variables = 0

    @staticmethod
    def from_matrix(A, maximize):
        operator = '<=' if maximize else '>='

        constraints = []
        for row in A[:-1]:
            constraints.append(LinearConstraint(*row[:-1], operator, row[-1]))
        objective = LinearCombination(*A[-1, :-1])

        return LinearProgram(constraints, objective, maximize)

    @property
    def is_canonical(self):
        "Checks if this linear program is in canonical form."
        if self.maximize:
            return all(c.is_less_than for c in self.constraints)
        else:
            return all(c.is_greater_than for c in self.constraints)

    @property
    def is_standard(self):
        "Checks if this linear program is in standard form."
        return all(c.is_equality for c in self.constraints)

    @property
    def matrix_form(self):
        coefs = []
        for constraint in self.constraints:
            coefs.append(constraint.vector_form)

        if self.is_standard:
            coefs.append(self.objective.vector_form)
        else:
            coefs.append([*self.objective.coefs, 1])
        #print([*self.objective.coefs, 1])
        return np.vstack(coefs)

    @property
    def num_decision_variables(self):
        return len(self.constraints[0]) - self.num_free_variables

    def standardize(self):
        num_free_variables = 0
        for constraint in self.constraints:
            if constraint.is_equality:
                continue

            constraint.make_equality()

            for c in self.constraints:
                if c == constraint:
                    continue

                c.add_free_variable()

            num_free_variables += 1

        coefs = [-c for c in self.objective.coefs]
        zeros = [0] * num_free_variables
        self.objective = LinearConstraint(*coefs, *zeros, 1, '=', 0)
        for c in self.constraints:
            c.add_free_variable()

        self.num_free_variables = num_free_variables + 1

    @property
    def dual(self):
        return self.from_matrix(self.matrix_form.T, not self.maximize)

    def solve(self):
        assert self.is_canonical, 'Linear program must be in canonical form!'

        if self.maximize:
            self._solve_primal()
        else:
            self._solve_dual()

    def _solve_primal(self):
        if not self.is_standard:
            self.standardize()
        A = self.matrix_form
        for column, value in enumerate(A[-1]):
            print(column, value)

        while True:
            pivot_column = -1
            for column, value in enumerate(A[-1]):
                if value < 0:
                    pivot_column = column
                    break

            # Optimum reached
            if pivot_column == -1:
                break

            pivot_row = -1
            min_ratio = np.inf
            for row, value in enumerate(A[:-1, -1]):
                elem = A[row, pivot_column]
                if elem <= 0:
                    continue

                ratio = value / elem
                if ratio < min_ratio:
                    pivot_row = row
                    min_ratio = ratio

            # Optimum is unbounded
            if pivot_row == -1:
                print('optim = inf')
                return

            self._pivot(A, pivot_row, pivot_column)

        print(f'optim = {A[-1][-1]}')
        print(len(self.constraints[0]))
        print(self.num_free_variables)
        for index in range(self.num_decision_variables):
            print(f'x_{index + 1} = {A[-1, index]}')

    def _solve_dual(self):
        #print(self)
        orig = self.dual
        #print(orig)
        orig.standardize()
        #print(orig)
        A = orig.matrix_form
        A = np.delete(A, -2, axis=1)
        A[-1] *= -1
        print(A)

        while True:
            pivot_row = -1
            for row, value in enumerate(A[:, -1]):
                if value < 0:
                    pivot_row = row
                    break

            # optimum reached
            if pivot_row == -1:
                break

            pivot_column = -1
            min_ratio = np.inf
            for column, value in enumerate(A[-1, :-1]):
                elem = A[pivot_row][column]
                if elem >= 0:
                    continue

                ratio = value / elem
                if ratio < min_ratio:
                    pivot_column = column
                    min_ratio = ratio

            # Optimum is unbounded
            if pivot_column == -1:
                print('optim = -inf')
                return

            self._pivot(A, pivot_row, pivot_column)

        print(f'optim = {A[-1][-1]}')
        for index in range(self.num_decision_variables):
            print(f'y_{index + 1} = {-A[-1, self.num_decision_variables + index]}')

    @staticmethod
    def _pivot(matrix, pivot_row, pivot_column):
        "Performs a step of Gaussian reduction for the given pivot."
        pivot = matrix[pivot_row, pivot_column]
        matrix[pivot_row] /= pivot

        for row in range(len(matrix)):
            if row == pivot_row:
                continue

            matrix[row] -= matrix[row, pivot_column] * matrix[pivot_row]

    def __repr__(self):
        buffer = ''
        if self.maximize:
            buffer += 'max'
        else:
            buffer += 'min'

        buffer += f' {self.objective}\n'

        for constraint in self.constraints:
            buffer += f'{constraint}\n'
        return buffer

prog = LinearProgram([
    LinearConstraint(2, 1, '<=', 3),
    LinearConstraint(1, 2, '<=', 9),
], LinearCombination(8, 8))

#print(prog)
#prog.solve()

prog = LinearProgram([
    LinearConstraint(-1, -1, 2, '<=', -3),
    LinearConstraint(-4, -2, 1, '<=', -4),
    LinearConstraint(1, 1, -4, '<=', 2),
], LinearCombination(-4, -2, -1))

pd = prog.dual

#print(pd)
pd.solve()

#print(LinearCombination(-4, -2, -1))
#print(LinearCombination(8, 8))

#print(LinearConstraint(-1, -1, 2, '<=', -3).vector_form)

'''constraints = [
        [2., 1., 80., '>='],
        [1., 2., 60., '>=']
    ]
type = "min"

target = [600.,500.]'''
'''constraints = [
        [2., 1., 3., '<='],
        [1., 2., 9., '<=']
    ]
type = "max"

target = [8.,8.]'''
'''constraints = [
        [2., 1., 8., '>='],
        [1., 2., 8., '>=']
    ]
type = "max"

target = [3.,9.]'''
'''constraints = [
        [4., 8., 12., '<='],
        [2., 1., 3., '<=']
    ]
type = "max"

target = [2.,3.]'''

'''constraints = [
        [2., 3., 2., 1000., '<='],
        [1., 1., 2., 800., '<=']
    ]
type = "max"

target = np.array([7,8,10])'''