import math
import copy


class Matrix:

    def size_m(self):
        return self.dimensions()[0]

    def size_n(self):
        return self.dimensions()[1]

    def dimensions(self):
        return len(self.matrix_rows), len(self.matrix_rows[0])

    def calculate_vector_norm(self):
        if self.dimensions()[1] == 1:
            norm = 0
            for element in self.matrix_rows:
                norm += element[0] ** 2
            return math.sqrt(norm)
        else:
            print("Matrix is not a vector!")
            exit()

    def __add__(self, second_matrix):
        if self.dimensions() == second_matrix.dimensions():
            tmp_mat = Matrix([])
            tmp_rows = self.matrix_rows
            for M in range(0, self.size_m()):
                for N in range(0, self.size_n()):
                    tmp_rows[M][N] += second_matrix.matrix_rows[M][N]
            tmp_mat.matrix_rows = tmp_rows
            return tmp_mat
        else:
            print("Matrix dimensions don't match!")
            exit()

    def __sub__(self, second_matrix):
        if self.dimensions() == second_matrix.dimensions():
            tmp_mat = Matrix([])
            tmp_rows = self.matrix_rows
            for M in range(0, self.size_m()):
                for N in range(0, self.size_n()):
                    tmp_rows[M][N] -= second_matrix.matrix_rows[M][N]
            tmp_mat.matrix_rows = tmp_rows
            return tmp_mat
        else:
            print("Matrix dimensions don't match!")
            exit()

    def to_string(self):
        text = ""
        for row in self.matrix_rows:
            for element in row:
                text += str(element)+" "
            text += "\n"
        return text

    def __mul__(self, second_matrix):
        if self.size_n() == second_matrix.size_m():
            mat = Matrix().uniform_matrix(self.size_m(), second_matrix.size_n(), 0)
            for rows1 in range(0, self.size_m()):
                for columns2 in range(0, second_matrix.size_n()):
                    dot_prod = 0
                    for i in range(0, self.size_n()):
                        dot_prod += self.matrix_rows[rows1][i] * second_matrix.matrix_rows[i][columns2]
                    mat.matrix_rows[rows1][columns2] = dot_prod
            return mat

        else:
            print("Matrix dimensions don't match!")
            exit()

    def __init__(self, *rows):
        self.matrix_rows = []
        for row in rows:
            self.matrix_rows.append(row)

    @staticmethod
    def uniform_matrix(m, n, fill):
        mat = Matrix()
        for i_m in range(0, m):
            row = [fill] * n
            mat.matrix_rows.append(row)
        return mat

    @staticmethod
    def identity_matrix(M, N):
        mat_ones = Matrix().uniform_matrix(M, N, 0)
        for m in range(mat_ones.size_m()):
            for n in range(mat_ones.size_n()):
                if m == n:
                    mat_ones.matrix_rows[m][n] = 1
        return mat_ones

    @staticmethod
    def vector(n, fill):
        mat = Matrix().uniform_matrix(n, 1, fill)
        return mat

    @staticmethod
    def get_upper_and_lower_matrix(matrix):
        mat_U = matrix
        mat_L = Matrix().identity_matrix(matrix.size_m(), matrix.size_n())
        for k in range(matrix.size_m()-1):
            for j in range(k+1, matrix.size_m()):
                mat_L.matrix_rows[j][k] = mat_U.matrix_rows[j][k] / mat_U.matrix_rows[k][k]
                for x in range(k, matrix.size_m()):
                    mat_U.matrix_rows[j][x] = mat_U.matrix_rows[j][x] - mat_L.matrix_rows[j][k]*mat_U.matrix_rows[k][x]
        return mat_U, mat_L
