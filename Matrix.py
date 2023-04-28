class Matrix:

    def dimensions(self):
        return len(self.matrix_rows), len(self.matrix_rows[0])

    def __add__(self, second_matrix):
        if self.dimensions() == second_matrix.dimensions():
            tmp_mat = Matrix([])
            tmp_rows = self.matrix_rows
            for M in range(len(tmp_rows)):
                for N in range(len(tmp_rows[0])):
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
            for M in range(len(tmp_rows)):
                for N in range(len(tmp_rows[0])):
                    tmp_rows[M][N] -= second_matrix.matrix_rows[M][N]
            tmp_mat.matrix_rows = tmp_rows
            return tmp_mat
        else:
            print("Matrix dimensions don't match!")
            exit()

    # def __mul__(self, second_matrix):
    #     if self.dimensions()[0] == second_matrix.dimensions()[1]:
    #
    #     else:
    #         print("matrix dimensions don't match!")
    #         return Matrix([])

    def __init__(self, *rows):
        self.matrix_rows = []
        for row in rows:
            self.matrix_rows.append(row)
