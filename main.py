from Matrix import *
import math


def create_A_matrix(N, a1, a2, a3):
    A = Matrix().uniform_matrix(N, N, 0)
    for i in range(0, N):
        for j in range(0, N):
            if i == j:
                A.matrix_rows[i][j] = a1
                if i > 0:
                    A.matrix_rows[i-1][j] = a2
                if i < N - 1:
                    A.matrix_rows[i+1][j] = a2
                if i > 1:
                    A.matrix_rows[i - 2][j] = a3
                if i < N - 2:
                    A.matrix_rows[i + 2][j] = a3
    return A


def Jacobi(A, b, x, res_error):
    x_new = Matrix().vector(x.size_m(), 0)
    MAX_JACOBI_ITER = 1000

    for r in range(MAX_JACOBI_ITER):
        for i in range(x.size_m()):
            iter_sum = 0
            for j in range(A.size_n()):
                if j != i:
                    iter_sum += x.matrix_rows[j][0] * A.matrix_rows[i][j]

            x_new.matrix_rows[i][0] = (b.matrix_rows[i][0] - iter_sum) / A.matrix_rows[i][i]
        x = x_new
        res = A * x - b
        if res.calculate_vector_norm() <= res_error:
            print(r)
            return x


def Gauss_Seidler(A, b, x, res_error):
    x_new = Matrix().vector(x.size_m(), 0)
    MAX_GAUSS_ITER = 1000

    for r in range(MAX_GAUSS_ITER):
        for i in range(x.size_m()):
            iter_sum = 0
            for j in range(A.size_n()):

                if j < i:
                    iter_sum += x_new.matrix_rows[j][0] * A.matrix_rows[i][j]
                if j > i:
                    iter_sum += x.matrix_rows[j][0] * A.matrix_rows[i][j]

            x_new.matrix_rows[i][0] = (b.matrix_rows[i][0] - iter_sum) / A.matrix_rows[i][i]
        x = x_new
        res = A * x - b
        if res.calculate_vector_norm() <= res_error:
            print(r)
            return x

def zadA(f, e, c, d):
    a1 = 5 + e
    a2 = -1
    a3 = a2
    N = 900 + 10 * c + d
    A = create_A_matrix(N, a1, a2, a3)
    b = Matrix().vector(N, 0)
    for n in range(0, N):
        b.matrix_rows[n][0] = math.sin(n * (f + 1))

    x = Matrix().vector(N, 0)
    res_error = 10 ** -9
    x_jacobi = Jacobi(A, b, x, res_error)
    x_gauss = Gauss_Seidler(A, b, x, res_error)

if __name__ == '__main__':
    zadA(8, 8, 1, 3)
