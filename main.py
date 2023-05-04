from Matrix import *
import math
import matplotlib.pyplot as plt
import time


def generate_plot(x_axis, y_axis, title, x_title, y_title, filename):
    if filename != "":
        plt.plot(x_axis, y_axis)
        plt.ylabel(y_title)
        plt.yscale('log')
        plt.xlabel(x_title)
        plt.title(title)
        plt.savefig(filename + ".png")
        plt.show()


def calc_res_error_norm(A, x, b):
    res = A * x - b
    return res.calculate_vector_norm()


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


def create_b_matrix(N, f):
    b = Matrix().vector(N, 0)
    for n in range(0, N):
        b.matrix_rows[n][0] = math.sin(n * (f + 1))
    return b


def Jacobi(A, b, res_error, plot_filename):
    x = Matrix().vector(A.size_n(), 0)
    x_new = Matrix().vector(x.size_m(), 0)
    MAX_JACOBI_ITER = 100

    plot_iters = []
    plot_res_error_norm = []

    start = time.time()
    for r in range(MAX_JACOBI_ITER):
        plot_iters.append(r)

        for i in range(x.size_m()):
            iter_sum = 0
            for j in range(i):
                iter_sum += (x.matrix_rows[j][0] * A.matrix_rows[i][j])
            for j in range(i + 1, A.size_n()):
                iter_sum += (x.matrix_rows[j][0] * A.matrix_rows[i][j])

            x_new.matrix_rows[i][0] = (b.matrix_rows[i][0] - iter_sum) / A.matrix_rows[i][i]

        for i in range(x.size_m()):
            x.matrix_rows[i][0] = x_new.matrix_rows[i][0]

        curr_res_error = calc_res_error_norm(A, x, b)
        plot_res_error_norm.append(curr_res_error)
        if curr_res_error <= res_error:
            end = time.time()
            generate_plot(plot_iters, plot_res_error_norm, "wykres norm bledow rezydualnych dla metody Jacobiego",
                          "nr iteracji", "norma bledu rezydualnego", plot_filename)

            print("Metoda Jacobiego: ")
            print("czas trwania [s]: ", (end - start))
            print("norma bledu rezydualnego: ", curr_res_error)
            print("wymagana ilosc iteracji: ", r)
            return end-start
    end = time.time()
    generate_plot(plot_iters, plot_res_error_norm, "wykres norm bledow rezydualnych dla metody Jacobiego",
                  "nr iteracji", "norma bledu rezydualnego", plot_filename)
    return end-start


def Gauss_Seidel(A, b, res_error, plot_filename):
    x = Matrix().vector(A.size_n(), 0)
    x_new = Matrix().vector(x.size_m(), 0)
    MAX_GAUSS_ITER = 100

    plot_iters = []
    plot_res_error_norm = []

    start = time.time()
    for r in range(MAX_GAUSS_ITER):
        plot_iters.append(r)

        for i in range(x.size_m()):
            iter_sum = 0
            for j in range(i):
                iter_sum += (x_new.matrix_rows[j][0] * A.matrix_rows[i][j])
            for j in range(i+1, A.size_n()):
                iter_sum += (x.matrix_rows[j][0] * A.matrix_rows[i][j])

            x_new.matrix_rows[i][0] = (b.matrix_rows[i][0] - iter_sum) / A.matrix_rows[i][i]
        for i in range(x.size_m()):
            x.matrix_rows[i][0] = x_new.matrix_rows[i][0]

        curr_res_error = calc_res_error_norm(A, x, b)
        plot_res_error_norm.append(curr_res_error)
        if curr_res_error <= res_error:
            end = time.time()
            generate_plot(plot_iters, plot_res_error_norm, "wykres norm bledow rezydualnych dla metody Gaussa-Seidela",
                          "nr iteracji", "norma bledu rezydualnego", plot_filename)

            print("Metoda Gaussa-Seidla: ")
            print("czas trwania [s]: ", (end - start))
            print("norma bledu rezydualnego: ", curr_res_error)
            print("wymagana ilosc iteracji: ", r)
            return end-start
    end = time.time()
    generate_plot(plot_iters, plot_res_error_norm, "wykres norm bledow rezydualnych dla metody Gaussa-Seidela",
                  "nr iteracji", "norma bledu rezydualnego", plot_filename)
    return end-start


def back_substitution(upper_triangular, right_side_vector):
    x = Matrix.vector(upper_triangular.size_m(), 0)
    n = upper_triangular.size_m()
    for i in range(n-1, -1, -1):
        subSum = 0
        for k in range(n-1, i, -1):
            subSum += x.matrix_rows[k][0]*upper_triangular.matrix_rows[i][k]
        x.matrix_rows[i][0] = (right_side_vector.matrix_rows[i][0] - subSum)/upper_triangular.matrix_rows[i][i]
    return x


def forward_substitution(lower_triangular, right_side_vector):
    y = Matrix.vector(lower_triangular.size_n(), 0)
    for i in range(1, y.size_m()):
        subSum = 0
        for j in range(0, i):
            subSum += (lower_triangular.matrix_rows[i][j] * y.matrix_rows[j][0])
        y.matrix_rows[i][0] = right_side_vector.matrix_rows[i][0] - subSum
    return y


def zadE(f, e):
    a1 = 5 + e
    a2 = -1
    a3 = a2
    N_arr = [100, 500, 1000, 2000]
    res_error = 10 ** -9
    y_axis_gauss = []
    y_axis_jacobi = []
    y_axis_LU = []
    for N in N_arr:
        A = create_A_matrix(N, a1, a2, a3)
        A_tmp = create_A_matrix(N, a1, a2, a3)
        b = create_b_matrix(N, f)
        y_axis_jacobi.append(Jacobi(A, b, res_error, ""))
        y_axis_gauss.append(Gauss_Seidel(A, b, res_error, ""))
        y_axis_LU.append(LU_decomposition(A, A_tmp, b))

    plt.plot(N_arr, y_axis_jacobi, label="Jacobi")
    plt.plot(N_arr, y_axis_gauss, label="Gauss-Seidel")
    plt.plot(N_arr, y_axis_LU, label="LU")
    plt.ylabel("Czas dzialania [s]")
    plt.yscale('log')
    plt.xlabel("Numer iteracji")
    plt.title("Porownanie czasu dzialania poszczegolnych algorytmow")
    plt.legend()
    plt.savefig("compared.png")
    plt.show()


def LU_decomposition(A, A_copy, b):
    start = time.time()
    U, L = Matrix().get_upper_and_lower_matrix(A_copy)
    y = forward_substitution(L, b)
    x = back_substitution(U, y)
    end = time.time()
    print("Metoda faktoryzacji LU: ")
    print("czas trwania [s]: ", (end - start))
    print("norma bledu rezydualnego: ", calc_res_error_norm(A,x,b))
    return end - start


def zadD(f, c, d):
    a1 = 3
    a2 = -1
    a3 = a2
    N = 900 + 10 * c + d
    A = create_A_matrix(N, a1, a2, a3)
    A_tmp = create_A_matrix(N, a1, a2, a3)
    b = create_b_matrix(N, f)
    LU_decomposition(A, A_tmp, b)


def zadABC(f, e, c, d, plot_jacobi_filename, plot_gauss_filename):
    a1 = 5 + e
    a2 = -1
    a3 = a2
    N = 900 + 10 * c + d
    A = create_A_matrix(N, a1, a2, a3)
    b = create_b_matrix(N, f)

    res_error = 10 ** -9
    Jacobi(A, b, res_error, plot_jacobi_filename)
    Gauss_Seidel(A, b, res_error, plot_gauss_filename)


if __name__ == '__main__':
    #zad A i B
    zadABC(7, 5, 8, 9, "Jacobi B", "Gauss B")

    #zad C - e equals -2 so a1 will be 3
    zadABC(7, -2, 8, 9, "Jacobi C", "Gauss C")

    #zad D
    zadD(7, 8, 9)

    #zad E
    zadE(7, 5)
