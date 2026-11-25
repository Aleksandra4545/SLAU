#include <iostream>
#include <fstream>
#include <array>
#include <type_traits>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

// Предварительное объявление шаблонного класса
template <typename T, size_t N>
class LinearSystem;

// Предварительное объявление дружественной функции
template <typename T, size_t N>
void scaleEquation(LinearSystem<T, N>& system,
                   size_t eqIndex, T factor);

// 1. ОСНОВНОЙ ШАБЛОННЫЙ КЛАСС

template <typename T, size_t N>
class LinearSystem {
    static_assert(is_floating_point_v<T>,
                  "Only floating point types are allowed for LinearSystem");

private:
    T A[N][N];       // Матрица коэффициентов
    T B[N];          // Вектор свободных членов

    static inline size_t totalSystemsSolved = 0;

public:

    // Конструктор
    LinearSystem() {
        for (size_t i = 0; i < N; i++) {
            B[i] = 0;

            for (size_t j = 0; j < N; j++) {
                A[i][j] = 0;
            }
        }
    }

    // Получение коэффициента матрицы
    T getA(size_t i, size_t j) const {
        return A[i][j];
    }

    // Получение свободного члена
    T getB(size_t i) const {
        return B[i];
    }

    // Изменение коэффициента
    void setA(size_t i, size_t j, T value) {
        A[i][j] = value;
    }

    // Изменение свободного члена
    void setB(size_t i, T value) {
        B[i] = value;
    }

    // Увеличение счетчика решенных систем
    static void incrementCounter() {
        totalSystemsSolved++;
    }

    // Получение количества решенных систем
    static size_t getTotalSystemsSolved() {
        return totalSystemsSolved;
    }

    // Дружественная шаблонная функция
    template <typename U, size_t M>
    friend void scaleEquation(LinearSystem<U, M>& system,
                              size_t eqIndex, U factor);
};

template <>
class LinearSystem<double, 3> {

private:
    double A[3][3];
    double B[3];

    static inline size_t totalSystemsSolved = 0;

public:

    LinearSystem() {
        for (int i = 0; i < 3; i++) {
            B[i] = 0;

            for (int j = 0; j < 3; j++) {
                A[i][j] = 0;
            }
        }
    }

    double getA(size_t i, size_t j) const {
        return A[i][j];
    }

    double getB(size_t i) const {
        return B[i];
    }

    void setA(size_t i, size_t j, double value) {
        A[i][j] = value;
    }

    void setB(size_t i, double value) {
        B[i] = value;
    }

    static void incrementCounter() {
        totalSystemsSolved++;
    }

    static size_t getTotalSystemsSolved() {
        return totalSystemsSolved;
    }

    // Дружественная функция
    template <typename U, size_t M>
    friend void scaleEquation(LinearSystem<U, M>& system,
                              size_t eqIndex, U factor);
};

template <typename T, size_t N>
void scaleEquation(LinearSystem<T, N>& system,
                   size_t eqIndex, T factor) {

    if (eqIndex >= N) {
        throw out_of_range("Incorrect equation index!");
    }

    // Умножаем все коэффициенты выбранного уравнения
    for (size_t j = 0; j < N; j++) {
        system.A[eqIndex][j] *= factor;
    }

    // Умножаем свободный член
    system.B[eqIndex] *= factor;
}

template <typename T, size_t N>
void printSystem(const LinearSystem<T, N>& system) {

    cout << "System of equations:\n";

    for (size_t i = 0; i < N; i++) {

        for (size_t j = 0; j < N; j++) {

            cout << setw(9) << system.getA(i, j)
                << " * x" << j + 1;
            if (j < N - 1) { 
                cout << " + ";

            }
        }

        cout << " = " << system.getB(i) << '\n';
    }
}

double determinant(double a[3][3]) {

    return
        a[0][0] * (a[1][1] * a[2][2] -
                   a[1][2] * a[2][1])

      - a[0][1] * (a[1][0] * a[2][2] -
                   a[1][2] * a[2][0])

      + a[0][2] * (a[1][0] * a[2][1] -
                   a[1][1] * a[2][0]);
}

