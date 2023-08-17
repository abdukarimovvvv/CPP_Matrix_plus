#include "s21_matrix_oop.h"

#include <cmath>
#include <cstring>
#include <iostream>

void S21Matrix::Create(int rows, int cols) {
  rows_ = rows;  // Устанавливаем количество строк матрицы
  cols_ = cols;  // Устанавливаем количество столбцов матрицы

  matrix_ = new double *[rows_]();  // Выделяем память для указателей на строки
                                    // матрицы (динамический массив указателей)
  matrix_[0] =
      new double[rows_ * cols_]();  // Выделяем память для хранения элементов
                                    // матрицы в виде одномерного массива

  // Устанавливаем указатели каждой строки на соответствующую область в
  // одномерном массиве
  for (size_t i = 1; i < (size_t)rows_; ++i) {
    matrix_[i] = matrix_[i - 1] + cols_;
  }
}

S21Matrix::S21Matrix() {
  Create(1, 1);  // Создаем объект S21Matrix с одной строкой и одним столбцом,
                 // используя функцию Create()
}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows < 1 || cols < 1) {
    throw std::out_of_range(
        "Incorrect input, matrices should have cols and rows");
  }
  Create(rows, cols);  // Создаем объект S21Matrix с заданным количеством строк
                       // и столбцов, используя функцию Create()
}

S21Matrix::S21Matrix(const S21Matrix &other) {
  rows_ = other.rows_;
  cols_ = other.cols_;

  Create(rows_, cols_);  // Создаем объект S21Matrix с таким же количеством
                         // строк и столбцов, как в объекте other

  for (size_t i = 0; i < (size_t)rows_; ++i) {
    for (size_t j = 0; j < (size_t)cols_; ++j) {
      matrix_[i][j] = other.matrix_[i][j];  // Копируем элементы матрицы из
                                            // объекта other в текущий объект
    }
  }
}

S21Matrix::S21Matrix(S21Matrix &&other) noexcept {
  matrix_ = other.matrix_;  // Перемещаем указатель на матрицу из объекта other
                            // в текущий объект
  rows_ = other.rows_;  // Перемещаем количество строк из объекта other в
                        // текущий объект
  cols_ = other.cols_;  // Перемещаем количество столбцов из объекта other в
                        // текущий объект

  other.matrix_ =
      nullptr;  // Устанавливаем указатель на матрицу в объекте other в nullptr,
                // чтобы избежать двойного удаления памяти
  other.cols_ = 0;  // Устанавливаем количество столбцов в объекте other в 0
  other.rows_ = 0;  // Устанавливаем количество строк в объекте other в 0
}

S21Matrix::~S21Matrix() {
  if (matrix_ != nullptr) {
    delete[] matrix_[0];  // Освобождаем память для первого ряда матрицы
    delete[] matrix_;  // Освобождаем память для массива указателей на ряды
                       // матрицы
    matrix_ = nullptr;  // Устанавливаем указатель на матрицу в nullptr, чтобы
                        // избежать ошибок при дальнейшем использовании
    rows_ = 0;  // Устанавливаем количество строк в 0
    cols_ = 0;  // Устанавливаем количество столбцов в 0
  }
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  if (&other != this) {  // Проверяем, не является ли присваиваемый объект
                         // текущим объектом
    if (rows_ != other.rows_ ||
        cols_ != other.cols_) {  // Проверяем, отличаются ли размеры матриц
      delete[] matrix_[0];  // Освобождаем память, выделенную для строк
      delete[] matrix_;  // Освобождаем память, выделенную для указателей на
                         // строки
      Create(other.rows_, other.cols_);  // Выделяем новую память для матрицы с
                                         // размерами другого объекта
    }

    rows_ = other.rows_;  // Присваиваем новые значения для количества строк
    cols_ = other.cols_;  // Присваиваем новые значения для количества столбцов

    for (size_t i = 0; i < (size_t)rows_;
         ++i) {  // Копируем элементы из другой матрицы в текущую матрицу
      for (size_t j = 0; j < (size_t)cols_; ++j) {
        matrix_[i][j] = other.matrix_[i][j];
      }
    }
  }
  return *this;  // Возвращаем ссылку на текущий объект для поддержки цепочечных
                 // присваиваний
}

S21Matrix &S21Matrix::operator=(S21Matrix &&other) noexcept {
  if (&other !=
      this) {  // Проверяем, не является ли перемещаемый объект текущим объектом
    delete[] matrix_[0];  // Освобождаем память, выделенную для строк в текущем
                          // объекте
    delete[] matrix_;  // Освобождаем память, выделенную для указателей на
                       // строки в текущем объекте

    rows_ = other.rows_;  // Присваиваем новые значения для количества строк из
                          // перемещаемого объекта
    cols_ = other.cols_;  // Присваиваем новые значения для количества столбцов
                          // из перемещаемого объекта
    matrix_ = other.matrix_;  // Перемещаем указатель на матрицу из
                              // перемещаемого объекта в текущий объект

    other.matrix_ = nullptr;  // Устанавливаем указатель на матрицу в
                              // перемещаемом объекте в nullptr
    other.rows_ =
        0;  // Устанавливаем количество строк в перемещаемом объекте в 0
    other.cols_ =
        0;  // Устанавливаем количество столбцов в перемещаемом объекте в 0
  }
  return *this;  // Возвращаем ссылку на текущий объект для поддержки цепочечных
                 // присваиваний
}

double &S21Matrix::operator()(int row, int col) {
  // Проверяем, что индексы row и col находятся в допустимом диапазоне
  if (row >= rows_ || col >= cols_ || col < 0 || row < 0) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  // Возвращаем ссылку на элемент матрицы с указанными индексами row и col
  return matrix_[row][col];
}

double &S21Matrix::operator()(int row, int col) const {
  // Проверяем, что индексы row и col находятся в допустимом диапазоне
  if (row >= rows_ || col >= cols_ || col < 0 || row < 0) {
    throw std::out_of_range("Incorrect input, index is out of range");
  }
  // Возвращаем ссылку на константный элемент матрицы с указанными индексами row
  // и col
  return matrix_[row][col];
}

bool S21Matrix::operator==(const S21Matrix &other) {
  // Используется метод EqMatrix для проверки равенства матриц
  return EqMatrix(other);
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  // Используется метод SumMatrix для прибавления матрицы other к текущей
  // матрице
  SumMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  // Создается новая матрица result, копирующая текущую матрицу
  S21Matrix result(*this);
  // Используется метод SumMatrix для прибавления матрицы other к матрице result
  result.SumMatrix(other);
  return result;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  // Используется метод SubMatrix для вычитания матрицы other из текущей матрицы
  SubMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  // Создается новая матрица result, копирующая текущую матрицу
  S21Matrix result(*this);
  // Используется метод SubMatrix для вычитания матрицы other из матрицы result
  result.SubMatrix(other);
  return result;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  // Создается новая матрица result, копирующая текущую матрицу
  S21Matrix result(*this);
  // Используется метод MulMatrix для умножения матрицы result на матрицу other
  result.MulMatrix(other);
  return result;
}

S21Matrix operator*(double num, S21Matrix &matrix) {
  // Создается новая матрица result, копирующая матрицу matrix
  S21Matrix result(matrix);
  // Используется метод MulNumber для умножения матрицы result на число num
  result.MulNumber(num);
  return result;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  // Используется метод MulMatrix для умножения текущей матрицы на матрицу other
  MulMatrix(other);
  return *this;
}

S21Matrix S21Matrix::operator*(const double num) {
  // Создается новая матрица result, копирующая текущую матрицу
  S21Matrix result(*this);
  // Используется метод MulNumber для умножения матрицы result на число num
  result.MulNumber(num);
  return result;
}

S21Matrix &S21Matrix::operator*=(const double num) {
  // Используется метод MulNumber для умножения текущей матрицы на число num
  MulNumber(num);
  return *this;
}

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  bool result = true;

  // Проверяем, имеют ли матрицы одинаковые размеры
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    for (size_t i = 0; i < (size_t)rows_; ++i) {
      for (size_t j = 0; j < (size_t)cols_; ++j) {
        // Сравниваем элементы матриц с заданной точностью
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) > eps) {
          result = false;
        }
      }
    }
  } else {
    // Если размеры матриц не совпадают, считаем их не равными
    result = false;
  }

  return result;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  // Проверяем, что матрицы имеют одинаковый размер
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }

  // Прибавляем элементы другой матрицы к текущей матрице
  for (size_t i = 0; i < (size_t)rows_; ++i) {
    for (size_t j = 0; j < (size_t)cols_; ++j) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  // Проверяем, что матрицы имеют одинаковый размер
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::out_of_range(
        "Incorrect input, matrices should have the same size");
  }

  // Вычитаем элементы другой матрицы из текущей матрицы
  for (size_t i = 0; i < (size_t)rows_; ++i) {
    for (size_t j = 0; j < (size_t)cols_; ++j) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  // Умножаем каждый элемент матрицы на заданное число
  for (size_t i = 0; i < (size_t)rows_; ++i) {
    for (size_t j = 0; j < (size_t)cols_; ++j) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  // Проверяем размерности матриц для выполнения умножения
  if (cols_ != other.rows_) {
    throw std::out_of_range(
        "Incorrect input, number of columns in the first matrix should be "
        "equal to the number of rows in the second matrix");
  }

  // Создаем временную матрицу для хранения результирующей матрицы
  S21Matrix result(rows_, other.cols_);

  // Выполняем умножение матриц
  for (size_t i = 0; i < (size_t)rows_; ++i) {
    for (size_t j = 0; j < (size_t)other.cols_; ++j) {
      for (size_t s = 0; s < (size_t)other.rows_; ++s) {
        result.matrix_[i][j] += matrix_[i][s] * other.matrix_[s][j];
      }
    }
  }
  // Копируем результат в текущую матрицу
  *this = result;
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result(cols_,
                   rows_);  // Создаем новую матрицу с перевернутыми размерами

  for (size_t i = 0; i < (size_t)cols_; ++i) {
    for (size_t j = 0; j < (size_t)rows_; ++j) {
      // Заполняем новую матрицу значениями из текущей матрицы в обратном
      // порядке
      result.matrix_[i][j] = matrix_[j][i];
    }
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  S21Matrix result(rows_, cols_);  // Создание матрицы-результата с тем же
                                   // размером, что и исходная матрица
  if (rows_ != cols_) {
    throw std::invalid_argument(
        "The matrix is not square");  // Проверка, является ли исходная матрица
                                      // квадратной
  }
  for (size_t i = 0; i < (size_t)rows_; ++i) {
    for (size_t j = 0; j != (size_t)cols_; ++j) {
      // Создание минорной матрицы, исключая строку i и столбец j
      S21Matrix minor_matrix = Minor(i, j);
      // Расчет значения комплементарного элемента по формуле: (-1)^(i+j) *
      // определитель минорной матрицы
      result.matrix_[i][j] = pow((-1), i + j) * minor_matrix.Determinant();
    }
  }
  return result;  // Возвращение матрицы-результата, содержащей комплементарные
                  // элементы
}

double S21Matrix::Determinant() {
  // Проверка, является ли матрица квадратной
  if (rows_ != cols_) {
    throw std::invalid_argument("The matrix is not square");
  }
  double result = 0.0;

  // Если матрица состоит из одного элемента, возвращаем значение этого элемента
  if (rows_ == 1) {
    result = matrix_[0][0];
  }
  // Если матрица 2x2, вычисляем определитель по формуле
  else if (rows_ == 2) {
    result = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  }
  // Для более крупных матриц вычисляем определитель методом разложения по
  // минорам
  else {
    for (size_t j = 0; j < (size_t)cols_; ++j) {
      // Получаем минорную матрицу, исключая первую строку и столбец j
      S21Matrix minor_matrix = Minor(0, j);

      // Вычисляем вклад j-го элемента в определитель с учетом знака (-1) в
      // степени j
      result += matrix_[0][j] * pow(-1, j) * minor_matrix.Determinant();
    }
  }
  return result;
}

S21Matrix S21Matrix::Minor(int row, int col) {
  // Создание матрицы-минора с размерами на 1 меньше исходной матрицы
  S21Matrix result(rows_ - 1, cols_ - 1);

  // Итерация по строкам исходной матрицы
  for (size_t i = 0, min_i = 0; min_i < (size_t)result.rows_; ++min_i) {
    // Пропуск строки, соответствующей заданному параметру row
    if ((size_t)row == i) ++i;

    // Итерация по столбцам исходной матрицы
    for (size_t j = 0, min_j = 0; min_j < (size_t)result.cols_; ++min_j) {
      // Пропуск столбца, соответствующего заданному параметру col
      if ((size_t)col == j) ++j;

      // Копирование элемента из исходной матрицы в матрицу-минор
      result.matrix_[min_i][min_j] = matrix_[i][j];
      ++j;  // Увеличение счетчика столбца исходной матрицы
    }
    ++i;  // Увеличение счетчика строки исходной матрицы
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  // Вычисление определителя матрицы
  double det = Determinant();
  if (fabs(det) < eps) {
    throw std::invalid_argument("Matrix determinant is 0");
  }
  // Создание матрицы-результата с теми же размерами, что и исходная матрица
  S21Matrix result(rows_, cols_);

  // Если матрица имеет размер 1x1, обратная матрица вычисляется путем деления 1
  // на этот элемент
  if (rows_ == 1) {
    result.matrix_[0][0] = 1 / matrix_[0][0];
  } else {
    // Вычисление алгебраических дополнений исходной матрицы
    S21Matrix tmp = CalcComplements();
    // Транспонирование матрицы алгебраических дополнений
    result = tmp.Transpose();
    // Умножение матрицы алгебраических дополнений на обратное значение
    // определителя
    result.MulNumber(1 / det);
  }

  return result;
}

void S21Matrix::SetRows(const int rows) {
  if (rows < 1) {
    throw std::out_of_range(
        "Incorrect input, matrices should have cols and rows");
  }

  // Создание временной матрицы с новым количеством строк и тем же количеством
  // столбцов
  S21Matrix result(rows, cols_);

  // Копирование элементов из исходной матрицы во временную матрицу
  for (size_t i = 0; i < (size_t)rows; ++i) {
    for (size_t j = 0; j < (size_t)cols_; ++j) {
      // Проверка, что индекс строки не превышает оригинальное количество строк
      if (i < (size_t)rows_) {
        result.matrix_[i][j] = matrix_[i][j];
      }
    }
  }
  // Присвоение временной матрицы текущему объекту
  *this = result;
}

void S21Matrix::SetCols(const int cols) {
  if (cols <= 0) {
    throw std::out_of_range(
        "Incorrect input, matrices should have cols and rows");
  }
  // Создание временной матрицы с тем же количеством строк и новым количеством
  // столбцов
  S21Matrix result(rows_, cols);

  // Копирование элементов из исходной матрицы во временную матрицу
  for (size_t i = 0; i < (size_t)rows_; ++i) {
    for (size_t j = 0; j < (size_t)cols; ++j) {
      // Проверка, что индекс столбца не превышает оригинальное количество
      // столбцов
      if (j < (size_t)cols_) {
        result.matrix_[i][j] = matrix_[i][j];
      }
    }
  }
  // Проверка, что индекс столбца не превышает оригинальное количество столбцов
  *this = result;
}