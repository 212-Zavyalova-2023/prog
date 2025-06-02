
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

struct Point { // структура для хранения координат точки
  int x, y;
};

struct Color { // структура для хранения rgb цвета
  uint8_t r, g, b;
};



class Server { //Класс, который отвечает за логирование сервера (Contrl)
  std::ofstream server_log; // Серверный лог
  bool logging_enabled =
      false; // Нужно для конфига, чтобы включать или выключать логи
  std::string log_filename = "server_log.txt"; // Имя файла лога (по умолчанию)
public:
  Server() {}

  ~Server() {
    if (server_log.is_open()) {
      server_log.close();
    }
  }

  void enable_logging(const std::string &filename) { // Функция для включения
                                                     // логирования сервера
    logging_enabled = true;
    log_filename = filename;
    if (server_log.is_open()) {
      server_log.close();
    }
    server_log.open(log_filename, std::ios::app);
  }

  void disable_logging() { // Функция для отключения логирования сервера
    logging_enabled = false;
    if (server_log.is_open()) {
      server_log.close();
    }
  }

  void log(const std::string &message) { // Функция-шаблон лога
    if (logging_enabled && server_log.is_open()) {
      server_log << current_time() << " - " << message << std::endl;
    }
  }

  std::string current_time() { // Как будет отображаться время в логах
    std::time_t t = std::time(nullptr);
    std::tm *now = std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(now, "%d.%m.%Y %H:%M:%S");
    return oss.str();
  }
};



class Client { // Класс, овтечающий за логирвоание клиента (Interface)
  std::ofstream client_log; // Клиентский лог
  bool logging_enabled =
      false; // //Нужно для конфига, чтобы включать или выключать логи
  std::string log_filename = "client_log.txt"; // Имя файла лога (по умолчанию)

public:
  Client() {}

  ~Client() {
    if (client_log.is_open()) {
      client_log.close();
    }
  }

  void enable_logging(const std::string &filename) { // Функция для включения
                                                     // логирования клиента
    logging_enabled = true;
    log_filename = filename;
    if (client_log.is_open()) {
      client_log.close();
    }
    client_log.open(log_filename, std::ios::app);
  }

  void disable_logging() { // Функция для выключения логирования клиента
    logging_enabled = false;
    if (client_log.is_open()) {
      client_log.close();
    }
  }

  void log(const std::string &message) { // Функция-шаблон лога
    if (logging_enabled && client_log.is_open()) {
      client_log << current_time() << " - " << message << std::endl;
    }
  }

  std::string current_time() { // Как будет отоборажаться время в логах
    std::time_t t = std::time(nullptr);
    std::tm *now = std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(now, "%d.%m.%Y %H:%M:%S");
    return oss.str();
  }
};



class Component { // Класс, хранящий все компоненты в нужном нам срезе
public:
  std::vector<std::vector<double>> comp;
  Component(const std::vector<std::vector<double>> &inpcomp) : comp(inpcomp) {}

  Component(int A, int B) { comp.resize(A, std::vector<double>(B, 255)); }
};

class Gauss { // Класс, хранящий Гауссы
public:
  double h, x0, y0, sigma_x, sigma_y;
  Gauss(double h, double x0, double y0, double sigma_x, double sigma_y)
      : h(h), x0(x0), y0(y0), sigma_x(sigma_x), sigma_y(sigma_y) {}

  double calculate(int x, int y) const { //Функция вычисления Гаусса
    double dx = x - x0;
    double dy = y - y0;
    double denom_x = 2 * sigma_x * sigma_x;
    double denom_y = 2 * sigma_y * sigma_y;

    if ((int)denom_x == 0 || (int)denom_y == 0) {
      std::cerr << "Error: sigma_x or sigma_y iz zero!" << std::endl;
      return 0;
    }

    double exponent = -((dx * dx) / denom_x + (dy * dy) / denom_y);

    if (exponent < -700) {
      return 0;
    }

    else if (exponent > 700) {
      return std::numeric_limits<double>::infinity();
    }

    return h * exp(exponent);
  }
};



class Field { //Класс поле
public:
  int length, width;
  std::vector<std::vector<double>> matrix;
  std::vector<std::vector<double>> weight_matrix; // Матрица сумм весов
  std::vector<std::vector<Color>> colors; //Матрица цветов для bmp (если нужно)

  Field(int l, int w) : length(l), width(w) {
    matrix.resize(length, std::vector<double>(width, 0));
    weight_matrix.resize(length, std::vector<double>(width, 0));
  }

  void apply_gauss(const Gauss &g) { //Добавляет гаусс на поле

    double value;

    for (int i = 0; i < length; ++i) {
      for (int j = 0; j < width; ++j) {
        value = g.calculate(i, j);
        matrix[i][j] += value * g.h; // Взвешенное добавление значения
        weight_matrix[i][j] += g.h; // Вес в точке (на основе амплитуд)
      }
    }
  }

  void normalize(void) {
    for (int i = 0; i < length; ++i) {
      for (int j = 0; j < width; ++j) {
        if (weight_matrix[i][j] > 0) {
          matrix[i][j] /= weight_matrix[i][j];
        }
      }
    }
  }

  void
  save_to_gnuplot(const std::string &filename) { // Функция сохранения гнуплота
    std::ofstream bmp_file(filename);
    for (int i = 0; i < length; ++i) {
      for (int j = 0; j < width; ++j) {
        bmp_file << i << " " << j << " " << matrix[i][j] << std::endl;
      }
      bmp_file << std::endl;
    }
    bmp_file.close();
  }

  void bmp_write(const std::string &filename, int k) { // Запись BMP
    const int BMP_HEADER_SIZE = 54;
    const int PIXEL_SIZE = 3;
    int file_size = BMP_HEADER_SIZE + PIXEL_SIZE * length * width;

    unsigned char bmp_header[BMP_HEADER_SIZE] = {0};

    // Заголовок BMP файла
    bmp_header[0] = 'B';
    bmp_header[1] = 'M';
    bmp_header[2] = file_size & 0xFF;
    bmp_header[3] = (file_size >> 8) & 0xFF;
    bmp_header[4] = (file_size >> 16) & 0xFF;
    bmp_header[5] = (file_size >> 24) & 0xFF;
    bmp_header[10] = BMP_HEADER_SIZE;

    // Заголовок DIB
    bmp_header[14] = 40; // Размер заголовка DIB
    bmp_header[18] = width & 0xFF;
    bmp_header[19] = (width >> 8) & 0xFF;
    bmp_header[20] = (width >> 16) & 0xFF;
    bmp_header[21] = (width >> 24) & 0xFF;
    bmp_header[22] = length & 0xFF;
    bmp_header[23] = (length >> 8) & 0xFF;
    bmp_header[24] = (length >> 16) & 0xFF;
    bmp_header[25] = (length >> 24) & 0xFF;
    bmp_header[26] = 1;  // Число цветовых плоскостей
    bmp_header[28] = 24; // Количество бит на пиксель

    std::ofstream bmp_file(filename, std::ios::binary);
    bmp_file.write(reinterpret_cast<char *>(bmp_header), BMP_HEADER_SIZE);

    // Записываем пиксели (матрицу)
    if (k == 1) { //Для ч\б картинки
      for (int i = length - 1; i >= 0; --i) {
        for (int j = 0; j < width; ++j) {
          int value = static_cast<int>(
              matrix[i][j] * 100); // Умножаю на коэффициент 100, чтобы
                                   // отображалось красиво и ярко
          unsigned char pixel[3] = {
              static_cast<unsigned char>(std::min(std::max(value, 0), 255)),
              static_cast<unsigned char>(std::min(std::max(value, 0), 255)),
              static_cast<unsigned char>(std::min(std::max(value, 0), 255))};
          bmp_file.write(reinterpret_cast<char *>(pixel), PIXEL_SIZE);
        }
      }
    }

    else { //для разукраски класстеров
      for (int i = length - 1; i >= 0; --i) {
        for (int j = 0; j < width; ++j) {
          unsigned char pixel[3] = {static_cast<unsigned char>(colors[i][j].r),
                                    static_cast<unsigned char>(colors[i][j].g),
                                    static_cast<unsigned char>(colors[i][j].b)};
          bmp_file.write(reinterpret_cast<char *>(pixel), PIXEL_SIZE);
        }
      }
    }

    bmp_file.close();
  }

  std::pair<int, int> bmp_read(
      const std::string
          &filename) { // Чтение BMP. Узнаем размер поля и инициализируем его.

    std::ifstream bmp_file(filename, std::ios::binary);

    if (!bmp_file) {
      std::cerr << "Failed to open BMP file: " << filename << std::endl;
      return {0, 0};
    }

    unsigned char header[54];
    bmp_file.read(reinterpret_cast<char *>(header), 54);

    if (header[0] != 'B' || header[1] != 'M') {
      std::cerr << "Invalid BMP file:" << filename << std::endl;
      return {0, 0};
    }

    int width = header[18] | (header[19] << 8) | (header[20] << 16) |
                (header[21] << 24);
    int height = header[22] | (header[23] << 8) | (header[24] << 16) |
                 (header[25] << 24);

    return {height, width};
  }

  void load_data(std::ifstream &bmp_file, int length,
                 int width) { // Чтение BMP. Выгружаем информацию в матрицу.
    for (int i = length - 1; i >= 0; --i) {
      for (int j = 0; j < width; j++) {
        unsigned char color = bmp_file.get();
        bmp_file.get();
        bmp_file.get();
        matrix[i][j] = color;
      }
      bmp_file.ignore((4 - (width * 3) % 4) % 4);
    }
  }
};



class Srez { // Класс со срезом (bin) и алгоритмом волна (wave)
  int length, width;
  std::vector<std::vector<double>> CopyField;
  std::vector<Component> components;
  int count = 0;

public:
  void bin(int h, Field &f) { // Функция для получения разреза на высоте h
    CopyField.resize(f.matrix.size(),
                     std::vector<double>(f.matrix[0].size(), 0));
    length = f.matrix.size();
    width = f.matrix[0].size();

    for (int i = length - 1; i >= 0; --i) {
      for (int j = 0; j < width; ++j) {
        if (f.matrix[i][j] < h) {
          CopyField[i][j] = 255;
        } else {
          CopyField[i][j] = 0;
        }
      }
    }
    Field pole(f.length, f.width);
    pole.matrix = CopyField;
    pole.bmp_write("bin.bmp", 1);
  }

  void wave(int n) { // Функция запуска алгоритма wave и записи соотв. компонент
                     // в вектор
    Component Componenta(length, width);
    int c = 0; // Тут будет записано значение count

    for (int i = length - 1; i >= 0; --i) {
      for (int j = 0; j < width; ++j) {
        c = inc(Componenta.comp, i, j, 1);

        if (c > n) {
          components.emplace_back(Componenta);
          Componenta = Component(length, width);
        }

        else if (c > 0 && c < n) {
          Componenta = Component(length, width);
        }

        count = 0;
      }
    }

    for (int i = 0; i < (int)components.size(); i++) {
      Field compole(length, width);
      compole.matrix = components[i].comp;
      compole.bmp_write("Comp" + std::to_string(i + 1) + ".bmp", 1);
    }
  }

  int inc(std::vector<std::vector<double>> &component, int x, int y,
          int k) { // Сама функция wave
    if (x < 1 || y < 1 || x > (int)component.size() - 1 ||
        y > (int)component[0].size() - 1 || (int)CopyField[x][y] == 255)
      return -1;

    else {
      CopyField[x][y] = 255;
      // count = count < k + 1 ? k + 1 : count;
      count++;
      component[x][y] = 0;
      inc(component, x + 1, y, k + 1);
      inc(component, x - 1, y, k + 1);
      inc(component, x, y + 1, k + 1);
      inc(component, x, y - 1, k + 1);
      inc(component, x - 1, y - 1, k + 1);
      inc(component, x + 1, y + 1, k + 1);
      inc(component, x + 1, y - 1, k + 1);
      inc(component, x - 1, y + 1, k + 1);
    }
    return count;
  }

  std::vector<Color> generateColor(
      int k) { //Генерируем свой цвет для каждого получившегося класстера
    std::vector<Color> colors;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, 255);

    for (int i = 0; i < k; i++) {
      colors.push_back({static_cast<uint8_t>(dis(gen)),
                        static_cast<uint8_t>(dis(gen)),
                        static_cast<uint8_t>(dis(gen))});
    }

    return colors;
  }

  void kMeans(int k, int p) { // Алгоритм kMeands

    //Первый этап: стандартный kMeans для поиска кластеров
    std::vector<Point> points;

    for (const auto &comp : components) {
      for (int i = 0; i < length; i++) {
        for (int j = 0; j < width; j++) {
          if (comp.comp[i][j] == 0) {
            points.push_back({i, j}); // Вытаскиваем точки из разреза
          }
        }
      }
    }

    if (points.empty()) { // Проверка, смогли ли что-то вытащить
      std::cerr << "No points available for clustering!" << std::endl;
      return;
    }

    std::vector<Point> centroids;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis_x(0, length - 1);
    std::uniform_int_distribution<> dis_y(0, width - 1);

    for (int i = 0; i < k; i++) {
      centroids.push_back(
          {dis_x(gen), dis_y(gen)}); // Рандомно выбираем k центроид
    }

    bool changed = true;
    std::vector<int> labels(points.size(), -1);

    while (changed) {
      changed = false;

      //Шаг 1: назначаем каждую точку ближайшему центроиду
      for (size_t i = 0; i < points.size(); i++) {
        double minDist = std::numeric_limits<double>::max();
        int cluster = -1;

        for (int j = 0; j < k;
             j++) { //Считаем расстояние от точки до каждого класстера, выбираем
                    //среди них минимальное
          double dist = std::pow(points[i].x - centroids[j].x, 2) +
                        std::pow(points[i].y - centroids[j].y, 2);
          if (dist < minDist) {
            minDist = dist;
            cluster = j;
          }
        }

        if (labels[i] !=
            cluster) { //Проверяем, поменяла хоть одна точка класстер
          labels[i] = cluster; // Если поменяла, то меняем ее класстер
          changed = true;
        }
      }

      //Шаг 2: пересчитываем центроиды
      std::vector<Point> newCentroids(k, {0, 0});
      std::vector<int> counts(k, 0);

      for (size_t i = 0; i < points.size();
           i++) { // Для каждой центроиды складываем все координаты точек,
                  // вычисляем количество точек в них
        newCentroids[labels[i]].x += points[i].x;
        newCentroids[labels[i]].y += points[i].y;
        counts[labels[i]]++;
      }

      for (int j = 0; j < k; j++) {
        if (counts[j] > 0) {
          newCentroids[j].x /= counts[j];
          newCentroids[j].y /= counts[j];
        }
      }

      centroids = newCentroids;
    }

    std::vector<Color> clusterColors =
        generateColor(k); //Генерируем цвета для кластеров
    std::vector<std::vector<Color>> clusterImage(
        length, std::vector<Color>(
                    width, {255, 255,
                            255})); //Заополняем матрицу цветов каждого пикселя
    std::vector<std::vector<Point>> clusteredPoints(k);

    for (size_t i = 0; i < points.size(); i++) { //Каждой точке -- свой цвет
      int cluster = labels[i];
      clusterImage[points[i].x][points[i].y] = clusterColors[cluster];
      clusteredPoints[cluster].push_back(
          points[i]); //Добвавляем каждую точку в свой кластер
    }

    // Второй этап: запуск kMeans для каждого из k кластеров с числом центроид p
    std::vector<std::vector<Point>> subCentroids(k);

    for (int clusterIdx = 0; clusterIdx < k;
         clusterIdx++) { //Для каждого кластера выполняем алгоритм kMeans
      if (clusteredPoints[clusterIdx].empty())
        continue;

      std::vector<Point> &clusterPoints = clusteredPoints[clusterIdx];
      std::vector<Point> subCenters(p);

      //Инициализируем p случайных центроид
      for (int i = 0; i < p; i++) {
        subCenters[i] = clusterPoints[std::uniform_int_distribution<>(
            0, clusterPoints.size() - 1)(gen)];
      }

      bool subChanged = true;
      std::vector<int> subLabels(clusterPoints.size(), -1);

      while (subChanged) {
        subChanged = false;

        //Назначение точек
        for (size_t i = 0; i < clusterPoints.size(); i++) {
          double minDist = std::numeric_limits<double>::max();
          int subCluster = -1;

          for (int j = 0; j < p; j++) {
            double dist = std::pow(clusterPoints[i].x - subCenters[j].x, 2) +
                          std::pow(clusterPoints[i].y - subCenters[j].y, 2);
            if (dist < minDist) {
              minDist = dist;
              subCluster = j;
            }
          }

          if (subLabels[i] != subCluster) {
            subLabels[i] = subCluster;
            subChanged = true;
          }
        }

        //Пересчет субцентроид
        std::vector<Point> newSubCenters(p, {0, 0});
        std::vector<int> subCounts(p, 0);

        for (size_t i = 0; i < clusterPoints.size(); i++) {
          newSubCenters[subLabels[i]].x += clusterPoints[i].x;
          newSubCenters[subLabels[i]].y += clusterPoints[i].y;
          subCounts[subLabels[i]]++;
        }

        for (int j = 0; j < p; j++) {
          if (subCounts[j] != 0) {
            newSubCenters[j].x /= subCounts[j];
            newSubCenters[j].y /= subCounts[j];
          }
        }

        subCenters = newSubCenters;
      }
      subCentroids[clusterIdx] = subCenters;

      //Визуализация субцентроид
      for (const auto &subCenter : subCenters) {
        clusterImage[subCenter.x][subCenter.y] = {0, 0,
                                                  0}; //Цвет для центров тяжести
      }
    }

    Field pole(length, width);
    pole.colors = clusterImage;
    pole.bmp_write("clusters.bmp", 2);
  }

  // void kMeans(int k, int p){
  //     if(p <= 0){
  //         std::cerr << "Parametr p must be greater than 0!" << std::endl;
  //         return;
  //     }

  //     std::vector<Point> points;

  //     for(const auto& comp : components){
  //         for(int i = 0; i < length; i++){
  //             for(int j = 0; j < width; j++){
  //                 if(comp.comp[i][j] == 0){
  //                     points.push_back({i, j}); //Извлекаем точки из разреза
  //                 }
  //             }
  //         }
  //     }

  //     if(points.empty()){
  //         std::cerr << "No points avalable for clustering!" << std:: endl;
  //         return;
  //     }

  //     //Инициализируем p центроид для каждого класса
  //     std::vector<std::vector<Point>> cores (k, std::vector<Point>(p));
  //     std::random_device rd; //Высокорандомное число
  //     std::mt19937 gen(rd()); //Псевдорандом на основе истинно рандомного
  //     первого числа std::uniform_int_distribution<> dis (0, points.size() -
  //     1);

  //     std::set<int> usedIndexes; //Множество для проверки, выбран ли уже
  //     данный индекс

  //     for(int i = 0; i < k; i++){
  //         for(int j = 0; j < p; j++){
  //             int idx;
  //             //Пытаемся выбрать уникальный индекс

  //             do{
  //                 idx = dis(gen);
  //             } while(usedIndexes.find(idx) != usedIndexes.end()); //
  //             Проверка на уникальность выбранного индекса

  //             usedIndexes.insert(idx); //Запоминаем индекс как использованый
  //             cores[i][j] = points[idx]; //Присваиваем центроид
  //         }
  //     }

  //     bool changed = true;
  //     std::vector<int> labels(points.size(), -1); //Метки для каждой точки (к
  //     каким кластерам они относятся)

  //     while(changed){
  //         changed = false;

  //         //Шаг 1: Назначаем каждую точку одному из k классов, учитывая p
  //         центроид for(size_t i = 0; i < points.size(); i++){
  //             double minDist = std::numeric_limits<double>::max();
  //             int cluster = -1;

  //             for(int j = 0; j < k; j++){
  //                 for(int l = 0; l < p; l++){
  //                     double dist = std::pow(points[i].x - cores[j][l].x, 2)
  //                     +
  //                                   std::pow(points[i].y - cores[j][l].y, 2);
  //                     if(dist < minDist){
  //                         minDist = dist;
  //                         cluster = j;
  //                     }
  //                 }
  //             }

  //             if(labels[i] != cluster){ //Проверка, изменила ли точка свой
  //             кластер
  //                 labels[i] = cluster;
  //                 changed = true;
  //             }
  //         }

  //         //Шаг 2: Пересчет центрид
  //         for(int j = 0; j < k; j++){
  //             std::vector<Point> clusterPoints; // Точки, которые принадлежат
  //             данному класстеру

  //             for(size_t i = 0; i < points.size(); i++){ // Если точка
  //             принадлежит кластеру j, то добавляем ее в clusterPoints
  //                 if(labels[i] == j){
  //                     clusterPoints.push_back(points[i]);
  //                 }
  //             }

  //             if(clusterPoints.empty()){
  //                 std::cerr << "Cluster" << j << "has no points" <<
  //                 std::endl; continue;
  //             }

  //             //Далее применяем kMeans внутри класса для p центроид ядра
  //             std::vector<Point> newCentroids(p, {0, 0}); //Центроиды класса
  //             (p штук) std::vector<int> counts(p, 0); //Количество точек в
  //             каждом из p кластеров

  //             for(const auto& pt : clusterPoints){ //Для каждой точки из
  //             кластера...
  //                 double minDist = std::numeric_limits<double>::max();
  //                 int centroidIndex = -1;

  //                 //Находим ближайшую центроиду внутри класса
  //                 for(int l = 0; l < p; l++){
  //                     double dist = std::pow(pt.x - cores[j][l].x, 2) +
  //                                   std::pow(pt.y - cores[j][l].y, 2); //
  //                                   Считаем минимальное расстояние до каждой
  //                                   центроиды класстера
  //                     if(dist < minDist){
  //                         minDist = dist;
  //                         centroidIndex = l;
  //                     }
  //                 }

  //                 //Добавляем к центроиде координаты точки, которая ей
  //                 принадлежит newCentroids[centroidIndex].x += pt.x;
  //                 newCentroids[centroidIndex].y += pt.y;
  //                 counts[centroidIndex]++;
  //             }

  //             //Считаем центры тяжести кластеров
  //             for(int l = 0; l < p; l++){
  //                 if(counts[l] != 0){ // Если в данном кластере есть хотя бы
  //                 одна точка
  //                     newCentroids[l].x /= counts[l];
  //                     newCentroids[l].y /= counts[l];
  //                 }
  //             }
  //             cores[j] = newCentroids; // Обновляем ядро класса
  //         }
  //     }

  //     //Генерация цветов для кластеров
  //     std::vector<Color> clusterColors = generateColor(k);
  //     std::vector<std::vector<Color>> clusterImage(length,
  //     std::vector<Color>(width, {255, 255, 255})); Color brightRed = {0, 0,
  //     255}; //Центроиды будут подсвечены ярко-красным

  //     for(size_t i = 0; i < points.size(); i++){ //Каждой точке свой цвет в
  //     зависимости от кластера
  //         int cluster = labels[i];
  //         clusterImage[points[i].x][points[i].y] = clusterColors[cluster];
  //     }

  //     //Подсвечиваем центроиды
  //     for(int j = 0; j < k; j ++){
  //         for(int l = 0; l < p; l++){
  //             int x = cores[j][l].x;
  //             int y = cores[j][l].y;

  //             if(x >= 0 && x < length && y>= 0 && y < width){
  //                 clusterImage[x][y] = brightRed;
  //             }
  //         }
  //     }

  //     Field pole(length, width);
  //     pole.colors = clusterImage;
  //     pole.bmp_write("cluster2.bmp", 2); //Сохраняем результат
  // }
};



class Control { // Класс Control выполняет роль диспетчера
  std::vector<Gauss> gausses;
  Field *field = nullptr;
  Srez srez;

public:
  Server server;

  void wave_cntrl(int n) { // Алгоритм wave

    if (!field) {
      std::cerr << "Error! Field not initialized!" << std::endl;
      server.log("Control received 'bin', but field was not initialized.");
      return;
    }

    srez.wave(n);
    server.log("Control completed the wave program");
  }

  void init(int length, int width) { // Инициализация поля
    if (field) {
      server.log("Control received 'init', but field was already initialized.");
      return;
    }
    field = new Field(length, width); // new возвращает указатель
    server.log("Field created with size " + std::to_string(length) + "x" +
               std::to_string(width));
    server.log("Control passed 'init' command to Field");
  }

  void bin_cntrl(int h) { // Разрез
    if (!field) {
      std::cerr << "Error! Field not initialized!" << std::endl;
      server.log("Control received 'bin', but field was not initialized.");
      return;
    }
    srez.bin(h, *field);
    server.log("A section on height " + std::to_string(h) + " is obtained");
  }

  void kMeans_cntrl(int k, int p) { // Алгоритм kMeans
    if (!field) {
      std::cerr << "Error! Field not initialized!" << std::endl;
      server.log("Control received 'KMEANS', but field was not initialized.");
      return;
    }
    server.log("Algotitgm KMEANS for the k =" + std::to_string(k) +
               " is launched");
    srez.kMeans(k, p);
    server.log("Algorithm KMEANS for the k =" + std::to_string(k) +
               " is completed");
  }

  void add_gauss(double h, double x0, double y0, double sigma_x,
                 double sigma_y) { // Добавление гаусса в список

    if (!field) {
      std::cerr << "Error! Field not initialized!" << std::endl;
      server.log("Control received 'g', but field was not initialized.");
      return;
    }

    gausses.emplace_back(h, x0, y0, sigma_x, sigma_y);

    // Логирование факта передачи команды Gauss
    server.log(
        "Control passed Gauss command to Field with (h=" + std::to_string(h) +
        ", x0=" + std::to_string(x0) + ", y0=" + std::to_string(y0) +
        ", sigma_x=" + std::to_string(sigma_x) +
        ", sigma_y=" + std::to_string(sigma_y) + ")");
  }

  void generate() { //Функция, которая добавляет список гауссов
    if (!field) {
      std::cerr << "Error! Field not initialized!" << std::endl;
      server.log("Control received 'generate', but field was not initialized.");
      return;
    }

    if (gausses.empty()) {
      server.log("No gausses to apply in 'generate'");
      return;
    }

    server.log("Control passed 'generate' command to Field");

    for (size_t i = 0; i < gausses.size(); ++i) {
      server.log("Control is applying Gauss #" + std::to_string(i + 1) +
                 " to Field");
      field->apply_gauss(gausses[i]);
      server.log("Gauss #" + std::to_string(i + 1) + "applied");
    }

    field->normalize();
    gausses.clear();
    server.log("Gauss completed applying all Gausses");
  }

  void gnuplot(const std::string &filename) { // Гнуплот
    if (!field) {
      std::cerr << "Error! Field not initialized!" << std::endl;
      server.log("Control received 'gnuplot', but field was not initialized.");
      return;
    }

    field->save_to_gnuplot("gnuplot_" + filename + ".txt");
    server.log("Field saved to Gnuplot file: " + filename);

    std::ofstream gp_file("gnuplot_comm.txt");
    gp_file << "set view 60,30\n"; // угол в градусах по вертикальной и по
                                   // горизонтальной оси
    gp_file << "set palette defined (0 \"blue\", 1 \"red\")\n";
    gp_file << "set pm3d at s\n";
    gp_file << "splot 'gnuplot_" << filename << ".txt' with lines\n";
    gp_file << "pause -1";
    gp_file.close();
    server.log("Field generated Gnuplot file: " + filename);
  }

  void bmp_write_cntrl(const std::string &filename, int k) { //Запись BMP
    if (!field) {
      std::cerr << "Error! Field not initialized!" << std::endl;
      server.log("Control received 'save bmp', but field was not initialized.");
      return;
    }
    field->bmp_write(filename, k);
    server.log("Field saved to BMP file:" + filename);
  }

  void bmp_read_cntrl(const std::string &filename) { // Чтение BMP
    if (!field) {
      std::cerr << "Error! Field not initialized!" << std::endl;
      server.log("Control received 'read bmp', but field was not initialized.");
      return;
    }

    std::pair<int, int> dimensions = field->bmp_read(filename);
    int new_length = dimensions.first;
    int new_width = dimensions.second;

    if (new_length == 0 || new_width == 0) {
      std::cerr << "Error reading BMP file. No changes made" << std::endl;
      server.log("Error reading BMP file");
    }

    delete field;
    field = new Field(new_length, new_width);

    std::ifstream bmp_file(filename, std::ios::binary);
    bmp_file.ignore(54);
    field->load_data(bmp_file, new_length, new_width);
    bmp_file.close();

    server.log("Field loaded bmp file:" + filename);
  }
};



class Interface { // Клиент (интерфейс)
  Control control;
  Client client;

public:
  std::string trim(const std::string &str) { // Функция, чтобы убирать пробелы в
                                             // начале и конце строки
    size_t first = str.find_first_not_of(' ');
    if (first == std::string::npos) {
      return "";
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
  }

  int run() { // функция, которая непосредственно и является интерфейсом
    std::string filename;
    std::cout << "Write name of the config" << std::endl;
    std::cin >> filename;
    std::cin.ignore(
        std::numeric_limits<std::streamsize>::max(),
        '\n'); // Очещаем символ новой строки для корректной работы Batch-pause
    std::ifstream config(filename);
    bool gaussi = false; // Индикатор, нужен для гаусса по умолчанию в конфиге
    bool batch_pause = false; // Индикатор для Batch_pause
    int n = 0; // Количество пикселей, которые мы считаем шумом

    if (!config.is_open()) {
      std::cerr << "Failed to open file:" << filename << std::endl;
      return 1;
    }

    std::cout << "Succesfully opened file:" << filename << std::endl;
    std::string line;
    bool noize = false;

    while (!noize && std::getline(config, line)) {
      std::istringstream iss(line);
      std::string key, value;

      if (std::getline(iss, key, '=') && std::getline(iss, value)) {
        key = trim(key);
        value = trim(value);

        if (key == "Noize") {
          std::istringstream ss(value);
          ss >> n;
          noize = true;
        }
      }
    }

    config.clear(); //Сбрасываем флаги ошибок
    config.seekg(0, std::ios::beg); // На начало конфига

    while (std::getline(config, line)) {
      std::istringstream iss(line);
      std::string key, value;

      if (std::getline(iss, key, '=') && std::getline(iss, value)) {
        key = trim(key);
        value = trim(value);

        if (key == "Log_Server") {
          if (value == "ON") {
            std::string log_filename = "server_log.txt";
            std::getline(config, line);
            std::istringstream iss_filename(line);
            if (std::getline(iss_filename, key, '=') &&
                std::getline(iss_filename, value)) {
              key = trim(key);
              value = trim(value);
              if (key == "Log_Server_filename") {
                log_filename = value;
              }
            }
            control.server.enable_logging(log_filename);
            std::cout << "Server client enabled with file: " << log_filename
                      << std::endl;
          }
        }

        else if (key == "Log_Client") {

          if (value == "ON") {
            std::string log_filename = "client_log.txt";
            std::getline(config, line);
            std::istringstream iss_filename(line);

            if (std::getline(iss_filename, key, '=') &&
                std::getline(iss_filename, value)) {
              key = trim(key);
              value = trim(value);
              if (key == "Log_Client_filename") {
                log_filename = value;
              }
            }

            client.enable_logging(log_filename);
            std::cout << "Server logging enabled with file: " << log_filename
                      << std::endl;

            if (!noize) {
              client.log("The noise value was not found, so we assumed that "
                         "there is no noise");
            }

            else {
              client.log("The noise value is set to" + std::to_string(n));
            }
          }
        }

        else if (key == "Batch_pause") {
          if (value == "ON") {
            batch_pause = true;
          } else {
            batch_pause = false;
          }
        }

        else if (key == "Batch_file") {
          if (value == "ON") {

            std::string Batch_filename = "Batch_file.txt";
            std::getline(config, line);
            std::istringstream iss_filename(line);

            if (std::getline(iss_filename, key, '=') &&
                std::getline(iss_filename, value)) {
              key = trim(key);
              value = trim(value);
              if (key == "Batch_filename") {
                Batch_filename = value;
              }
            }

            std::ifstream batch(Batch_filename);

            if (!batch.is_open()) {
              std::cerr << "Failed to open file:" << Batch_filename
                        << std::endl;
              return 1;
            }

            std::string batch_line, batch_key, batch_value;
            std::getline(batch, batch_line);
            std::istringstream iss_batch(batch_line);
            iss_batch >> batch_key;

            if (batch_key != "INI") {
              std::getline(config, line);
              iss.clear();
              iss.str(line);
              std::getline(iss, key, '=');
              std::getline(iss, value);
              key = trim(key);
              value = trim(value);

              if (key == "Field_DEF") {
                int length, width;
                std::istringstream iss_Field(value);
                iss_Field >> length >> width;
                client.log("The interface has set the default value of the "
                           "field: length=" +
                           std::to_string(length) +
                           ", width=" + std::to_string(width));
                control.init(length, width);
              }

              else {
                client.log("Field initialization error");
                return -1;
              }

              batch.clear(); // Сбрасываем флаги ошибок (например, EOF)
              batch.seekg(0, std::ios::beg);
            }

            else {
              batch.clear(); // Сбрасываем флаги ошибок (например, EOF)
              batch.seekg(0, std::ios::beg);
            }

            while (std::getline(batch, batch_line)) {

              if (batch_pause == true) {
                std::cout << "Нажмите Enter, чтобы продолжить.." << std::endl;
                std::cin.get();
              }
              std::istringstream iss_batch(batch_line);

              if (iss_batch >> batch_key) {
                std::getline(iss_batch, batch_value);
                std::istringstream iss_info(batch_value);

                if (batch_key == "INI") {
                  int length, width;
                  iss_info >> length >> width;
                  if (batch_pause == true) {
                    std::cout << "Interface received 'INI' command with "
                                 "parameters: length=" +
                                     std::to_string(length) +
                                     ", width=" + std::to_string(width)
                              << std::endl;
                  }
                  client.log("Interface received 'INI' command with "
                             "parameters: length=" +
                             std::to_string(length) +
                             ", width=" + std::to_string(width));
                  control.init(length, width);

                } else if (batch_key == "G") {
                  double h, x0, y0, sigma_x, sigma_y;
                  iss_info >> h >> x0 >> y0 >> sigma_x >> sigma_y;
                  if (batch_pause == true) {
                    std::cout << "Interface received 'g' command with "
                                 "parameters: h=" +
                                     std::to_string(h) +
                                     ", x0=" + std::to_string(x0) +
                                     ", y0=" + std::to_string(y0) +
                                     ", sigma_x=" + std::to_string(sigma_x) +
                                     ", sigma_y=" + std::to_string(sigma_y)
                              << std::endl;
                  }
                  client.log(
                      "Interface received 'g' command with parameters: h=" +
                      std::to_string(h) + ", x0=" + std::to_string(x0) +
                      ", y0=" + std::to_string(y0) +
                      ", sigma_x=" + std::to_string(sigma_x) +
                      ", sigma_y=" + std::to_string(sigma_y));
                  control.add_gauss(h, x0, y0, sigma_x, sigma_y);
                  gaussi = true;

                } else if (batch_key == "GEN") {
                  if (gaussi == false) {
                    std::streampos saved_position = config.tellg();
                    while (std::getline(config, line)) {
                      iss.str(line);
                      iss.clear();
                      if (std::getline(iss, key, '=') &&
                          std::getline(iss, value)) {
                        key = trim(key);
                        value = trim(value);
                        if (key == "Gauss_DEF") {
                          double h, x0, y0, sigma_x, sigma_y;
                          iss.clear();
                          iss.str(value);
                          iss >> h >> x0 >> y0 >> sigma_x >> sigma_y;
                          if (batch_pause == true) {
                            std::cout
                                << "The default Gauss value is set: h=" +
                                       std::to_string(h) +
                                       ", x0=" + std::to_string(x0) +
                                       ", y0=" + std::to_string(y0) +
                                       ", sigma_x=" + std::to_string(sigma_x) +
                                       ", sigma_y=" + std::to_string(sigma_y)
                                << std::endl;
                          }
                          client.log("The default Gauss value is set: h=" +
                                     std::to_string(h) +
                                     ", x0=" + std::to_string(x0) +
                                     ", y0=" + std::to_string(y0) +
                                     ", sigma_x=" + std::to_string(sigma_x) +
                                     ", sigma_y=" + std::to_string(sigma_y));
                          control.add_gauss(h, x0, y0, sigma_x, sigma_y);
                          gaussi = true;
                        }
                      }
                    }
                    config.seekg(saved_position);
                  }

                  if (batch_pause == true) {
                    std::cout << "Interface received 'generate' command."
                              << std::endl;
                  }
                  client.log("Interface received 'generate' command.");
                  control.generate();

                } else if (batch_key == "GNU") {
                  std::string filename;
                  iss_info >> filename;

                  if (batch_pause == true) {
                    std::cout << "Interface received 'gnuplot' command with "
                                 "filename: " +
                                     filename
                              << std::endl;
                  }

                  client.log(
                      "Interface received 'gnuplot' command with filename: " +
                      filename);
                  control.gnuplot(filename);

                } else if (batch_key == "BMP") {
                  std::string filename;
                  iss_info >> filename;
                  if (batch_pause == true) {
                    std::cout << "Interface received 'write bmp' command with "
                                 "filename:" +
                                     filename
                              << std::endl;
                  }
                  client.log(
                      "Interface received 'write bmp' command with filename:" +
                      filename);
                  control.bmp_write_cntrl(filename, 1);

                } else if (batch_key == "RBMP") {
                  std::string filename;
                  iss_info >> filename;
                  if (batch_pause == true) {
                    std::cout << "Interface received 'read bmp' command with "
                                 "filename: " +
                                     filename
                              << std::endl;
                  }
                  client.log(
                      "Interface received 'read bmp' command with filename: " +
                      filename);
                  control.bmp_read_cntrl(filename);
                  gaussi = true;

                } else if (batch_key == "BIN") {
                  int h;
                  iss_info >> h;

                  if (batch_pause == true) {
                    std::cout << "Interface received 'BIN' command with h = " +
                                     std::to_string(h)
                              << std::endl;
                  }

                  client.log("Interface received 'BIN' command with h = " +
                             std::to_string(h));
                  control.bin_cntrl(h);

                } else if (batch_key == "WAVE") {

                  if (batch_pause == true) {
                    std::cout << "Interface received 'WAVE' command with pixel "
                                 "noise:" +
                                     std::to_string(n)
                              << std::endl;
                  }

                  client.log(
                      "Interface received 'WAVE' command with pixel noise:" +
                      std::to_string(n));
                  control.wave_cntrl(n);
                } else if (batch_key == "KMEANS") {
                  int k, p;
                  iss_info >> k >> p;

                  if (batch_pause == true) {
                    std::cout
                        << "Interface received 'KMEANDS' coommand with k =" +
                               std::to_string(k)
                        << std::endl;
                  }

                  client.log("Interface received 'KMEANDS' coommand with k =" +
                             std::to_string(k));
                  control.kMeans_cntrl(k, p);
                }
              }
            }
          }
        }
      }
    }

    // else if (input_type == "2") {
    //     std::string command;
    //     while (true) {
    //         std::cin >> command;

    //         if (command == "init") {
    //             int length, width;
    //             std::cin >> length >> width;
    //             client.log("Interface received 'init' command with
    //             parameters: length=" + std::to_string(length) + ", width=" +
    //             std::to_string(width)); control.init(length, width);

    //         } else if (command == "g") {
    //             double h, x0, y0, sigma_x, sigma_y;
    //             std::cin >> h >> x0 >> y0 >> sigma_x >> sigma_y;
    //             client.log("Interface received 'g' command with parameters:
    //             h=" + std::to_string(h) + ", x0=" + std::to_string(x0) + ",
    //             y0=" + std::to_string(y0) +
    //                 ", sigma_x=" + std::to_string(sigma_x) + ", sigma_y=" +
    //                 std::to_string(sigma_y));
    //             control.add_gauss(h, x0, y0, sigma_x, sigma_y);

    //         } else if (command == "generate") {
    //             client.log("Interface received 'generate' command.");
    //             control.generate();

    //         } else if (command == "gnuplot") {
    //             std::string filename;
    //             std::cin >> filename;
    //             client.log("Interface received 'gnuplot' command with
    //             filename: " + filename); control.gnuplot(filename);

    //         } else if (command == "save") {
    //             std::string filetype;
    //             std::cin >> filetype;
    //             if (filetype == "bmp") {
    //                 client.log("Interface received 'save bmp' command.");
    //                 control.bmp_write_cntrl();
    //             }
    //         } else if (command == "exit") {
    //             break;
    //         }
    //     }
    // }
    return 0;
  }
};


int main() {

  Interface interface;
  interface.run();
  return 0;
}
