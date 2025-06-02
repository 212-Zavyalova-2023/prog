//#include "config.txt"
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

//config file read it and what is it
class Config {
public:
  int fieldWidth;
  int fieldHeight;
  double defaultX, defaultY, defaultSx, defaultSy, defaultH;
  std::string logFileNameForControl;
  std::string logFileNameForInterface;
  bool loggingInterfaceEnabled;
  bool loggingControlEnabled;

  Config(const std::string &filename) {
    std::ifstream configFile(filename);
    if (!configFile.is_open()) {
      std::cerr << "Failed to open config file." << std::endl;
      return;
    }

    std::string key;
    while (configFile >> key) {
      if (key == "fieldWidth")
        configFile >> fieldWidth;
      else if (key == "fieldHeight")
        configFile >> fieldHeight;
      else if (key == "defaultX")
        configFile >> defaultX;
      else if (key == "defaultY")
        configFile >> defaultY;
      else if (key == "defaultSx")
        configFile >> defaultSx;
      else if (key == "defaultSy")
        configFile >> defaultSy;
      else if (key == "defaultH")
        configFile >> defaultH;
      else if (key == "logFileNameForControl")
        configFile >> logFileNameForControl;
      else if (key == "logFileNameForInterface")
        configFile >> logFileNameForInterface;
      else if (key == "loggingInterfaceEnabled")
        configFile >> std::boolalpha >> loggingInterfaceEnabled;
      else if (key == "loggingControlEnabled")
        configFile >> std::boolalpha >> loggingControlEnabled;
    }

    configFile.close();
  }
};
//logging messages
class Logger {
private:
  std::ofstream logFile;

public:
  Logger(const std::string &fileName) {
    logFile.open(fileName, std::ios::out | std::ios::trunc);
    if (!logFile.is_open()) {
      // logFile.open(fileName, std::ios::out | std::ios::app);
      std::cerr << "File opening error: " << fileName << std::endl;
    }
  }


  Logger(Logger &&other) noexcept : logFile(std::move(other.logFile)) {}


  Logger &operator=(Logger &&other) noexcept {
    if (this != &other) {
      if (logFile.is_open()) {
        logFile.close();
      }
      logFile = std::move(other.logFile);
    }
    return *this;
  }

  ~Logger() {
    if (logFile.is_open()) {
      logFile.close();
    }
  }

  void logMessage(const std::string &message, bool b) {
    if (logFile.is_open() && b) {
      auto now = std::chrono::system_clock::now();
      std::time_t now_c = std::chrono::system_clock::to_time_t(now);
      std::stringstream timeStamp;
      timeStamp << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S");
      logFile << "[" << timeStamp.str() << "] " << message << std::endl;
    }
  }
};

//hills
class Gaus {
public:
  double h, x0, y0, sx, sy;

  Gaus(double h, double x0, double y0, double sx, double sy)
      : h(h), x0(x0), y0(y0), sx(sx), sy(sy) {}
};

class Pole {
public:
  std::vector<std::vector<double>> field;

  Pole(int A, int B) { field.resize(A, std::vector<double>(B, 0)); }
};

class Component {
public:
  std::vector<std::vector<double>> componenta;
  Component(const std::vector<std::vector<double>> &inputComponenta)
      : componenta(inputComponenta) {}

  Component(int A, int B) { componenta.resize(A, std::vector<double>(B, 0)); }
};

class Control {
private:
  std::vector<Gaus> gaussi;
  Config config;
  Logger logger;
  bool b = true;
  std::vector<std::vector<double>> CopyPole;
  int count = 0;

  int incrementAndCollect(std::vector<std::vector<double>> &componenta, int x,
                          int y, int i) {

    if (x < 1 || y < 1 || x > (int)componenta[0].size() - 2 ||
        y > (int)componenta.size() - 2 || CopyPole[y][x] < 250)
      return -1;

    if (CopyPole[y][x] >= 255 && CopyPole[y][x] <= 255) {
      CopyPole[y][x] = 0; // � � �
      count = count < i + 1 ? i + 1 : count;
      componenta[y][x] = 255; // � � � Componenta
      incrementAndCollect(componenta, x + 1, y, i + 1);
      incrementAndCollect(componenta, x - 1, y, i + 1);
      incrementAndCollect(componenta, x, y + 1, i + 1);
      incrementAndCollect(componenta, x, y - 1, i + 1);
    }

    return count;
  }

public:
  std::vector<Component> componenti;
  Pole *p; // Pointer to Pole
  Control() : config("config.txt"), logger("control_log.txt") {}

  Control(const std::string &configFileName)
      : config(configFileName), logger(""), p(nullptr) {
    if (!config.logFileNameForControl.empty()) {
      logger = Logger(config.logFileNameForControl);
    }

    if (config.loggingControlEnabled) {

      logger.logMessage("Logging Control is enabled.", b);
      std::cout << "Logging Control is enabled." << std::endl;
    } else {
      logger.logMessage("Logging Control is disabled.", b);
      std::cout << "Logging Control is disabled." << std::endl;
      b = false;
    }
  }

  ~Control() {
    delete p; // Free memory
  }

  void addgauss(double h, double x0, double y0, double sigma_x,
                double sigma_y) {
    logger.logMessage("Gauss addition will begin", b);
    gaussi.emplace_back(h, x0, y0, sigma_x, sigma_y);
    logger.logMessage("Added gauss", b);
  }

  void init(int A, int B) {
    delete p; // Delete previous pointer
    logger.logMessage("Field generation will begin", b);
    p = new Pole(A, B);
    logger.logMessage("Added field", b);
    CopyPole = p->field; //� �
  }

  void generate() {
    double value;
    for (const auto &g : gaussi) {
      for (long unsigned int x = 0; x < p->field[0].size(); ++x) {
        for (long unsigned int y = 0; y < p->field.size(); ++y) {
          value =
              g.h *
              exp(-((pow((x - g.x0) / g.sx, 2)) + (pow((y - g.y0) / g.sy, 2))) /
                  2);
          p->field[y][x] += value;
        }
      }
    }
  }

  void gnuplot() {
    int rows = p->field.size();
    int cols = p->field[0].size();
    // Open a pipe to gnuplot
    FILE *gnuplotPipe = popen("gnuplot -p", "w");
    if (!gnuplotPipe) {
      std::cerr << "Could not open pipe to gnuplot." << std::endl;
      logger.logMessage("Could not open pipe to gnuplot.", b);
      return;
    }
    // Send gnuplot commands
    fprintf(gnuplotPipe, "set contour base\n");
    fprintf(gnuplotPipe, "set view 60,30\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", cols - 1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", rows - 1);
    fprintf(gnuplotPipe, "set terminal png\n");
    fprintf(gnuplotPipe, "set output 'landscape.png'\n");
    fprintf(gnuplotPipe, "splot '-' with lines\n");
    // Write data directly to gnuplot
    for (int y = 0; y < rows; ++y) {
      for (int x = 0; x < cols; ++x) {
        fprintf(gnuplotPipe, "%d %d %f\n", x, y, p->field[y][x]);
      }
      fprintf(gnuplotPipe, "\n"); // Newline to separate rows
    }
    fprintf(gnuplotPipe, "EOF\n"); // Close the pipe
    pclose(gnuplotPipe);
  }
  /*void gnuplot() {
    int m = 0, k = 0; //Позиция
    gnuplot.open("Gnuplot.txt");
    for (double i = 0.0; i <= sizex; i += 0.3) {
      for (double j = 0.0; j <= sizey; j += 0.3) {
        gnuplot << i << " " << j << " " << Points[m][k] << "\n";
        k += 1;
      }
      gnuplot << "\n";
      k = 0;
      m += 1;
    }
  }*/

  void bmp_write(const std::vector<std::vector<double>> &pixelMatrix,
                 const std::string &filename) {
    int width = pixelMatrix[0].size();
    int height = pixelMatrix.size();
    int padding = (4 - (width * 3) % 4) % 4; // Padding for alignment to 4 bytes
    std::ofstream bmpFile(filename, std::ios::binary);
    if (!bmpFile) {
      std::cerr << "Failed to create BMP file." << std::endl;
      logger.logMessage("Failed to create BMP file.", b);
      return;
    }
    // Write BMP header
    unsigned char bmpHeader[54] = {
        'B',  'M',        // Identifier
        0,    0,    0, 0, // Size of file (will be set later)
        0,    0,    0, 0, // Reserved
        54,   0,    0, 0, // Header size
        40,   0,    0, 0, // Info header size
        0,    0,    0, 0, // Width (will be set later)
        0,    0,    0, 0, // Height (will be set later)
        1,    0,          // Number of color planes
        24,   0,          // Bits per pixel
        0,    0,    0, 0, // Compression
        0,    0,    0, 0, // Image size (will be set later)
        0x13, 0x0B, 0, 0, // Horizontal resolution
        0x13, 0x0B, 0, 0, // Vertical resolution
        0,    0,    0, 0, // Number of colors in palette
        0,    0,    0, 0  // Important colors
    };
    // Set width and height in header
    bmpHeader[18] = (width & 0xFF);
    bmpHeader[19] = (width >> 8) & 0xFF;
    bmpHeader[20] = (width >> 16) & 0xFF;
    bmpHeader[21] = (width >> 24) & 0xFF;
    bmpHeader[22] = (height & 0xFF);
    bmpHeader[23] = (height >> 8) & 0xFF;
    bmpHeader[24] = (height >> 16) & 0xFF;
    bmpHeader[25] = (height >> 24) & 0xFF;
    // Write header
    bmpFile.write(reinterpret_cast<char *>(bmpHeader), 54);
    // Write pixel data
    for (int y = height - 1; y >= 0; --y) { // BMP stores pixels bottom-to-top
      for (int x = 0; x < width; ++x) {
        unsigned char color =
            255 - static_cast<unsigned char>(pixelMatrix[y][x]); // Color
        bmpFile.put(color);                                      // B
        bmpFile.put(color);                                      // G
        bmpFile.put(color);                                      // R
      }
      // Add padding
      for (int p = 0; p < padding; ++p) {
        bmpFile.put(0);
      }
    }
    bmpFile.close();
  }

  void bmp_read(const std::string &filename) {
    std::ifstream bmpFile(filename, std::ios::binary);
    if (!bmpFile) {
      std::cerr << "Failed to open BMP file." << std::endl;
      logger.logMessage("Failed to open BMP file.", b);
      return;
    }

    unsigned char header[54];
    bmpFile.read(reinterpret_cast<char *>(header), 54);

    int width = header[18] | (header[19] << 8) | (header[20] << 16) |
                (header[21] << 24);
    int height = header[22] | (header[23] << 8) | (header[24] << 16) |
                 (header[25] << 24);

    init(height, width);


    for (int y = height - 1; y >= 0; --y) {
      for (int x = 0; x < width; ++x) {
        unsigned char color = bmpFile.get();
        bmpFile.get();
        bmpFile.get();
        double value = 255 - color;
        p->field[y][x] = value;
      }
      bmpFile.ignore((4 - (width * 3) % 4) % 4);
    }
    bmpFile.close();
  }

  void bin(int slise) {
    for (const auto &g : gaussi) {
      (void)g;
      for (int x = 0; x < (int)p->field[0].size(); ++x) {
        for (int y = 0; y < (int)p->field.size(); ++y) {
          CopyPole[y][x] = p->field[y][x] > slise ? 255 : 0;
        }
      }
    }

    bmp_write(CopyPole, "slise.bmp");
    logger.logMessage("Created BMP file.", b);
  }

  void wave() {

    Component Componenta(p->field.size(), p->field[0].size());

    for (int y = 0; y < (int)p->field.size(); ++y) {
      for (int x = 0; x < (int)p->field[y].size(); ++x) {
        if (CopyPole[y][x] <= 255 && CopyPole[y][x] >= 255) {
          count = 0;
          if (incrementAndCollect(Componenta.componenta, x, y, 1) > 6) {
            componenti.emplace_back(Componenta);
          }
        }
      }
    }

    logger.logMessage("Wave used, amount component = " +
                          std::to_string(componenti.size()),
                      b);
  }
};

class Interface {
private:
  bool b = true;
  Logger logger;

public:
  Config config;
  Control c;
  Interface() : logger("loginterface.txt"), config("config.txt") {}

  Interface(const std::string &configFileName)
      : logger(""), config(configFileName), c(configFileName) {
    if (!config.logFileNameForInterface.empty()) {
      logger = Logger(config.logFileNameForInterface);
    }

    if (config.loggingInterfaceEnabled) {

      logger.logMessage("Logging Interface is enabled.", b);
      std::cout << "Logging Interface is enabled." << std::endl;
    } else {
      logger.logMessage("Logging Interface is disabled.", b);
      std::cout << "Logging Interface is disabled." << std::endl;
      b = false;
    }
  }

  void print() {
    double x, y, sx, sy, h;
    std::string s;
    bool a;
    std::string filename;
    std::ifstream file;
    std::cout
        << "Hello, dear user, this program builds Gaussians.\nEnter commands "
           "from a text file (PRESS 0) or from the keyboard (PRESS 1)?\n"
        << std::endl;
    std::cin >> a;
    logger.logMessage("User chose input method: " + std::to_string(a), b);
    if (a == 0) {
      std::cout << "You will enter commands from a text file.\nEnter filename:"
                << std::endl;
      std::cin >> filename;
      logger.logMessage("Reading commands from file: " + filename, b);
      file.open(filename);
      if (!file) {
        std::cout << "File not found" << std::endl;
        logger.logMessage("Error: File not found.", b);
        return;
      }
    } else {
      std::cout << "You will enter commands from the keyboard" << std::endl;
      logger.logMessage("User chose to input commands from the keyboard.", b);
    }
    if (a == 0) {
      int n = 0;
      while (file >> s) {
        logger.logMessage("Received command: " + s, b);
        if (s == "init") {

          if (n != 0) {
            std::cout << "The init command has already been called.\nError\n";
            logger.logMessage("Error: Multiple init commands.", b);
            return;
          }
          n = 1;

          int A = config.fieldWidth;
          int B = config.fieldHeight;

          logger.logMessage("Initializing field with size: " +
                                std::to_string(A) + " x " + std::to_string(B),
                            b);
          c.init(A, B);
          logger.logMessage("Field initialized.", b);
        } else if (n != 1) {
          std::cout << "The init command was not use.\nError\n";
          logger.logMessage("Error: The init command was not use.", b);
          return;
        } else if (s == "g") {


          file >> x >> y >> sx >> sy >> h;


          if (file.fail()) {

            if (file.eof()) {
              x = config.defaultX;
              y = config.defaultY;
              sx = config.defaultSx;
              sy = config.defaultSy;
              h = config.defaultH;
            } else {

              file.clear();


              if (!(file >> x))
                x = config.defaultX;
              if (!(file >> y))
                y = config.defaultY;
              if (!(file >> sx))
                sx = config.defaultSx;
              if (!(file >> sy))
                sy = config.defaultSy;
              if (!(file >> h))
                h = config.defaultH;
            }
          }

          logger.logMessage(
              "Adding Gaussian: x=" + std::to_string(x) +
                  ", y=" + std::to_string(y) + ", sx=" + std::to_string(sx) +
                  ", sy=" + std::to_string(sy) + ", h=" + std::to_string(h),
              b);
          c.addgauss(h, x, y, sx, sy);
        } else if (s == "generate") {
          c.generate();
          logger.logMessage("Generated values in the field.", b);
        } else if (s == "gnuplot") {
          c.gnuplot();
          logger.logMessage("Called gnuplot.", b);
        } else if (s == "bmp_write") {
          c.bmp_write(c.p->field, "output.bmp");
          logger.logMessage("Created BMP file.", b);
        } else if (s == "bmp_read") {
          file >> filename;
          c.bmp_read(filename);
          logger.logMessage("Read BMP file: " + filename, b);
        } else if (s == "bin") {
          int slice;
          file >> slice;
          c.bin(slice);
          logger.logMessage("Slice applied: slice=" + std::to_string(slice), b);
          logger.logMessage("Wave will be used", b);
          c.wave();
          std::cout << "Amount component = " << c.componenti.size()
                    << std::endl;
        }
      }
    } else {
      int n = 0;
      while (true) {
        std::cout << "Enter command (init, g, generate, gnuplot, bmp_write, "
                     "bmp_read, bin, end):";
        std::cin >> s;
        std::cout << "\n";
        logger.logMessage("Received command: " + s, b);
        if (s == "init") {

          if (n != 0) {
            std::cout << "The init command has already started.\nError\n";
            logger.logMessage("Error: Multiple init commands.", b);
            return;
          }
          n = 1;

          int A = config.fieldWidth;
          int B = config.fieldHeight;

          logger.logMessage("Initializing field with size: " +
                                std::to_string(A) + " x " + std::to_string(B),
                            b);
          c.init(A, B);
          logger.logMessage("Field initialized.", b);
        }
        if (n != 1) {
          std::cout << "The init command was not use.\nError\n";
          logger.logMessage("Error: The init command was not use.", b);
          return;
        }

        if (s == "g") {
          std::string input;
          std::getline(std::cin, input);

          std::istringstream inputStream(input);


          x = config.defaultX;
          y = config.defaultY;
          sx = config.defaultSx;
          sy = config.defaultSy;
          h = config.defaultH;


          if (!(inputStream >> x)) {
            std::cout << "The default value for x is used: " << config.defaultX
                      << std::endl;
          }
          if (!(inputStream >> y)) {
            std::cout << "The default value for y is used: " << config.defaultY
                      << std::endl;
          }
          if (!(inputStream >> sx)) {
            std::cout << "The default value for sx is used: "
                      << config.defaultSx << std::endl;
          }
          if (!(inputStream >> sy)) {
            std::cout << "The default value for sy is used: "
                      << config.defaultSy << std::endl;
          }
          if (!(inputStream >> h)) {
            std::cout << "The default value for h is used: " << config.defaultH
                      << std::endl;
          }

          logger.logMessage(
              "Adding Gaussian: x=" + std::to_string(x) +
                  ", y=" + std::to_string(y) + ", sx=" + std::to_string(sx) +
                  ", sy=" + std::to_string(sy) + ", h=" + std::to_string(h),
              b);
          c.addgauss(h, x, y, sx, sy);
        }
        if (s == "generate") {
          c.generate();
          std::cout << "Generated values in the field." << std::endl;
          logger.logMessage("Generated values in the field.", b);
        }
        if (s == "gnuplot") {
          c.gnuplot();
          std::cout << "Called gnuplot." << std::endl;
          logger.logMessage("Called gnuplot.", b);
        }
        if (s == "bmp_write") {
          c.bmp_write(c.p->field, "output.bmp");
          std::cout << "Created BMP file." << std::endl;
          logger.logMessage("Created BMP file.", b);
        }
        if (s == "bmp_read") {
          std::cout << "Enter the filename to read:" << std::endl;
          std::cin >> filename;
          c.bmp_read(filename);
          std::cout << "Read BMP file: " + filename << std::endl;
          logger.logMessage("Read BMP file: " + filename, b);
        }
        if (s == "end") {
          std::cout << "Ending the program" << std::endl;
          logger.logMessage("Ending the program.", b);
          break;
        }
        if (s == "bin") {
          int slice;
          std::cout << "Enter slice level:" << std::endl;
          std::cin >> slice;
          c.bin(slice);
          logger.logMessage("Slice applied: slice=" + std::to_string(slice), b);
          std::cout << "Slice applied: slice=" << slice << std::endl;
          logger.logMessage("Wave will be used", b);
          c.wave();
          std::cout << "Amount component = " << c.componenti.size()
                    << std::endl;
        }
      }

      if (file.is_open()) {
        file.close();
        logger.logMessage("Closed input file.", b);
      }
    }
  }
};

int main() {

  std::string configFileName = "config.txt";


  Interface i(configFileName);


  i.print();

  return 0;
}
