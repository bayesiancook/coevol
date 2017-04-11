#ifndef CHRONO_H
#define CHRONO_H

#include <sys/time.h>
#include <chrono>
#include <ctime>
#include <iostream>
#include <sstream>

class Chrono {
  public:
    void Reset();
    void Start();
    void Stop();

    inline int operator++() { return N++; }

    inline double GetTime() { return TotalTime; }

    inline double GetTimePerCount() { return TotalTime / N; }

    inline int GetCount() { return N; }

  private:
    // this is in milli seconds
    double sec1;
    double sec2;
    double milli1;
    double milli2;
    double TotalTime;
    int N;
};

class MeasureTime : public std::stringstream {
  public:
    MeasureTime() : stopped(false) { start(); }

    void start() {
        counter = std::chrono::high_resolution_clock::now();
        stopped = false;
    }

    void stop() {
        stopped = true;
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - counter);
    }

    template <int i>
    void print(std::string message) {
        if (!stopped) {
            stop();
        }
        std::string left(2 * i, ' ');
        std::cout << left << "* " << message << str() << "Time: " << duration.count() << "ms."
                  << std::endl;
        str("");
        start();
    }

    template <int i>
    void print() {
        print<i>("");
    }

  private:
    std::chrono::time_point<std::chrono::high_resolution_clock> counter;
    std::chrono::milliseconds duration;
    bool stopped;
    std::stringstream ss;
};

#endif  // CHRONO_H
