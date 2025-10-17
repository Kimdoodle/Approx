#pragma once
#include <chrono>
#include <iostream>
using namespace std::chrono;


time_point<high_resolution_clock> cur_time();
void calculate_time(time_point<high_resolution_clock> time1, time_point<high_resolution_clock> time2);