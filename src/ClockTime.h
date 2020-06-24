#ifndef __ClockTime__h
#define __ClockTime__h

#include <bits/stdc++.h> 
#include <chrono> 

/************************************************************
 * This class is used to control time spent by the program.						
************************************************************/
class ClockTime {
    
private:
    std::chrono::_V2::system_clock::time_point START;

public:
    static std::chrono::_V2::system_clock::time_point getTimeNow() { return std::chrono::high_resolution_clock::now(); }
    
    ClockTime(std::chrono::_V2::system_clock::time_point t) : START(t) {};
    std::chrono::_V2::system_clock::time_point getStart() const { return START; }
    double getTimeInSecFromStart() const { return 1e-9*(std::chrono::duration_cast<std::chrono::nanoseconds>(getTimeNow() - getStart()).count()); }

    void setStart(const std::chrono::_V2::system_clock::time_point &t){ START = t; }
};

#endif