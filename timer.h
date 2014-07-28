#ifndef TIMER_H
#define TIMER_H

#include <vector>
#include <unordered_map>
#include <boost/timer/timer.hpp>

using boost::timer::cpu_timer;
using namespace std;

class Timer
{
    public:
        Timer(bool isOn);
        cpu_timer* TPointer(string _tName); 
        
        void StopAll();
        void StartAll();
        void ResumeAll();
        
        void stop(string ttype);
        void start(string ttype);
        int resume(string ttype);
 
    private:
        bool state;
        vector <string>   types;
        unordered_map <string,cpu_timer*> TList;
}; 
#endif
