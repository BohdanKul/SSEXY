#include "timer.h"


using namespace std;
using boost::timer::cpu_timer;

Timer::Timer(bool isOn)
{
    state = false;
    if  (isOn){
        types = vector<string> {"Total","RandomUpdate","LinksConstruct","DeterUpdate","DetPathTracing","DetLoopsConnect"};
        for (auto type=types.begin(); type!=types.end(); type++){
            TList[*type] = new cpu_timer();
            TList[*type]->stop();
        }
        state = true;
    }
}

cpu_timer* Timer::TPointer(string _tName)
{
    return TList[_tName];
}

void Timer::stop(string ttype)
{
    if (state)
        TPointer(ttype)->stop();
}

void Timer::start(string ttype)
{
    if (state)
        TPointer(ttype)->start();
}

int Timer::resume(string ttype)
{
    //if (state){
        //cout << ttype << endl;
        //TPointer(ttype)->resume();
    //}
}



void Timer::StopAll()
{
    if (state)
        for (auto type=types.begin(); type!=types.end(); type++){
            TPointer(*type)->stop();
        }

}

void Timer::StartAll()
{
    if (state)
        for (auto type=types.begin(); type!=types.end(); type++){
            TPointer(*type)->start();
        }

}
void Timer::ResumeAll()
{
    if (state)
        for (auto type=types.begin(); type!=types.end(); type++){
            TPointer(*type)->resume();
        }

}

//std::cout << "Total: " << boost::timer::format(timerTotal.elapsed(),5) << std::endl;
