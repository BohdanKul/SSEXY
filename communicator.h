#include <fstream>
#include <vector>
#include <unordered_map>

using namespace std;

class Communicator
{
    public:
        Communicator(int Nx, int Ny, float T);
        fstream* stream(string _fileName); 
        string dataName;
        string outDir;
        vector <string>   types;
        vector <fstream> files;   

    private:
        unsigned long id;
        void GenerateId();
        unordered_map <string,fstream*> mFStreams;
}; 
