#ifndef RANDOMBASE_H
#define RANDOMBASE_H

#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

using namespace std;

class RandomBase {
    public:
        //Define data-types of random objects 
        typedef boost::mt19937 t_eng;                 //Mersenne twister
        typedef boost::uniform_real<float>   t_uReal;   //Uniform real distribution
        typedef boost::uniform_int<uint64_t> t_uInt;    //Uniform longeger distribution

        //Random generation objects
        t_eng   eng;                 //Mersenne twister
        t_uReal uReal;               //Uniform real distribution
        t_uInt  uInt;                //Uniform distribution 

        boost::variate_generator<t_eng&, t_uInt > uRandInt; //Generator of integers
        boost::variate_generator<t_eng&, t_uReal> uRand;    //Generator of real number
    
        RandomBase(long seed);
};
#endif
