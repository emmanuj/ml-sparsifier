#ifndef ETIMER_H
#define ETIMER_H

#include <sys/times.h>
#include <iostream>

class ETimer{

private:
    struct tms st_cpu;
    struct tms en_cpu;
    clock_t st_time;
    clock_t en_time;

public:

    ETimer():st_cpu(), en_cpu(), st_time(), en_time(){}

    void start(){
    	st_time = times(&st_cpu);
    }
    void stop(const std::string desc){
        en_time = times(&en_cpu);
        std::cout<< desc << " : "<<(double)(en_time - st_time)/100<<" secs."<<std::endl;
    }
};


#endif // ETIMER_H