#ifndef ETIMER_H
#define ETIMER_H


#include <chrono>
#include <iostream>

class ETimer{

private:

    std::chrono::high_resolution_clock::time_point t1,t2;

public:

    ETimer():
        t1(), t2(){}

    void start(){
    	t1 = std::chrono::high_resolution_clock::now();
    }
    void stop(const std::string desc){
        t2 = 	std::chrono::high_resolution_clock::now();
        double duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
        std::cout << desc<<" : " << duration << " microseconds " << std::endl;
    }
};

#endif // ETIMER_H
