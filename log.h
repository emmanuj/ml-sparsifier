#ifndef LOG_H
#define LOG_H
#include "OptionParser.h"
#include <iostream>
#include <string>
//simple logger class
class Log{
public:
	static Log& log(optparse::Values& opt){
		static Log logger(opt);
		return logger;
	}
	//info
	void d(std::string& msg){
		if((bool)options.get("debug")){
			std::cout<<msg<<std::endl;
		}
	}
	//verbose
	void i(std::string& msg){
		if((bool)options.get("log-info")){
			std::cout<<msg<<std::endl;
		}
	}

	//error
	void e(std::string& msg){
		std::cout<<msg<<std::endl;
	}
	
private:
	optparse::Values& options;
	Log(optparse::Values& op):options(op){}
	Log();
	Log(const Log&);
	Log& operator=(const Log&);
	~Log();

};

#endif
