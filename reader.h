#ifndef READER_H
#define READER_H

#include "graph.h"
#include <iostream>

class Reader
{
    public:
        virtual Graph& read(std::ifstream& in, long numnodes) = 0;
};
#endif
