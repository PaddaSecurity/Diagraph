

#ifndef ROADMAPREADER_HPP
#define ROADMAPREADER_HPP

#include "RoadMap.hpp"
#include "InputReader.hpp"



class RoadMapReader 
{
public:
    // readRoadMap() reads a RoadMap from the given InputReader.  The
    // RoadMap is expected to be described in the format given in the
    // project write-up.
    RoadMap readRoadMap(InputReader& in);
};



#endif

