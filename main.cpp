

#include <iostream>
#include "RoadMapReader.hpp"
#include "RoadMapWriter.hpp"
#include "TripReader.hpp"
#include <map>



int main()
{

    InputReader myReader  = InputReader(std::cin);
	RoadMapReader roadMapReader;
	RoadMap roadMap = roadMapReader.readRoadMap(myReader);




    return 0;
}

