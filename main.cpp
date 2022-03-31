// main.cpp
//
// ICS 46 Winter 2022
// Project #5: Rock and Roll Stops the Traffic
//
// This is the program's main() function, which is the entry point for your
// console user interface.

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

