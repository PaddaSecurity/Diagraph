
#ifndef INPUTREADER_HPP
#define INPUTREADER_HPP

#include <istream>
#include <string>



class InputReader
{
public:
    // Initializes an InputReader so that it reads from the given input
    // stream.  For example, pass std::cin as a parameter to the constructor
    // if you want to read input from std::cin.
    InputReader(std::istream& in);

    // readLine() reads a line of input from the input stream associated
    // with this InputReader, skipping non-meaningful lines.
    std::string readLine();

    // readLineInt() reads a line of input from the input stream associated
    // with this InputReader, assuming that the line of input contains an
    // integer value (e.g., "7").
    int readIntLine();

private:
    std::istream& in_;
};



inline InputReader::InputReader(std::istream& in)
    : in_{in}
{
}



#endif

