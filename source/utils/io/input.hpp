#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

template<typename T>
std::vector<T> readInput(const std::string &fileName)
{
	std::vector<T> vector;
	std::fstream file;
	file.open(fileName, std::ios::in);
	if(!file.good())
    {
		std::cerr << "Error in reading from " << fileName << std::endl;
		exit(EXIT_FAILURE);
	}
	else
    {
		while(true)
        {
			if(file.eof())
			{
				break;
			}
			std::string str;
			std::getline(file, str);
			if(str.length() > 0)
			{
				std::stringstream strstream(str);
				T x;
				strstream >> x;
				vector.push_back(x);
			}
        }
	}
	file.close();
	return vector;
}