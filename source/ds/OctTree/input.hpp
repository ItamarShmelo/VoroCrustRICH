/*

#ifndef _INPUT_HPP
#define _INPUT_HPP

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

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
            T x;
			std::string str;
            std::getline(file, str);
            std::stringstream sstr(str);
            sstr >> x;
            vector.push_back(x);
        }
	}
	file.close();
	return vector;
}

#endif // _INPUT_HPP
*/