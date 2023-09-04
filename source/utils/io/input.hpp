#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <filesystem>

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

template<typename T>
std::vector<T> readDirectory(const std::string &dirName)
{
	std::vector<T> result;
	for(const auto &file : std::filesystem::directory_iterator(dirName))
	{
		if(file.is_regular_file())
		{
			std::vector<T> tempRes = readInput<T>(file.path());
			result.insert(result.end(), tempRes.begin(), tempRes.end());
		}
	}
	return result;
}

template<typename T>
std::vector<T> readListFiles(const std::vector<std::string> &listOfFiles)
{
	std::vector<T> result;
	for(const auto &fileName : listOfFiles)
	{
		std::vector<T> tempRes = readInput<T>(fileName);
		result.insert(result.end(), tempRes.begin(), tempRes.end());
	}
	return result;
}