#pragma once
/**
* Stores all the command line options
*/
class options
{
public:
	std::string folderName,
	            inputFile,
	            domFile,
	            systemFile;
	int  cudaDevice;
	bool writeData;
	options()
	{
		folderName = "new";
		domFile = "domains/cavity.dom";
	}
};