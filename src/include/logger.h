/*
*  Copyright (C) 2012 by Anush Krishnan, Simon Layton, Lorena Barba
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*/

#pragma once

#include <sys/time.h>
#include <map>
#include <string>
#include <iostream>

typedef std::map<std::string, double>            Event;
typedef std::map<std::string, double>::iterator  E_iter;

class Logger {
private:
	std::ofstream file;
	std::ofstream stepFile;
	std::ofstream legendFile;
	Event         tic;
	Event         timer;    // times for total
	Event         timeStep; // times for 1 step
	Event         memory;

	double get_time() // Timer function
	{
		struct timeval tv;                        // Time value
		gettimeofday(&tv, NULL);                  // Get time of day in seconds and microseconds
		return double(tv.tv_sec+tv.tv_usec*1e-6); // Combine seconds and microseconds and return
	}

public:
	bool printNow;

	Logger(){}
	
	Logger(std::string folder)
	{
		file.open((folder + "/time").c_str());
		stepFile.open((folder + "/profiling").c_str());
		legendFile.open((folder + "/profiling_legend").c_str());
		printNow = false;
	}
	~Logger()
	{
		writeLegend();
		file.close();
		stepFile.close();
	}

	void startTimer(std::string event)
	{
		tic[event] = get_time();
	}

	void stopTimer(std::string event, bool print=false)
	{
		double toc = get_time();
		timeStep[event] += toc - tic[event];
		timer[event] += toc - tic[event];
		if(print) 
			std::cout << event << " : " << timer[event] << std::endl;
	}

	void eraseTimer(std::string event)
	{
		timer.erase(event);
	}

	void resetTimer()
	{
		for( E_iter E=timer.begin(); E!=timer.end(); ++E )
			E->second = 0;
	}
  
	void resetTimeStep()
	{
		for (E_iter E=timeStep.begin(); E!=timeStep.end(); ++E )
			E->second = 0;
	}

	void allocMemory(std::string event, double bytes)
	{
		memory[event] += bytes;
	}

	void freeMemory(std::string event, double bytes)
	{
		memory[event] -= bytes;
	}

	void printTime(std::string event)
	{
		std::cout << event << " : " << timer[event] << std::endl;
	}

	void printMemory(std::string event)
	{
		std::cout << event << " : " << memory[event] << std::endl;
	}

	void printAllTime()
	{
		double totalTime = 0.0;
		std::cout << std::endl;
		for( E_iter E=timer.begin(); E!=timer.end(); ++E )
		{
			totalTime += E->second;
			std::cout << std::setw(24) << E->first << std::setw(13) << std::fixed \
			          << std::setprecision(4) << E->second << std::endl;
		}
		std::cout << "-------------------------------------" << std::endl;
		std::cout << std::setw(24) << "TOTAL" << std::setw(13) << std::fixed \
			          << std::setprecision(4) << totalTime << std::endl;
	}
  
	void writeLegend()
	{
		for( E_iter E=timer.begin(); E!=timer.end(); ++E )
		{
			legendFile <<  E->first << std::endl;
		}
	}

	void writeTime()
	{
		for( E_iter E=timer.begin(); E!=timer.end(); ++E )
			file <<  E->first << " " << E->second << std::endl;
	}
  
	void writeTimeStep(int n)
	{
		// std::cout << "writeTimeStep called" << std::endl;
		stepFile << n << "\t";
		for ( E_iter E=timeStep.begin(); E!=timeStep.end(); ++E )
			stepFile << E->second << "\t"; //  << std::endl;
		stepFile << std::endl;
	}

};
