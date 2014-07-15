/***************************************************************************//**
* \file logger.h
* \author Krishnan, A. (anush@bu.edu)
* \brief Definition of the class \c Logger
*/

#pragma once

#include <sys/time.h>
#include <map>
#include <string>
#include <iostream>

/**
* \typedef
* \brief Map from a \c std::string to a \c double
*/
typedef std::map<std::string, double>            Event;

/**
* \typedef
* \brief Iterator over a STL map
*/
typedef std::map<std::string, double>::iterator  E_iter;

/**
* \class Logger
* \brief Track the time spent to achieve a certain task
*/
class Logger
{
private:
	std::ofstream file;			///< file to write the time spent
	std::ofstream stepFile;		///< file to write the time spent for 1 step
	std::ofstream legendFile;	///< file to write the name of the events
	Event         tic;			///< contain the starting time of events
	Event         timer;		///< contain the total time spent of events during the simulation
	Event         timeStep; 	///< contain the tiem spent times of events for 1 step
	Event         memory;		///< contain info about the memomry allocated for events

	/**
	* \brief Timer function
	*
	* \return the time of the day
	*/
	double get_time()
	{
		struct timeval tv;                        // time value
		gettimeofday(&tv, NULL);                  // get time of day in seconds and microseconds
		return double(tv.tv_sec+tv.tv_usec*1e-6); // combine seconds and microseconds and return
	}

public:
	bool printNow;	///< do you want to profile timing on the go?

	/**
	* \brief Constructor of the class \c Logger
	*/
	Logger(){}
	
	/**
	* \brief Overloaded constructor of the class \c Logger
	*
	* \param folder the case folder
	*/
	Logger(std::string folder)
	{
		file.open((folder + "/time").c_str());
		stepFile.open((folder + "/profiling").c_str());
		legendFile.open((folder + "/profiling_legend").c_str());
		printNow = false;
	}

	/**
	* \brief Destructor of the class \c Logger
	*/
	~Logger()
	{
		writeLegend();
		file.close();
		stepFile.close();
	}

	/**
	* \brief Start the timer
	*
	* \param event a string to describe the event
	*/
	void startTimer(std::string event)
	{
		tic[event] = get_time();
	}

	/**
	* \brief Stop the timer
	*
	* \param event a string that describes the event
	* \param print do you want to print the time?
	*/
	void stopTimer(std::string event, bool print=false)
	{
		double toc = get_time();
		timeStep[event] += toc - tic[event];
		timer[event] += toc - tic[event];
		if(print) 
			std::cout << event << " : " << timer[event] << std::endl;
	}

	/**
	* \brief Erase the timer
	*
	* \param event a string that describes the event
	*/
	void eraseTimer(std::string event)
	{
		timer.erase(event);
	}

	/**
	* \brief Reset the timer
	*/
	void resetTimer()
	{
		for( E_iter E=timer.begin(); E!=timer.end(); ++E )
			E->second = 0;
	}

	/**
	* \brief Reset the time-step timer
	*/
	void resetTimeStep()
	{
		for (E_iter E=timeStep.begin(); E!=timeStep.end(); ++E )
			E->second = 0;
	}

	/**
	* \brief Allocate memory
	*
	* \param event a string that describes the event
	* \param bytes the number of bytes
	*/
	void allocMemory(std::string event, double bytes)
	{
		memory[event] += bytes;
	}

	/**
	* \brief Free memory
	*
	* \param event a string that describes the event
	* \param bytes the number of bytes
	*/
	void freeMemory(std::string event, double bytes)
	{
		memory[event] -= bytes;
	}

	/**
	* \brief Print the time spent to achieve a task
	*
	* \param event a string that describes the event
	*/
	void printTime(std::string event)
	{
		std::cout << event << " : " << timer[event] << std::endl;
	}

	/**
	* \brief Print the memory allocated for given event
	*
	* \param event a string that describes the event
	*/
	void printMemory(std::string event)
	{
		std::cout << event << " : " << memory[event] << std::endl;
	}

	/**
	* \brief Print time spent of every event and the total time
	*/
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

	/**
	* \brief Write the events into a file
	*/
	void writeLegend()
	{
		for( E_iter E=timer.begin(); E!=timer.end(); ++E )
		{
			legendFile <<  E->first << std::endl;
		}
	}

	/**
	* \brief Write time spent into a file
	*/
	void writeTime()
	{
		for( E_iter E=timer.begin(); E!=timer.end(); ++E )
			file <<  E->first << " " << E->second << std::endl;
	}
  
	/**
	* \brief Write time spent for an iteration into a file
	*
	* \param n the iteration number
	*/
	void writeTimeStep(int n)
	{
		stepFile << n << "\t";
		for ( E_iter E=timeStep.begin(); E!=timeStep.end(); ++E )
			stepFile << E->second << "\t";
		stepFile << std::endl;
	}
};
