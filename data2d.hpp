#ifndef _DATA2D_H
#define _DATA2D_H

#include <vector>
#include <iostream>
#include <fstream>


typedef std::vector <double> record_t;
typedef std::vector <record_t> data_t;

class Data2d 
{
public:
	Data2d(void);
	Data2d(const unsigned nrecords,const unsigned nvalues);
	~Data2d(void);
	record_t & operator[](std::size_t idx){return (*this)[idx];}
	void readfile(std::string & infilename);
	void transpose(void);
	inline void setval(const unsigned i, const unsigned j, double val){data[i][j] = val;}
	void copyrow(const unsigned i,std::vector<double> &row);
	void copycol(const unsigned j,std::vector<double> &col);

private:
	bool transposecopy;
	record_t record;
	data_t data;
};
#endif
