#include "data2d.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

//-----------------------------------------------------------------------------
// Let's overload the stream input operator to read a list of CSV fields (which a CSV record).
// Remember, a record is a list of doubles separated by commas ','.
std::istream& operator >> ( std::istream& ins, record_t& record )
{
  record.clear();
  std::string line;
  getline( ins, line );
  std::istringstream iss( line ); // I think istringstream is an in-stringstream, we should try using ostringstream to get outfiles
  double value;
  while (iss >> value){
        record.push_back(value);
  }
  return ins;
}

//-----------------------------------------------------------------------------
// Let's likewise overload the stream input operator to read a list of CSV records.
// This time it is a little easier, just because we only need to worry about reading
// records, and not fields.
std::istream& operator >> ( std::istream& ins, data_t& data )
{
  data.clear();
  record_t record;
  while (ins >> record)
    {
    data.push_back( record );
    }
  return ins;
}

Data2d::Data2d(void)
{

}

Data2d::Data2d(const unsigned nrecords,const unsigned nvalues){
	data.resize(nrecords);
	for (unsigned i=0;i<data.size();i++){
		data[i].resize(nvalues);
	}
}

Data2d::~Data2d(void){
}



void Data2d::readfile(std::string & infilename)
{
	std::ifstream infile;
	infile.open(infilename.c_str(),std::ios::in);
        infile >> data; // because of nested overloaded operator >> from the typedef lines. 
        if (!infile.eof()){
                std::cerr << "nested overloaded >> must have failed" << std::endl;
                return;
        }
        infile.close();
}
