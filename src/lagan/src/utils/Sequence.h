// Sequence.h
// ----------
// Class file to hold a sequence object.

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <string>
#include "SafeVector.h"

using namespace std;

class Sequence {

 private:

  // Read header of MFA/XMFA file.
  bool readHeader (ifstream &infile, bool &isXMFA){
    string header;
    
    while (true){
      
      // check to make sure that the there is more data in the file
      if (infile.fail() || infile.eof()) return false;
      
      // get new header line
      getline (infile, header);
      
      // check that header line is not empty
      if (header.length() != 0) break;
    }
    
    // check for appropriate header
    if (header[0] != '>') return false;
    
    // attempt to read XMFA format
    isXMFA = true;
    char buffer[1024];
    int numread = sscanf (header.c_str(), ">%d:%d-%d %c %s", &id, &startCoord, &endCoord, &direction, buffer);
    
    // if basic requirements for XMFA not met, then MFA file
    if (numread < 4){
      comment = header.substr(1);
      isXMFA = false;
    }
    
    // basic requirements for XMFA met, no comments
    else if (numread < 5)
      comment = "";
    
    // otherwise full XMFA format
    else
      comment = buffer;
    
    return true;
  }
  
 protected:

  SafeVector<char> data;     // character data for the sequence
  bool isValid;              // is the sequence valid?
  int length;                // length of the sequence
  int id;                    // sequence ID (for XMFA)
  int startCoord;            // sequence position of first character
  int endCoord;              // sequence position of last character
  char direction;            // + or -
  string comment;            // comments                             

 public:

  Sequence (){
    isValid = true;
    length = 1;
    data.resize (1, ' ');
    startCoord = 1; endCoord = 1;
    direction = '+';
  }

  // Constructor.  Reads in a sequence from the input file.
  Sequence (ifstream &infile){

    bool isXMFA = true;
    
    // sequence starts out not valid
    isValid = false;
    
    // check to make sure that the header is read first
    if (readHeader (infile, isXMFA)){
      
      // put in a dummy character to fill the zero position
      data.push_back ('@');
      
      // read in character data
      char ch;
      
      // loop until no more character data or end of sequence found
      while (infile.get(ch)){
	
	// check to make sure that the end of a section is not reached
	if (ch == '>' || ch == '='){
	  infile.unget();
	  break;
	}
	
	// check for white space
	if (ch == ' ' || ch == '\f' || ch == '\n' || ch == '\r' || ch == '\t' || ch == '\v') continue;
	
	// convert lowercase letters to uppercase
	if (ch >= 'a' && ch <= 'z') ch = ch - 'a' + 'A';
	
	// check that characters are letters OR contig breaks OR gaps
	assert ((ch >= 'A' && ch <= 'Z') || ch == '.' || ch == '-');
	
	
	// add character to list
	data.push_back (ch);
      }
      
      // check to see if any data was read
      if (data.size() > 1){
	
	// if so, the sequence is valid, and compute the length
	isValid = true;
	length = data.size() - 1;
	
	// if the sequence is not originally XMFA
	if (!isXMFA){
	  
	  // assign it some temporary values for XMFA format
	  id = 0;
	  startCoord = 1;
	  endCoord = length;
	  direction = '+';
	}
      }
    }
    
    // some sanity checks
    if (isValid){
      assert (id >= 0);
      assert (startCoord >= 0);
      assert (endCoord >= 0);
      assert (startCoord <= endCoord);
      assert (direction == '+' || direction == '-');
      assert (length > 0);
    }
  }

  // Constructor.  Gets sequence from array data.
  Sequence (SafeVector<char> data, string comment) : data(data), comment(comment) {
    length = data.size() - 1;
    id = 0;
    startCoord = 1;
    endCoord = length;
    direction = '+';
    isValid = true;
    comment = "";

    assert (length > 0);
  }

  SafeVector<char> getData (){
    SafeVector<char> temp;
    for (int i = 1; i <= length; i++) temp.push_back (data[i]);
    return temp;
  }

  const string getComment () const {
    return comment;
  }

  void setLength (int num){
    if (num > length){
      length = num;
      endCoord = length;
      data.resize(length+1, ' ');
    }
  }

  SafeVector<char>::iterator getIterator (){
    return data.begin();
  }

  const char operator[] (int index) const {
    assert (index >= 1 && index <= length);
    return data[index];
  }

  // Used to check for sequence validity after construction.
  const bool fail () const { return !isValid; }

  // Return the length of the sequence.
  const int getLength () const { assert (isValid); return length; }
  const char getStrand () const { assert (isValid); return direction; }
  
  const int getStartCoord () const { assert (isValid); return startCoord; }
  const int getEndCoord () const { assert (isValid); return endCoord; }

  // Print XMFA header only.
  void writeXMFAHeader (ostream &outfile) const {
    assert (isValid);
    outfile << '>' << id << ':' << startCoord << '-' << endCoord << ' ' << direction << ' ' << comment << endl;
  }

  // Return sequence ID.
  const int getID () const { assert (isValid); return id; }

  // Set sequence ID.
  void setID (int id) { assert (isValid); this->id = id; }

  // Writes sequence to XMFA format.
  void writeToXMFA (ostream &outfile, int numColumns) const {

    assert (isValid);
    
    // print XMFA header
    outfile << ">" << comment << endl;
    //  outfile << '>' << id << ':' << startCoord << '-' << endCoord << ' ' << direction << ' ' << comment << endl;
    
    // print character data
    for (int i = 1; i <= length; ++i){
      outfile << data[i];      
      if (i % numColumns == 0) outfile << endl;
    }
    if (length % numColumns != 0) outfile << endl;
  }
};

#endif
