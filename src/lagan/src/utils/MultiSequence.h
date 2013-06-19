// MultiSequence.h
// ---------------
// Multiple sequence class

#ifndef MULTISEQUENCE_H
#define MULTISEQUENCE_H

#include <vector>
#include <string>
#include <fstream>
#include <stdio.h>
#include "Sequence.h"
#include "SafeVector.h"

using namespace std;

class MultiSequence {
 private:
  SafeVector<Sequence> sequences;     // sequences
  SafeVector<char> cache;
  bool cacheEnabled;
    
 public:

  MultiSequence (): cacheEnabled (false) {}

  void buildCache (){
    assert (!cacheEnabled);
    cacheEnabled = true;

    int length = sequences[0].getLength();
    int numSeqs = getNumSeqs();

    cache.resize ((length + 1) * numSeqs, (char) 0);
    for (int i = 0; i < numSeqs; i++){
      Sequence &seq = (*this)[i];
      cache[i] = '@';      
      for (int j = 1; j <= length; j++){
	cache[j * numSeqs + i] = seq[j];
      }
    }
  }

  // return letter cache for fast processing
  SafeVector<char>::iterator getCache (){
    assert (cacheEnabled);
    return cache.begin();
  }
  
  // add a sequence to the alignment
  void addSequence (Sequence &sequence){
    sequences.push_back (sequence);
  }

  // Read in all of the Sequences in an MFA file and append them to the
  // existing MultiSequence object.
  void addRawFromMFA (const string& filename){
    
    // open up file for reading
    ifstream infile (filename.c_str());
    
    // check for error
    assert (!infile.fail());
    
    // add only sequences that check out ok
    while (true){
      Sequence seq (infile);
      if (seq.fail()) break;
      sequences.push_back (seq);
    }
    
    // close up the input file
    infile.close();
  }

  // Read in all of the Sequences in an MFA file and append them to the
  // existing MultiSequence object.
  void addRawFromMFA (ifstream &infile){
    
    // check for error
    assert (!infile.fail());
    
    // add only sequences that check out ok
    while (true){
      Sequence seq (infile);
      if (seq.fail()) break;
      sequences.push_back (seq);
    }
  }

  // Writes sequences to outfile in XMFA format.
  void writeToXMFA (ostream &outfile, int numColumns) const {
    for (int i = 0; i < (int) sequences.size(); ++i){
      sequences[i].writeToXMFA (outfile, numColumns);
    }
  }

  // Returns a sequence.
  Sequence& operator[] (int index){

    // error checking on bounds
    assert (index >= 0 && index < (int) sequences.size());
    
    // return the correct sequence
    return sequences[index];
  }

  // Returns a sequence.
  const Sequence& operator[] (int index) const {
    
    // error checking on bounds
    assert (index >= 0 && index < (int) sequences.size());
    
    // return the correct sequence
    return sequences[index];
  }

  // Returns number of sequences.
  const int getNumSeqs() const {
    return sequences.size();
  }
};

#endif
