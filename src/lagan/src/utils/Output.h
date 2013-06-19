#ifndef OUTPUT_H
#define OUTPUT_H

// print reversed string in MFA format
void printMFA (ostream &outfile, SafeVector<char> &data, string comment, int numColumns){

  int charsWritten = 0;

  outfile << ">" << comment << endl;
  for (int i = 0; i < (int) data.size(); i++){
    outfile << data[i];
    charsWritten++;
    if (charsWritten % numColumns == 0) outfile << endl;
  }
  
  if (charsWritten % numColumns != 0) outfile << endl;
}


#endif
