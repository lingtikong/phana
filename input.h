#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>
#include <string>
#include <vector>
#include <string.h>

using namespace std;

class UserInput {
public:
   UserInput(int);
   ~UserInput();

   void read_stdin(char *);

private:
   FILE *fp;

};
#endif
