#ifndef INPUT_H
#define INPUT_H

#include <stdio.h>

class UserInput {
public:
   UserInput(int);
   ~UserInput();

   void read_stdin(char *);

private:
   FILE *fp;

};
#endif
