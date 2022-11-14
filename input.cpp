#include "input.h"
#include "global.h"

/* -------------------------------------------------------------------
 * Constructor. If flag = 1, output user inputs as stdin.log
 * ---------------------------------------------------------------- */ 
UserInput::UserInput(int flag)
{
   if (flag) fp = fopen("script.inp", "w");
   return;
}

/* -------------------------------------------------------------------
 * Deconstructor. Output user inputs as required and clear workspace.
 * ---------------------------------------------------------------- */ 
UserInput::~UserInput()
{
   if (fp) fclose(fp);
   fp = NULL;
}

/* -------------------------------------------------------------------
 * Read stdin and keep a record of it.
 * ---------------------------------------------------------------- */ 
void UserInput::read_stdin(char *str)
{
   fgets(str, MAXLINE, stdin);
   if (fp) fprintf(fp, "%s", str);

   return;
}
/* ---------------------------------------------------------------- */
