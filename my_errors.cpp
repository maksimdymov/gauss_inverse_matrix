#include "my_errors.h"

//Prints error message by code of error and returns this code
int
PrintErrorMsgByCode (int code, char* msg)
{
  switch (code)
  {
    case SUCCESS:
      break;

    case MAIN_ARGS_ERROR:
      fprintf (stderr, "Usage: %s <n> <m> <r> <s> [file]\n", msg);
      break;

    case FILE_OPENING_ERROR:
      fprintf (stderr, "Cannot open file %s\n", msg);
      break;

    case FILE_READING_ERROR:
      fprintf (stderr, "Cannot read from file %s\n", msg);
      break;

    case FILE_PRINTING_ERROR:
      fprintf (stderr, "Cannot print in file %s\n", msg);
      break;

    case ALLOCATE_MEMORY_ERROR:
      fprintf (stderr, "Cannot allocate %s size\n", msg);
      break;

    case END_OF_FILE:
      fprintf (stderr, "Reached end of file %s, but matrix wasn't completely filled\n", msg);
      break;

    case TOO_MANY_INFO:
      fprintf (stderr, "Too many information in file %s\n", msg);
      break;

    default:
      fprintf (stderr, "Sorry, something went wrong...\nError code: %d\nError message: %s", code, msg);
      break;
  }
  return code;
}