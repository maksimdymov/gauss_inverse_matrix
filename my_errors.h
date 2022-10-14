#ifndef ERRORS_H_
#define ERRORS_H_

#include <cstdio>

int PrintErrorMsgByCode (int code, char* msg);

enum ReturnValues
{
  TOO_MANY_INFO = -7,
  END_OF_FILE,
  ALLOCATE_MEMORY_ERROR,
  FILE_PRINTING_ERROR,
  FILE_READING_ERROR,
  FILE_OPENING_ERROR,
  MAIN_ARGS_ERROR,
  SUCCESS,
  CANNOT_SOLVE,
};

#endif