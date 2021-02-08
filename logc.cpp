#include <stddef.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdbool.h>

//int verbosity = 10;
#include "logc.h"

void _error(const char* msg) {
  _logc(ERROR, msg);
}

void error(const char* fmt, ...) {
  char buf[MAX_LOG_MSG_SIZE];
  va_list vl;
  va_start(vl, fmt);
  vsnprintf(buf, sizeof(buf), fmt, vl);
  va_end(vl);
  _error(buf);
}

void _timing(const char* msg) {
  _logc(TIME, msg);
}

void timing(const char* fmt, ...) {
    char buf[MAX_LOG_MSG_SIZE];
    va_list vl;
    va_start(vl, fmt);
    vsnprintf(buf, sizeof(buf), fmt, vl);
    va_end(vl);
    _timing(buf);
}

void _info(const char* msg) {
  _logc(INFO, msg);
}

void info(const char* fmt, ...) {
    char buf[MAX_LOG_MSG_SIZE];
    va_list vl;
    va_start(vl, fmt);
    vsnprintf(buf, sizeof(buf), fmt, vl);
    va_end(vl);
    _info(buf);
  
}


void dbgcollect(const char* fmt, ...) {
    char buf[MAX_LOG_MSG_SIZE];
    va_list vl;
    va_start(vl, fmt);
    vsnprintf(buf, sizeof(buf), fmt, vl);
    va_end(vl);
    _logc(COLLECT, buf);
}


void dbggrid(const char* fmt, ...) {
    char buf[MAX_LOG_MSG_SIZE];
    va_list vl;
    va_start(vl, fmt);
    vsnprintf(buf, sizeof(buf), fmt, vl);
    va_end(vl);
    _logc(GRID, buf);
}


void dbgconnect(const char* fmt, ...) {
    char buf[MAX_LOG_MSG_SIZE];
    va_list vl;
    va_start(vl, fmt);
    vsnprintf(buf, sizeof(buf), fmt, vl);
    va_end(vl);
    _logc(CONNECT, buf);
}


void dbgfit(const char* fmt, ...) {
    char buf[MAX_LOG_MSG_SIZE];
    va_list vl;
    va_start(vl, fmt);
    vsnprintf(buf, sizeof(buf), fmt, vl);
    va_end(vl);
    _logc(FIT, buf);
}

void dbgmerge(const char* fmt, ...) {
    char buf[MAX_LOG_MSG_SIZE];
    va_list vl;
    va_start(vl, fmt);
    vsnprintf(buf, sizeof(buf), fmt, vl);
    va_end(vl);
    _logc(MERGE, buf);
}

void dbgtrkz(const char* fmt, ...) {
    char buf[MAX_LOG_MSG_SIZE];
    va_list vl;
    va_start(vl, fmt);
    vsnprintf(buf, sizeof(buf), fmt, vl);
    va_end(vl);
    _logc(TRKZ, buf);
}

void dbgtrkerror(const char* fmt, ...) {
    char buf[MAX_LOG_MSG_SIZE];
    va_list vl;
    va_start(vl, fmt);
    vsnprintf(buf, sizeof(buf), fmt, vl);
    va_end(vl);
    _logc(TRKERROR, buf);
}

void _logc(int level, const char* msg) {

  if (verbosity.test(level)) {
    char *ts = timestamp();
  
    switch (level) {
    case ERROR:
      fprintf(stderr, WHT"%s " RED "ERROR " WHT "%s\n" RESET, ts, msg); fflush(stderr); break;
    case TIME:
      fprintf(stdout, WHT"%s " YEL "TIME  " WHT "%s\n" RESET, ts, msg); fflush(stdout); break;
    case INFO:
     fprintf(stdout, WHT"%s "  GRN "INFO  " WHT "%s\n" RESET, ts, msg); fflush(stdout); break;
    case COLLECT:
      fprintf(stdout, WHT"%s " CYN "COLLECT " WHT "%s\n" RESET, ts, msg); fflush(stdout); break;
    case GRID:
      fprintf(stdout, WHT"%s " CYN "GRID " WHT "%s\n" RESET, ts, msg); fflush(stdout); break;
    case CONNECT:
      fprintf(stdout, WHT"%s " BLU "CONNECT " WHT "%s\n" RESET, ts, msg); fflush(stdout); break;
    case FIT:
      fprintf(stdout, WHT"%s " BLU "FIT " WHT "%s\n" RESET, ts, msg); fflush(stdout); break;
    case MERGE:
      fprintf(stdout, WHT"%s " BLU "MERGE " WHT "%s\n" RESET, ts, msg); fflush(stdout); break;
    case TRKERROR:
      fprintf(stdout, WHT"%s " MAG "TRKERROR " WHT "%s\n" RESET, ts, msg); fflush(stdout); break;
    case TRKZ:
      fprintf(stdout, WHT"%s " BLU "Z COORD " WHT "%s\n" RESET, ts, msg); fflush(stdout); break;
    default:
      fprintf(stderr, WHT"%s "  RED "INVALID LEVEL " WHT "%s\n" RESET, ts, msg); fflush(stderr); break;
    }

    free(ts);
  }
}

void logc(int level, const char* fmt, ...) {
  char buf[MAX_LOG_MSG_SIZE];
  va_list vl;
  va_start(vl, fmt);
  vsnprintf(buf, sizeof(buf), fmt, vl);
  va_end(vl);
  _logc(level, buf);
}

void set_verbosity(bool v[10]) {
  /* to upper case */
  for (int i = 0; i < 10; ++i) {
    verbosity.set(i, v[i]);
  }

  // (std::string("0000000011"));
  /* if (equals(v, "OFF")) {
    verbosity = OFF;
  } else if (equals(v, "ERROR")) {
    verbosity = ERROR;
  } else if (equals(v, "WARN")) {
    verbosity = WARN;
  } else if (equals(v, "INFO")) {
    verbosity = INFO;
  } else if (equals(v, "DEBUG")) {
    verbosity = DEBUG;
  } else if (equals(v, "TRACE")) {
    verbosity = TRACE;
  } else if (equals(v, "TIMING")) {
    verbosity = TIMING;
  } else if (equals(v, "ALL")) {
    verbosity = ALL;
  } else {
    _error("No valid verbosity level supplied");
    }*/
}

bool equals(const char *a, const char *b) {
  return (strcmp(a, b) == 0);
}

char *timestamp(void) {
  time_t rawtime;
  struct tm *info;
  char *buf = (char *) malloc(80 * sizeof(char));
  time( &rawtime );
  info = localtime( &rawtime );

  strftime(buf, 80, "%H:%M:%S", info);
  return buf;
}
