#include "logc.h"
#include <stddef.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <stdbool.h>
/* verbosity levels are defined in logc.h */
int verbosity = 10;


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

void _warn(const char* msg) {
  _logc(WARN, msg);
}

void warn(const char* fmt, ...) {
  char buf[MAX_LOG_MSG_SIZE];
  va_list vl;
  va_start(vl, fmt);
  vsnprintf(buf, sizeof(buf), fmt, vl);
  va_end(vl);
  _warn(buf);
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

void _debug(const char* msg) {
  _logc(DEBUG, msg);
}

void debug(const char* fmt, ...) {
  char buf[MAX_LOG_MSG_SIZE];
  va_list vl;
  va_start(vl, fmt);
  vsnprintf(buf, sizeof(buf), fmt, vl);
  va_end(vl);
  _debug(buf);
}

void _trace(const char* msg) {
  _logc(TRACE, msg);
}

void trace(const char* fmt, ...) {
  char buf[MAX_LOG_MSG_SIZE];
  va_list vl;
  va_start(vl, fmt);
  vsnprintf(buf, sizeof(buf), fmt, vl);
  va_end(vl);
  _trace(buf);
}

void _timing(const char* msg) {
  _logc(TIMING, msg);
}

void timing(const char* fmt, ...) {
    char buf[MAX_LOG_MSG_SIZE];
    va_list vl;
    va_start(vl, fmt);
    vsnprintf(buf, sizeof(buf), fmt, vl);
    va_end(vl);
    _timing(buf);
  
}

void _logc(int level, const char* msg) {
  if (verbosity >= level) {
    char *ts = timestamp();
    switch (level) {
    case ERROR:
      fprintf(stderr, CYN"%s " YEL "%d " RED "ERROR " WHT "%s\n" RESET, ts, 0, msg); fflush(stderr); break;
    case WARN:
      fprintf(stdout, CYN"%s " YEL "%d " YEL "WARN  " WHT "%s\n" RESET, ts, 0, msg); fflush(stdout); break;
    case TIMING:
      fprintf(stdout, CYN"%s " YEL "%d " BLU "TIME  " WHT "%s\n" RESET, ts, 0, msg); fflush(stdout); break;
    case INFO:
     fprintf(stdout, CYN"%s " YEL "%d " GRN "INFO  " WHT "%s\n" RESET, ts, 0, msg); fflush(stdout); break;
     //fprintf(stdout, "%s\n", msg); fflush(stdout); break;
    case DEBUG:
      fprintf(stdout, CYN"%s " YEL "%d " CYN "DEBUG " WHT "%s\n" RESET, ts, 0, msg); fflush(stdout); break;
    case TRACE:
      fprintf(stdout, CYN"%s " YEL "%d " BLU "TRACE " WHT "%s\n" RESET, ts, 0, msg); fflush(stdout); break;
    default:
      fprintf(stderr, CYN"%s " YEL "%d " RED "INVALID LEVEL " WHT "%s\n" RESET, ts, 0, msg); fflush(stderr); break;
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

void set_verbosity(char *v) {
  /* to upper case */
  for (char *p = v; *p != '\0'; ++p) {
    *p = toupper(*p);
  }

  if (equals(v, "OFF")) {
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
  }
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

  strftime(buf, 80, "%d-%m-%Y %H:%M:%S", info);
  return buf;
}
