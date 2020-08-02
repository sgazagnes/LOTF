#ifndef LOGC_H
#define LOGC_H


#define MAX_LOG_MSG_SIZE 255

#define RED  "\x1B[91m"
#define GRN  "\x1B[92m"
#define YEL  "\x1B[93m"
#define BLU  "\x1B[94m"
#define MAG  "\x1B[95m"
#define CYN  "\x1B[96m"
#define WHT  "\x1B[97m"
#define RESET "\033[0m"

/* verbosity levels */
enum {
  OFF = INT_MIN, ERROR = 0, WARN, TIMING, INFO, DEBUG, TRACE, ALL = INT_MAX
};

void _error(const char* msg);
void error(const char* fmt, ...);
void _warn(const char* msg);
void warn(const char* fmt, ...);
void _info(const char* msg);
void info(const char* fmt, ...);
void _debug(const char* msg);
void debug(const char* fmt, ...);
void _trace(const char* msg);
void trace(const char* fmt, ...);
void _timing(const char* msg);
void timing(const char* fmt, ...);

void _logc(int level, const char* msg);
void logc(int level, const char* fmt, ...);

void set_verbosity(char *v);
bool equals(const char *a, const char *b);
char* timestamp(void);
#endif
