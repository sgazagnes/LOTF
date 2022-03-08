#ifndef LOGC_H
#define LOGC_H

#include <bitset>


#define MAX_LOG_MSG_SIZE 255

#define RED  "\x1B[91m"
#define GRN  "\x1B[92m"
#define YEL  "\x1B[93m"
#define BLU  "\x1B[94m"
#define MAG  "\x1B[95m"
#define CYN  "\x1B[96m"
#define WHT  "\x1B[97m"
#define RESET "\033[0m"

std::bitset<10> verbosity; // 0000000000, 1111111111

/* verbosity levels */
enum {
      ERROR_ = 0, TIME, INFO_, COLLECT, GRID, CONNECT, FIT, MERGE, TRKZ, TRKERROR 
};

void _error(const char* msg);
void error(const char* fmt, ...);
void _timing(const char* msg);
void timing(const char* fmt, ...);
void _info(const char* msg);
void info(const char* fmt, ...);
void dbgcollect(const char* fmt, ...);
void dbggrid(const char* fmt, ...);
void dbgconnect(const char* fmt, ...);
void dbgfit(const char* fmt, ...);
void dbgmerge(const char* fmt, ...);
void dbgtrkerror(const char* fmt, ...);
void dbgtrkz(const char* fmt, ...);

void _logc(int level, const char* msg);
void logc(int level, const char* fmt, ...);

void set_verbosity(bool v[10]);
bool equals(const char *a, const char *b);
char* timestamp(void);
#endif
