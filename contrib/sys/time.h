/*
 * MSVC compatibility shim for <sys/time.h>
 *
 * Provides gettimeofday() and struct timeval for Windows builds.
 * On POSIX systems this header is never reached because the real
 * <sys/time.h> is found first.
 */

#ifndef COMPAT_SYS_TIME_H
#define COMPAT_SYS_TIME_H

#ifdef _MSC_VER

#ifndef NOMINMAX
#define NOMINMAX  /* Prevent windows.h from defining min/max macros */
#endif
#include <winsock2.h>  /* struct timeval */
#include <windows.h>

static __inline int gettimeofday(struct timeval *tp, void *tzp) {
    (void)tzp;
    FILETIME ft;
    unsigned __int64 tmp;
    GetSystemTimeAsFileTime(&ft);
    tmp = ((unsigned __int64)ft.dwHighDateTime << 32) | ft.dwLowDateTime;
    /* FILETIME epoch is 1601-01-01; convert to Unix epoch 1970-01-01 */
    tmp -= 116444736000000000ULL;
    tp->tv_sec  = (long)(tmp / 10000000ULL);
    tp->tv_usec = (long)((tmp % 10000000ULL) / 10);
    return 0;
}

#else
/* Non-MSVC: just include the real header */
#include_next <sys/time.h>
#endif

#endif /* COMPAT_SYS_TIME_H */
