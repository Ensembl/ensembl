
CC = gcc
CFLAGS = -c -g -O
IFLAGS = -I/usr/local/lib/glib/include


libsom.a : timetolease.o simpleobjectmanager.o memwatch.o
	ar -ru libsom.a $?


.c.o :
	$(CC) $(CFLAGS) $(IFLAGS) $<
