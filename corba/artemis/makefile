
IFLAGS = -I/usr/local/lib/glib/include -I/usr/local/include -I../popt -I. -I../simpleobjectman

LFLAGS = -L/usr/lib/mysql -Llib -L../timetolease

test-artseq-server : test-artseq-server.o artemis-mysql-impl.o artemis-exon-impl.o artemis-skels.o artemis-common.o
	$(CC) -o test-artseq-server test-artseq-server.o artemis-mysql-impl.o artemis-skels.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm -lresolv -lmysqlclient -lsom


test-artentry-server : test-artentry-server.o artemis-mysql-impl.o artemis-exon-impl.o artemis-qual-helper.o artemis-skels.o artemis-common.o
	$(CC) -o test-artentry-server test-artentry-server.o artemis-mysql-impl.o artemis-qual-helper.o artemis-exon-impl.o artemis-skels.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm  -lresolv -lmysqlclient -lsom


ensembl-artemis-server : ensembl-artemis-server.o artemis-mysql-impl.o artemis-exon-impl.o artemis-db-mysql.o artemis-qual-helper.o artemis-skels.o artemis-common.o
	$(CC) -o ensembl-artemis-server ensembl-artemis-server.o artemis-mysql-impl.o artemis-db-mysql.o artemis-qual-helper.o artemis-exon-impl.o artemis-skels.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm  -lresolv -lmysqlclient ../popt/.libs/libpopt.al -lsom

efs : ensembl-artemis-server.o artemis-mysql-impl.o artemis-exon-impl.o artemis-db-mysql.o artemis-qual-helper.o artemis-skels.o artemis-common.o
	$(CC) -o ensembl-artemis-server ensembl-artemis-server.o artemis-mysql-impl.o artemis-db-mysql.o artemis-qual-helper.o artemis-exon-impl.o artemis-skels.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm  -lresolv -lmysqlclient ../popt/.libs/libpopt.al -lsom -lefence

ensembl-artemis-client : ensembl-artemis-client.o artemis-stubs.o artemis-common.o
	$(CC) -o ensembl-artemis-client ensembl-artemis-client.o artemis-stubs.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm  -lresolv ../popt/.libs/libpopt.al -lsom


test-artseq-client : test-artseq-client.o artemis-stubs.o artemis-common.o
	$(CC) -o test-artseq-client test-artseq-client.o artemis-stubs.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm  -lresolv -lsom

test-artentry-client : test-artentry-client.o artemis-stubs.o artemis-common.o
	$(CC) -o test-artentry-client test-artentry-client.o artemis-stubs.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm -lresolv -lsom

.c.o :
	cc -c -g $(IFLAGS) $<






