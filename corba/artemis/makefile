
IFLAGS = -I/usr/local/lib/glib/include -I/usr/local/include -I/home/birney/downloads/ORBit-0.4.95/popt

LFLAGS = -L/usr/lib/mysql 

test-artseq-server : test-artseq-server.o artemis-mysql-impl.o artemis-exon-impl.o artemis-skels.o artemis-common.o
	$(CC) -o test-artseq-server test-artseq-server.o artemis-mysql-impl.o artemis-skels.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm -lnsl -lresolv -lmysqlclient 

test-artentry-server : test-artentry-server.o artemis-mysql-impl.o artemis-exon-impl.o artemis-qual-helper.o artemis-skels.o artemis-common.o
	$(CC) -o test-artentry-server test-artentry-server.o artemis-mysql-impl.o artemis-qual-helper.o artemis-exon-impl.o artemis-skels.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm -lnsl -lresolv -lmysqlclient -lefence


ensembl-artemis-server : ensembl-artemis-server.o artemis-mysql-impl.o artemis-exon-impl.o artemis-db-mysql.o artemis-qual-helper.o artemis-skels.o artemis-common.o
	$(CC) -o ensembl-artemis-server ensembl-artemis-server.o artemis-mysql-impl.o artemis-db-mysql.o artemis-qual-helper.o artemis-exon-impl.o artemis-skels.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm -lnsl -lresolv -lmysqlclient ../popt/.libs/libpopt.al

ensembl-artemis-client : ensembl-artemis-client.o artemis-stubs.o artemis-common.o
	$(CC) -o ensembl-artemis-client ensembl-artemis-client.o artemis-stubs.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm -lnsl -lresolv ../popt/.libs/libpopt.al


test-artseq-client : test-artseq-client.o artemis-stubs.o artemis-common.o
	$(CC) -o test-artseq-client test-artseq-client.o artemis-stubs.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm -lnsl -lresolv 

test-artentry-client : test-artentry-client.o artemis-stubs.o artemis-common.o
	$(CC) -o test-artentry-client test-artentry-client.o artemis-stubs.o artemis-common.o $(LFLAGS) -lIIOP -lORBit -lORBitutil -lglib -lm -lnsl -lresolv 

.c.o :
	cc -c -g $(IFLAGS) $<


