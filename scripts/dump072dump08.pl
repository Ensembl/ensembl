#!/usr/local/bin/perl



while( <> ) {
    s/(insert into contig .*)\)/$1,'international-dummy-id'\)/i;
    s/(insert into exon .*)\)/$1,1\)/i;
    s/(insert into feature .*)\)/$1,0.0,0\)/i;
    print;
}

