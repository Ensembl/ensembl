#!/usr/local/bin/perl -w
# $Header$
# 
# Script to create a tab-separated file (to be loaded into the
# gene_description table) from a file sp-descriptions.dat
# (swiss-prot/sptrembl -> description mapping, from running
# sp-descriptions.pl on a swissprot/trembl flatfile) and from the
# mapping.dat (ensg -> sp/trembl accno mapping, obtained from script
# sp-mapping-dump.sh). 


## list of regexps tested, in order of increasing preference, when
## deciding which description to use (used by sub compare_desc):

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my @word_order = qw(unknown hypothetical putative novel probable [0-9]{3} kDa fragment cdna protein);

my $Usage = "Usage: $0 sp-descriptions.dat mapping.dat > gene-descriptions.dat\n";

die $Usage if @ARGV == 0;

my $spdesc = shift;
#my $map = shift;

open SPDESC,$spdesc || die "$spdesc:$!";

my %sp_desc;

while (<SPDESC>) {
    chomp;
    my ($db, $acc, $desc) = split(/\t/);
    $sp_desc{"$db:$acc"}=$desc;
}

close SPDESC || die "$!";

my $host = 'ecs1d.sanger.ac.uk';
my $dbname = 'homo_sapiens_core_4_28';
my $dbuser = 'ensro';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor (-host => $host,
					     -user => $dbuser,
					     -dbname => $dbname );

my $sth = $db->prepare("SELECT tsc.translation_id, tsc.gene_id, xdb.db_name, x.dbprimary_id, ix.query_identity, ix.target_identity FROM transcript tsc, objectXref ox, Xref x, externalDB xdb, identityXref ix WHERE tsc.translation_id = ox.ensembl_id AND ox.xrefId = x.xrefId AND x.externalDBId = xdb.externalDBId AND xdb.db_name in ('SWISSPROT', 'SPTREMBL') AND ox.objectxrefId = ix.objectxrefId order by tsc.gene_id asc, xdb.db_name desc, x.dbprimary_id asc");

unless ($sth->execute()) {
  $db->throw("Failed execution of a select query");
}

my %gene_desc;

LINE:

while (my ($ensp, $ensg, $db, $acc, $qy_percid, $tg_percid) = $sth->fetchrow_array()) {

  my $desc;
  if ( defined($gene_desc{$ensg}) ) {
    my ($prevdb, $prev_desc, $prev_acc, $prev_qy_percid, $prev_tg_percid)  = @{$gene_desc{$ensg}}; 
#    if ($prev_qy_percid < $qy_percid) { 
#      $desc = $sp_desc{"$db:$acc"};
#      $gene_desc{$ensg} = [ $db, $desc, $acc, $qy_percid, $tg_percid ]; # kick out the SPTREMBL desc.
#      next LINE;
#    }
#    next LINE;
    if ($prevdb  eq 'SWISSPROT') {   
      next LINE;                  # nothing to change
    }
    if ($db  eq 'SWISSPROT') { 
      $desc = $sp_desc{"$db:$acc"};
      $gene_desc{$ensg} = [ $db, $desc, $acc, $qy_percid, $tg_percid ]; # kick out the SPTREMBL desc.
      next LINE;
    }
    
    if ($db  eq 'SPTREMBL' &&  $prevdb eq $db ) {   
      $desc = $sp_desc{"$db:$acc"}; 
      if ( &compare_desc($prev_desc, $desc) < 0 ) {
	# new desc is better
	# warn "new better: $desc (old was: $prev_desc)\n";
	$gene_desc{$ensg} = [ $db, $desc, $acc, $qy_percid, $tg_percid ];
	next LINE;
      } else {
	# warn "old better: $prev_desc (new is: $desc)\n";
	next LINE;
      }
      die "should not reach this point: only know SWISS-PROT and SPTREMBL";
    }
  } else {
    $desc = $sp_desc{"$db:$acc"};
    $gene_desc{$ensg} = [ $db, $desc, $acc, $qy_percid, $tg_percid ];
  }
}

#  now dump the stuff to stdout.

foreach my $ensg ( keys %gene_desc )  { 

  my ($db,$desc,$acc,$qy_percid,$tg_percid) = @{$gene_desc{$ensg}};

  ### final cleanup 
  ### get rid of the Rik mess:
  $_ = $desc;
  if (s/[0-9A-Z]{10}Rik protein[ \.]//g) {
    warn "throwing away: $desc\n";
  }
  s/^\s*\(Fragment\)\.?\s*$//g;
  s/^\s*\(\s*\)\s*$//g;
  ### add more as appropriate

#  print STDOUT "$ensg\t $_ [Source:$db;Acc:$acc,%qy:$qy_percid,\%tg:$tg_percid]\n" if $_; #  =~ /[a-z]/;???
  print STDOUT "$ensg\t $_ [Source:$db;Acc:$acc]\n" if $_; #  =~ /[a-z]/;???
}

#### following taken from ensembl-external/scripts/family-input.pl

sub compare_desc {
    my ($a, $b) = @_; 

    my ($am, $bm);
    foreach my $w (@word_order) {
        $am = ($a =~ /$w/i)?1:0;
        $bm = ($b =~ /$w/i)?1:0;

        if ($am  != $bm ) {             # ie, one matches, other doesn't
            if ( $am == 1 ) {           # first one worse than second 
                return -1;
            } else { 
                return 1; 
            }
        }
    }
    # still look same; base result on length: longer is better
    return length($a) <=> length($b);
}                                       # compare_desc
