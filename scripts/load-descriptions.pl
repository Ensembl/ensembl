#!/usr/local/bin/perl 

# $Id$ 

# script to find gene descriptions from swissprot, Tribe, or interpro,
# then load them into the gene_description table

# This needs Manu's corrections to the Swissprot parser (otherwise,
# multiple AC-lines will be wrong; second AC line overwrites the first
# one! E.g, make /work1/birney/mongin/src/bioperl-live part of PERLLIB)

my $Usage = "Usage: $0 -E connectstring -F connectstring swissprot-file ";

my $max_desc_len = 255;                 # max val for MySQL strings
my $sptr_prefix ="";
my $tribe_prefix = "[Family] ";
my $interpro_prefix = "[InterPro] ";
my $unknown_string = "unknown";

my $max = 100;

use DBI;
use Bio::SeqIO;
use strict;
use Getopt::Std;

use vars qw($opt_E $opt_F);
my $opts = 'E:F:';
getopts($opts) || die $Usage; 

my $spfile=shift;
die "no such file: $spfile" unless -s $spfile;
# e.g., /work1/birney/mongin/mapping_dev/data/data/total_10_00_sp.txt

my $ensdb_connect_string = 
  ($opt_E || 'database=ensembl080;host=ensrv5;user=ensadmin');

my $famdb_connect_string = 
  ($opt_F || 'database=family080;host=ecs1b;user=root');

#     my ($host, $user, $db) = ('ensrv5', 'ensro', 'ensembl080');
#     my ($host, $user, $db) = ('ensrv3', 'ensro', 'simon_oct07');

my $ensdb = db_connect($ensdb_connect_string) || die "couldn't connect";
my $famdb = db_connect($famdb_connect_string) || die "couldn't connect";

# use global vars, so we can use bind, should be faster;

my $sp_acc_q = "";
my $tribe_q = "";

sub open_cursors {
  $sp_acc_q = 
  "SELECT external_id 
   FROM genedblink x
   WHERE gene_id = ?
     AND external_db LIKE 'SP%'  
     ORDER BY external_db ASC, external_id DESC";
  $sp_acc_q = $ensdb->prepare($sp_acc_q) || die $ensdb->errstr;


  $tribe_q = 
  "SELECT description
   FROM family f, family_members fm
   WHERE f.internal_id = fm.family 
     AND fm.db_name = 'ENSEMBLGENE'
     AND fm.db_id = ?";
  $tribe_q = $famdb->prepare($tribe_q) || die $famdb->errstr;
  
}

my %sp_desc;                                 # ac->description hash

&readsp($spfile);
&open_cursors();
&main();

$ensdb->disconnect;
$famdb->disconnect;

sub getgenes {
# return qw(ENSG00000000351 ENSG00000000419 ENSG00000000457
# ENSG00000000460
#           ENSG00000000893);

# warn "test driving, doing just $max genes";
    my $q= "SELECT id from gene";
    $q = $ensdb->prepare($q) || die $ensdb->errstr;

    my @out;

my $n = 0;
    $q->execute || die $ensdb->errstr;
    while( my $rowhash = $q->fetchrow_hashref) {
        push(@out,$rowhash->{'id'});
# last if $n++ > $max; 
    }
    @out;
}                                       # getgenes

sub main {

    my @genes = getgenes;

    my $insert_q = 
      "INSERT INTO gene_description(gene_id, description)
       VALUES(?,?)";
    $insert_q = $ensdb->prepare($insert_q) || die  $ensdb->errstr;

    my ($desc, $pep, $score, $newscore);

    my $n =0;
    foreach my $gene ( @genes ) {
        $desc = find_desc($gene);
# warn $gene, $desc;
# last if $n++ > $max; 
        $insert_q->execute($gene, $desc) || die $ensdb->errstr;
    }
}                                       # main

sub find_sp_acc {
    my ($gene) = @_;

    $sp_acc_q->execute($gene) || die $ensdb->errstr;
    my ($pep) = $sp_acc_q->fetchrow;

warn "gene=$gene, spacc=$pep";

    return $pep;                        # may be empty
}

sub find_desc {
    my ($gene) = @_;
    
   ### First try SwissProt:
    my $sp_acc=find_sp_acc($gene);

    if ($sp_acc ne '') { 
        if (!defined($sp_desc{$sp_acc})) { 
            return "bogus";
        } else { 
warn "HOOOOOOOOrrraayay*******\n";
            return  $sptr_prefix .  $sp_desc{$sp_acc}; 
        }
    };

    $tribe_q->execute($gene) || die $famdb->errstr;
    my ($desc) = $tribe_q->fetchrow();

    if ($desc ne '' and $desc ne 'UNKNOWN') { 
warn 'tribe';
        return $tribe_prefix . $desc; 
    }
warn "unknown";
    return $unknown_string;
    ### forget interpro for now, try later if there's time.

}                                       # find_desc

#    ### try InterPro. This is hairy, since things have to come from
#    ### protein_feature in a round about way. This may be multiple too.
#    my @pepids = ();
#    foreach my $tsc ($gene->each_Transcript) {
#        push @pepids, $tsc->translation->id;
#    }
# 
#    ###
#    my @analysis_dbs = qw(Pfam-protein PRINTS EMOTIF);
#    
#    $q = 
#      "SELECT id.description
#       FROM  protein_feature pf,
#             analysis a,
#             interpro ix,
#             interpro_description id
#       WHERE
#             pf.translation in (".join(',', @pepids). ")
#         AND a.db in (".join(',', @analysis_dbs).") 
#         AND pf.analysis = a.id
#         AND pf.hid = ix.id
#         AND ix.interpro_ac = id.interpro_ac
#       ORDER BY pf.score DESC
#       LIMIT 5";
# 
#    $q = $self->_db_obj->prepare($q) || die $self->_db_obj->errstr;
#    
#    while( my ($d) = $q->fetchrow ) {
#        $desc .= "$d. ";
#    }   
# 
#     
#     if (
#    
# }

## slightly more general function, which takes string that looks like
## "database=foo;host=bar;user=jsmith;passwd=secret", connects to mysql
## and return the handle
sub db_connect { 
    my ($dbcs) = @_;

    my %keyvals= split('[=;]', $dbcs);
    my $user=$keyvals{'user'};
    my $paw=$keyvals{'password'};
#    $dbcs =~ s/user=[^;]+;?//g;
#    $dbcs =~ s/password=[^;]+;?//g;
# (mysql doesn't seem to mind the extra user/passwd values, leave them)

    my $dsn = "DBI:mysql:$dbcs";

    my $dbh=DBI->connect($dsn, $user, $paw) ||
      die "couldn't connect using dsn $dsn, user $user, password $paw:" 
         . $DBI::errstr;
    $dbh;
}


# fill big ac->description hash
sub readsp {
    my ($spfile) = @_;    
# warn "test driving, truncating after $max or so prtiens";
    
    my $sp  = Bio::SeqIO->new('-file' => $spfile,
                          '-format' => 'swiss') || die "couldn't open";
    # kindly cut-n-pasted from Manu's AltaVista dump script
    my $n=0;
    while ( my $seq = $sp->next_seq() ) {
#        my @refs;
#        my @comments;
        
        my $ac =  $seq->accession;
# warn $ac;         
#        my $ke = $seq->keywords;
        
        my $desc = $seq->desc;

        if (length($desc) > $max_desc_len) {
            warn "Description for $ac longer than $max_desc_len; truncating it\n";
            $desc = substr($desc, 0, $max_desc_len);
        }

        $sp_desc{$ac}=$desc;

#        $insert_q->execute($ac, $desc) || die $ensdb->errstr;
# last if $n++ > $max; 

#        foreach my $ref ( $seq->annotation->each_Reference ) {
#            push (@refs,$ref->comment());
#        }
#         foreach my $comment ( $seq->annotation->each_Comment ) {
#             push (@comments,$comment->text)
#         }
#        
#        my $tot_sp = "$ke|$desc|@refs|@comments";
#        
#        if (!defined $sp{$ac}) {
#            $sp{$ac} = [];
#        }
#        
#        $sp{$ac} = $tot_sp;
    }                                   # while
    undef;
}                                       # readsp

