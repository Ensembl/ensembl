#!/usr/local/bin/perl
# 
# $Id$
# 
#

=head1 NAME

 sattelite_dbdump_bychr

=head1 SYNOPSIS

  This script generates a dump of an EnsEMBL satellite database for
  particular chromosome. Useful to create a small but fully functional
  ensembl installation, e.g. for a laptop. It needs access to an ensembl-lite
  database (for golden path etc.)

  (1) Needs to be called within a new directory where you want
     all the files to be written

  (2) with a user that is allowed to use mysqldump

  (3) needs to be run on the host that runs the daemon

  (4) Usage: 

       satellite_dbdump_bychr  -<dbtype> <dbinstance>
  
     e.g
  
       satellite_dbdump_bychr  -disease homo_sapiens_disease_110

     Known types are: family disease maps expression est # snp 

=head1 DESCRIPTION

This script generates a full dump of one or several EnsEMBL sattelite
database for a particular chromosome. Useful to create a small but fully
functional EnsEMBL db (e.g. laptop mini-mirror) 

Based on make_dbdumpk_bychr (which should be used for the core, embl and
EST database. embl is a problem still, however.

=cut

;

use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

my $workdir = `pwd`; chomp($workdir);
my $host = "localhost";
my $port   = '';
my $litedb = ''; # 'homo_sapiens_lite_110'; # force user to provide it
my $dbuser = 'ensadmin';
my $dbpass = undef;
my $module = 'Bio::EnsEMBL::DBSQL::DBAdaptor';
my $chr = 'chr21';                      # smaller than chr22
# my $lim;
my $mysql = 'mysql'; 
my $mysqldump = 'mysqldump'; # in $PATH we trust
#                  /mysql/current/bin/mysqldump

# satellites:
my $famdb;
my $diseaseb;
my $mapsdb;
my $expressiondb;
my $snpdb;
my $embldb;
my $estdb;
# end of satellites

&GetOptions( 
            'port:n'     => \$port,
            'litedb:s'   => \$litedb,
            'dbuser:s'   => \$dbuser,
            'dbpass:s'   => \$dbpass,
            'module:s'   => \$module,
            'chr:s'      => \$chr,
            'workdir:s'  => \$workdir,
            'limit:n'    => \$lim,
            'family:s' => \$famdb,
            'disease:s' => \$diseasedb,
            'maps:s' => \$mapsdb,
            'expression:s' => \$expressiondb,
            'snp:s' => \$snpdb,
            'embl:s' => \$embldb,
            'est:s' => \$estdb,
           );

die "need a litedb; use -litedb something " unless $litedb;
die "chromosome names should start with 'chr'" unless $chr =~ /^chr/;
my $pass_arg=""; $pass_arg="-p$dbpass" if $dbpass;

my $limit;
if ($lim) {
    $limit = "limit $lim";
}

&dump_family($famdb);
&dump_disease($diseasedb);
&dump_maps($mapsdb);
&dump_expression($expressiondb);
&dump_snp($snpdb);
&dump_embl($embldb);

&dump_est($estdb);

sub dump_family { 
    my ($satdb) = @_;
    return unless $satdb;

    dump_schema($satdb);

    my $sql;
    $sql = "
SELECT distinct f.* 
FROM $satdb.family f, $litedb.gene g
WHERE g.chr_name = '$chr'
  and g.family = f.id
  $limit
";
    dump_data($sql, $satdb, 'family' );

    $sql = "
SELECT distinct fm.* 
FROM $satdb.family_members fm, $satdb.family f, $litedb.gene g
WHERE g.chr_name = '$chr'
  and g.family = f.id
  and f.internal_id  = fm.family
  $limit
";
    dump_data($sql, $satdb, 'family_members' );
}                                       # family

sub dump_disease {
    my ($satdb) = @_;
    return unless $satdb;

    dump_schema($satdb);

# may need an ALTER TABLE gene ADD KEY(gene_symbol);
    my $sql;
    $sql = "
SELECT dg.*
FROM  $satdb.gene dg, 
      $litedb.gene lg, 
      $litedb.gene_xref lgx
WHERE lg.chr_name = '$chr' 
  AND lg.gene = lgx.gene 
  AND lgx.display_id = dg.gene_symbol
";
    dump_data($sql, $satdb, 'gene' );

    $sql = "
SELECT dd.*
FROM  $satdb.gene dg, 
      $satdb.disease dd,
      $litedb.gene lg, 
      $litedb.gene_xref lgx
WHERE lg.chr_name = '$chr' 
  AND lg.gene = lgx.gene 
  AND lgx.display_id = dg.gene_symbol
  AND dd.id = dg.id;
";
    dump_data($sql, $satdb, 'disease' );

# here's the sql to restrict the disease_index_*list, but they're so small
# it's really not worth the trouble. Left here in case anyone is interested
#     $sql = "
# SELECT ddl.*
# FROM  $satdb.gene dg, 
#       $satdb.disease_index_doclist ddl,
#       $litedb.gene lg, 
#       $litedb.gene_xref lgx
# WHERE lg.chr_name = '$chr' 
#   AND lg.gene = lgx.gene 
#   AND lgx.display_id = dg.gene_symbol
#   AND ddl.id  = dg.id
# ";

    foreach my $w ( qw(doc stop vector word) ) {
        my $table = "disease_index_${w}list";
        $sql = "select distinct * from $satdb.$table";
        dump_data($sql, $satdb, $table );
    }
}                                       # disease

sub dump_maps {
    my ($satdb) = @_;
    return unless $satdb;

    warn "ignoring non-RHdb markers !\n";
    dump_schema($satdb);

    my $chr_short = $chr;
    $chr_short =~ s/^chr//;

    my $sql;

    # the simple ones having a chromosome column:
    foreach my $table ( qw(ChromosomeBands CytogeneticMap RHMaps Fpc_Contig)) {
        $sql = "
SELECT * FROM $satdb.$table WHERE chromosome = '$chr_short'
";
        dump_data($sql, $satdb, $table );
    }

    $sql = "SELECT * FROM $satdb.Map";  # 4 rows
    dump_data($sql, $satdb, $table );

    # less simple ones that can both use the RHMaps table
    foreach my $table ( qw(Marker MarkerSynonym) ) {              
        $sql = "
SELECT distinct t.* 
FROM $satdb.$table t,
     $satdb.RHMaps r
WHERE t.marker=r.marker 
  AND r.chromosome = '$chr_short'
";
        dump_data($sql, $satdb, $table );
    }    

    # this one needs a join 
    $sql="
SELECT distinct cl.*
FROM $satdb.Fpc_Clone cl,
     $satdb.Fpc_Contig cg
WHERE cg.chromosome = '$chr_short'
  AND cl.contig_id = cg.contig_id
";
    dump_data($sql, $satdb, 'Fpc_Clone' );
}                                       # maps

sub dump_expression  {
    my ($satdb) = @_;
    return unless $satdb;

    warn "ignoring any non-ENSG aliases";
    my $dumpdir = "$workdir/$satdb";
    dump_schema($satdb);

    # small ones:
    foreach $table ( qw(key_word lib_key library source ) ) {
        $sql = "select distinct * from $satdb.$table";
        dump_data($sql, $satdb, $table);
    }
    # frequency                            ;
    # seqtag                               ;
    # seqtag_alias                         ;
    $sql = "
SELECT distinct sa.*
FROM $satdb.seqtag_alias sa, 
     $litedb.gene lg
WHERE sa.db_name = 'ensgene'
  AND sa.external_name =lg.name
  AND lg.chr_name = '$chr'
";
    dump_data($sql, $satdb, 'seqtag_alias');

    $sql = "
SELECT distinct  st.*
FROM  $satdb.seqtag st,
      $satdb.seqtag_alias sa, 
      $litedb.gene lg
WHERE sa.db_name = 'ensgene'
  AND sa.external_name =lg.name
  AND lg.chr_name = '$chr'
  AND st.seqtag_id = sa.seqtag_id
";
    dump_data($sql, $satdb, 'seqtag');

    $sql = "
SELECT distinct f.*
FROM  $satdb.frequency f,
      $satdb.seqtag_alias sa, 
      $litedb.gene lg
WHERE sa.db_name = 'ensgene'
  AND sa.external_name =lg.name
  AND lg.chr_name = '$chr'
  AND f.seqtag_id = sa.seqtag_id
";
    dump_data($sql, $satdb, 'frequency');
}                                       # expression

sub dump_snp  {
    my ($satdb) = @_;
    return unless $satdb;

    warn "ignoring any non-ENSG aliases";
    my $dumpdir = "$workdir/$satdb";
    dump_schema($satdb);

    my @small_ones = qw(Assay ContigHit Locus  Pop Resource Submitter);
    foreach my $table ( @small_ones ) { 
        $sql = "select distinct * from $satdb.$table";
        dump_data($sql, $satdb, $table);
    }

    #  RefSNP:
    $sql = "
SELECT distinct rs.*
FROM   $satdb.RefSNP rs, 
       $litedb.gene_snp lgs,
       $litedb.gene lg
WHERE  lg.chr_name = '$chr'
 AND   lg.gene = lgs.gene
 AND   lgs.refsnpid = rs.id 
";
    dump_data($sql, $satdb, 'RefSNP');
    
    #  SubSNP
    $sql = "
SELECT distinct ss.*
FROM   $satdb.SubSNP ss, 
       $litedb.gene_snp lgs,
       $litedb.gene lg
WHERE  lg.chr_name = '$chr'
 AND   lg.gene = lgs.gene
 AND   lgs.refsnpid = ss.refsnpid 
";


# (or should the last bit be ``lgs.refsnpid = ss.id'') ? 
    dump_data($sql, $satdb, 'SubSNP');

#  Hit        
    $sql = "
SELECT distinct h.*
FROM   $satdb.Hit h, 
       $litedb.gene_snp lgs,
       $litedb.gene lg
WHERE  lg.chr_name = '$chr'
 AND   lg.gene = lgs.gene
 AND   lgs.refsnpid = h.refsnpid 
";
    dump_data($sql, $satdb, 'Hit');

# ignore these (says Heikki):
#  Freq  
#  GPHit      
#  SubPop     

}                                       # snp


sub dump_embl  {
    my ($satdb) = @_;
    return unless $satdb;

    dump_schema($satdb);
    warn "This may take a while...\n";

    my @small_ones = qw(externalDB);
    foreach my $table ( @small_ones ) { 
        $sql = "select distinct * from $satdb.$table";
        dump_data($sql, $satdb, $table);
    }

    my $sql;

    $sql="
SELECT distinct cl.*
  FROM $litedb.gene lg,
       $satdb.clone cl
 WHERE lg.chr_name = '$chr'
   and lg.contig like concat(cl.id, '%') 
";
    dump_data($sql, $satdb, 'clone');

    $sql="
SELECT distinct ctg.*
  FROM $litedb.gene lg,
       $satdb.clone cl, 
       $satdb.contig ctg
 WHERE lg.chr_name = '$chr'
   and lg.contig like concat(cl.id, '%') 
   and ctg.clone = cl.internal_id
";
    dump_data($sql, $satdb, 'contig');
    $sql ="
SELECT distinct dna.*
  FROM $litedb.gene lg,
       $satdb.clone cl, 
       $satdb.contig ctg,
       $satdb.dna dna
 WHERE lg.chr_name = '$chr'
   and lg.contig like concat(cl.id, '%') 
   and ctg.clone = cl.internal_id
   and ctg.dna = dna.id 
";
     dump_data($sql, $satdb, 'dna');
    
    $sql="
SELECT DISTINCT e.*
  FROM $litedb.gene lg,
       $satdb.clone cl, 
       $satdb.contig ctg,
       $satdb.exon e
 WHERE lg.chr_name = '$chr'
   AND lg.contig like concat(cl.id, '%') 
   AND ctg.clone = cl.internal_id
   AND ctg.internal_id = e.contig
";
     dump_data($sql, $satdb, 'exon');

    $sql="
SELECT distinct et.*
  FROM $litedb.gene lg,
       $satdb.clone cl, 
       $satdb.contig ctg,
       $satdb.exon e,
       $satdb.exon_transcript et
 WHERE lg.chr_name = '$chr'
   and lg.contig like concat(cl.id, '%') 
   and ctg.clone = cl.internal_id
   and ctg.internal_id = e.contig
   and e.id = et.exon
";
    dump_data($sql, $satdb, 'exon_transcript');

    $sql="
SELECT distinct tsc.*
  FROM $litedb.gene lg,
       $satdb.clone cl, 
       $satdb.contig ctg,
       $satdb.exon e,
       $satdb.exon_transcript et,
       $satdb.transcript tsc
 WHERE lg.chr_name = '$chr'
   and lg.contig like concat(cl.id, '%') 
   and ctg.clone = cl.internal_id
   and ctg.internal_id = e.contig
   and e.id = et.exon
   and et.transcript=tsc.id
";
    dump_data($sql, $satdb, 'transcript');

    $sql="
SELECT distinct g.*
  FROM $litedb.gene lg,
       $satdb.clone cl, 
       $satdb.contig ctg,
       $satdb.exon e,
       $satdb.exon_transcript et,
       $satdb.transcript tsc,
       $satdb.gene  g
 WHERE lg.chr_name = '$chr'
   and lg.contig like concat(cl.id, '%') 
   and ctg.clone = cl.internal_id
   and ctg.internal_id = e.contig
   and e.id = et.exon
   and et.transcript=tsc.id
   and tsc.gene=g.id
";
    dump_data($sql, $satdb, 'gene');

    $sql="
SELECT distinct gt.*
  FROM $litedb.gene lg,
       $satdb.clone cl, 
       $satdb.contig ctg,
       $satdb.exon e,
       $satdb.exon_transcript et,
       $satdb.transcript tsc,
       $satdb.gene  g,
       $satdb.genetype  gt
 WHERE lg.chr_name = '$chr'
   and lg.contig like concat(cl.id, '%') 
   and ctg.clone = cl.internal_id
   and ctg.internal_id = e.contig
   and e.id = et.exon
   and et.transcript=tsc.id
   and tsc.gene=g.id
   and g.id = gt.gene_id
";
    dump_data($sql, $satdb, 'genetype');

    $sql="
SELECT distinct tl.*
  FROM $litedb.gene lg,
       $satdb.clone cl, 
       $satdb.contig ctg,
       $satdb.exon e,
       $satdb.exon_transcript et,
       $satdb.transcript tsc,
       $satdb.translation tl
 WHERE lg.chr_name = '$chr'
   and lg.contig like concat(cl.id, '%') 
   and ctg.clone = cl.internal_id
   and ctg.internal_id = e.contig
   and e.id = et.exon
   and et.transcript=tsc.id
   and tsc.translation = tl.id
";
    dump_data($sql, $satdb, 'translation');

    $sql="
SELECT distinct ox.*
  FROM $litedb.gene lg,
       $satdb.clone cl, 
       $satdb.contig ctg,
       $satdb.exon e,
       $satdb.exon_transcript et,
       $satdb.transcript tsc,
       $satdb.translation tl,
       $satdb.objectXref ox
 WHERE lg.chr_name = '$chr'
   and lg.contig like concat(cl.id, '%') 
   and ctg.clone = cl.internal_id
   and ctg.internal_id = e.contig
   and e.id = et.exon
   and et.transcript=tsc.id
   and tsc.translation = tl.id
   and tl.id = ox.ensembl_id
   and ox.ensembl_object_type = 'Translation'
";
    dump_data($sql, $satdb, 'objectXref');

    $sql="
SELECT distinct x.*
  FROM $litedb.gene lg,
       $satdb.clone cl, 
       $satdb.contig ctg,
       $satdb.exon e,
       $satdb.exon_transcript et,
       $satdb.transcript tsc,
       $satdb.translation tl,
       $satdb.objectXref ox,
       $satdb.Xref  x
 WHERE lg.chr_name = '$chr'
   and lg.contig like concat(cl.id, '%') 
   and ctg.clone = cl.internal_id
   and ctg.internal_id = e.contig
   and e.id = et.exon
   and et.transcript=tsc.id
   and tsc.translation = tl.id
   and tl.id = ox.ensembl_id
   and ox.ensembl_object_type = 'Translation'
   and ox.xrefId = x.xrefId
   and x.xrefId = 7036
";
    dump_data($sql, $satdb, 'Xref');
    return;
}                                       # embl


sub dump_est  {
    my ($satdb) = @_;
    return unless $satdb;
    warn "no written, doing nohting";
    return undef;
}                                       # est

sub dump_schema {
    my ($satdb) = @_;

    my $destdir = "$workdir/$satdb";
    my $destfile = "$satdb.sql";

    unless (-d $destdir) {
        mkdir $destdir, 0755 || die "mkdir $destdir: $!";
    }

    my $d = "$destdir/$destfile";

    warn "Dumping database schema of $satdb to $d\n";
    die "$d exists" if -s $d ;
    $command = "$mysqldump -u $dbuser $pass_arg -d $satdb > $d ";
    if ( system($command) ) {
        die "Error: ``$command'' ended with exit status $?";
    }
}

sub dump_data {
    my($sql, $satdb, $tablename) = @_;
    my ($destdir) = "$workdir/$satdb";
    my ($datfile)=  "$tablename.txt";

    unless (-d $destdir) {
        mkdir $destdir, 0755 || die "mkdir $destdir: $!";
    }
    
    $sql =~ s/\s+/ /g;
    
    my $cmd = "echo \"$sql\" | $mysql -q --batch -u $dbuser -p$dbpass $litedb > $destdir/$datfile";
    # warn "dumping: $cmd\n"; too verbose
    warn "dumping $tablename ...\n";

    if ( system($cmd) ) { 
        die "``$cmd'' exited with exit-status $?";
    }
}

## stuff below is not used (yet), since everything is done by plain SQL

## This comes from family-input.pl, and should at one point be put somewhere
## more central (the ones in EnsEMBL load modules etc. that are not relevant)
## Takes string that looks like
## "database=foo;host=bar;user=jsmith;passwd=secret", connects to mysql
## and return the handle
sub db_connect { 
    my ($dbcs) = @_;

    my %keyvals= split('[=;]', $dbcs);
    my $user=$keyvals{'user'};
    my $paw=$keyvals{'pass'};
#    $dbcs =~ s/user=[^;]+;?//g;
#    $dbcs =~ s/password=[^;]+;?//g;
# (mysql doesn't seem to mind the extra user/passwd values, leave them)

    my $dsn = "DBI:mysql:$dbcs";

    my $dbh=DBI->connect($dsn, $user, $paw) ||
      die "couldn't connect using dsn $dsn, user $user, password $paw:" 
         . $DBI::errstr;
    $dbh->{RaiseError}++;
    $dbh;
}                                       # db_connect

sub unique {
    
    my @unique;
    my %seen = ();
    foreach my $item (@_) {
	push(@unique,$item) unless $seen{$item}++;
    }
    return @unique;
}

sub get_inlist {
    my $string_flag = shift (@_);
    my $string;
    foreach my $element (@_) {
	if ($string_flag) {
	    $string .= "'$element',";
	}
	else {
	    $string .= "$element,";
	}
    }
    $string =~ s/,$//;
    return "($string)";
} 
