#!/usr/local/bin/perl
# 
# $Id$
# 
#

die "Not yet ready for the MAIN branch; use Tbranch-ensembl-110!"

use Carp;

=head1 NAME

 dump_by_chr.pl

=head1 SYNOPSIS

  This script generates a dump of an EnsEMBL database + satellites for
  a particular chromosome. Useful to create a small but fully functional
  ensembl installation, e.g. for a laptop. It needs access to an ensembl-lite
  database (for golden path etc.). It is used for ensembl-pruned. It replaces 
  the  former make_dbdumpk_bychr.pl and satellite_dbdump_bychr.pl. 

  It will not overwrite previous dumps, so if one or more dumps failed,
  correct the bug, delete the files produced (if any), and simply run
  again. It will complain a bit, that is all. 

  (1) Needs to be called within a new directory where you want
     all the files to be written

  (2) with a user that is allowed to use mysqldump

  (3) may need to be run on the host that runs the daemon

  (4) Usage: 

       dump_by_chr.pl -litedb <litedbinstance>  -<dbtype> <dbinstance> 
           [ -<anotherdb> <anotherdb instance>  ... ]  [ -chr <a-chromosome> ] 
           -h host -u user -p password
     e.g
  
       dump_by_chr.pl  -disease homo_sapiens_disease_110 \n
            + host,user,passwd args

     Alternatively, dump all things in one go, relying on the consistent
     naming scheme, as 

       dump_by_chr.pl -litedb <litedbinstance> \
           -all 'homo_sapiens_\%s_110' [ -chr <a-chromosome> ]  \
            + host,user,passwd args

     In this case, the %s will be replaced by all the known database
     types, and all will be dumped. Is the theory.

     Known types are: (see source code)
  
   (5) When changing something in this script or dumping from a different
       database, run first with the -check option, to see if all indices
       are defined. Without it, the dump might take very long!

=head1 DESCRIPTION

This script generates a full dump of one or several EnsEMBL database (core
or satellite) for a particular chromosome. Useful to create a small but
fully functional EnsEMBL db (e.g. laptop mini-mirror). It is used for
ensembl-pruned.

Based on make_dbdumpk_bychr and satellite_dbdump_bychr (which it replaces).

=cut

;

use strict;
use Carp;
use Bio::EnsEMBL::DBLoader;
use Getopt::Long;

my $workdir = `pwd`; chomp($workdir);
my $host = "localhost";
my $port   = '';

my $user = 'ensadmin';
my $pass = undef;
# my $module = 'Bio::EnsEMBL::DBSQL::DBAdaptor';

my $litedb = ''; # 'homo_sapiens_lite_110'; # force user to provide it
my $dumplite = 0;                       # we always specify lite, but it
                                        # doesn't always need dumping

my $chr = 'chr21';                      # smaller than chr22

my $static_golden_path_type = 'UCSC';
warn "using a hard-coded \$static_golden_path_type being $static_golden_path_type! \n";

my $template = undef;                   # when dumping all

my $check_sql;
my $lim;

my $mysql = 'mysql'; 
# my $mysqldump = 'mysqldump'; # in $PATH we trust.
# No we don't; the above is version 8.10, which is too strict. Use an
# older version:
my $mysqldump = '/nfs/croc/michele/mysql/bin/mysqldump';
#                  /mysql/current/bin/mysqldump

# the db's that have been implemented:
my $coredb;
my $famdb;
my $diseasedb;
my $mapsdb;
my $expressiondb;
my $snpdb;
my $embldb;
my $estdb;
my $mousedb;
# end of db's

&GetOptions( 
#             'port:n'     => \$port,
            'litedb:s'   => \$litedb, # for referencing only
            'dumplite'  ,  => \$dumplite,  # specify if lite should dumped
            'h:s'   => \$host,
            'u:s'   => \$user,
            'p:s'   => \$pass,
#            'module:s'   => \$module,
            'chr:s'      => \$chr,
            'workdir:s'  => \$workdir,
            'limit:n'    => \$lim,      # no more than so many rows (debugging)
            'core:s' => \$coredb,
            'family:s' => \$famdb,
            'disease:s' => \$diseasedb,
            'maps:s' => \$mapsdb,
            'expression:s' => \$expressiondb,
            'snp:s' => \$snpdb,
            'embl:s' => \$embldb,
            'est:s' => \$estdb,
            'mouse:s' => \$mousedb,
            'all:s' => \$template,  #dump all known ones, using $template
            'check' => \$check_sql, # does not really dump
           );

die "need a litedb; use -litedb something " unless $litedb;
die "chromosome names should start with 'chr'" unless $chr =~ /^chr/;
my $pass_arg=""; $pass_arg="-p$pass" if $pass;

my $limit;
if ($lim) {
    $limit = "limit $lim";
}

if ($template) {                        # 
    # dump all, using the naming template
    my $db;

    &dump_lite( template_fill($template, 'lite'));
    &dump_core( template_fill($template, 'core'));
    &dump_family( template_fill($template, 'family'));
    &dump_disease( template_fill($template, 'disease'));
    &dump_maps( template_fill($template, 'maps'));
    &dump_expression( template_fill($template, 'expression'));
    &dump_snp( template_fill($template, 'snp'));
    &dump_embl( template_fill($template, 'embl'));
    &dump_est( template_fill($template, 'est'));
    &dump_mouse( template_fill($template, 'mouse'));

} else { 
    # dumping just the odd one. 
    &dump_lite($litedb) if $dumplite;
    &dump_core($coredb);
    &dump_family($famdb);
    &dump_disease($diseasedb);
    &dump_maps($mapsdb);
    &dump_expression($expressiondb);
    &dump_snp($snpdb);
    &dump_embl($embldb);
    &dump_est($estdb);
    &dump_mouse($mousedb);
}

warn "done\n";
exit(0);

sub template_fill {
    my ($tpl,  $actual) = @_;
    confess "template should look like 'prefix\%ssuffix'\n"
      unless $template =~ /%s/;
    my $result = $tpl;
    $result =~ s/%s/$actual/;
    $result;
}

sub dump_core { 
### replacement for make_dump. First commit on branch-ensembl-110, cause
### that's stable. 

    my ($satdb) = @_;
    return unless $satdb;

    dump_schema($satdb);
    warn "This may take a while...\n";

    my $sql;

    # first the small ones in full
    my @tables = qw(analysis analysisprocess chromosome externalDB meta
                       species interpro interpro_description);
    
    foreach my $table ( @tables ) { 
        $sql = "
SELECT distinct complete_table.* 
FROM $satdb.$table complete_table
";
        dump_data($sql, $satdb, $table);
    }

    # those that have a chr_name:
    @tables = qw(karyotype contig_landmarkMarker static_golden_path);
    foreach my $table ( @tables ) { 
        $sql = "
SELECT distinct t.*
FROM $satdb.$table t 
WHERE chr_name = '$chr'"
          ;
        dump_data($sql, $satdb, $table);
    }

    # ... or chrname (sigh);

    @tables = qw(map_density);
    foreach my $table ( @tables ) { 
        $sql = "
SELECT distinct t.*
FROM $satdb.$table t 
WHERE chrname = '$chr'";
        dump_data($sql, $satdb, $table);
    }

    $sql="
SELECT distinct ctg.*
  FROM $satdb.static_golden_path sgp,
       $satdb.contig ctg
 WHERE sgp.chr_name = '$chr'
   AND sgp.type = '$static_golden_path_type'
   AND sgp.raw_id = ctg.internal_id
";
    dump_data($sql, $satdb, 'contig');

    $sql="
SELECT distinct cl.*
  FROM $satdb.static_golden_path sgp,
       $satdb.contig ctg,
       $satdb.clone cl
 WHERE sgp.chr_name = '$chr'
   AND sgp.type = '$static_golden_path_type'
   AND sgp.raw_id = ctg.internal_id
   AND ctg.clone = cl.internal_id   
";
    dump_data($sql, $satdb, 'clone');

    $sql="
SELECT distinct d.* 
  FROM $satdb.static_golden_path sgp,
       $satdb.contig ctg,
       $satdb.dna d
 WHERE sgp.chr_name = '$chr'
   AND sgp.type = '$static_golden_path_type'
   AND sgp.raw_id = ctg.internal_id
   AND ctg.dna = d.id   
";
    dump_data($sql, $satdb, 'dna');

    foreach my $table ( qw(feature repeat_feature) ) { 
        $sql="
SELECT distinct f.*
  FROM $satdb.static_golden_path sgp,
       $satdb.$table f
 WHERE sgp.chr_name = '$chr'
   AND sgp.type = '$static_golden_path_type'
   AND sgp.raw_id = f.contig
";
        dump_data($sql, $satdb, $table);
    }

    $sql="
SELECT distinct fsf.* 
  FROM $satdb.static_golden_path sgp,
       $satdb.feature f, 
       $satdb.fset_feature fsf
 WHERE sgp.chr_name = '$chr'
   AND sgp.type = '$static_golden_path_type'
   AND sgp.raw_id = f.contig
   AND fsf.feature = f.id
";
    dump_data($sql, $satdb, 'fset_feature');

    $sql="
SELECT distinct fs.*
  FROM $satdb.static_golden_path sgp,
       $satdb.feature f, 
       $satdb.fset fs,
       $satdb.fset_feature fsf
 WHERE sgp.chr_name = '$chr'
   AND sgp.type = '$static_golden_path_type'
   AND sgp.raw_id = f.contig
   AND fsf.feature = f.id
   AND fsf.fset = fs.id
";
    dump_data($sql, $satdb, 'fset');

    $sql="
SELECT distinct e.*
  FROM $satdb.static_golden_path sgp,
       $satdb.exon e
 WHERE sgp.chr_name = '$chr'
   AND sgp.type = '$static_golden_path_type'
   AND sgp.raw_id = e.contig
";
    dump_data($sql, $satdb, 'exon');

    $sql="
SELECT distinct et.*
  FROM $litedb.gene_exon ge,
       $satdb.exon_transcript et
 WHERE ge.chr_name = '$chr'
   AND ge.exon = et.exon
";
    dump_data($sql, $satdb, 'exon_transcript');

    $sql="
SELECT distinct tsc.*
  FROM $litedb.gene_exon ge,
       $satdb.exon_transcript et,
       $satdb.transcript tsc
 WHERE ge.chr_name = '$chr'
       ge.exon = et.exon
   and et.transcript = tsc.id
";
    dump_data($sql, $satdb, 'transcript');

    $sql="
SELECT distinct g.*
  FROM $litedb.gene lg,
       $satdb.gene g
 WHERE lg.chr_name = '$chr' 
   AND lg.name = g.id
";
    dump_data($sql, $satdb, 'gene');

    $sql="
SELECT distinct trl.*
  FROM $litedb.gene lg,
       $litedb.gene_prot lgp,
       $satdb.translation trl
 WHERE lg.chr_name = '$chr'
   AND lg.gene = lgp.gene
   AND lgp.translation = trl.id
";
    dump_data($sql, $satdb, 'translation');

    $sql="
SELECT distinct pf.*
  FROM $litedb.gene lg,
       $litedb.gene_prot lgp,
       $satdb.protein_feature pf
 WHERE lg.chr_name = '$chr'
   AND lg.gene = lgp.gene
   AND lgp.translation = pf.translation
";
    dump_data($sql, $satdb, 'protein_feature');

    $sql="
SELECT distinct gt.* 
FROM   $satdb.genetype gt, 
       $litedb.gene lg
WHERE  lg.chr_name = '$chr'
  AND  lg.name = gt.gene_id
";
    dump_data($sql, $satdb, 'genetype');

    $sql="
SELECT distinct gd.* 
FROM   $satdb.gene_description gd, 
       $litedb.gene lg
WHERE  lg.chr_name = '$chr'
  AND  lg.name = gd.gene_id
";
    dump_data($sql, $satdb, 'gene_description');

    $sql="
SELECT distinct sf.*
FROM   $litedb.gene_exon ge, 
       $satdb.supporting_feature sf
WHERE  ge.chr_name = '$chr'
  AND  sf.exon  = ge.exon
";
    dump_data($sql, $satdb, 'supporting_feature');

    $sql="
SELECT distinct ox.* 
FROM   $satdb.objectXref ox, 
       $litedb.gene_prot lgp, 
       $litedb.gene lg
WHERE  lg.chr_name = '$chr'
  AND  lg.gene = lgp.gene
  AND  ox.ensembl_id = lgp.translation
";
    dump_data($sql, $satdb, 'objectXref');

    $sql="
SELECT distinct x.* 
FROM   $satdb.Xref x,
       $satdb.objectXref ox, 
       $litedb.gene_prot lgp, 
       $litedb.gene lg
WHERE  lg.chr_name = '$chr'
  AND  lg.gene = lgp.gene
  AND  ox.ensembl_id = lgp.translation
  AND  x.xrefId = ox.xrefId
";
    dump_data($sql, $satdb, 'Xref');
    
    $sql="
SELECT distinct xs.* 
FROM   $satdb.externalSynonym xs,
       $satdb.objectXref ox, 
       $litedb.gene_prot lgp, 
       $litedb.gene lg
WHERE  lg.chr_name = '$chr'
  AND  lg.gene = lgp.gene
  AND  ox.ensembl_id = lgp.translation
  AND  xs.xrefId = ox.xrefId
";
    dump_data($sql, $satdb, 'externalSynonym');

    $sql="
SELECT distinct ix.* 
FROM  $satdb.identityXref ix,
      $satdb.objectXref ox, 
      $litedb.gene_prot lgp, 
      $litedb.gene lg
WHERE lg.chr_name = '$chr'
 AND  lg.gene = lgp.gene
 AND  ox.ensembl_id = lgp.translation
 AND  ix.objectXrefId = ox.objectXrefId
";
    dump_data($sql, $satdb, 'identityXref');

}                                       # dump_core

sub dump_lite {
# needed for a few things
    my ($satdb) = @_;
    return unless $satdb;
    
    dump_schema($satdb);

    my $sql;

    # tables that have chr_name:
    foreach my $table ( qw(gene gene_exon karyotype location snp) ) {
        $sql="
SELECT distinct t.*
FROM   $litedb.$table t
WHERE  t.chr_name = '$chr' 
";
        dump_data($sql, $satdb, $table );
    }


    # tables to be linked to the gene table:
    foreach my $table ( qw(gene_disease gene_prot gene_snp gene_xref) ) {
        $sql="
SELECT distinct t.*
FROM   $litedb.$table t, $litedb.gene g
WHERE  g.chr_name = '$chr' 
  AND  t.gene = g.gene
";
        dump_data($sql, $satdb, $table );
    }
}                                       # dump_lite

sub dump_family { 
    my ($satdb) = @_;
    return unless $satdb;

    dump_schema($satdb);

    my $sql;
    $sql = "
SELECT distinct f.* 
FROM $satdb.family f, $litedb.gene g
WHERE g.chr_name = '$chr'
  AND g.family = f.id
  $limit
";
    dump_data($sql, $satdb, 'family' );

    $sql = "
SELECT distinct fm.* 
FROM $satdb.family_members fm, $satdb.family f, $litedb.gene g
WHERE g.chr_name = '$chr'
  AND g.family = f.id
  AND f.internal_id  = fm.family
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
SELECT distinct dg.*
FROM  $satdb.gene dg, 
      $litedb.gene lg, 
      $litedb.gene_xref lgx
WHERE lg.chr_name = '$chr' 
  AND lg.gene = lgx.gene 
  AND lgx.display_id = dg.gene_symbol
";
    dump_data($sql, $satdb, 'gene' );

    $sql = "
SELECT distinct dd.*
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
        $sql = "
SELECT distinct complete_table.* 
FROM $satdb.$table complete_table
";
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
SELECT t.* 
FROM $satdb.$table t
WHERE chromosome = '$chr_short'
";
        dump_data($sql, $satdb, $table );
    }

    $sql = "
SELECT complete_table.* 
FROM $satdb.Map complete_table
";  # 4 rows
    dump_data($sql, $satdb, 'Map' );

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
    dump_schema($satdb);

    my $sql;
    # small ones:
    foreach my $table ( qw(key_word lib_key library source ) ) {
        $sql = "
SELECT distinct complete_table.* 
FROM $satdb.$table complete_table
";
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
    dump_schema($satdb);

    my $sql;

    my @small_ones = qw(Assay ContigHit Locus  Pop Resource Submitter);
    foreach my $table ( @small_ones ) { 
        $sql = "
SELECT distinct complete_table.* 
FROM $satdb.$table complete_table
";
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

# GPHit (i.e., those on the golden path); denormalized table used by contigview
    $sql="
SELECT distinct h.*
FROM   $satdb.Hit h, 
       $litedb.gene_snp lgs,
       $litedb.gene lg
WHERE  lg.chr_name = '$chr'
 AND   lg.gene = lgs.gene
 AND   lgs.refsnpid = h.refsnpid 
";
    dump_data($sql, $satdb, 'GPHit');


# ignore these (says Heikki):
#  Freq  
#  SubPop     

}                                       # snp


sub dump_embl  {
    my ($satdb) = @_;
    return unless $satdb;

    dump_schema($satdb);
    warn "This may take a while...\n";

    my $sql;

    my @small_ones = qw(externalDB);
    foreach my $table ( @small_ones ) { 
        $sql = "
SELECT DISTINCT complete_table.*
FROM $satdb.$table complete_table
";
        dump_data($sql, $satdb, $table);
    }

    $sql="
SELECT distinct cl.*
  FROM $litedb.gene lg,
       $satdb.clone cl
 WHERE lg.chr_name = '$chr'
   AND lg.contig like concat(cl.id, '%') 
";
    dump_data($sql, $satdb, 'clone');

    $sql="
SELECT distinct ctg.*
  FROM $litedb.gene lg,
       $satdb.clone cl, 
       $satdb.contig ctg
 WHERE lg.chr_name = '$chr'
   AND lg.contig like concat(cl.id, '%') 
   AND ctg.clone = cl.internal_id
";
    dump_data($sql, $satdb, 'contig');

    $sql ="
SELECT distinct dna.*
  FROM $litedb.gene lg,
       $satdb.clone cl, 
       $satdb.contig ctg,
       $satdb.dna dna
 WHERE lg.chr_name = '$chr'
   AND lg.contig like concat(cl.id, '%') 
   AND ctg.clone = cl.internal_id
   AND ctg.dna = dna.id 
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
   AND lg.contig like concat(cl.id, '%') 
   AND ctg.clone = cl.internal_id
   AND ctg.internal_id = e.contig
   AND e.id = et.exon
   AND e.sticky_rank = 1
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
   AND lg.contig like concat(cl.id, '%') 
   AND ctg.clone = cl.internal_id
   AND ctg.internal_id = e.contig
   AND e.id = et.exon
   AND e.sticky_rank = 1
   AND et.transcript=tsc.id
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
   AND lg.contig like concat(cl.id, '%') 
   AND ctg.clone = cl.internal_id
   AND ctg.internal_id = e.contig
   AND e.id = et.exon
   AND e.sticky_rank = 1
   AND et.transcript=tsc.id
   AND tsc.gene=g.id
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
   AND lg.contig like concat(cl.id, '%') 
   AND ctg.clone = cl.internal_id
   AND ctg.internal_id = e.contig
   AND e.id = et.exon
   AND e.sticky_rank = 1
   AND et.transcript=tsc.id
   AND tsc.gene=g.id
   AND g.id = gt.gene_id
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
   AND lg.contig like concat(cl.id, '%') 
   AND ctg.clone = cl.internal_id
   AND ctg.internal_id = e.contig
   AND e.id = et.exon
   AND e.sticky_rank = 1
   AND et.transcript=tsc.id
   AND tsc.translation = tl.id
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
   AND lg.contig like concat(cl.id, '%') 
   AND ctg.clone = cl.internal_id
   AND ctg.internal_id = e.contig
   AND e.id = et.exon
   AND e.sticky_rank = 1
   AND et.transcript=tsc.id
   AND tsc.translation = tl.id
   AND tl.id = ox.ensembl_id
   AND ox.ensembl_object_type = 'Translation'
";
    dump_data($sql, $satdb, 'objectXref');

    # this one is a total dog ... I need the LIKE statement, given the
    # strange embl id mapping ...  I suppose I should create a temp table,
    # index it, then do it on that. Left as an exercise :-(
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
   AND lg.contig like concat(cl.id, '%') 
   AND cl.internal_id = ctg.clone 
   AND ctg.internal_id = e.contig
   AND e.sticky_rank = 1
   AND e.id = et.exon
   AND et.transcript=tsc.id
   AND tsc.translation = tl.id
   AND tl.id = ox.ensembl_id
   AND ox.ensembl_object_type = 'Translation'
   AND ox.xrefId = x.xrefId
";

    warn "**** Not dumping Xref, too slow; rewrite using a temporary table";
    # then just re-run
#    dump_data($sql, $satdb, 'Xref');
    return;
}                                       # dump_embl


sub dump_est  {
    my ($satdb) = @_;
    return unless $satdb;
    dump_schema($satdb);
    my $sql;

    # small tables:
    foreach my $table ( qw(analysis analysisprocess) )  {
        $sql="
SELECT complete_table.* 
FROM $satdb.$table complete_table
";
        dump_data($sql, $satdb, $table);
    }

    $sql="
SELECT distinct c.*
FROM   $satdb.contig  c, $litedb.gene g
WHERE  g.chr_name = '$chr' 
  AND  g.contig = c.id
";
    dump_data($sql, $satdb, 'contig');

    $sql="
SELECT distinct f.*
FROM   $satdb.contig  c, 
       $satdb.feature f,
       $litedb.gene g
WHERE  g.chr_name = '$chr' 
  AND  g.contig = c.id
  AND  f.contig = c.internal_id
";    
    dump_data($sql, $satdb, 'feature');
    return undef;
}                                       # dump_est

sub dump_mouse  {
    my ($satdb)=@_;
    return unless $satdb;
    &dump_est($satdb);
}

sub dump_schema {
    my ($satdb) = @_;
    
    my $destdir = "$workdir/$satdb";
    my $destfile = "$satdb.sql";

    warn "Dumping/checking schema of $satdb\n";

    unless (-d $destdir) {
        mkdir $destdir, 0755 || die "mkdir $destdir: $!";
    }

    my $cmd = "$mysqldump -u $user -h $host $pass_arg -d $satdb";

    if ($check_sql) {
        # just see if it's there
        my $out = `$cmd`;

        if ( !$out || length($out) < 100  or $?  ) { 
            warn "``$cmd'' exited with no or little  output and exit-status $?\n";
            return;
        } 
        warn "schema dumps OK\n";
        return;
    }        

    # else: really dump it 
    my $d = "$destdir/$destfile";
    $cmd .= " > $d";
    warn "Dumping schema to $d\n";
    if  ( -s $d  )  { 
        warn "$d exists; not dumping it";
        return; 
    }
    
    if ( system($cmd) ) {
        die "Error: ``$cmd'' ended with exit status $?";
    }
}                                       # dump_schema

sub dump_data {
    my($sql, $satdb, $tablename) = @_;
    my ($destdir) = "$workdir/$satdb";
    my ($datfile)=  "$tablename.txt";

    if ($check_sql) { 
        check_sql($satdb, $sql, $tablename);
        return;
    }

    unless (-d $destdir) {
        mkdir $destdir, 0755 || warn "mkdir $destdir: $!";
    }

    my $outfile="$destdir/$datfile";    ;
    if (-s $outfile) {
        warn "File $outfile exists and not empty; not dumping this!";
        return;
    }

    $sql =~ s/\s+/ /g;
    
    my $cmd = "echo \"$sql\" | $mysql -N -q --batch -h $host -u $user $pass_arg $litedb > $outfile";
    # warn "dumping: $cmd\n"; too verbose
    warn "dumping $tablename ...\n";

    if ( system($cmd) ) { 
        die "``$cmd'' exited with exit-status $?";
    }
}                                       # dump_data

sub check_sql {
    my($satdb, $sql, $tablename)=@_;
    $sql =~ s/\s+/ /g;
    $sql =~ s/\bdistinct\b/ /gi;        # would last too long
    
    warn "checking $tablename ...\n";

    # see if all indices are there:
    my $cmd = "echo \"explain $sql\" | $mysql -N -q --batch -h $host -u $user $pass_arg $litedb";
    my $out = `$cmd`;

    if (!$out or $? ) {
        warn "ERROR: ``$cmd'' exited without output and exit-status $?\n";
        return;
    }

    if ( $out =~ /\bALL\b/ ) {
        carp "$sql: missing/not using index:\n$out\n";
    } else { 
        warn "query for $tablename OK\n";
    } 

    # now see if the query is not empty (due to a e.g. a wrong join
    # condition):
    $cmd = "echo \"$sql limit 1\" | $mysql -N -q --batch -h $host -u $user $pass_arg $litedb ";
    $out=`$cmd`;
    
    if ( !$out or length($out) < 8 or $?  ) { 
        warn "``$cmd'' exited with no or little output and exit-status $?\n";
    } 

    return;
}                                       # check_sql

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
