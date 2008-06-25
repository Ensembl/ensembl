#!/usr/local/bin/perl

# Dump variation information to an XML file for indexing by the EBI's search engine.
#
# To copy files to the EBI so that they can be picked up:
# scp *.xml.gz glenn@puffin.ebi.ac.uk:xml/
#
# Email eb-eye@ebi.ac.uk after copying so the files can be indexed.

package ebi_search_dump;

use strict;
use DBI;
use Carp;
use Getopt::Long;
use IO::Zlib;

use Data::Dumper;
use HTML::Entities;

my (
    $host,    $user,        $pass,   $port,     $species, $ind,
    $release, $max_entries, $nogzip, $parallel, $dir,     $inifile
);

my %rHash = map { $_ } @ARGV;
if ( $inifile = $rHash{'-inifile'} ) {
    my $icontent = `cat $inifile`;
    warn $icontent;
    eval $icontent;
}

GetOptions(
    "host=s",        \$host,        "port=i",    \$port,
    "user=s",        \$user,        "pass=s",    \$pass,
    "species=s",     \$species,     "release=s", \$release,
    "index=s",       \$ind,         "nogzip!",   \$nogzip,
    "max_entries=i", \$max_entries, "parallel",  \$parallel,
    "dir=s",         \$dir,         "help",      \&usage,
    "inifile=s",     \$inifile,
);

$species ||= 'ALL';
$ind     ||= 'ALL';
$dir     ||= ".";
$release ||= 'LATEST';

usage() and exit unless ( $host && $port && $user );

my $entry_count;
my $global_start_time = time;
my $total             = 0;

my $fh;
## HACK 1 - if the INDEX is set to all grab all dumper methods...
my @indexes = split ',', $ind;
@indexes = map { /dump(\w+)/ ? $1 : () } keys %ebi_search_dump::
  if $ind eq 'ALL';

warn Dumper \@indexes;

my $dbHash = get_databases();
warn Dumper $dbHash;

#warn Dumper $dbcHash;

foreach my $species ( sort keys %$dbHash ) {
    foreach my $index (@indexes) {
        my $function = "dump$index";
        no strict "refs";
        &$function( $species, $dbHash->{$species} );
    }
}

print_time($global_start_time);
warn " Dumped $total entries ...\n";

# -------------------------------------------------------------------------------

sub text_month {

    my $m = shift;

    my @months = qw[JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC];

    return $months[$m];

}

# -------------------------------------------------------------------------------

sub print_time {

    my $start = shift;

    my $t = time - $start;
    my $s = $t % 60;
    $t = ( $t - $s ) / 60;
    my $m = $t % 60;
    $t = ( $t - $m ) / 60;
    my $h = $t % 60;

    print "Time taken: " . $h . "h " . $m . "m " . $s . "s\n";

}

#------------------------------------------------------------------------------------------------
sub usage {
    print <<EOF; exit(0);

Usage: perl $0 <options>

  -host         Database host to connect to. Defaults to ens-staging.

  -port         Database port to connect to. Defaults to 3306.

  -species      Species name. Defaults to ALL.

  -index        Index to create. Defaults to ALL.

  -release      Release of the database to dump. Defaults to 'latest'.

  -user         Database username. Defaults to ensro.

  -pass         Password for user.

  -dir          Directory to write output to. Defaults to /lustre/scratch1/ensembl/gp1/xml.

  -nogzip       Don't compress output as it's written.

  -max_entries  Only dump this many entries for testing.

  -help         This message.

  -inifile      First take the arguments from this file. Then overwrite with what is provided in the command line

EOF

}

sub get_databases {

    my ( $dbHash, $dbcHash );
    my $dsn = "DBI:mysql:host=$host";
    $dsn .= ";port=$port" if ($port);

    my $db = DBI->connect( $dsn, $user, $pass );

    #    warn "DSN: $dsn";

    my @dbnames =
      map { $_->[0] } @{ $db->selectall_arrayref("show databases") };

    $db->disconnect();

    my $latest_release = 0;
    my ( $db_species, $db_release, $db_type );
    my $compara_hash;
    for my $dbname (@dbnames) {
        if ( ( $db_species, $db_type, $db_release ) =
            $dbname =~ /^([a-z]+_[a-z]+)_([a-z]+)_(\d+)_\w+$/ )
        {

            next if ( $species ne 'ALL' ) && ( $db_species ne $species );
            $latest_release = $db_release if ( $db_release > $latest_release );
            $dbHash->{$db_species}->{$db_type}->{$db_release} = $dbname;

        }
        if ( ($db_release) = $dbname =~ /ensembl_compara_(\d+)/ ) {

            #N.B Re:COMAPARA for now using
            #ensembl_compara_VERSION. Others will follow
            $compara_hash->{$db_release} = $dbname;
        }

    }

    map { $dbHash->{$_}->{'compara'} = $compara_hash } keys %$dbHash;
    $release = $latest_release if ( $release eq 'LATEST' );

    return $dbHash;

}

sub footer {
    my ($ecount) = @_;
    p("</entries>");
    p("<entry_count>$ecount</entry_count>");
    p("</database>");

    print "Dumped $ecount entries\n";
    if ($nogzip) {
        close(FILE) or die $!;
    }
    else {
        $fh->close();
    }
    $total += $ecount;
}

sub header {
    my ( $dbname, $dbspecies, $dbtype ) = @_;

    p("<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>");
    p("<!DOCTYPE database [ <!ENTITY auml \"&#228;\">]>");
    p("<database>");
    p("<name>$dbname</name>");
    p("<description>Ensembl $dbspecies $dbtype database</description>");
    p("<release>$release</release>");
    p("");
    p("<entries>");
}

sub p {
    my ($str) = @_;

    # TODO - encoding
    $str .= "\n";
    if ($nogzip) {
        print FILE $str or die "Can't write to file ", $!;
    }
    else {
        print $fh $str or die "Can't write string: $str";
    }
}

sub format_date {
    my $t = shift;

    my ( $y, $m, $d, $ss, $mm, $hh ) = ( localtime($t) )[ 5, 4, 3, 0, 1, 2 ];
    $y += 1900;
    $d = "0" . $d if ( $d < 10 );
    my $mm = text_month($m);
    return "$d-$mm-$y";
}

sub format_datetime {
    my $t = shift;

    my ( $y, $m, $d, $ss, $mm, $hh ) = ( localtime($t) )[ 5, 4, 3, 0, 1, 2 ];
    $y += 1900;
    $d = "0" . $d if ( $d < 10 );
    my $ms = text_month($m);
    return sprintf "$d-$ms-$y %02d:%02d:%02d", $hh, $mm, $ss;
}


sub dumpSNP {
    my ( $dbspecies, $conf ) = @_;

    #    warn Dumper $conf;

    warn "\n", '*' x 20, "\n";


    my $dbname = $conf->{variation}->{$release} or next;
    my $file = "$dir/SNP_$dbname.xml";
    $file .= ".gz" unless $nogzip;
    my $start_time = time;
    warn "Dumping $dbname to $file ... ", format_datetime($start_time),
      "\n";
    unless ($nogzip) {
	$fh = new IO::Zlib;
	$fh->open( "$file", "wb9" )
	  || die "Can't open compressed stream to $file: $!";
    }
    else {
	open( FILE, ">$file" ) || die "Can't open $file: $!";
    }

    header( $dbname, $dbspecies, $dbname );
    my $dsn = "DBI:mysql:host=$host";
    $dsn .= ";port=$port" if ($port);
    my $ecount;
    my $dbh = DBI->connect( "$dsn:$dbname", $user, $pass )
      or die "DBI::error";



    my $source_hash = $dbh->selectall_hashref(qq{SELECT source_id, name FROM source} ,[qw(source_id)]);



    my $sth = $dbh->prepare("select vf.variation_name, vf.source_id, group_concat(vs.source_id, ' ',vs.name) 
from variation_feature vf 
left join variation_synonym vs on vf.variation_id = vs.variation_id 
 group by vf.variation_id");


    $sth->execute() or die "Error:", $DBI::errstr;

    my $xml;
    my @variation_synonym;


    my $rows = [];
    while (my $row = (shift(@$rows) || 
		      shift( @{ $rows= $sth->fetchall_arrayref(undef, 10_000) || []  } ))
	  ){

	my @synonyms = split /,/, @$row->[2];
	my $snp_source = $source_hash->{$row->[1]}->{name};

	my $description =
        "A $snp_source SNP with "
      . scalar @synonyms
      . ' synonym'
      . (  @synonyms > 1 | @synonyms < 1 ? 's '   : ' ' )
      . ( @synonyms > 0 ? "( " . (join "",  map{  map{  $source_hash->{$_->[0]}->{name} , ':', $_->[1] , ' ' } [split]  } @synonyms ) . ")" : '' );

	$xml =   qq{<entry id="$row->[0]">
  <description>$description</description>
  <additional_fields>
    <field name="species">$species</field>
    <field name="featuretype">SNP</field>
  <cross_references>
    <ref dbname="$source_hash->{$row->[1]}->{name}" dbkey="$row->[0]"/>
  </cross_references>
</entry>};

	p($xml);


    }

    footer($sth->rows);
    print_time($start_time);

}


sub dumpGenomicAlignment {
    my ( $dbspecies, $conf ) = @_;

    #    warn Dumper $conf;

    warn "\n", '*' x 20, "\n";
    my %tables = (
        'dna_align_feature' => [ 'DnaAlignFeature', 'DNA alignment feature' ],
        'protein_align_feature' =>
          [ 'ProteinAlignFeature', 'Protein alignment feature' ]
    );
    my $ecount;

    foreach my $db (  'core', 'cdna', 'otherfeatures' ) {

        my $ecount  = 0;
        my $DB_NAME = $conf->{$db}->{$release} or next;
        my $file    = "$dir/GenomicAlignment_$DB_NAME.xml";
        $file .= ".gz" unless $nogzip;
        my $start_time = time;
        warn "Dumping $DB_NAME to $file ... ", format_datetime($start_time),
          "\n";

        unless ($nogzip) {
            $fh = new IO::Zlib;
            $fh->open( "$file", "wb9" )
              || die "Can't open compressed stream to $file: ", $!;
        }
        else {
            open( FILE, ">$file" ) || die "Can't open $file: ", $!;
        }

        header( $DB_NAME, $dbspecies, $db );
        my $dsn = "DBI:mysql:host=$host";
        $dsn .= ";port=$port" if ($port);

        my $dbh = DBI->connect( "$dsn:$DB_NAME", $user, $pass )
          or die "DBI::error";
        foreach my $table ( keys %tables ) {
            my $source = $tables{$table}[0];

            #	    $source .= ";db=$db" unless $db eq 'core';\

            my $type = $tables{$table}[1];

            # Due to the sheer number of features - generating this
            # dump causes temp table to be written on the mysqld
            # filesys. /tmp space can (and has) be easily filled.
            # Have split the feature fetching to happen by an
            # analysis_id at a time.  The changes below gave a x3
            # speed up on XML dumping as compared to fetching all at
            # features in one query once.

            # make a lookup for the analysis display labels.
            my $logic_name_lookup = $dbh->selectall_hashref(
                "select a.analysis_id, a.logic_name
                                                    from $DB_NAME.analysis as a
                                                ", [qw(analysis_id)]
            ) or die $DBI::Err;

            my $display_label_lookup = $dbh->selectall_hashref(
                "select ad.analysis_id, ad.display_label
                                                    from $DB_NAME.analysis_description as ad
                                                ", [qw(analysis_id)]
            ) or die $DBI::Err;

            my $sth = $dbh->prepare(
                "select t.analysis_id,t.hit_name, count(*) as hits
                                         from $DB_NAME.$table as t where t.analysis_id = ?
                                        group by t.analysis_id, t.hit_name"
            ) or die $DBI::Err;

            foreach my $ana_id (
                @{$dbh->selectall_arrayref("select distinct distinct(analysis_id) from $DB_NAME.$table")
                }
              )
            {

                my $adesc =
                  ( $display_label_lookup->{ $ana_id->[0] }->{display_label}
                      || $logic_name_lookup->{ $ana_id->[0] }->{logic_name} );

                $sth->execute( $ana_id->[0] ) or die $DBI::Err;

                my $rows = [];    # cache for batches of rows
                while (
                    my $row = (
                        shift(@$rows) ||    # get row from cache,
                                            # or reload cache:
                          shift(
                            @{
                                $rows =
                                  $sth->fetchall_arrayref( undef, 100_000 )
                                  || []
                              }
                          )
                    )
                  )
                {

                    my $hid   = $row->[1];
                    my $count = $row->[2];

                    my $xml = qq{
  <entry id="$hid">
    <name>$type:$adesc:$hid</name>
    <description>$adesc $hid hits the genome in $count locations.</description>
    <additional_fields>
        <field name="species">$dbspecies</field>
        <field name="analysis">$adesc</field>
        <field name="featuretype">GenomicAlignment</field>
        <field name="source">$source</field>
        <field name="db">$db</field>
    </additional_fields>
  </entry>};
                    p($xml);

                }

                $ecount += $sth->rows;

            }

        }
        print_time($start_time);
        footer($ecount);
    }
}

sub dumpQTL {
    my ( $dbspecies, $conf ) = @_;

    # print Dumper($conf);
    my $xml_data;
    $xml_data->{species} = $dbspecies;
    my $db     = 'core';
    my $dbname = $conf->{$db}->{$release} or next;
    my $file   = "$dir/QTL_$dbname.xml";
    $file .= ".gz" unless $nogzip;
    my $start_time = time;
    warn "Dumping $dbname to $file ... ", format_datetime($start_time), "\n";

    unless ($nogzip) {
        $fh = new IO::Zlib;
        $fh->open( "$file", "wb9" )
          || die("Can't open compressed stream to $file: $!");
    }
    else {
        open( FILE, ">$file" ) || die "Can't open $file: $!";
    }

    header( $dbname, $dbspecies, $db );
    my $dsn = "DBI:mysql:host=$host";
    $dsn .= ";port=$port" if ($port);
    my $ecount;
    my $dbh = DBI->connect( "$dsn:$dbname", $user, $pass ) or die "DBI::error";

    my $sth = $dbh->prepare(
        "select c.name as chr, qf.seq_region_start, qf.seq_region_end,
        a.logic_name as analysis, q.qtl_id,
        q.trait, qs.source_database, qs.source_primary_id,
        fms1.source as fm1_source, fms1.name as fm1_name,
        fms2.source as fm2_source, fms2.name as fm2_name,
        pms.source  as pm_source,  pms.name  as pm_name
   from ((((((seq_region as c, qtl_feature as qf, qtl_synonym as qs,
        analysis as a, qtl as q) left join marker as fm1 on
        fm1.marker_id = q.flank_marker_id_1) left join marker_synonym as fms1 on
        fm1.display_marker_synonym_id = fms1.marker_synonym_id) left join marker as fm2 on
        fm2.marker_id = q.flank_marker_id_2) left join marker_synonym as fms2 on
        fm2.display_marker_synonym_id = fms2.marker_synonym_id) left join marker as pm on
        pm.marker_id = q.peak_marker_id) left join marker_synonym as pms on
        pm.display_marker_synonym_id = pms.marker_synonym_id 
  where c.seq_region_id = qf.seq_region_id and qs.qtl_id = q.qtl_id and 
        qf.analysis_id = a.analysis_id and qf.qtl_id = q.qtl_id
  "
    );
    $sth->execute();
    my $desc    = '';
    my $old_qtl = 0;
    my $old_ID  = '';
    my $old_pos = '';
    my $counter = make_counter(0);
    while ( my $T = $sth->fetchrow_hashref() ) {

        if ( $T->{qtl_id} eq $old_qtl ) {

            #	    $IDS  .= " $T->{source_primary_id}";
            $desc .= " $T->{source_database}:$T->{source_primary_id}";

            $xml_data->{cross_ref}->{ $T->{source_database} } =
              $T->{source_primary_id};

            #	    $xml_data->{source_primary_id} = $T->{source_primary_id};
            #	    print Dumper($T);
        }
        else {
            $xml_data->{pm_name} = $T->{pm_name};
            $old_pos =
                "$T->{chr}:"
              . ( $T->{seq_region_start} - 1e4 ) . '-'
              . ( $T->{seq_region_end} + 1e4 );
            $desc = "QTL exhibiting '$T->{trait}' has ";
            my $f2 = $T->{pm_name} ? 1 : 0;
            if ( $T->{fm1_name} || $T->{fm2_name} ) {
                my $f1 = ( $T->{fm1_name} ) && ( $T->{fm2_name} ) ? 1 : 0;

                $desc .=
                    'flanking marker'
                  . ( $f1 ? 's ' : ' ' )
                  . $T->{fm1_name}
                  . ( $f1 ? ' and ' : '' )
                  . $T->{fm2_name}
                  . ( $f2 ? '; ' : '' );
                $xml_data->{f1} = $T->{fm1_name};
                $xml_data->{f2} = $T->{fm2_name};

            }
            if ($f2) {
                $desc .= "peak marker $T->{pm_name};";
                $xml_data->{pm} = $T->{pm_name};
            }
            $desc .=
              " and names: $T->{source_database}:$T->{source_primary_id}";

            #my $sd = $T->{source_database};
            $xml_data->{description} = $desc;
            $xml_data->{cross_ref}->{ $T->{source_database} } =
              $T->{source_primary_id};
            $old_qtl = $T->{qtl_id};
            $xml_data->{pos} = $old_pos;
            if ( $xml_data->{pm_name} ) {
                p( QTLXML( $xml_data, $counter ) );
            }
        }
    }

    $xml_data->{description} = $desc;

    if ( $xml_data->{pm_name} ) {
        p( QTLXML( $xml_data, $counter ) );
    }
    $dbh->disconnect();
    footer( $counter->() );
}

sub QTLXML {
    my ( $xml_data, $counter ) = @_;

    my $xml = qq{
<entry id="$xml_data->{pm_name}">   
  <name>$xml_data->{pm_name}</name>
  <description>$xml_data->{description}</description>
  <additional_fields>
    <field name="species">$xml_data->{species}</field>
    <field name="flanking">$xml_data->{f1}</field>
    <field name="flanking">$xml_data->{f2}</field>
    <field name="peak_marker">$xml_data->{pm}</field>
    <field name="pos">$xml_data->{pos}</field>
    <field name="featuretype">QTL</field>
  </additional_fields>
  <cross_references> };

    foreach ( keys( %{ $xml_data->{cross_ref} } ) ) {
        $xml .=
          qq{\n    <ref dbname="$_" dbkey="$xml_data->{cross_ref}->{$_}"/>};
    }

    $xml .= qq{\n  </cross_references>
</entry>};
    $counter->();
    return $xml;

}

sub dumpMarker {
    my ( $dbspecies, $conf ) = @_;

    #    my $xml_data;
    #   $xml_data->{species} = $dbspecies;

    my $db     = 'core';
    my $dbname = $conf->{$db}->{$release} or next;
    my $file   = "$dir/Marker_$dbname.xml";
    $file .= ".gz" unless $nogzip;
    my $start_time = time;
    warn "Dumping $dbname to $file ... ", format_datetime($start_time), "\n";

    unless ($nogzip) {
        $fh = new IO::Zlib;
        $fh->open( "$file", "wb9" )
          || die("Can't open compressed stream to $file: $!");
    }
    else {
        open( FILE, ">$file" ) || die "Can't open $file: $!";
    }
    header( $dbname, $dbspecies, $db );
    my $dsn = "DBI:mysql:host=$host";
    $dsn .= ";port=$port" if ($port);
    my $ecount;
    my $dbh = DBI->connect( "$dsn:$dbname", $user, $pass ) or die "DBI::error";

#   my $sth = $dbh->prepare(q{
#   SELECT @rownum := @rownum+1 AS rownum,  ms2.name as marker, ms1.name
#     from (select @rownum := 0) r, (marker_synonym as ms1, marker as m) left join
#          marker_synonym as ms2 on ms2.marker_synonym_id = m.display_marker_synonym_id
#   where ms1.marker_id = m.marker_id
#  order by m.marker_id}
#  );
#  $sth->execute( );
    my $sth = $dbh->prepare(
        q{ SELECT  ms2.name as marker, ms1.name
       from  (marker_synonym as ms1, marker as m) left join
            marker_synonym as ms2 on ms2.marker_synonym_id = m.display_marker_synonym_id
      where ms1.marker_id = m.marker_id
      order by m.marker_id}
    );

    my $data = $dbh->selectall_hashref( $sth, [ 'marker', 'name' ] );
    foreach my $marker ( keys %$data ) {
        p markerXML( $marker, $data, $dbspecies );
    }

    footer( scalar keys(%$data) );

}

sub markerXML {
    my ( $marker, $xml_data, $species ) = @_;

    my $xml;

    my @keys = keys %{ $xml_data->{$marker} };
    my $desc =
        'A marker with '
      . scalar @keys
      . ' synonym'
      . ( scalar @keys > 1 ? 's ' : ' ' ) . '('
      . join( " ", @keys ) . ')';
    $xml = qq{ 
<entry id="$marker">   
  <name>$marker</name>
   <description>$desc</description>
   <additional_fields>};

    foreach (@keys) {
        $xml .= qq{
      <field name="synonym">$_</field>}

    }
    $xml .= qq{
     <field name="species">$species</field>
    <field name="featuretype">marker</field>
  </additional_fields>
</entry>};

    return $xml;

}

sub dumpOligoProbe {
    my ( $dbspecies, $conf ) = @_;

    my $db     = 'core';
    my $dbname = $conf->{$db}->{$release} or next;
    my $file   = "$dir/OligoProbe_$dbname.xml";
    $file .= ".gz" unless $nogzip;
    my $start_time = time;
    warn "Dumping $dbname to $file ... ", format_datetime($start_time), "\n";

    unless ($nogzip) {
        $fh = new IO::Zlib;
        $fh->open( "$file", "wb9" )
          || die("Can't open compressed stream to $file: $!");
    }
    else {
        open( FILE, ">$file" ) || die "Can't open $file: $!";
    }
    header( $dbname, $dbspecies, $db );
    my $dsn = "DBI:mysql:host=$host";
    $dsn .= ";port=$port" if ($port);
    my $ecount;
    my $dbh = DBI->connect( "$dsn:$dbname", $user, $pass ) or die "DBI::error";

    my $sth = $dbh->prepare(
        "select p.probeset, count(*) as hits, a.type
       from oligo_probe as p, oligo_feature as f, oligo_array as a
      where p.oligo_probe_id = f.oligo_probe_id and p.oligo_array_id = a.oligo_array_id
      group by p.probeset"
    );
    $sth->execute();

    my $count = 0;
    while ( my $data = $sth->fetchrow_arrayref ) {
        p OligoProbeXML( $data, $species );

#      print  join "\t", ("$dbspecies\t")."$type Oligo Probe set",
#     $hid, "/$species/featureview?type=OligoProbe;id=$hid",
#    $hid, qq($type oligo probeset $hid hits the genome in $count locations.\n);
    }

    footer( $sth->rows );

}

sub OligoProbeXML {
    my ( $xml_data, $species ) = @_;

    my $desc =
qq{$xml_data->[0], $xml_data->[2] oligo probeset $xml_data->[0] hits the genome in $xml_data->[1] locations.};

    return qq{ 
<entry id="$xml_data->[0]"> 
  <name>$xml_data->[0]</name>
   <description>$desc</description>
   <additional_fields>
      <field name="type">$xml_data->[2]</field>
     <field name="species">$species</field>
    <field name="featuretype">OligoProbe</field>
  </additional_fields>
</entry>};

}

sub dumpDomain {
    my ( $dbspecies, $conf ) = @_;

    my $db = 'core';
    my $dbname = $conf->{$db}->{$release} or next;

    my $file = "$dir/Domain_$dbname.xml";
    $file .= ".gz" unless $nogzip;
    my $start_time = time;
    warn "Dumping $dbname to $file ... ", format_datetime($start_time), "\n";

    unless ($nogzip) {
        $fh = new IO::Zlib;
        $fh->open( "$file", "wb9" )
          || die("Can't open compressed stream to $file: $!");
    }
    else {
        open( FILE, ">$file" ) || die "Can't open $file: $!";
    }
    header( $dbname, $dbspecies, $db );
    my $dsn = "DBI:mysql:host=$host";
    $dsn .= ";port=$port" if ($port);

    my $dbh = DBI->connect( "$dsn:$dbname", $user, $pass ) or die "DBI::error";
    my $sth = $dbh->prepare(
        "select dbprimary_acc, i.id, x.description
       from xref as x, interpro as i
      where x.dbprimary_acc = i.interpro_ac
      order by x.dbprimary_acc"
    );
    $sth->execute();
    my $old_acc     = '';
    my $IDS         = '';
    my $description = '';
    my $counter     = 0;
    my ( $acc, $id, $desc, $old_desc, $ecount );

    while ( ( $acc, $id, $desc ) = $sth->fetchrow_array() ) {
        if ( $acc eq $old_acc ) {

            #    $IDS         .= " $id";
            push @$IDS, $id;
            $description .= ", $id";
            $counter++;
        }
        else {
            if ( $old_acc ne '' ) {
                p domainLineXML(
                    $dbspecies, $old_acc,     $IDS,
                    $old_desc,  $description, $counter
                );
                $ecount++;
            }

            $description = $id;
            $IDS         = undef;
            $IDS->[0]    = $id;
            $old_acc     = $acc;
            $old_desc    = $desc;
            $counter     = 1;
        }

    }
    if ( $old_acc ne '' ) {
        p domainLineXML( $dbspecies, $old_acc, $IDS, $old_desc, $description,
            $counter );
        $ecount++;
    }

    footer($ecount);

}

sub domainLineXML {
    my ( $species, $did, $IDS, $desc, $description, $counter ) = @_;

    my $xml = qq{
<entry id="$did"> 
  <name>$did</name>
   <description>InterPro domain $did [$desc] has $counter associated external database identifiers: $description.</description>
   <additional_fields>
      <field name="type">Interpro domain</field>};

    map {
        $xml .= qq{
      <field name="external_id">$_</field>}
    } @$IDS;

    $xml .= qq{
      <field name="species">$species</field>
      <field name="featuretype">Domain</field>
  </additional_fields>
</entry>};
    return $xml;

}

sub dumpFamily {
    my ( $dbspecies, $conf ) = @_;

    my $db = 'core';
    my $dbname = $conf->{$db}->{$release} or next;

    my $file = "$dir/Family_$dbname.xml";
    $file .= ".gz" unless $nogzip;
    my $start_time = time;
    warn "Dumping $dbname to $file ... ", format_datetime($start_time), "\n";

    unless ($nogzip) {
        $fh = new IO::Zlib;
        $fh->open( "$file", "wb9" )
          || die("Can't open compressed stream to $file: $!");
    }
    else {
        open( FILE, ">$file" ) || die "Can't open $file: $!";
    }
    header( $dbname, $dbspecies, $db );
    my $dsn = "DBI:mysql:host=$host";
    $dsn .= ";port=$port" if ($port);
    my $ecount;
    my $dbh = DBI->connect( "$dsn:$dbname", $user, $pass ) or die "DBI::error";

    my $FAMDB = $conf->{'compara'}->{$release};

    my $CORE  = $conf->{'core'}->{$release};
    my $t_sth = $dbh->prepare(
qq{select meta_value from $CORE.meta where meta_key='species.taxonomy_id'}
    );
    $t_sth->execute;
    my $taxon_id = ( $t_sth->fetchrow );

    return unless $taxon_id;

    $dbh->do("SET SESSION group_concat_max_len = 100000");
    my $sth = $dbh->prepare(
qq{ select f.stable_id as fid , f.description, group_concat(m.stable_id, unhex('1D') ,m.source_name) as IDS 
from $FAMDB.family as f, $FAMDB.family_member as fm, $FAMDB.member as m 
 where fm.family_id = f.family_id and fm.member_id = m.member_id and m.taxon_id = $taxon_id  group by fid}
    );
    $sth->execute;
    foreach my $xml_data ( @{ $sth->fetchall_arrayref( {} ) } ) {

        my @bits = split /,/, delete $xml_data->{IDS};
        map { push @{ $xml_data->{IDS} }, [ split /\x1D/ ] } @bits;
        $xml_data->{species} = $dbspecies;
        p familyLineXML($xml_data);

    }

    footer( $sth->rows );

}

sub familyLineXML {
    my ( $xml_data, $counter ) = @_;

    my $description_line =
        "Ensembl protein family $xml_data->{fid} "
      . "[$xml_data->{description}] "
      . ( join( " ", map { $_->[0] } @{ $xml_data->{IDS} }[ 1 .. 8 ] ) )
      . ' ...' . " has "
      . scalar @{ $xml_data->{IDS} }
      . " members.";

    my $xml = qq{ 
<entry id="$xml_data->{fid}"> 
  <name>$xml_data->{fid}</name>
   <description>$description_line</description>
   <cross_references>} .

      (
        join "",
        (
            map {
                qq{
     <ref dbname="$1" dbkey="$_->[0]"/>} if $_->[1] =~ /(Uniprot|ENSEMBL).*/
              } @{ $xml_data->{IDS} }
        )
      )
    .

      qq{
  </cross_references>
  <additional_fields>
     <field name="species">$xml_data->{species}</field>
    <field name="featuretype">Family</field>
  </additional_fields>
</entry>};

    return $xml;

}

sub dumpSequence {
    my ( $dbspecies, $conf ) = @_;

    #    my $sanger = sanger_project_names( $conf );
    my $sanger = 'SANGER STUFF';
    my %config = (
        "homo_sapiens" => [
            [
                'Clone',
                'tilepath, cloneset_1mb, cloneset_30k, cloneset_32k',
'name,well_name,clone_name,synonym,embl_acc,sanger_project,alt_well_name,bacend_well_name'
            ],
            [ 'NT Contig',     'ntctgs', 'name' ],
            [ 'Encode region', 'encode', 'name,synonym,description' ],
        ],
        "mus_musculus" => [
            [
                'BAC',
                'cloneset_0_5mb,cloneset_1mb,bac_map,tilingpath_cloneset',
                'embl_acc,name,clone_name,well_name,synonym,alt_embl_acc'
            ],
            [ 'Fosmid',      'fosmid_map', 'name,clone_name' ],
            [ 'Supercontig', 'superctgs',  'name' ],
        ],
        "anopheles_gambiae" => [
            [ 'BAC',      'bacs',       'name,synonym,clone_name' ],
            [ 'BAC band', 'bacs_bands', 'name,synonym,clone_name' ],
        ],
        "gallus_gallus" => [
            [ 'BAC', 'bac_map', 'name,synonym,clone_name' ],
            [
                'BAC ends',                'bacends',
                'name,synonym,clone_name', 'otherfeatures'
            ]
        ]
    );

    my $dbname = $conf->{'core'}->{$release} or next;

    my $file = "$dir/Sequence_$dbname.xml";
    $file .= ".gz" unless $nogzip;
    my $start_time = time;
    warn "Dumping $dbname to $file ... ", format_datetime($start_time), "\n";

    unless ($nogzip) {
        $fh = new IO::Zlib;
        $fh->open( "$file", "wb9" )
          || die("Can't open compressed stream to $file: $!");
    }
    else {
        open( FILE, ">$file" ) || die "Can't open $file: $!";
    }
    header( $dbname, $dbspecies, 'core' );
    my $dsn = "DBI:mysql:host=$host";
    $dsn .= ";port=$port" if ($port);
    my $ecount;
    my $dbh = DBI->connect( "$dsn:$dbname", $user, $pass ) or die "DBI::error";

    my $COREDB = $dbname;
    my $ESTDB  = $conf->{otherfeatures}->{$release};

    my @types = @{ $config{$dbspecies} || [] };
    my $ecounter;
    foreach my $arrayref (@types) {

        my ( $TYPE, $mapsets, $annotationtypes, $DB ) = @$arrayref;

        my $DB = $DB eq 'otherfeatures' ? $ESTDB : $COREDB;
        my @temp = split( ',', $mapsets );
        my @mapsets;
        foreach my $X (@temp) {
            my $ID = $dbh->selectrow_array(
                "select misc_set_id from $DB.misc_set where code = ?",
                {}, $X );
            push @mapsets, $ID if ($ID);
        }

        next unless @mapsets;
        @temp = split( ',', $annotationtypes );
        my @mapannotationtypes;
        foreach my $X (@temp) {
            my $ID = $dbh->selectrow_array(
                "select attrib_type_id from $DB.attrib_type where code = ?",
                {}, $X );
            push @mapannotationtypes, $ID if ($ID);
        }
        next unless @mapannotationtypes;
        my $Z       = " ma.value";
        my $MAPSETS = join ',', @mapsets;
        my $sth     = $dbh->prepare(
            "select mf.misc_feature_id, sr.name,
              ma.value, mf.seq_region_end-mf.seq_region_start+1 as len, 
              at.code
         from $DB.misc_feature_misc_set as ms, 
              $DB.misc_feature as mf,
              seq_region   as sr,
              $DB.misc_attrib  as ma,
              $DB.attrib_type  as at 
        where mf.seq_region_id = sr.seq_region_id and mf.misc_feature_id = ms.misc_feature_id and ms.misc_set_id in ($MAPSETS) and
              mf.misc_feature_id = ma.misc_feature_id and ma.attrib_type_id = at.attrib_type_id
        order by mf.misc_feature_id, at.code"
        );
        $sth->execute();
        my ( $oldtype, $old_ID, $oldchr, $emblaccs, $oldlen, $synonyms, $NAME );

        while ( my ( $ID, $chr, $val, $len, $type ) = $sth->fetchrow_array() ) {

            if ( $ID == $old_ID ) {
                $NAME = $val
                  if $type eq 'well_name'
                      || $type eq 'clone_name'
                      || $type eq 'name'
                      || $type eq 'non_ref';
                $NAME = $val if !$NAME && $type eq 'embl_acc';
                $NAME = $val if !$NAME && $type eq 'synonym';
                $NAME = $val if !$NAME && $type eq 'sanger_project';
                push @{$emblaccs}, $val if $val;
            }
            else {
                p seqLineXML(
                    $dbspecies, $TYPE,   $NAME, $oldchr,
                    $emblaccs,  $oldlen, $sanger
                ) if $old_ID;
                $NAME     = undef;
                $emblaccs = undef;
                $NAME     = $val
                  if $type eq 'well_name'
                      || $type eq 'clone_name'
                      || $type eq 'name'
                      || $type eq 'non_ref';
                $NAME = $val if !$NAME && $type eq 'embl_acc';
                $NAME = $val if !$NAME && $type eq 'synonym';
                $NAME = $val if !$NAME && $type eq 'sanger_project';
                $emblaccs->[0] = $val;
                ( $old_ID, $oldchr, $oldlen ) = ( $ID, $chr, $len );
                $ecounter += 1;
            }
        }
        p seqLineXML( $dbspecies, $TYPE, $NAME, $oldchr, $emblaccs, $oldlen,
            $sanger )
          if $old_ID;
    }

    footer($ecounter);

#   my $sth = $conf->{'dbh'}->prepare(
#     "select c.name, c.length, cs.name
#        from seq_region as c, coord_system as cs
#       where c.coord_system_id = cs.coord_system_id" );
#   $sth->execute();
#   while( my($name,$length,$type) = $sth->fetchrow_array() ) {
#     my $extra_IDS = ''; mysql $extra_desc = '';
#     if( %{$sanger->{$name}||{}} ) {
#       $extra_IDS  = join ' ', '',sort keys %{$sanger->{$name}};
#       $extra_desc = " and corresponds to the following Sanger projects: ".join( ', ',sort keys %{$sanger->{$name}});
#     }
#     print_time O join "\t",
#       (INC_SPECIES?"$conf->{'species'}\t":"").ucfirst($type),       $name,
#       ($type eq 'chromosome' && length( $name ) < 5) ?
#         "/$conf->{'species'}/mapview?chr=$name" :
#         ($length > 0.5e6 ? "/$conf->{'species'}/cytoview?region=$name" :
#               "/$conf->{'species'}/contigview?region=$name" ),
#       "$name$extra_IDS", "$name isnull a @{[ucfirst($type)]} (of length $length)$extra_desc\n";
#   }
}

sub seqLineXML {
    my ( $species, $type, $name, $chr, $val, $len, $sanger ) = @_;

    pop @$val;
    my $description = "$type $name is mapped to Chromosome $chr" .

      (
        @$val > 0
        ? ' and has ' 
          . @$val 
          . " EMBL accession"
          . (
            @$val > 1
            ? 's'
            : ''
          )
          . "/synonym"
          . (
            @$val > 1
            ? 's '
            : ' '
          )
          . "@$val"
        : ''
      )

      . " and length $len bps\n";

    my $xml = qq{ 
 <entry id="$name"> 
   <name>$name</name>
    <description>$description</description>
    <cross_references>}

      . (
        join "",
        (
            map {
                qq{
      <ref dbname="EMBL" dbkey="$_"/>}
              } @$val
        )
      )

      . qq{</cross_references>
    <additional_fields>
      <field name="species">$species</field>
      <field name="type">$type</field>
      <field name="chromosome">$chr</field>
      <field name="length">$len</field>
    <field name="featuretype">Genomic</field>
   </additional_fields>
 </entry>};

    return $xml;

}

sub dumpGene {
    my ( $dbspecies, $conf ) = @_;

    foreach my $DB ( 'core', 'otherfeatures', 'vega' ) {
        my $counter = make_counter(0);

        my $DBNAME = $conf->{$DB}->{$release}
          or warn "$dbspecies $DB $release: no database not found";
        next unless $DBNAME;
        print "START... $DB";
        my $file = "$dir/Gene_$DBNAME.xml";
        $file .= ".gz" unless $nogzip;
        my $start_time = time;

        unless ($nogzip) {
            $fh = new IO::Zlib;
            $fh->open( "$file", "wb9" )
              || die("Can't open compressed stream to $file: $!");
        }
        else {
            open( FILE, ">$file" ) || die "Can't open $file: $!";
        }
        header( $DBNAME, $dbspecies, $DB );
        my $dsn = "DBI:mysql:host=$host";
        $dsn .= ";port=$port" if ($port);

        warn "Dumping $DBNAME to $file ... ", format_datetime($start_time),
          "\n";
        my $extra = $DB ne 'core' ? ";db=$DB" : '';

        my $dbh = DBI->connect( "$dsn:$DBNAME", $user, $pass )
          or die "DBI::error";

        my %xrefs      = ();
        my %xrefs_desc = ();
        my %disp_xrefs = ();
        foreach my $type (qw(Gene Transcript Translation)) {
            my $T = $dbh->selectall_arrayref(
                "select ox.ensembl_id,
                x.display_label, x.dbprimary_acc, ed.db_display_name, es.synonym, x.description
           from ($DBNAME.object_xref as ox, $DBNAME.xref as x, $DBNAME.external_db as ed) left join $DBNAME.external_synonym as es on es.xref_id = x.xref_id
          where ox.ensembl_object_type = '$type' and ox.xref_id = x.xref_id and x.external_db_id = ed.external_db_id"
            );
            foreach (@$T) {

                $xrefs{$type}{ $_->[0] }{ $_->[3] }{ $_->[1] } = 1 if $_->[1];
                $xrefs{$type}{ $_->[0] }{ $_->[3] }{ $_->[2] } = 1 if $_->[2];
                $xrefs{$type}{ $_->[0] }{ $_->[3] }{ $_->[4] } = 1 if $_->[4];
                $xrefs_desc{$type}{ $_->[0] }{ $_->[5] }       = 1 if $_->[5];
            }

            warn "XREF $type query...";
        }

        my %exons = ();
        my $T     = $dbh->selectall_arrayref(
            "select distinct t.gene_id, esi.stable_id
         from transcript as t, exon_transcript as et, exon_stable_id as esi
        where t.transcript_id = et.transcript_id and et.exon_id = esi.exon_id"
        );
        foreach (@$T) {
            $exons{ $_->[0] }{ $_->[1] } = 1;
        }
        my $gene_info = $dbh->selectall_arrayref( "
        select gsi.gene_id, tsi.transcript_id, trsi.translation_id,
             gsi.stable_id as gsid, tsi.stable_id as tsid, trsi.stable_id as trsid,
             g.description, ed.db_display_name, x.dbprimary_acc,x.display_label, ad.display_label, ad.description, g.source, g.status, g.biotype
        from (((( $DBNAME.gene_stable_id as gsi, $DBNAME.gene as g,
             $DBNAME.transcript_stable_id as tsi,
             $DBNAME.analysis_description as ad,
             $DBNAME.transcript as t) left join
             $DBNAME.translation as tr on t.transcript_id = tr.transcript_id) left join
             $DBNAME.translation_stable_id as trsi on tr.translation_id = trsi.translation_id) left join
             $DBNAME.xref as x on g.display_xref_id = x.xref_id) left join
             $DBNAME.external_db as ed on ed.external_db_id = x.external_db_id
       where t.gene_id = gsi.gene_id and t.transcript_id = tsi.transcript_id and t.gene_id = g.gene_id 
             and g.analysis_id = ad.analysis_id
       order by gsi.stable_id, tsi.stable_id;
    " );
        warn "Gene query...";

        my %hash = map { $_->[0] } @$gene_info;
        my $ecount = scalar keys %hash, "\n\n";

        my %old;

        foreach my $row (@$gene_info) {

            # g = gene_id, t = transcript_id , tr = translation_id ,
            # gs = gene_stable_id, ts = transcript_stable_id , trs =
            # translation_stable_id, d = description,
            # ddb= external_db_dispay_name,
            # dpa = xref_primary_accession,
            # dn = xref display_label,
            # a = analysis_description
            # display label,
            # ad = analysis description descripion,
            #s = gene.source, st = gene.status, bt = gene.biotype

            my (
                $gene_id,                            $transcript_id,
                $translation_id,                     $gene_stable_id,
                $transcript_stable_id,               $translation_stable_id,
                $gene_description,                   $extdb_db_display_name,
                $xref_primary_acc,                   $xref_display_label,
                $analysis_description_display_label, $analysis_description,
                $gene_source,                        $gene_status,
                $gene_biotype
            ) = @$row;
            if ( $old{'gene_id'} != $gene_id ) {
                if ( $old{'gene_id'} ) {
                    p &geneLineXML( $dbspecies, \%old, $counter );

                }
                %old = (
                    'gene_id'                => $gene_id,
                    'gene_stable_id'         => $gene_stable_id,
                    'description'            => $gene_description,
                    'translation_stable_ids' => {
                        $translation_stable_id ? ( $translation_stable_id => 1 )
                        : ()
                    },
                    'transcript_stable_ids' => {
                        $transcript_stable_id ? ( $transcript_stable_id => 1 )
                        : ()
                    },
                    'exons'                => {},
                    'external_identifiers' => {},
                    'alt'                  => $xref_display_label
                    ? "($extdb_db_display_name: $xref_display_label)"
                    : "(novel gene)",
                    'ana_desc_label' => $analysis_description_display_label,
                    'ad'             => $analysis_description,
                    'source'         => ucfirst($gene_source),
                    'st'             => $gene_status,
                    'biotype'        => $gene_biotype
                );
                $old{'source'} =~ s/base/Base/;
                $old{'exons'} = $exons{$gene_id};
                foreach my $K ( keys %{ $exons{$gene_id} } ) {
                    $old{'i'}{$K} = 1;
                }

                foreach my $db ( keys %{ $xrefs{'Gene'}{$gene_id} || {} } ) {
                    foreach my $K ( keys %{ $xrefs{'Gene'}{$gene_id}{$db} } ) {
                        $old{'external_identifiers'}{$db}{$K} = 1;

                    }
                }
                foreach my $db (
                    keys %{ $xrefs{'Transcript'}{$transcript_id} || {} } )
                {
                    foreach my $K (
                        keys %{ $xrefs{'Transcript'}{$transcript_id}{$db} } )
                    {
                        $old{'external_identifiers'}{$db}{$K} = 1;

                    }
                }
                foreach my $db (
                    keys %{ $xrefs{'Translation'}{$translation_id} || {} } )
                {
                    foreach my $K (
                        keys %{ $xrefs{'Translation'}{$translation_id}{$db} } )
                    {
                        $old{'external_identifiers'}{$db}{$K} = 1;
                    }
                }

            }
            else {
                $old{'transcript_stable_ids'}{$transcript_stable_id}   = 1;
                $old{'translation_stable_ids'}{$translation_stable_id} = 1;

                foreach my $db (
                    keys %{ $xrefs{'Transcript'}{$transcript_id} || {} } )
                {
                    foreach my $K (
                        keys %{ $xrefs{'Transcript'}{$transcript_id}{$db} } )
                    {
                        $old{'external_identifiers'}{$db}{$K} = 1;
                    }
                }
                foreach my $db (
                    keys %{ $xrefs{'Translation'}{$translation_id} || {} } )
                {
                    foreach my $K (
                        keys %{ $xrefs{'Translation'}{$translation_id}{$db} } )
                    {
                        $old{'external_identifiers'}{$db}{$K} = 1;

                    }
                }
            }
        }

        p geneLineXML( $dbspecies, \%old, $counter );

        footer( $counter->() );
        warn "FINISHED...... genes $DB ...";

    }

}

sub geneLineXML {
    my ( $species, $xml_data, $counter ) = @_;

    return warn "gene id not set" if $xml_data->{'gene_stable_id'} eq '';

    my $gene_id     = $xml_data->{'gene_stable_id'};
    my $altid       = $xml_data->{'alt'} or die "altid not set";
    my $transcripts = $xml_data->{'transcript_stable_ids'}
      or die "transcripts not set";
    my $peptides = $xml_data->{'translation_stable_ids'}
      or die "peptides not set";
    my $exons = $xml_data->{'exons'} or die "exons not set";
    my $external_identifiers = $xml_data->{'external_identifiers'}
      or die "external_identifiers not set";
    my $description = $xml_data->{'description'};
    my $type        = $xml_data->{'source'} . ' ' . $xml_data->{'biotype'}
      or die "problem setting type";

    my $exon_count       = scalar keys %$exons;
    my $transcript_count = scalar keys %$transcripts;

    my $xml = qq{ 
 <entry id="$gene_id"> 
   <name>$gene_id $altid</name>
    <description>$description</description>
    <cross_references>};

    foreach my $ext_db ( keys %$external_identifiers ) {
        foreach my $dbkey ( keys %{ $external_identifiers->{$ext_db} } ) {
            $ext_db =~ s/^Affymx.*/ensembl/;
            $ext_db =~ s/^Agilent.*/ensembl/;
            $ext_db =~ s/^Illumina.*/ensembl/;
            $ext_db =~ s/^GE Healthcare\/Amersham Codelink WGA/ensembl/;
            $ext_db =~ s/^Havana.*/Vega/;
            $ext_db =~ s/^Uniprot.*/Uniprot/i;

            $xml .= qq{<ref dbname="$ext_db" dbkey="$dbkey"/>
};
        }
    }

    $xml .= qq{
    </cross_references>
    <additional_fields>
      <field name="species">$species</field>
      <field name="featuretype">$type</field>
      <field name="transcript_count">$transcript_count</field> }

      . (
        join "",
        (
            map {
                qq{
      <field name="transcript">$_</field>}
              } keys %$transcripts
        )
      )

      . qq{  <field name="exon_count">$exon_count</field> }

      . (
        join "",
        (
            map {
                qq{
      <field name="exon">$_</field>}
              } keys %$exons
        )
      )

      . (
        join "",
        (
            map {
                qq{
      <field name="peptide">$_</field>}
              } keys %$peptides
        )
      )

      . qq{
   </additional_fields>
 </entry>};
    $counter->();
    return $xml;

}

sub dumpUnmappedFeatures {
    my ( $dbspecies, $conf ) = @_;

    my $db     = 'core';
    my $COREDB = $conf->{$db}->{$release} or next;
    my $file   = "$dir/UnmappedFeature_$COREDB.xml";
    $file .= ".gz" unless $nogzip;
    my $start_time = time;
    warn "Dumping $COREDB to $file ... ", format_datetime($start_time), "\n";

    unless ($nogzip) {
        $fh = new IO::Zlib;
        $fh->open( "$file", "wb9" )
          || die("Can't open compressed stream to $file: $!");
    }
    else {
        open( FILE, ">$file" ) || die "Can't open $file: $!";
    }
    header( $COREDB, $dbspecies, $db );
    my $dsn = "DBI:mysql:host=$host";
    $dsn .= ";port=$port" if ($port);
    my $ecount;
    my $dbh = DBI->connect( "$dsn:$COREDB", $user, $pass ) or die "DBI::error";

    my %unmapped_queries = (
        'None' => qq(
      select a.logic_name, e.db_display_name,
             uo.identifier, ur.summary_description,
             'Not mapped'
        from $COREDB.analysis as a, $COREDB.external_db as e, $COREDB.unmapped_object as uo,
             $COREDB.unmapped_reason as ur
       where a.analysis_id = uo.analysis_id and 
             uo.external_db_id = e.external_db_id and
             uo.unmapped_reason_id = ur.unmapped_reason_id and
               uo.ensembl_id = 0 
    ),
        'Transcript' => qq(
      select a.logic_name, e.db_display_name,
             uo.identifier, ur.summary_description,
             concat( 'Transcript: ', tsi.stable_id, '; Gene: ',gsi.stable_id )
        from $COREDB.analysis as a, $COREDB.external_db as e, $COREDB.unmapped_object as uo,
             $COREDB.unmapped_reason as ur, $COREDB.transcript_stable_id as tsi,
             $COREDB.transcript as t, $COREDB.gene_stable_id as gsi
       where a.analysis_id = uo.analysis_id and 
             uo.external_db_id = e.external_db_id and
             uo.unmapped_reason_id = ur.unmapped_reason_id and
             uo.ensembl_id = t.transcript_id and
             uo.ensembl_object_type = 'Transcript' and
             t.transcript_id = tsi.transcript_id and
             t.gene_id       = gsi.gene_id 
    ),
        'Translation' => qq(
      select a.logic_name, e.db_display_name, uo.identifier, ur.summary_description,
             concat( 'Translation: ',trsi.stable_id,'; Transcript: ', tsi.stable_id, '; Gene: ',gsi.stable_id )
        from $COREDB.analysis as a, $COREDB.external_db as e, $COREDB.unmapped_object as uo,
             $COREDB.unmapped_reason as ur, $COREDB.transcript_stable_id as tsi,
             $COREDB.translation as tr, $COREDB.translation_stable_id as trsi,
             $COREDB.transcript as t, $COREDB.gene_stable_id as gsi
       where a.analysis_id = uo.analysis_id and 
             uo.external_db_id = e.external_db_id and
             uo.unmapped_reason_id = ur.unmapped_reason_id and
             uo.ensembl_id = tr.translation_id and 
             tr.transcript_id = t.transcript_id and
             trsi.translation_id = tr.translation_id and
             uo.ensembl_object_type = 'Translation' and
             t.transcript_id = tsi.transcript_id and
             t.gene_id       = gsi.gene_id 
    )
    );
    my $entry_count = 0;
    foreach my $type ( keys %unmapped_queries ) {
        my $SQL = $unmapped_queries{$type};
        my $sth = $dbh->prepare($SQL);
        $sth->execute;
        while ( my $T = $sth->fetchrow_arrayref() ) {

  #            print join "\t", ("$species\t") . qq(Unmapped feature),
  #             "$T->[1] $T->[2]",
  #            "$dbspecies/featureview?type=Gene;id=$T->[2]", "$T->[2] $T->[4]",
  #           "$T->[3]; $T->[4]\n";
            p unmappedFeatureXML( $T, $dbspecies )

        }
        $entry_count += $sth->rows

    }

    footer($entry_count);

}

sub unmappedFeatureXML {
    my ( $xml_data, $dbspecies ) = @_;

    return qq{ 
 <entry id="$xml_data->[2]"> 
   <name>$xml_data->[1] $xml_data->[2]</name>
    <description>$xml_data->[3]; $xml_data->[4]</description>
    <additional_fields>
      <field name="species">$dbspecies</field>
      <field name="featuretype">UnmappedFeature</field>
    </additional_fields>
 </entry>};

}

sub dumpUnmappedGenes {
    my ( $dbspecies, $conf ) = @_;

    my $db = 'core';
    my $dbname = $conf->{$db}->{$release} or next;

    my $file = "$dir/UnmappedGene_$dbname.xml";
    $file .= ".gz" unless $nogzip;
    my $start_time = time;
    warn "Dumping $dbname to $file ... ", format_datetime($start_time), "\n";

    unless ($nogzip) {
        $fh = new IO::Zlib;
        $fh->open( "$file", "wb9" )
          || die("Can't open compressed stream to $file: $!");
    }
    else {
        open( FILE, ">$file" ) || die "Can't open $file: $!";
    }
    header( $dbname, $dbspecies, $db );
    my $dsn = "DBI:mysql:host=$host";
    $dsn .= ";port=$port" if ($port);
    my $ecount;
    my $dbh = DBI->connect( "$dsn:$dbname", $user, $pass ) or die "DBI::error";

    my $COREDB = $conf->{$db}->{$release};

    my %current_stable_ids = ();
    foreach my $type (qw(gene transcript translation)) {
        $current_stable_ids{$type} = {
            map { @$_ } @{
                $dbh->selectall_arrayref(
                    "select stable_id,1 from $COREDB." . $type . "_stable_id"
                )
              }
        };
    }
    my $species = $dbspecies;
    my $sth     = $dbh->prepare(
        qq(
    select sie.type, sie.old_stable_id, if(isnull(sie.new_stable_id),'NULL',sie.new_stable_id), ms.old_release*1.0 as X, ms.new_release*1.0 as Y
      from $COREDB.mapping_session as ms, $COREDB.stable_id_event as sie
     where ms.mapping_session_id = sie.mapping_session_id and ( old_stable_id != new_stable_id or isnull(new_stable_id) )
     order by Y desc, X desc 
  )
    );

    $sth->execute();
    my %mapping = ();
    while ( my ( $type, $osi, $nsi ) = $sth->fetchrow_array() ) {
        next
          if $current_stable_ids{$type}{ $osi
              };    ## Don't need to cope with current IDS already searchable...
        $mapping{$type}{$osi}{$nsi} = 1;
        if ( $mapping{$type}{$nsi} ) {
            foreach ( keys %{ $mapping{$type}{$nsi} } ) {
                $mapping{$type}{$osi}{$_} = 1;
            }
        }
    }

    foreach my $type ( keys %mapping ) {
        $ecount += scalar keys %{ $mapping{$type} }, '  ';

        foreach my $osi ( keys %{ $mapping{$type} } ) {

            my @current_sis    = ();
            my @deprecated_sis = ();
            foreach ( keys %{ $mapping{$type}{$osi} } ) {
                if ( $current_stable_ids{$_} ) {
                    push @current_sis, $_;
                }
                elsif ( $_ ne 'NULL' ) {
                    push @deprecated_sis, $_;
                }
            }
            if (@current_sis) {

                my $description =
qq{$type $osi is no longer in the Ensembl database but it has been mapped to the following current identifiers: @current_sis}
                  . (
                    @deprecated_sis
                    ? "; and the following deprecated identifiers: @deprecated_sis"
                    : ''
                  );
                p unmappedGeneXML( $osi, $dbspecies, $description, lc($type) );

            }
            elsif (@deprecated_sis) {

                my $description =
qq($type $osi is no longer in the Ensembl database but it has been mapped to the following identifiers: @deprecated_sis);
                p unmappedGeneXML( $osi, $dbspecies, $description, lc($type) );
            }
            else {

                my $description =
qq($type $osi is no longer in the Ensembl database and it has not been mapped to any newer identifiers);
                p unmappedGeneXML( $osi, $dbspecies, $description, lc($type) );
            }
        }
    }

    footer($ecount);
}

sub unmappedGeneXML {
    my ( $id, $dbspecies, $description, $type ) = @_;

    return qq{
 <entry id="$id">
   <name>$id</name>
    <description>$description</description>
    <additional_fields>
      <field name="species">$dbspecies</field>
      <field name="featuretype">Unmapped $type</field>
    </additional_fields>
 </entry>};

}

sub make_counter {
    my $start = shift;
    return sub { $start++ }
}

