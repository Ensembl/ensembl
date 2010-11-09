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
$port    ||= 3306;
usage() and exit unless ( $host && $port && $user );

my $entry_count;
my $global_start_time = time;
my $total             = 0;
my $FAMILY_DUMPED;

my $fh;
## HACK 1 - if the INDEX is set to all grab all dumper methods...
my @indexes = split ',', $ind;
@indexes = map { /dump(\w+)/ ? $1 : () } keys %ebi_search_dump::
  if $ind eq 'ALL';

#warn Dumper \@indexes;

my $dbHash = get_databases();

#warn Dumper $dbcHash;

foreach my $species ( sort keys %$dbHash ) {

    my $conf = $dbHash->{$species};
    foreach my $index (@indexes) {

        # we don't dump compara anymore
        next if $index =~ /Family/;
        my $function = "dump$index";
        no strict "refs";

        $species =~ s/_/ /g;
        if ( $index eq 'Gene' ) {
            &$function( ucfirst($species), $conf );
            print $function, "\n";
        }    #elsif ($index eq 'Family' && !$FAMILY_DUMPED) {
             # &dumpFamily($conf);

        #}

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

sub dumpGene {

    my ( $dbspecies, $conf ) = @_;

    my $orth_species = {
        'homo_sapiens'             => 'ensembl_ortholog',
        'mus_musculus'             => 'ensembl_ortholog',
        'drosophila_melanogaster'  => 'ensemblgenomes_ortholog',
        'caenorhabditis_elegans'   => 'ensemblgenomes_ortholog',
        'saccharomyces_cerevisiae' => 'ensemblgenomes_ortholog'
    };
    my $compara_sth;
    my $want_species_orthologs;

    my $COMPARA_DB_NAME = $conf->{compara}->{$release}
      || die "can't get compara dbname";
    my $ortholog_lookup;
    my $orth_target_species;
    ( $orth_target_species = lcfirst($dbspecies) ) =~ s/\s/_/;
    if ( $want_species_orthologs =
        delete( $orth_species->{$orth_target_species} ) )
    {

        print
"Fetching Orthogs [$orth_target_species]\n-------------------------\n";

        #die Dumper($conf);
        my $compara_dsn = "DBI:mysql:host=$host";
        $compara_dsn .= ";port=$port" if ($port);
        my $compara_db_name = $conf->{compara}->{$release}
          || die "can't get compara dbname";

        my $compara_dbh =
          DBI->connect( "$compara_dsn:$compara_db_name", $user, $pass )
          or die "DBI::error";

        #                  $compara_dbh->{TraceLevel} = "1|SQL|test";

        my @wanted_ortholog_species = keys %$orth_species;

        #  change perl array value interpolation
        local $" = qq{","};
        my $orth_species_string = qq/@wanted_ortholog_species/;

        # Get orthologs
        my $orthologs_sth = $compara_dbh->prepare(
            qq{SELECT
 m1.stable_id , m2.stable_id, gdb2.name
FROM 
 genome_db gdb1  JOIN member m1 USING (genome_db_id)
 JOIN homology_member hm1 USING (member_id)
 JOIN homology h USING (homology_id)
 JOIN homology_member hm2 USING (homology_id)
 JOIN member m2 ON (hm2.member_id = m2.member_id)
 JOIN genome_db gdb2 on (m2.genome_db_id = gdb2.genome_db_id)
WHERE
 gdb1.name = "$orth_target_species" 
 AND m2.source_name = "ENSEMBLGENE"
 AND gdb2.name IN ("$orth_species_string")
 AND h.description in ("ortholog_one2one", "apparent_ortholog_one2one",
   "ortholog_one2many", "ortholog_many2many")}
        );

        $orthologs_sth->execute;
        my $rows = [];    # cache for batches of rows
        while (
            my $row = (
                shift(@$rows) ||    # get row from cache, or reload cache:
                  shift(
                    @{
                        $rows =
                          $orthologs_sth->fetchall_arrayref( undef, 10_000 )
                          || []
                      }
                  )
            )
          )
        {
            push @{ $ortholog_lookup->{ $row->[0] } },
              [ $row->[1], $orth_species->{ $row->[2] } ];

        }
        print "Done Fetching Orthogs\n-------------------------\n";
    }

    #     foreach my $DB ( 'core', 'otherfeatures', 'vega' ) {
    foreach my $DB ( 'core', 'vega' ) {
        my $counter = make_counter(0);
        my $SNPDB   = eval { $conf->{variation}->{$release} };
        my $DBNAME  = $conf->{$DB}->{$release}
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

        # SNP query
        my $snp_sth = eval {
            $dbh->prepare(
"select distinct(vf.variation_name) from $SNPDB.transcript_variation as tv, $SNPDB.variation_feature as vf where vf.variation_feature_id = tv.variation_feature_id and tv.transcript_stable_id in(?)"
            );
        };

        my $haplotypes = $dbh->selectall_hashref(
            "select gene_id from gene g, assembly_exception ae where
g.seq_region_id=ae.seq_region_id and ae.exc_type='HAP'", [qw(gene_id)]
        ) or die $dbh->errstr;

        my $taxon_id = $dbh->selectrow_arrayref(
            "select meta_value from meta where meta_key='species.taxonomy_id'");

        my %xrefs      = ();
        my %xrefs_desc = ();
        my %disp_xrefs = ();
        foreach my $type (qw(Gene Transcript Translation)) {
            my $T = $dbh->selectall_arrayref(
                "select ox.ensembl_id,
                x.display_label, x.dbprimary_acc, ed.db_name, es.synonym, x.description
           from ($DBNAME.object_xref as ox, $DBNAME.xref as x, $DBNAME.external_db as ed) left join $DBNAME.external_synonym as es on es.xref_id = x.xref_id
          where ox.ensembl_object_type = '$type' and ox.xref_id = x.xref_id and x.external_db_id = ed.external_db_id"
            );
            print "XREF $type query...\n---------------------\n";
            foreach (@$T) {

                $xrefs{$type}{ $_->[0] }{ $_->[3] }{ $_->[1] } = 1 if $_->[1];
                $xrefs{$type}{ $_->[0] }{ $_->[3] }{ $_->[2] } = 1 if $_->[2];
                $xrefs{$type}{ $_->[0] }{ $_->[3] . "_synonym" }{ $_->[4] } = 1
                  if $_->[4];
                $xrefs_desc{$type}{ $_->[0] }{ $_->[5] } = 1 if $_->[5];

            }

            print "Done XREF $type query...\n---------------------\n";
        }

        # my %exons = ();
        # my $T     = $dbh->selectall_arrayref(
        #     "select distinct t.gene_id, esi.stable_id
        #  from transcript as t, exon_transcript as et, exon_stable_id as esi
        # where t.transcript_id = et.transcript_id and et.exon_id = esi.exon_id"
        # );
        # foreach (@$T) {
        #     $exons{ $_->[0] }{ $_->[1] } = 1;
        # }

        print "Get Genes query...\n---------------------\n";

        my %exons = ();
        my $get_genes_sth    = $dbh->prepare(
            "select distinct t.gene_id, esi.stable_id
         from transcript as t, exon_transcript as et, exon_stable_id as esi
        where t.transcript_id = et.transcript_id and et.exon_id = esi.exon_id"
        );


	$get_genes_sth->execute;
        my $gene_rows = [];    # cache for batches of rows

	while (
            my $row = (
                shift(@$gene_rows) ||    # get row from cache, or reload cache:
                  shift(
                    @{
                        $gene_rows =
                          $get_genes_sth->fetchall_arrayref( undef, 10_000 )
                          || []
                      }
                  )
            )
          )
        {
            push @{ $ortholog_lookup->{ $row->[0] } },
              [ $row->[1], $orth_species->{ $row->[2] } ];

            $exons{ $row->[0] }{ $row->[1] } = 1;

        }
    

        my $gene_info = $dbh->selectall_arrayref( "
        select gsi.gene_id, tsi.transcript_id, trsi.translation_id,
             gsi.stable_id as gsid, tsi.stable_id as tsid, trsi.stable_id as trsid,
             g.description, ed.db_name, x.dbprimary_acc,x.display_label, ad.display_label, ad.description, g.source, g.status, g.biotype
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

        print "Done Get Genes query...\n---------------------\n";

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

                    if ( $SNPDB && $DB eq 'core' ) {
                        my @transcript_stable_ids =
                          keys %{ $old{transcript_stable_ids} };
                        $snp_sth->execute("@transcript_stable_ids");
                        $old{snps} = $snp_sth->fetchall_arrayref;
                    }

                    if ($want_species_orthologs) {
                        $old{orthologs} =
                          $ortholog_lookup->{ $old{'gene_stable_id'} };
                    }

                    p geneLineXML( $dbspecies, \%old, $counter );

                }
                %old = (
                    'gene_id'   => $gene_id,
                    'haplotype' => $haplotypes->{$gene_id} ? 'haplotype'
                    : 'reference',
                    'gene_stable_id'         => $gene_stable_id,
                    'description'            => $gene_description,
                    'taxon_id'               => $taxon_id->[0],
                    'translation_stable_ids' => {
                        $translation_stable_id ? ( $translation_stable_id => 1 )
                        : ()
                    },
                    'transcript_stable_ids' => {
                        $transcript_stable_id ? ( $transcript_stable_id => 1 )
                        : ()
                    },
                    'transcript_ids' => {
                        $transcript_id ? ( $transcript_id => 1 )
                        : ()
                    },
                    'exons'                => {},
                    'external_identifiers' => {},
                    'alt'                  => $xref_display_label
                    ? "($extdb_db_display_name: $xref_display_label)"
                    : "(novel gene)",
                    'gene_name' => $xref_display_label ? $xref_display_label
                    : $gene_stable_id,
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
                $old{'transcript_ids'}{$transcript_id}                 = 1;
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

        if ( $SNPDB && $DB eq 'core' ) {
            my @transcript_stable_ids = keys %{ $old{transcript_stable_ids} };
            $snp_sth->execute("@transcript_stable_ids");
            $old{snps} = $snp_sth->fetchall_arrayref;
        }
        if ($want_species_orthologs) {
            $old{orthologs} = $ortholog_lookup->{ $old{'gene_stable_id'} };

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

    my $snps      = $xml_data->{'snps'};
    my $orthologs = $xml_data->{'orthologs'};

    my $peptides = $xml_data->{'translation_stable_ids'}
      or die "peptides not set";
    my $exons = $xml_data->{'exons'} or die "exons not set";
    my $external_identifiers = $xml_data->{'external_identifiers'}
      or die "external_identifiers not set";
    my $description = $xml_data->{'description'};
    my $gene_name   = $xml_data->{'gene_name'};
    my $type        = $xml_data->{'source'} . ' ' . $xml_data->{'biotype'}
      or die "problem setting type";
    my $haplotype        = $xml_data->{'haplotype'};
    my $taxon_id         = $xml_data->{'taxon_id'};
    my $exon_count       = scalar keys %$exons;
    my $transcript_count = scalar keys %$transcripts;
    $description =~ s/</&lt;/g;
    $description =~ s/>/&gt;/g;
    $description =~ s/'/&apos;/g;
    $description =~ s/&/&amp;/g;

    $gene_name =~ s/</&lt;/g;
    $gene_name =~ s/>/&gt;/g;
    $gene_name =~ s/'/&apos;/g;
    $gene_name =~ s/&/&amp;/g;

    $gene_id =~ s/</&lt;/g;
    $gene_id =~ s/>/&gt;/g;

    $altid =~ s/</&lt;/g;
    $altid =~ s/>/&gt;/g;

    my $xml = qq{
 <entry id="$gene_id">
   <name>$gene_id $altid</name>
    <description>$description</description>};

    my $synonyms = "";
    my $unique_synonyms;
    my $cross_references = qq{
       <cross_references>};
    $cross_references .= qq{
         <ref dbname="ncbi_taxonomy_id" dbkey="$taxon_id"/>};

    # for some types of xref, merge the subtypes into the larger type
    # e.g. Uniprot/SWISSPROT and Uniprot/TREMBL become just Uniprot
    # synonyms are stored as additional fields rather than cross references
    foreach my $ext_db_name ( keys %$external_identifiers ) {

        if ( $ext_db_name =~
            /(Uniprot|GO|Interpro|Medline|Sequence_Publications|EMBL)/ )
        {

            my $matched_db_name = $1;

            # synonyms
            if ( $ext_db_name =~ /_synonym/ ) {

                foreach
                  my $ed_key ( keys %{ $external_identifiers->{$ext_db_name} } )
                {

                    #		$unique_synonyms->{$ed_key} = 1;
                    $synonyms .= qq{ 
             <field name="${matched_db_name}_synonym">$ed_key</field>};
                }

            }
            else {    # non-synonyms

                map {
                    $cross_references .= qq{
         <ref dbname="$matched_db_name" dbkey="$_"/>};
                  } keys %{ $external_identifiers->{$ext_db_name} }

            }

        }
        else {

            foreach my $key ( keys %{ $external_identifiers->{$ext_db_name} } )
            {

                $key         =~ s/</&lt;/g;
                $key         =~ s/>/&gt;/g;
                $key         =~ s/&/&amp;/g;
                $ext_db_name =~ s/^Ens.*/ENSEMBL/;

                if ( $ext_db_name =~ /_synonym/ ) {
                    $unique_synonyms->{$key} = 1;
                    $synonyms .= qq{
        <field name="$ext_db_name">$key</field>};

                }
                else {

                    $cross_references .= qq{
        <ref dbname="$ext_db_name" dbkey="$key"/>};

                }
            }

        }
    }

    $cross_references .= (
        join "",
        (
            map {
                qq{
      <ref dbname="ensemblvariation" dbkey="$_->[0]"/>}
              } @$snps
        )
    );

    $cross_references .= (
        join "",
        (
            map {
                qq{
      <ref dbname="$_->[1]" dbkey="$_->[0]"/>}
              } @$orthologs
        )
    );

    $cross_references .= qq{
</cross_references>};

    map {
        $synonyms .=
          qq{     
      <field name="gene_synonym">$_</field> }
    } keys %$unique_synonyms;

    my $additional_fields .= qq{
    <additional_fields>
      <field name="species">$species</field>
      <field name="featuretype">Gene</field>
      <field name="source">$type</field>
      <field name="transcript_count">$transcript_count</field>
      <field name="gene_name">$gene_name</field>
      <field name="haplotype">$haplotype</field>}
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
      . $synonyms

      . qq{
   </additional_fields>
};

    $counter->();
    return $xml . $cross_references . $additional_fields . '</entry>';

}

sub dumpGenomicAlignment {
    warn "in dump Genomic";
    my ( $dbspecies, $conf ) = @_;

    #    warn Dumper $conf;

    warn "\n", '*' x 20, "\n";
    my %tables = (
        'dna_align_feature' => [ 'DnaAlignFeature', 'DNA alignment feature' ],
        'protein_align_feature' =>
          [ 'ProteinAlignFeature', 'Protein alignment feature' ]
    );
    my $ecount;

    foreach my $db ( 'core', 'cdna', 'otherfeatures' ) {

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

            # Due to the sheer number of features - generating this
            # dump causes temp table to be written on the mysqld
            # filesys. /tmp space can be easily filled.
            # Have split the feature fetching to happen by an
            # analysis_id at a time.  The changes below gave a x3
            # speed up on XML dumping as compared to fetching all
            # features in one query at once.

            # make a lookup for the analysis display labels.
            my $type              = $tables{$table}[1];
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
                @{
                    $dbh->selectall_arrayref(
"select distinct distinct(analysis_id) from $DB_NAME.$table"
                    )
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
    <additional_fields>
        <field name="species">$dbspecies</field>
        <field name="featuretype">$source</field>
        <field name="db">$db</field>
        <field name="genome_hits">$count</field>
        <field name="adesc">$adesc</field>
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

sub dumpMarker {
    warn "in dump MArker";
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
          || die("can't open compressed stream to $file: $!");
    }
    else {
        open( file, ">$file" ) || die "can't open $file: $!";
    }
    header( $dbname, $dbspecies, $db );
    my $dsn = "dbi:mysql:host=$host";
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

    $desc =~ s/</&lt;/g;
    $desc =~ s/>/&gt;/g;

    $xml = qq{
<entry id="$marker">
   <additional_fields>};

    foreach (@keys) {
        s/</&lt;/g;
        s/>/&gt;/g;

        $xml .= qq{
      <field name="synonym">$_</field>}

    }
    $xml .= qq{
     <field name="species">$species</field>
    <field name="featuretype">Marker</field>
  </additional_fields>
</entry>};

    return $xml;

}

sub dumpOligoProbe {
    warn "in dump Oligo";
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

    while ( my $rowcache = $sth->fetchall_arrayref( undef, '10_000' ) ) {
        my $xml;
        while ( my $data = shift( @{$rowcache} ) ) {
            $xml .= OligoProbeXML( $data, $dbspecies );

        }
        p $xml;
    }

    footer( $sth->rows );

}

sub OligoProbeXML {
    my ( $xml_data, $dbspecies ) = @_;

#     my $desc =qq{$xml_data->[0], $xml_data->[2] oligo probeset $xml_data->[0] hits the genome in $xml_data->[1] locations.};

    return qq{
<entry id="$xml_data->[0]">
   <additional_fields>
     <field name="type">$xml_data->[2]</field>
     <field name="species">$dbspecies</field>
     <field name="featuretype">OligoProbe</field>
     <field name="genome_hits">$xml_data->[1]</field>
   </additional_fields>
</entry>};

}

sub dumpQTL {
    warn "in dumpQTL";

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

sub dumpSequence {
    warn "in dump Sequence";
    my ( $dbspecies, $conf ) = @_;

    #    my $sanger = sanger_project_names( $conf );
    my $sanger = 'SANGER STUFF';
    my %config = (
        "Homo sapiens" => [
            [
                'Clone',
                'tilepath, cloneset_1mb, cloneset_30k, cloneset_32k',
'name,well_name,clone_name,synonym,embl_acc,sanger_project,alt_well_name,bacend_well_name'
            ],
            [ 'NT Contig',     'ntctgs', 'name' ],
            [ 'Encode region', 'encode', 'name,synonym,description' ],
        ],
        "Mus musculus" => [
            [
                'BAC',
                'cloneset_0_5mb,cloneset_1mb,bac_map,tilingpath_cloneset',
                'embl_acc,name,clone_name,well_name,synonym,alt_embl_acc'
            ],
            [ 'Fosmid',      'fosmid_map', 'name,clone_name' ],
            [ 'Supercontig', 'superctgs',  'name' ],
        ],
        "Anopheles gambiae" => [
            [ 'BAC',      'bacs',       'name,synonym,clone_name' ],
            [ 'BAC band', 'bacs_bands', 'name,synonym,clone_name' ],
        ],
        "Gallus gallus" => [
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

    #     my $description = "$type $name is mapped to Chromosome $chr" .

    #       (
    #         @$val > 0
    #         ? ' and has '
    #           . @$val
    #           . " EMBL accession"
    #           . (
    #             @$val > 1
    #             ? 's'
    #             : ''
    #           )
    #           . "/synonym"
    #           . (
    #             @$val > 1
    #             ? 's '
    #             : ' '
    #           )
    #           . "@$val"
    #         : ''
    #       )

    #       . " and length $len bps\n";

    my $xml = qq{
 <entry id="$name">
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

sub dumpSNP {
    my ( $dbspecies, $conf ) = @_;

    #    warn Dumper $conf;

    warn "\n", '*' x 20, "\n";

    my $COREDB = my $dbname = $conf->{'core'}->{$release};

    my $dbname = $conf->{variation}->{$release} or next;
    my $file = "$dir/SNP_$dbname.xml";
    $file .= ".gz" unless $nogzip;
    my $start_time = time;
    warn "Dumping $dbname to $file ... ", format_datetime($start_time), "\n";
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
    my $source_hash =
      $dbh->selectall_hashref( qq{SELECT source_id, name FROM source},
        [qw(source_id)] );

#     my $tid_to_gene = $dbh->selectall_hashref(qq{select  t.transcript_id, gsi.stable_id from $COREDB.gene as g, $COREDB.gene_stable_id as gsi, $COREDB.transcript as t where gsi.gene_id = g.gene_id and t.gene_id = g.gene_id limit 10;},[qw(transcript_id)]);

#     my $sth = $dbh->prepare("select vf.variation_name, vf.source_id, group_concat(vs.source_id, ' ',vs.name), vf.variation_feature_id,vf.variation_id from variation_feature vf , transcript_variation tv
# left join variation_synonym vs on vf.variation_id = vs.variation_id where tv.transcript_variation_id = vf.variation_feature_id group by vf.variation_id");

    my $sth = $dbh->prepare(
"select vf.variation_name, vf.source_id, group_concat(vs.source_id, ' ',vs.name), vf.consequence_type from variation_feature vf left join variation_synonym vs on vf.variation_id = vs.variation_id group by vf.variation_id"
    );

#     my $vfi2gene_sth = $dbh->prepare(qq{select distinct(gsi.stable_id) from $COREDB.gene as g, $COREDB.gene_stable_id as gsi, $COREDB.transcript as t where gsi.gene_id = g.gene_id and t.gene_id = g.gene_id and transcript_id in
# (select tv.transcript_id from transcript_variation tv , variation_feature vf where vf.variation_feature_id =tv.variation_feature_id and vf.variation_feature_id = ?)});

    $sth->execute() or die "Error:", $DBI::errstr;

    while ( my $rowcache = $sth->fetchall_arrayref( undef, 10_000 ) ) {

        my $xml;
        while ( my $row = shift( @{$rowcache} ) ) {

            #	    $vfi2gene_sth->execute($row->[3]);
            #	      my $gsi = $vfi2gene_sth->fetchall_arrayref;
            my $name       = $row->[0];
            my @synonyms   = split /,/, @$row->[2];
            my $snp_source = $source_hash->{ $row->[1] }->{name};

#	    my $description =
#	      "A $snp_source SNP with "
#		. scalar @synonyms
#		  . ' synonym'
#		    . (  @synonyms > 1 | @synonyms < 1 ? 's '   : ' ' )
#		      . ( @synonyms > 0 ? "( " . (join "",  map{  map{  $source_hash->{$_->[0]}->{name} , ':', $_->[1] , ' ' } [split]  } @synonyms ) . ")" : '' );

            $xml .= qq{<entry id="$name">
  <additional_fields>
    <field name="species">$dbspecies</field>
    <field name="featuretype">SNP</field>
    <field name="consequence">$row->[3]</field>};

            foreach my $syn (@synonyms) {
                my @syn_bits = split / /, $syn;
                $syn_bits[1] =~ s/:/ /;

                my $source = $source_hash->{ $syn_bits[0] }->{name};
                $xml .= qq{
     <field name="synonym">$syn_bits[1] [source; $source]</field>};
            }
            $xml .= qq{
  </additional_fields>
</entry>
};

        }

        p($xml);
    }

    footer( $sth->rows );
    print_time($start_time);

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
    <description>$description</description>
    <additional_fields>
      <field name="species">$dbspecies</field>
      <field name="featuretype">Unmapped$type</field>
    </additional_fields>
 </entry>};

}

sub make_counter {
    my $start = shift;
    return sub { $start++ }
}

sub FamilyDumped {
    my $is_dumped;
    return sub { $is_dumped }
}
