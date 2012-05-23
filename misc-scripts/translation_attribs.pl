#!/software/bin/perl

=head1 NAME

translation_attribs.pl - script to calculate peptide statistics and store 
                         them in translation_attrib table

=head1 SYNOPSIS

translation_attribs.pl [arguments]

Required arguments:

  --user=user                         username for the database

  --pass=pass                         password for database


Optional arguments:

  --pattern=pattern                   calculate translation attribs for databases matching pattern
                                      Note that this is a standard regular expression of the
                                      form '^[a-b].*core.*' for all core databases starting with a or b

  --binpath=PATH                      directory where the binary script to calculate 
                                      pepstats is stored (default: /software/pubseq/bin/emboss)

  --tmpfile=file                      file to store tmp results of pepstats (default=/tmp)

  --host=host                         server where the core databases are stored (default: ens-staging)

  --dbname=dbname                     if you want a single database to calculate the pepstats
                                      (all databases by default)

  --port=port                         port (default=3306)

  --help                              print help (this message)

=head1 DESCRIPTION

This script will calculate the peptide statistics for all core databases in the server 
and store them as a translation_attrib values

=head1 EXAMPLES

Calculate translation_attributes for all core databases in ens-staging 

  $ ./translation_attribs.pl --user ensadmin --pass password

Calculate translation_attributes for core databases starting with [a-c] in ens-staging (output LSF to PWD) 

  $ ./translation_attribs.pl --user ensadmin --pass password --pattern '^[a-c].*core_50.*'

Calculate translation_attribs for a single database in a ens-genomics1

  $ ./translation_attribs.pl  --host ens-genomics1 --user ensadmin --pass password --dbname my_core_db

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Daniel Rios <dani@ebi.ac.uk>, Ensembl core API team

=head1 CONTACT

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Attribute;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;
use DBI;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::CliHelper;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

# get the basic options for connecting to a database server
my $optsd = $cli_helper->get_dba_opts();
# add the print option
push( @{$optsd}, "binpath:s" );
push( @{$optsd}, "tmpdir:s" );
push( @{$optsd}, "verbose|v" );
# process the command line with the supplied options plus a help subroutine
my $opts = $cli_helper->process_args( $optsd, \&pod2usage );

$opts->{binpath} ||= '/software/pubseq/bin/emboss';
$opts->{tmpdir}  ||= '/tmp';
$opts->{port}    ||= '3306';
$opts->{host}    ||= 'ens-staging';
$opts->{user}    ||= 'ensro';
$opts->{pattern} ||= qr/.+_core_.+/;

my %PEPSTATS_CODES = ( 'Number of residues' => 'NumResidues',
		       'Molecular weight'   => 'MolecularWeight',
		       'Ave. residue weight' => 'AvgResWeight',
		       'Charge' => 'Charge',
		       'Isoelectric point' => 'IsoPoint'
		      );

my %MET_AND_STOP = ( 'Starts with methionine' => 'starts_met', 
		     'Contains stop codon' => 'has_stop_codon'
		    );



my %attributes_to_delete = (%PEPSTATS_CODES);
my $translation_attribs = {};
my $translation;
my $dbID;

# use the command line options to get an array of database details
for my $db_args ( @{ $cli_helper->get_dba_args_for_opts($opts) } ) {

	# use the args to create a DBA
	my $dba    = Bio::EnsEMBL::DBSQL::DBAdaptor->new( %{$db_args} );
	my $dbname = $dba->dbc()->dbname();
	print STDOUT "Processing species "
	  . $dba->species_id()
	  . " from database "
	  . $dba->dbc()->dbname()
	  . " on server "
	  . $dba->dbc()->host() . "\n";

	if ( $dbname =~ /(vega|otherfeatures)/ ) {
		my $other_dbname = $dbname;
		$other_dbname =~ s/$1/core/;

		$opts->{dbname} = $other_dbname;
		#for vega databases, add the core as the dna database
		my $core_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new( %{$db_args} );
		$dba->dnadb($core_db);
	}

	print "Removing attributes from database ", $dba->dbc->dbname, "\n";
	remove_old_attributes( $dba, \%attributes_to_delete);
	my $translationAdaptor = $dba->get_TranslationAdaptor();
	my $attributeAdaptor   = $dba->get_AttributeAdaptor();
	print
	  "Going to update translation_attribs in",
	  $dba->dbc->dbname, "\n";

#for all the translations in the database, run pepstats and update the translation_attrib table
	my @translation_ids = @{
		$dba->dbc()->sql_helper()->execute_simple(
-SQL=>"SELECT tl.translation_id FROM translation tl, transcript tr, seq_region s, coord_system c WHERE tl.transcript_id = tr.transcript_id AND tr.seq_region_id = s.seq_region_id AND s.coord_system_id = c.coord_system_id AND c.species_id = ? order by tl.translation_id",
			-PARAMS=>[$dba->species_id()] ) };
	my $translations = {};
	my $tmpfile      = $opts->{tmpdir} . "/$$.pep";
	open( TMP, "> $tmpfile" ) || warn "PEPSTAT: $!";
	print "Retrieving translations\n";
	for my $dbID (@translation_ids) {

		#foreach translation, retrieve object
		my $translation = $translationAdaptor->fetch_by_dbID($dbID);
		if ( $opts->{verbose} ) {
			print "Dumping translation dbID, $dbID...\n";
		}
		$translations->{$dbID} = $translation;
		my $peptide_seq = $translation->seq();
		if ( !( $peptide_seq =~ m/[BZX]/ig ) ) {
			if ( $peptide_seq !~ /\n$/ ) { $peptide_seq .= "\n" }
			$peptide_seq =~ s/\*$//;
			print TMP ">$dbID\n$peptide_seq";
		} else {
			print "Skipping translation dbID $dbID due to ambiguity codes...\n";
		}
	}
	close(TMP);
	print "Running pepstat\n";

	my $PEPSTATS = $opts->{binpath} . '/bin/pepstats';
	throw("pepstats binary at  $PEPSTATS cannot be executed")
	  if ( !-x $PEPSTATS );
	open( OUT, "$PEPSTATS -filter < $tmpfile 2>&1 |" )
	  || warn "PEPSTAT: $!";
	my @lines = <OUT>;
	close(OUT);
	unlink($tmpfile);
	my $attribs = {};
	my $tId;
	print "Parsing pepstat\n";

	foreach my $line (@lines) {
		if ( $line =~ /PEPSTATS of ([^ ]+)/ ) {
			$tId = $1;
		} elsif ( defined $tId ) {
			if ( $line =~ /^Molecular weight = (\S+)(\s+)Residues = (\d+).*/ ) {
				$attribs->{$tId}{'Number of residues'} = $3;
				$attribs->{$tId}{'Molecular weight'}   = $1;
			} elsif ( $line =~
/^Average(\s+)(\S+)(\s+)(\S+)(\s+)=(\s+)(\S+)(\s+)(\S+)(\s+)=(\s+)(\S+)/ )
			{
				$attribs->{$tId}{'Ave. residue weight'} = $7;
				$attribs->{$tId}{'Charge'}              = $12;
			} elsif ( $line =~ /^Isoelectric(\s+)(\S+)(\s+)=(\s+)(\S+)/ ) {
				$attribs->{$tId}{'Isoelectric point'} = $5;
			} elsif ( $line =~ /FATAL/ ) {
				print STDERR "pepstats: $line\n";
			}
		}
	}
	for my $id ( keys %$attribs ) {
		my $translation = $translations->{$id};
		my $aas         = $attribs->{$id};
		if ( $opts->{verbose} ) {
			print "Storing attribs for translation dbID, $id...\n";
		}
		store_translation_attrib( $attributeAdaptor, $translation, $aas );
	}
} ## end for my $db_args ( @{ $cli_helper...})

#will remove any entries in the translation_attrib table for the attributes, if any
#this method will try to remove the old starts_met and has_stop_codon attributes, if present
#this is to allow to be run on old databases, but removing the not used attributes
sub remove_old_attributes {
	my $dba        = shift;
	my $attributes = shift;
		print "removing all translation attributes for db, "
		  . $dba->{_dbc}->{_dbname} . "\n";
		foreach my $value ( values %{$attributes} ) {
		    my $sql = "delete ta FROM translation_attrib ta, attrib_type at, translation tl, transcript tr, seq_region s, coord_system c  WHERE at.attrib_type_id = ta.attrib_type_id AND ta.translation_id=tl.translation_id and  tl.transcript_id = tr.transcript_id AND tr.seq_region_id = s.seq_region_id AND s.coord_system_id = c.coord_system_id AND c.species_id = ? and at.code=?";
			$dba->dbc()->sql_helper()->execute_update(-SQL=>$sql,-PARAMS=>[$dba->species_id(),$value]);
		}
		return;

} ## end sub remove_old_attributes

sub store_translation_attrib {
	my ( $attributeAdaptor, $translation, $attribs ) = @_;
	my @attribs;
	for my $key ( keys %$attribs ) {
		my $value = $attribs->{$key};
		push @attribs,
		  Bio::EnsEMBL::Attribute->new( '-code'  => $PEPSTATS_CODES{$key},
										'-name'  => $key,
										'-value' => $value );
	}
	$attributeAdaptor->store_on_Translation( $translation, \@attribs );
	return;
}

