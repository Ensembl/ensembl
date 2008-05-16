#!/software/bin/perl

=head1 NAME

translation_attribs.pl - script to calculate peptide statistics, if the first aminoacid is
                         methionine and there is a stop codon in the aminoacid sequence
                         and store them in translation_attrib table

=head1 SYNOPSIS

translation_attribs.pl [arguments]

Required arguments:

  --user=user                         username for the database

  --pass=pass                         password for database

Optional arguments:

  --binpath=PATH                      directory where the binary script to calculate 
                                      pepstats is stored (default: /software/pubseq/bin/emboss)

  --tmpfile=file                      file to store tmp results of pepstats (default=/tmp)

  --host=host                         server where the core databases are stored (default: ens-staging)

  --dbname=dbname                     if you want a single database to calculate the pepstats
                                      (all databases by default)

  --port=port                         port (default=3306)

  --pepstats_only                     when used, will only run the pepstats calculation

  --met_and_stop_only                 when used, will only run the methionine and stop codon calculation
  
  --help                              print help (this message)

=head1 DESCRIPTION

This script will calculate the peptide statistics, if the first aminoacid is methionine and
there is a stop codon  for all core databases in the server and store them as a 
translation_attrib values

=head1 EXAMPLES

Calculate translation_attributes for all databases in ens-staging (default in the release process)

  $ ./translation_attribs.pl --user ensadmin --pass password

Calculate pepstats for a single database in a ens-genomics1

  $ ./translation_attribs.pl --binpath /ensemblweb/shared/bin/ --tmpfile file --host ens-genomics1 --user ensadmin --pass password --dbname my_core_db --pepstats_only

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

use Bio::EnsEMBL::Utils::Exception qw(throw);

##global variable containing all possible pepstats and the codes used

my %PEPSTATS_CODES = ( 'Number of residues' => 'NumResidues',
		       'Molecular weight'   => 'MolecularWeight',
		       'Ave. residue weight' => 'AvgResWeight',
		       'Charge' => 'Charge',
		       'Isoelectric point' => 'IsoPoint'
		      );

my %MET_AND_STOP = ( 'Initial methionine' => 'initial_met', 
		     'Contains stop codon' => 'has_stop_codon'
		    );


## Command line options

my $binpath = '/software/pubseq/bin/emboss'; 
my $tmpdir = '/tmp';
my $host = 'ens-staging';
my $dbname = undef;
my $user = undef;
my $pass = undef;
my $port = 3306;
my $help = undef;
my $pepstats_only = undef;
my $met_and_stop_only = undef;


GetOptions('binpath=s' => \$binpath,
	   'tmpdir=s' => \$tmpdir,
	   'host=s'    => \$host,
	   'dbname=s'  => \$dbname,
	   'user=s'    => \$user,
	   'pass=s'    => \$pass,
	   'port=s'    => \$port,
	   'pepstats_only' => \$pepstats_only,
	   'met_and_stop_only' => \$met_and_stop_only,
	   'help'    => \$help
	   );

pod2usage(1) if($help);
throw("--user argument required") if (!defined($user));
throw("--pass argument required") if (!defined($pass));

my $dbas;
#load registry with all databses when no database defined
if (!defined ($dbname)){
    Bio::EnsEMBL::Registry->load_registry_from_db(-host => $host,
						  -user => $user,
						  -pass => $pass,
						  -port => $port
						  );
      $dbas = Bio::EnsEMBL::Registry->get_all_DBAdaptors(-group=>'core'); #get all core adaptors for all species
  }
else{
#only get a single DBAdaptor, the one for the database specified
    my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-host => $host,
						  -user => $user,
						  -pass => $pass,
						  -port => $port,
						  -dbname => $dbname
						  );
    push @{$dbas},$dba;
}



my %attributes; #hash containing attributes to be stored and removed from the databse
my $run_both = undef;
if (!defined($pepstats_only) and !defined($met_and_stop_only)){
#this is the default option, calculate all attributes
    %attributes = (%PEPSTATS_CODES,%MET_AND_STOP);
    $run_both = 1;
}
if (defined($pepstats_only)){
    #will only store and remove the peptide stats
    %attributes = %PEPSTATS_CODES;
}
if (defined($met_and_stop_only)){
    #will only store and remove the methionince and stop codon
    %attributes = %MET_AND_STOP;
}

my $translation_attribs = {};
my $translation;
my $dbID;
#foreach of the species, calculate the pepstats
foreach my $dba (@{$dbas}){
    next if (defined $dbname and $dba->dbc->dbname ne $dbname);
    print "Removing attributes from database ", $dba->dbc->dbname,"\n";
    remove_old_attributes($dba,\%attributes);

    my $translationAdaptor = $dba->get_TranslationAdaptor();
    my $transcriptAdaptor = $dba->get_TranscriptAdaptor();
    my $attributeAdaptor = $dba->get_AttributeAdaptor();
    print "Going to update translation_attribs for ", $dba->dbc->dbname,"\n";
    #for all the translations in the database, run pepstats and update the translation_attrib table
    my $sth = $dba->dbc->prepare("SELECT translation_id from translation"); 
    $sth->execute();
    $sth->bind_columns(\$dbID);
    while($sth->fetch()){
	#foreach translation, retrieve object
	$translation = $translationAdaptor->fetch_by_dbID($dbID);
	#calculate pepstats
	get_pepstats($translation,$binpath,$tmpdir,$translation_attribs) if($run_both || $pepstats_only);
	#calculate methionine and stop codon
	get_met_and_stop($translation,$translation_attribs,$transcriptAdaptor) if($run_both || $met_and_stop_only);
	#and store results in database
	store_translation_attribs($attributeAdaptor,$translation_attribs,$translation,\%attributes);    	
	$translation_attribs = {};
    }        
}

#will remove any entries in the translation_attrib table for the attributes, if any
sub remove_old_attributes{
    my $dba = shift;
    my $attributes = shift;

    my $sth = $dba->dbc()->prepare("DELETE ta FROM translation_attrib ta, attrib_type at WHERE at.attrib_type_id = ta.attrib_type_id AND at.code = ?");
    #remove all possible entries in the translation_attrib table for the attributes
    foreach my $value (values %{$attributes}){
	$sth->execute($value);
    }
    $sth->finish;
}

#method that retrieves the pepstatistics for a translation

sub get_pepstats {
    my $translation = shift;
    my $binpath = shift;
    my $tmpdir = shift;
    my $translation_attribs = shift;

    
    my $peptide_seq ;
    eval { $peptide_seq = $translation->seq ; };
    return {} if ($@ || $peptide_seq =~ m/[BZX]/ig);
    if( $peptide_seq !~ /\n$/ ){ $peptide_seq .= "\n" }
    $peptide_seq =~ s/\*$//;
    
    my $tmpfile = $tmpdir."/$$.pep";
    open( TMP, "> $tmpfile" ) || warn "PEPSTAT: $!";
    print TMP "$peptide_seq";
    close(TMP);
    my $PEPSTATS = $binpath.'/bin/pepstats';
    open (OUT, "$PEPSTATS -filter < $tmpfile 2>&1 |") || warn "PEPSTAT: $!";
    my @lines = <OUT>;
    close(OUT);
    unlink($tmpfile);
    foreach my $line (@lines){
	if($line =~ /^Molecular weight = (\S+)(\s+)Residues = (\d+).*/){
	    $translation_attribs->{'Number of residues'} = $3 ;
	    $translation_attribs->{'Molecular weight'} = $1;
	}
	if($line =~ /^Average(\s+)(\S+)(\s+)(\S+)(\s+)=(\s+)(\S+)(\s+)(\S+)(\s+)=(\s+)(\S+)/){
	    $translation_attribs->{'Ave. residue weight'} = $7;
	    $translation_attribs->{'Charge'} = $12;
	}
	if($line =~ /^Isoelectric(\s+)(\S+)(\s+)=(\s+)(\S+)/){
	    $translation_attribs->{'Isoelectric point'} = $5;
	}
	if ($line =~ /FATAL/){            
	    print STDERR "pepstats: $line\n";
	    $translation_attribs = {};
	}
    }
}

#method to calculate if the first aminoacid is methionine and there is a stop codon
sub get_met_and_stop{
    my $translation = shift;
    my $translation_attribs = shift;
    my $transcriptAdaptor = shift;
    
    my $peptide_seq ;
    eval { $peptide_seq = $translation->seq ; };
    return {} if ($@ || $peptide_seq =~ m/[BZX]/ig);
    if( $peptide_seq !~ /\n$/ ){ $peptide_seq .= "\n" }

    #need to get transleatable_seq to find out if there is a stop codon
    my $transcript = $transcriptAdaptor->fetch_by_translation_id($translation->dbID);

    my $translateable_seq = $transcript->translateable_seq();

    if ($translateable_seq =~ /TAG$|TGA$|TAA$/i){
	#there is stop codon
	$translation_attribs->{'Contains stop codon'} = 1;
    }
    else{
	$translation_attribs->{'Contains stop codon'} = 0;
    }
    if ($peptide_seq =~ /^M/){
	$translation_attribs->{'Initial methionine'} = 1;
    }
    else{
	$translation_attribs->{'Initial methionine'} = 0;
    }

}

sub store_translation_attribs{
    my $attributeAdaptor = shift;
    my $translation_attribs = shift;
    my $translation = shift;
    my $attributes = shift;

    my $attribute;
    my @attributes;
    #each of the keys in the pepstats is an attribute for the translation
    foreach my $key (keys %{$translation_attribs}){

	$attribute = Bio::EnsEMBL::Attribute->new('-code' => $attributes->{$key},
						  '-name' => $key,
						  '-value' => $translation_attribs->{$key}
						  );
	push @attributes, $attribute;

    }
    $attributeAdaptor->store_on_Translation($translation,\@attributes);
}
