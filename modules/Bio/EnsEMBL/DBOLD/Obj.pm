#
# BioPerl module for DBOLD::Obj
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBOLD::Obj - Object representing an instance of an EnsEMBL DB

=head1 SYNOPSIS

    $db = Bio::EnsEMBL::DBOLD::Obj->new(
        -user   => 'root',
        -dbname => 'pog',
        -host   => 'caldy',
        -driver => 'mysql',
        );

    $clone  = $db->get_clone('X45667');

    $contig = $db->get_Contig("dJ52N12.02793");

    $gene   = $db->get_Gene('HG45501');

    

=head1 DESCRIPTION

This object represents a database that is implemented somehow (you shouldn\'t
care much as long as you can get the object). From the object you can pull
out other objects by their stable identifier, such as Clone (accession number),
Exons, Genes and Transcripts. The clone gives you a DB::Clone object, from
which you can pull out associated genes and features. 

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBOLD::Obj;

use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::Root::RootI;

use Bio::EnsEMBL::DB::ObjI;
use Bio::EnsEMBL::DBOLD::Gene_Obj;
use Bio::EnsEMBL::DBOLD::Feature_Obj;
use Bio::EnsEMBL::DBOLD::RawContig;
use Bio::EnsEMBL::DBOLD::GapContig;

use Bio::EnsEMBL::Clone;
use Bio::EnsEMBL::DBOLD::StaticGoldenPathAdaptor;
use Bio::EnsEMBL::DBOLD::KaryotypeBandAdaptor;
use Bio::EnsEMBL::DBOLD::AnalysisAdaptor;
use Bio::EnsEMBL::FeatureFactory;
use Bio::EnsEMBL::Chromosome;
use Bio::EnsEMBL::DBOLD::CloneAdaptor;

use DBI;

use Bio::EnsEMBL::DBOLD::SQL;
use Bio::EnsEMBL::DB::ObjI;

use Bio::EnsEMBL::DBOLD::DummyStatement;

@ISA = qw(Bio::EnsEMBL::DB::ObjI Bio::Root::RootI);

sub new {
  my($pkg, @args) = @_;

  my $self = bless {}, $pkg;
  $self->warn( "DBAdaptor is now copy of Obj. Please start using DBAdaptor instead of Obj" );

    my (
        $db,
        $mapdbname,
        $host,
        $driver,
        $user,
        $password,
        $debug,
        $perl,
        $perlonlysequences,
        $external,
        $port,
        $contig_overlap_source,
        $overlap_distance_cutoff,
        ) = $self->_rearrange([qw(
            DBNAME
	    MAPDBNAME
	    HOST
	    DRIVER
	    USER
	    PASS
	    DEBUG
	    PERLONLYFEATURES
	    PERLONLYSEQUENCES
	    EXTERNAL
	    PORT
            CONTIG_OVERLAP_SOURCE
            OVERLAP_DISTANCE_CUTOFF
	    )],@args);
    $db   || $self->throw("Database object must have a database name");
    $user || $self->throw("Database object must have a user");

    #
    # This needs to be rethought. We are caching sequences
    # here to allow multiple exons to be retrieved fine
    # And now more cache's. I think cache's might be a fact of life...
    # 

    $self->{'_contig_seq_cache'} = {};
    $self->{'_contig_seq_cnt'} = 0;
    $self->{'_lock_table_hash'} = {};
    $self->_analysis_cache({});
    $self->{'_external_ff'} = [];
  $self->static_golden_path_type('UCSC');

    if( $debug ) {
        $self->_debug($debug);
    } else {
        $self->_debug(0);
    }
    if( ! $driver ) {
        $driver = 'mysql';
    }
    if( ! $host ) {
        $host = 'localhost';
    }
    if ( ! $port ) {
        $port = 3306;
    }

    if( ! defined $perlonlysequences ) {
        $perlonlysequences = 0;
    }

    my $dsn = "DBI:$driver:database=$db;host=$host;port=$port";
	
    if( $debug && $debug > 10 ) {
        $self->_db_handle("dummy dbh handle in debug mode $debug");
    } else {

        my $dbh = DBI->connect("$dsn","$user",$password, {RaiseError => 1});

        $dbh || $self->throw("Could not connect to database $db user $user using [$dsn] as a locator");

        if( $self->_debug > 3 ) {
	    $self->warn("Using connection $dbh");
        }

        $self->_db_handle($dbh);
    }
    $self->username( $user );
    $self->host( $host );
    $self->dbname( $db );
    # following was added on branch; unclear if it is needed:
    $self->mapdbname( $mapdbname );

    if ($perl && $perl == 1) {
        $Bio::EnsEMBL::FeatureFactory::USE_PERL_ONLY = 1;
    }

    $self->perl_only_sequences($perlonlysequences);

    if( defined $external ){
        foreach my $external_f ( @{$external} ) {
	    $self->add_ExternalFeatureFactory($external_f);
        }
    }

    # Store info for connecting to a mapdb.
    {
      $mapdbname ||= 'maps';
      $self->{'_mapdb'} = {
          -DBNAME => $mapdbname,
          -HOST   => $host,
          -PORT   => $port,
          -DRIVER => $driver,
          -USER   => $user,
          -PASS   => $password,
          -ENSDB  => $db,
          };
    }
	#$self->store_mapdb_handle();
	
  # What source of contigoverlaps should we use?
  $self->contig_overlap_source($contig_overlap_source) if $contig_overlap_source;
  
    # What is the maximum distance allowed between contigs
    # in an overlap?
    $self->overlap_distance_cutoff($overlap_distance_cutoff) if $overlap_distance_cutoff;

    return $self; # success - we hope!
}
# only the get part of the 3 functions should be considered public

sub dbname {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_dbname} = $arg );
  $self->{_dbname};
}

 # make the cunningly stored map DB handle available to everybody
 # who wants to play with it...
 
sub store_mapdb_handle {
	my ($self, $arg ) = @_;
	my $dsn = "DBI:"       . $self->{'_mapdb'}->{'-DRIVER'}
			. ":database=" . $self->{'_mapdb'}->{'-DBNAME'}
			. ";host="     . $self->{'_mapdb'}->{'-HOST'}
			. ";port="     . $self->{'_mapdb'}->{'-PORT'}
		;
	$self->warn("Using mapdb DSN $dsn");
	my $mapdbh = DBI->connect("$dsn",$self->{'_mapdb'}->{'-USER'},$self->{'_mapdb'}->{'-PASS'}, {RaiseError => 1});
	$mapdbh || $self->throw("Could not connect to maps database!");
	$self->_mapdb_handle($mapdbh);
}

sub username {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_username} = $arg );
  $self->{_username};
}

sub host {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_host} = $arg );
  $self->{_host};
}


=head2 get_Feature_Obj

 Title   : get_Feature_Obj
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Feature_Obj {
    my ($self) = @_;
    
    my $update_obj = Bio::EnsEMBL::DBOLD::Feature_Obj->new($self);
 
    return $update_obj;
}

=head2 get_Gene_Obj

 Title   : get_Gene_Obj
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Gene_Obj {
    my ($self) = @_;
    
    my $update_obj = Bio::EnsEMBL::DBOLD::Gene_Obj->new($self);
 
    return $update_obj;
}

=head2 get_Protfeat_Adaptor

 Title   : get_Protfeat_Adaptor
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Protfeat_Adaptor {
    my ($self) = @_;
    
    my $update_obj = Bio::EnsEMBL::DBOLD::Protein_Feature_Adaptor->new($self);
 
    return $update_obj;
}


=head2 get_object_by_wildcard

 Title   : get_object_by_wildcard
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_object_by_wildcard{
   my ($self,$type,$string) = @_;

   print STDERR "Got type: $type and string: $string\n";
   my @ids;
   my $sth = $self->prepare("select id from $type where id like \'$string\'");
   print STDERR "mysql: select id from $type where id like \'$string\'\n";
   my $res = $sth->execute || $self->throw("Could not get any ids!");
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@ids,$rowhash->{'id'});
   }
   
   if ($type eq 'gene') {
       return $self->gene_Obj->get_array_supporting('without',@ids);
   }
   if ($type eq 'transcript') {
       my @trans;
       foreach my $id (@ids) {
	   push @trans, $self->gene_Obj->get_Transcript($id);
       }
       return @trans;
   }
   if ($type eq 'exon') {
       my @exons;
       foreach my $id (@ids) {
	   push @exons, $self->gene_Obj->get_Exon($id);
       }
       return @exons;
   }
   if ($type eq 'clone') {
       my @clones;
       foreach my $id (@ids) {
	   push @clones, $self->get_Cone($id);
       }
       return @clones;
   }
   else {
       $self->throw("Type $type not supported, only gene, transcript, exon and clone\n");
   }
   return;
}

=head2 write_Clone

 Title   : write_Clone
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_Clone {
    my ($self,$clone) = @_;

    my $clone_id = $clone->id;

    $clone || $self->throw("Trying to write a clone without a clone object!\n");
    if( !$clone->isa('Bio::EnsEMBL::DB::CloneI') ) {
	$self->throw("Clone must be a CloneI type, not a $clone");
    }
    
    my @sql;
    
    my $sth = $self->prepare('insert into clone (id, internal_id, version, embl_id, embl_version, htg_phase, created, modified, stored) values(?, ?, ?, ?, ?, ?, FROM_UNIXTIME(?), FROM_UNIXTIME(?), NOW())'); 
    my $rv = $sth->execute(
			   $clone_id,
			   "NULL",
			   $clone->version || "NULL",
			   $clone->embl_id || "NULL",
			   $clone->embl_version || "NULL",
			   $clone->htg_phase,
			   $clone->created,
			   $clone->modified
			   );
        
    $self->throw("Failed to insert clone $clone_id") unless $rv;
    $sth = $self->prepare("select last_insert_id()");
    my $res = $sth->execute;
    my $row = $sth->fetchrow_hashref;
    my $id  = $row->{'last_insert_id()'};
    #print(STDERR "Clone $clone_id - $id\n");
    
    foreach my $contig ( $clone->get_all_Contigs() ) {        
        $self->write_Contig($contig,$id);
    }
    
    foreach my $overlap ($clone->get_all_ContigOverlaps) {    
        $self->write_ContigOverlap($overlap);
    }
   
}

=head2 write_Contig

 Title   : write_Contig
 Usage   : $obj->write_Contig($contig,$clone)
 Function: Writes a contig and its dna into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Contig {
    my($self, $contig, $clone)  = @_;
       
    #Why do we have $clone if contig->cloneid is ok?
     
    $self->throw("$contig is not a Bio::EnsEMBL::DB::ContigI - cannot insert contig for clone $clone")
        unless $contig->isa('Bio::EnsEMBL::DB::ContigI');   
    my $dna = $contig->primary_seq  || $self->throw("No sequence in contig object");
    $dna->id                        || $self->throw("No contig id entered.");
    $clone                          || $self->throw("No clone entered.");
    
#   (defined($contig->species)    && $contig->species   ->isa("Bio::EnsEMBL::Species"))    || $self->throw("No species object defined");
#    (defined($contig->chromosome) && $contig->chromosome->isa("Bio::EnsEMBL::Chromosome")) 
#                                    || $self->throw("No chromosomeobject defined");
                                    
#   my $species_id    = $self->write_Species   ($contig->species);
#   my $chromosome_id = $self->write_Chromosome($contig->chromosome,$species_id);    
    my $contigid      = $contig->id;
    my $date          = $contig->seq_date;
    my $len           = $dna   ->length;
    my $seqstr        = $dna   ->seq;
    my $offset        = $contig->embl_offset();
    my $order         = $contig->embl_order();
    #my $chromosome_id = $contig->chromosome->get_db_id;
    my  $chromosome_id = 25;

    # Insert the sequence into the dna table
    $self->_insertSequence($seqstr, $date);
    
    my @sql;
    
    my $sth = $self->prepare("
        insert into contig(id, internal_id, dna, length, clone, offset, corder, chromosomeId ) 
        values(?, ?, LAST_INSERT_ID(), ?, ?, ?, ?, ?)
        "); 
        
    my $rv = $sth->execute(
        $contigid,
        'null',
        $len,
        $clone,
        $offset,
        $order,
        $chromosome_id    
        );  
          
    $self->throw("Failed to insert contig $contigid") unless $rv;
       
    
    $sth = $self->prepare("select last_insert_id()");
    my $res = $sth->execute;
    my $row = $sth->fetchrow_hashref;
    my $id  = $row->{'last_insert_id()'};

    # this is a nasty hack. We should have a cleaner way to do this.
    my @features = $contig->get_all_SeqFeatures;
    #print(STDERR "Contig $contigid - $id\n");
    $contig->internal_id($id);
    
    # write sequence features. We write all of them together as it
    # is more efficient
    my $feature_obj=Bio::EnsEMBL::DBOLD::Feature_Obj->new($self);
    $feature_obj->write($contig, @features);
    
    return 1;
}

=head2 _insertSequence

 Title   : _insertSequence
 Usage   : $obj->_insertSequence
 Function: Insert the dna sequence and date into the dna table.
 Example :
 Returns : 
 Args    : $sequence, $date


=cut

sub _insertSequence {
    my ($self, $sequence, $date) = @_;
    
    $sequence =~ tr/atgcn/ATGCN/;
    
    my $statement = $self->prepare("
        insert into dna(sequence,created) 
        values(?, FROM_UNIXTIME(?))
        "); 
        
    my $rv = $statement->execute($sequence, $date); 
    
    $self->throw("Failed to insert dna $sequence") unless $rv;    
}


=head2 write_Chromosome

 Title   : write_Chromosome
 Usage   : $obj->write_Chromosome
 Function: writes a chromosome into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Chromosome {
    my ($self,$chromosome,$species_id) = @_;

    $self->throw("No chromosome argument input") unless defined($chromosome);
    $self->throw("No species_id argument input") unless defined($species_id);

    if (!$chromosome->isa("Bio::EnsEMBL::Chromosome")) {
	$self->throw("[$chromosome] is not a Bio::EnsEMBL::Chromosome object");
    }

    my $query = "select chromosome_id " .
	        "from   chromosome " .
		"where  name       = '" . $chromosome->name . "' " .
		"and    species_id = "  . $species_id . 
		"and    id         = "  . $chromosome->id;

    my $sth = $self->prepare($query);
    my $res = $sth->execute;

    if ($sth->rows == 1) {
	my $rowhash       = $sth->fetchrow_hashref;
	my $chromosome_id = $rowhash->{chromosome_id};
	return $chromosome_id;
    } 

    $query =  "insert into chromosome(chromosome_id,name,id,species_id) " . 
	      "            values(null,'" . $chromosome->name . "'," . $chromosome->id . "," . $species_id . ")";
	
    
    $sth = $self->prepare($query);
    $res = $sth->execute;

    $sth = $self->prepare("select last_insert_id()");
    $res = $sth->execute;

    my $rowhash       = $sth->fetchrow_hashref;
    my $chromosome_id = $rowhash->{'last_insert_id()'};
   
    return $chromosome_id;
}

=head2 mapdb

    $obj->mapdb($mapdb);
    my $mapdb = $obj->mapdb;

Sets or gets a mapdb connection, which is a
C<Bio::EnsEMBL::Map::DBOLD::Obj> object.

If a mapdb connection doesn't exist, a new
connection is made using the information provided
to the C<_initialize> method.  This will produce
an exception if the
C<Bio::EnsEMBL::Map::DBOLD::Obj> module can't be
found.

=cut

sub mapdb {
    my( $self, $value ) = @_;
    
    if ($value) {
        $self->throw("$value is not a valid mapdb object")
            unless $value->isa("Bio::EnsEMBL::Map::DBOLD::Obj");
        $self->{'_mapdb'} = $value;
    }
    else {
        my $map = $self->{'_mapdb'}
            or $self->throw("No mapdb information");

        # If $map is just an unblessed hash (first time
        # mapdb is called), connect to the map database.
        if (ref($map) eq 'HASH') {
            require Bio::EnsEMBL::Map::DBOLD::Obj;
            $self->{'_mapdb'} = Bio::EnsEMBL::Map::DBOLD::Obj->new(%$map);
        }
    }
    return $self->{'_mapdb'};
}

# was added on branch; not clear if needed:
sub mapdbname {
  my ($self, $arg ) = @_;
  ( defined $arg ) &&
    ( $self->{_mapdbname} = $arg );
  $self->{_mapdbname};
}



=head2 write_Species

 Title   : write_Species
 Usage   : $obj->write_Species
 Function: writes a species object into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Species {
    my ($self,$species) = @_;

    if (!defined($species)) {
	$self->throw("No species argument input");
    }
    if (!$species->isa("Bio::EnsEMBL::Species")) {
	$self->throw("[$species] is not a Bio::EnsEMBL::Species object");
    }

    my $query = "select species_id " .
	        "from   species " .
		"where  nickname    = '" . $species->nickname . "' " . 
		"and    taxonomy_id = "  . $species->taxonomy_id;

    my $sth = $self->prepare($query);
    my $res = $sth->execute;

    if ($sth->rows == 1) {
	my $rowhash    = $sth->fetchrow_hashref;
	my $species_id = $rowhash->{species_id};
	return $species_id;
    } 

    $query =  "insert into species(species_id,nickname,taxonomy_id) " . 
	      "            values(null,'" . $species->nickname . "'," . $species->taxonomy_id . ")";
	
    
    $sth = $self->prepare($query);
    $res = $sth->execute;

    $sth = $self->prepare("select last_insert_id()");
    $res = $sth->execute;

    my $rowhash = $sth->fetchrow_hashref;
    my $species_id = $rowhash->{'last_insert_id()'};
   
    return $species_id;
}

=head2 write_ContigOverlap

 Title   : write_ContigOverlap
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_ContigOverlap {
    my ($self, $overlap) = @_;

    unless (defined($overlap)) {
	$self->throw("No overlap object");
    }

    unless (($overlap->isa("Bio::EnsEMBL::ContigOverlap"))) {
	$self->throw("[$overlap] is not a Bio::EnsEMBL::ContigOverlap");
    }

    my $contiga           = $overlap->contiga;
    my $contigb           = $overlap->contigb;
    my $contig_a_position = $overlap->positiona;
    my $contig_b_position = $overlap->positionb;
    my $overlap_type      = $overlap->overlap_type;

    # Firstly check that both contigs involved in the overlap are present in the db.
    # If they are new they may not have been inserted yet. In this case when they are inserted, 
    # this ContigOverlap should be found again and will be correctly inserted, hopefully!
    
    my $r_contig_a;
    my $r_contig_b;
    
    eval {
       $r_contig_a = $self->get_Contig($contiga->id);
       $r_contig_b = $self->get_Contig($contigb->id);
    };
    if( $@ ) {
       $self->warn("Cannot write Overlaps - contig is not in database - $@");
       return;
    }
  
    #print(STDERR "contiga "         . $contiga->id . "\t" . $contiga->internal_id . "\n");
    #print(STDERR "contigb "         . $contigb->id . "\t" . $contigb->internal_id . "\n");
    #print(STDERR "contigaposition " . $contig_a_position . "\n");
    #print(STDERR "contigbposition " . $contig_b_position . "\n");
    #print(STDERR "overlap type "    . $overlap_type . "\n");

    # Get dna_id's
    my( $dna_a_id, $dna_b_id ) = map $_->dna_id, ($r_contig_a, $r_contig_b);
    $self->throw("Can't get both dna_a_id '$dna_a_id' and dna_b_id '$dna_b_id'")
        unless $dna_a_id and $dna_b_id;
    $self->throw("dna_id's are the same: dna_a_id '$dna_a_id' and dna_b_id '$dna_b_id'")
        if $dna_a_id == $dna_b_id;

    my $source   = $overlap->source;
    my $distance = $overlap->distance;

    #print(STDERR "DNA ids are $dna_a_id : $dna_b_id\n");

    # Insert data into table
    my $insert = $self->prepare(q{
        INSERT contigoverlap(
            dna_a_id, dna_b_id
          , contig_a_position, contig_b_position
          , overlap_source
          , overlap_size, overlap_type)
        VALUES(?,?,?,?,?,?,?)
        });
    eval{
        $insert->execute(
            $dna_a_id, $dna_b_id
          , $contig_a_position, $contig_b_position
          , $source
          , $distance, $overlap_type);
    };
    
    # Check that insert succeeded
    if ($@) {
        if ($@ =~ /Duplicate/) {
            $self->warn("Insert failed: $@");
        } else {
            $self->throw($@);
        }
    } else {
        return 1;
    }
}

=head2 prepare

 Title   : prepare
 Usage   : $sth = $dbobj->prepare("select seq_start,seq_end from feature where analysis = \" \" ");
 Function: prepares a SQL statement on the DBI handle

           If the debug level is greater than 10, provides information into the
           DummyStatement object
 Example :
 Returns : A DBI statement handle object
 Args    : a SQL string


=cut

sub prepare {
   my ($self,$string) = @_;

   if( ! $string ) {
       $self->throw("Attempting to prepare an empty SQL query!");
   }
   if( !defined $self->_db_handle ) {
      $self->throw("Database object has lost its database handle! getting otta here!");
   }
      
   if ($self->diffdump) {
       my $fh=$self->diff_fh;
       open (FILE,">>$fh");
       if ($string =~ /insert|delete|replace/i) {
	   print FILE "$string\n";
       }
       
   }
   
   if( $self->_debug > 10 ) {
       print STDERR "Prepared statement $string\n";
       my $st = Bio::EnsEMBL::DBOLD::DummyStatement->new();
       $st->_fileh(\*STDERR);
       $st->_statement($string);
       return $st;
   }

   # should we try to verify the string?

   return $self->_db_handle->prepare($string);
}


=head2 add_ExternalFeatureFactory

 Title   : add_ExternalFeatureFactory
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_ExternalFeatureFactory{
   my ($self,$value) = @_;

   unless( ref $value && $value->isa('Bio::EnsEMBL::DB::ExternalFeatureFactoryI') ) {
       $self->throw("[$value] is not a Bio::EnsEMBL::DB::ExternalFeatureFactoryI but it should be!");
   }

   push(@{$self->{'_external_ff'}},$value);
}

=head2 _each_ExternalFeatureFactory

 Title   : _each_ExternalFeatureFactory
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _each_ExternalFeatureFactory{
   my ($self) = @_;

   return @{$self->{'_external_ff'}}
}


=head2 _analysis_cache

 Title   : _analysis_cache
 Usage   : $obj->_analysis_cache()
 Function: 
 Returns : reference to a hash
 Args    : newvalue (optional)


=cut

sub _analysis_cache{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_analysis_cache'} = $value;
    }
    return $obj->{'_analysis_cache'};

}

=head2 _contig_seq_cache

 Title   : _contig_seq_cache
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _contig_seq_cache{
   my ($self,$id,$seq) = @_;

   if( $seq ) {
       
       #
       # Every 100 hits, flush the cache
       #
       if( $self->{'_contig_seq_cnt'} > 100 ) {
	   $self->_flush_seq_cache;
	   $self->{'_contig_seq_cnt'} = 0;
       }

       $self->{'_contig_seq_cnt'}++;
       $self->{'_contig_seq_cache'}->{$id} = $seq;
   }

   return $self->{'_contig_seq_cache'}->{$id};
}

=head2 _flush_seq_cache

 Title   : _flush_seq_cache
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _flush_seq_cache{
   my ($self,@args) = @_;

   $self->{'_contig_seq_cache'} = {};

}

=head2 _debug

 Title   : _debug
 Usage   : $obj->_debug($newval)
 Function: 
 Example : 
 Returns : value of _debug
 Args    : newvalue (optional)


=cut

sub _debug{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_debug'} = $value;
    }
    return $self->{'_debug'};
    
}


=head2 _db_handle

 Title   : _db_handle
 Usage   : $obj->_db_handle($newval)
 Function: 
 Example : 
 Returns : value of _db_handle
 Args    : newvalue (optional)


=cut

sub _db_handle{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'_db_handle'} = $value;
    }
    return $self->{'_db_handle'};

}

=head2 _lock_tables

 Title   : _lock_tables
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _lock_tables{
   my ($self,@tables) = @_;
   
   my $state;
   foreach my $table ( @tables ) {
       if( $self->{'_lock_table_hash'}->{$table} == 1 ) {
	   $self->warn("$table already locked. Relock request ignored");
       } else {
	   if( $state ) { $state .= ","; } 
	   $state .= "$table write";
	   $self->{'_lock_table_hash'}->{$table} = 1;
       }
   }

   my $sth = $self->prepare("lock tables $state");
   my $rv = $sth->execute();
   $self->throw("Failed to lock tables $state") unless $rv;

}

=head2 _unlock_tables

 Title   : _unlock_tables
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub _unlock_tables{
   my ($self,@tables) = @_;

   my $sth = $self->prepare("unlock tables");
   my $rv  = $sth->execute();
   $self->throw("Failed to unlock tables") unless $rv;
   %{$self->{'_lock_table_hash'}} = ();
}


=head2 DESTROY

 Title   : DESTROY
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub DESTROY {
   my ($obj) = @_;

   #$obj->_unlock_tables();

   if( $obj->{'_db_handle'} ) {
       $obj->{'_db_handle'}->disconnect;
       $obj->{'_db_handle'} = undef;
   }
}

##################DEPRECATED METHODS######################
#                                                        #
#All the methods below are deprecated methods,           #
#only kept here to allow old scripts to work             #
#They all send a warning and call the new method instead #
#                                                        #
##########################################################

=head2 get_Gene

 Title   : get_Gene
 Usage   : $obj->get_Gene($geneid, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Example : $obj->get_Gene('ENSG00000009151','evidence')
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : gene id and supporting tag (if latter not specified, assumes without
	   Note that it is much faster to get genes without supp.evidence!


=cut

sub get_Gene {
   my ($self,$geneid, $supporting) = @_;

   $self->warn("Obj->get_Gene is a deprecated method!\nCalling gene_Obj->get instead!");

   return $self->gene_Obj->get($geneid,$supporting);
}

=head2 get_Gene_by_Transcript_id

 Title   : get_Gene_by_Transcript_id
 Usage   : $obj->get_Gene_by_Transcript_id($transid, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Example : $obj->get_Gene_by_Transcript_id('ENST00000009151')
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified, assumes without
	   Note that it is much faster to get genes without supp.evidence!


=cut

sub get_Gene_by_Transcript_id {
   my ($self,$tid, $supporting) = @_;

   $self->warn("Obj->get_Gene_by_Transcript_id is a deprecated method! 
Calling gene_Obj->get instead!");

   return $self->gene_Obj->get_Gene_by_Transcript_id($tid,$supporting);
}



=head2 get_Gene_by_DBLink

 Title   : get_Gene_by_DBLink
 Usage   : $obj->get_Gene_by_DBLink($ext_id, $supporting)
 Function: gets one gene out of the db with or without supporting evidence
 Example : $obj->get_Gene_by_DBLink( 'MC1R')
 Returns : gene object (with transcripts, exons and supp.evidence if wanted)
 Args    : transcript id and supporting tag (if latter not specified, assumes without
	   Note that it is much faster to get genes without supp.evidence!


=cut

sub get_Gene_by_DBLink {
   my ($self,$ext_id, $supporting) = @_;

  # $self->warn("Obj->get_Gene_by_DBLink is a deprecated method! Calling gene_Obj->get instead!");

   return $self->gene_Obj->get_Gene_by_DBLink($ext_id,$supporting);
}








=head2 get_Gene_array

 Title   : get_Gene_array
 Usage   :
 Function: old deprecated method, points to new method
           get_gene_array_supporting without asking for supp.evidence
 Example :
 Returns : 
 Args    :


=cut

sub get_Gene_array {
    my ($self,@geneid) = @_;

    $self->throw("Very deprecated method, should call methods with supporting evidence and from
gene_Obj!");
}

=head2 get_Gene_array_supporting

 Title   : get_Gene_array_supporting
 Usage   : $obj->get_Gene_array_supporting($supporting,@geneid)
 Function: Gets an array of genes, with transcripts and exons. If $supporting
           equal to 'evidence' the supporting evidence for each exon is also read
           from the supporting evidence table
 Example : $obj->get_Gene_array_supporting ('evidence',@geneid)
 Returns : an array of gene objects
 Args    : 'evidence' and gene id array

=cut

sub get_Gene_array_supporting {
    my ($self,$supporting,@geneid) = @_;

    $self->warn("Obj->get_Gene_array_supporting is a deprecated method!
Calling gene_Obj->get_array_supporting instead!");

    return $self->gene_Obj->get_array_supporting($supporting,@geneid);
}



=head2 get_Virtual_Contig
    
 Title   : get_Virtual_Contig
 Usage   : $obj->get_Virtual_Contig($transcript,$max_length)
 Function: Gets a Bio::EnsEMBL::DB::Virtual Contig object which 
           spans the whole sequence on which the given 
           Bio::EnsEMBL::Transcript object lies, as long 
           as its length does not exceed max_length. If max_length
           is exceeded, undef is returned instead.
 Example : $obj->get_Virtual_Contig($transcript,50000)
 Returns : VirtualContig Object (or undef)
 Args    : Bio::EnsEMBL::Transcript object and max_length int variable


=cut



sub get_Virtual_Contig_by_Transcript_id {
   my ($self,$tid, $maxlen) = @_;

#   $self->warn("Obj->get_Virtual_contig is a deprecated method! 
#Calling gene_Obj->get_Virtual_contig instead!");

   return my $vc =$self->gene_Obj->get_Virtual_Contig($tid,$maxlen);
}





=head2 get_Transcript_in_VC_coordinates
    
 Title   : get_Transcript_in_VC_coordinates
 Usage   : $obj->get_Transcript_in_VC_coordinates($transcript_id)
 Function: Gets a Bio::EnsEMBL::Transcript object in vc coordinates
 Example : $obj->get_Virtual_Contig($transcript_id)
 Returns : Bio::EnsEMBL::Transcript
 Args    : transcript id


=cut




sub get_Transcript_in_VC_coordinates{
   my ($self,$tid) = @_;

   $self->warn("Obj->get_Transcript_in_VC_coordinates is a deprecated method! 
Calling gene_Obj->get_Virtual_contig instead!");

   return my $transcript =$self->gene_Obj->get_Transcript_in_VC_coordinates($tid);
}
    
=head2 get_Protein_annseq

 Title   : get_Protein_annseq
 Usage   : get_Protein_annseq ($ENSP); 
 Function: Creates an annseq object for a particular peptide, storing the peptide
           sequence in $annseq->primary_seq, and adding all the protein features as generic
           Seqfeatures
 Example : 
 Returns : $annseq
 Args    : $ENSP


=cut

sub get_Protein_annseq{
    my ($self,$ENSP) = @_;

    $self->warn("Obj->get_Protein_annseq is a deprecated method! 
Calling Feature_Obj->get_Protein_annseq instead!");
    
    my $feature_obj=Bio::EnsEMBL::DBOLD::Feature_Obj->new($self);
    return $feature_obj->get_Protein_annseq($ENSP);
} 

=head2 get_Transcript
    
 Title   : get_Transcript
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut
    
sub get_Transcript{
    my ($self,$transid) = @_;
 
    $self->warn("Obj->get_Transcript is a deprecated method! 
Calling gene_Obj->get_Transcript instead!");

    return $self->gene_Obj->get_Transcript($transid);
}

=head2 get_Translation

 Title   : get_Translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Translation{
   my ($self,$translation_id) = @_;

   $self->warn("Obj->get_Translation is a deprecated method! 
Calling gene_Obj->get_Translation instead!");

   return $self->gene_Obj->get_Translation($translation_id);
}

=head2 get_Exon

 Title   : get_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Exon{
   my ($self,$exonid) = @_;

   $self->warn("Obj->get_Exon is a deprecated method! 
Calling gene_Obj->get_Exon instead!");

   return $self->gene_Obj->get_Exon($exonid);
}

=head2 get_all_Gene_id

 Title   : get_all_Gene_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Gene_id{
   my ($self) = @_;

   $self->warn("Obj->get_all_Gene_id is a deprecated method! 
Calling gene_Obj->get_all_Gene_id instead!");

   return $self->gene_Obj->get_all_Gene_id();
}



=head2 get_all_Transcript_id

 Title   : get_all_Transcript_id
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Transcript_id{
   my ($self) = @_;

   $self->warn("Obj->get_all_Transcript_id is a deprecated method! 
Calling gene_Obj->get_all_Gene_id instead!");

   return $self->gene_Obj->get_all_Transcript_id();
}

=head2 delete_Exon

 Title   : delete_Exon
 Usage   : $obj->delete_Exon($exon_id)
 Function: Deletes exon, including exon_transcript rows
 Example : $obj->delete_Exon(ENSE000034)
 Returns : nothing
 Args    : $exon_id


=cut

sub delete_Exon{
    my ($self,$exon_id) = @_;

    $self->warn("Obj->delete_Exon is a deprecated method
Calling gene_Obj->delete_Exon instead!");

    return $self->gene_Obj->delete_Exon($exon_id);
}

=head2 delete_Supporting_Evidence

 Title   : delete_Supporting_Evidence
 Usage   : $obj->delete_Supporting_Evidence($exon_id)
 Function: Deletes exon\'s supporting evidence entries
 Example : $obj->delete_Supporting_Evidence(ENSE000034)
 Returns : nothing
 Args    : $exon_id


=cut

sub delete_Supporting_Evidence {
    my ($self,$exon_id) = @_;
 
    $self->warn("Obj->delete_Supporting_Evidence is a deprecated method
Calling gene_Obj->delete_Supporting_Evidence instead!");

    return $self->gene_Obj->delete_Supporting_Evidence($exon_id);
}

=head2 delete_Features

 Title   : delete_Features
 Usage   :
 Function: deletes all features from a contig;
 Example :
 Returns : 
 Args    :


=cut

sub delete_Features {
    my ($self,$contig) = @_;

    $self->warn("Obj->delete_Features is a deprecated method! 
Calling Feature_Obj->delete instead!");

   my $feature_obj=Bio::EnsEMBL::DBOLD::Feature_Obj->new($self);
   return $feature_obj->delete($contig);
} 

=head2 delete_Gene

 Title   : delete_Gene
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub delete_Gene{
   my ($self,$geneid) = @_;

   $self->warn("Obj->delete_Gene is a deprecated method! 
Calling gene_Obj->delete instead!");

   return $self->gene_Obj->delete($geneid);
}

=head2 geneid_to_cloneid

 Title   : geneid_to_cloneid
 Usage   : @cloneid = $db->geneid_to_cloneid($geneid);
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub geneid_to_cloneid{
    my ($self,$geneid) = @_;
    
    $self->throw("Obj->geneid_to_cloneid is a deprecated method, called gene_Obj->each_cloneid instead!
All the gene, transcript, and exon methods are now to be found in gene_Obj");
    return $self->gene_Obj->each_cloneid($geneid);
}

=head2 write_Gene

 Title   : write_Gene
 Usage   : $obj->write_Gene($gene)
 Function: writes a particular gene into the database
           
 Example :
 Returns : 
 Args    :


=cut


sub write_Gene{
   my ($self,$gene) = @_;

  # $self->warn("Obj->write_Gene is a deprecated method! Calling gene_Obj->write instead!");

   return $self->gene_Obj->write($gene);
}

=head2 write_all_Protein_features

 Title   : write_all_Protein_features
 Usage   : $obj->write_all_Protein_features($ENSP)
 Function: writes all protein features of a particular peptide into the database          
 Example :
 Returns : 
 Args    :


=cut

sub write_all_Protein_features {
    my ($self,$prot_annseq,$ENSP) = @_;

    $self->warn("Obj->write_all_Protein_features is a deprecated method! 
Calling Feature_Obj->write_all_Protein_features instead!");
    
    my $feature_obj=Bio::EnsEMBL::DBOLD::Feature_Obj->new($self);
    return $feature_obj->write_all_Protein_features($prot_annseq,$ENSP);
} 


=head2 write_Protein_feature

 Title   : write_Protein_feature
 Usage   : $obj->write_Protein_feature($ENSP, $feature)
 Function: writes a protein feature object of a particular peptide into the database          
 Example :
 Returns : 
 Args    :


=cut

sub write_Protein_feature {
    my ($self,$ENSP,$feature) = @_;
 
    $self->warn("Obj->write_Protein_feature is a deprecated method! 
Calling Feature_Obj->write_Protein_feature instead!");
    
    my $feature_obj=Bio::EnsEMBL::DBOLD::Feature_Obj->new($self);
    return $feature_obj->write_Protein_feature($ENSP,$feature);
} 

=head2 write_Feature

 Title   : write_Feature
 Usage   : $obj->write_Feature($contig,@features)
 Function: Writes a feature on the genomic sequence of a contig into the database
 Example :
 Returns : nothing
 Args    : Bio::EnsEMBL::SeqFeatureI


=cut

sub write_Feature {
    my ($self,$contig,@features) = @_;

    $self->warn("Obj->write_Feature is a deprecated method! 
Calling Feature_Obj->write_Feature instead!");
    
    my $feature_obj=Bio::EnsEMBL::DBOLD::Feature_Obj->new($self);
    return $feature_obj->write($contig,@features);
} 

=head2 write_supporting_evidence

 Title   : write_supporting_evidence
 Usage   : $obj->write_supporting_evidence
 Function: Writes supporting evidence features to the database
 Example :
 Returns : nothing
 Args    : None


=cut

sub write_supporting_evidence {
    my ($self,$exon) = @_;

    $self->warn("Obj->write_supporting_evidence is a deprecated method!
Calling gene_Obj->write_supporting_evidence instead!");

    return $self->gene_Obj->write_supporting_evidence($exon);
}

=head2 get_supporting_evidence

 Title   : get_supporting_evidence
 Usage   : $obj->get_supporting_evidence
 Function: Writes supporting evidence features to the database
 Example :
 Returns : nothing
 Args    : array of exon objects, needed to know which exon to attach the evidence to


=cut

sub get_supporting_evidence {
    my ($self,@exons) = @_;

    $self->warn("Obj->get_supporting_evidence is a deprecated method! 
Calling gene_Obj->get_supporting_evidence instead!");

   return $self->gene_Obj->get_supporting_evidence(@exons);
}

=head2 write_Analysis

 Title   : write_Analysis
 Usage   : $obj->write_Analysis($anal)
 Function: Writes analysis details to the database
           Checks first whether this analysis entry already exists
 Example :
 Returns : int
 Args    : Bio::EnsEMBL::AnalysisI


=cut

sub write_Analysis {
    my ($self,$analysis) = @_;

    $self->warn("Obj->write_Analysis is a deprecated method! 
Calling Feature_Obj->write_Analysis instead!");
    
    my $feature_obj=Bio::EnsEMBL::DBOLD::Feature_Obj->new($self);
    return $feature_obj->write_Analysis($analysis);
} 
    
=head2 exists_Homol_Feature

 Title   : exists_Homol_Feature
 Usage   : $obj->exists_Homol_Feature($feature)
 Function: Tests whether this feature already exists in the database
 Example :
 Returns : nothing
 Args    : Bio::SeqFeature::Homol


=cut

sub exists_Homol_Feature {
    my ($self,$feature,$analysisid,$contig) = @_;

    $self->warn("Obj->exists_Homol_Feature is a deprecated method! 
Calling Feature_Obj->exists_Homol_Feature instead!");
    
    my $feature_obj=Bio::EnsEMBL::DBOLD::Feature_Obj->new($self);
    return $feature_obj->exists($feature,$analysisid,$contig);
} 
    
=head2 get_Analysis

 Title   : get_Analysis
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Analysis {
    my ($self,$id) = @_;

    $self->warn("Obj->get_Analysis is a deprecated method! 
Calling Feature_Obj->get_Analysis instead!");
    
    my $feature_obj=Bio::EnsEMBL::DBOLD::Feature_Obj->new($self);
    return $feature_obj->get_Analysis($id);
} 

=head2 exists_Analysis

 Title   : exists_Analysis
 Usage   : $obj->exists_Analysis($anal)
 Function: Tests whether this feature already exists in the database
 Example :
 Returns : Analysis id if the entry exists
 Args    : Bio::EnsEMBL::Analysis


=cut

sub exists_Analysis {
    my ($self,$analysis) = @_;
    
    $self->warn("Obj->exists_Analysis is a deprecated method! 
Calling Feature_Obj->exists_Analysis instead!");
    
    my $feature_obj=Bio::EnsEMBL::DBOLD::Feature_Obj->new($self);
    return $feature_obj->exists_Analysis($analysis);
} 
 
    
=head2 write_Transcript

 Title   : write_Transcript
 Usage   : $obj->write_Transcript($trans,$gene)
 Function: writes a particular transcript *but not the exons* into
           the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Transcript{
   my ($self,$trans,$gene) = @_;

   $self->warn("Obj->write_Transcript is a deprecated method! 
Calling gene_Obj->write_Transcript instead!");

   return $self->gene_Obj->write_Transcript($trans,$gene);
}

=head2 write_Translation

 Title   : write_Translation
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub write_Translation{
    my ($self,$translation) = @_;

    $self->warn("Obj->write_Translation is a deprecated method
Calling gene_Obj->write_Translation instead!");

    return $self->gene_Obj->write_Translation($translation);
}


=head2 write_Exon

 Title   : write_Exon
 Usage   : $obj->write_Exon($exon)
 Function: writes a particular exon into the database
 Example :
 Returns : 
 Args    :


=cut

sub write_Exon {
   my ($self,$exon) = @_;

   $self->warn("Obj->write_Exon is a deprecated method! 
Calling gene_Obj->write_Exon instead!");

   return $self->gene_Obj->write_Exon($exon);
}

=head2 get_Clone

 Title   : get_Clone
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Clone { 
    my ($self,$accession) = @_;

    my $ca= Bio::EnsEMBL::DBOLD::CloneAdaptor->new($self);

    return $ca->fetch_by_accession($accession);
}
  
=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Contig{
   my ($self,$id) = @_;

   #$self->warn("Obj->get_Contig is a deprecated method! 
#Calling Contig->fetch instead!");

   my $contig      = new Bio::EnsEMBL::DBOLD::RawContig ( -dbobj => $self,
							  -id    => $id,
							  -perlonlysequences => $self->perl_only_sequences(),
							  -userawcontigacc => !$self->perl_only_contigs );
   
   return $contig->fetch();
}

sub get_Contig_by_international_id{
   my ($self,$int_id) = @_;
   #$self->warn("Obj->get_Contig is a deprecated method! 
#Calling Contig->fetch instead!");
   my $sth=$self->prepare("select id from contig where international_id = '$int_id'");
   $sth->execute;
   my $row = $sth->fetchrow_hashref;
   my $id  = $row->{'id'};

   return $self->get_Contig($id);
}

=head2 get_Contigs_by_Chromosome

 Title   : get_Contig_by_Chromosome
 Usage   : @contigs = $dbobj->get_Contig_by_Chromosome( $chrObj );
 Function: retrieve contigs belonging to a certain chromosome from the
           database 
 Example :
 Returns : A list of Contig objects. Probably an empty list.
 Args    :


=cut

sub get_Contigs_by_Chromosome {
    my ($self,$chromosome ) = @_;
    
    $self->throw("Obj->get_Contigs_by_Chromosome is a deprecated method! 
Calling Contig->get_by_Chromosome instead!");
}

=head2 get_all_Clone_id

 Title   : get_all_Clone_id
 Usage   : @cloneid = $obj->get_all_Clone_id
 Function: returns all the valid (live) Clone ids in the database
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Clone_id{
   my ($self) = @_;
   my @out;

   my $sth = $self->prepare("select id from clone");
   my $res = $sth->execute;

   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'id'});
   }

   return @out;
}

=head2 get_all_fpcctg_ids

 Title   : get_all_fpcctg_ids
 Usage   : @cloneid = $obj->get_all_fpcctg_ids
 Function: returns all the valid FPC contigs from given golden path
 Example :
 Returns : 
 Args    : static golden path type (typically, 'UCSC')


=cut

sub get_all_fpcctg_ids {
   my ($self, $type) = @_;

   $self->throw("no static_gold_path given") unless defined $type;
   my @out;

   my $q= "SELECT DISTINCT fpcctg_name 
           FROM static_golden_path sgp
           WHERE type = '$type'";
   my $sth = $self->prepare($q) || $self->throw("can't prepare: $q");
   my $res = $sth->execute || $self->throw("can't prepare: $q");

   while( my ($id) = $sth->fetchrow_array) {
       push(@out, $id);
   }
   return @out;
}



=head2 perl_only_sequences

 Title   : perl_only_sequences
 Usage   : $obj->perl_only_sequences($newval)
 Function: 
 Returns : value of perl_only_sequences
 Args    : newvalue (optional)


=cut

sub perl_only_sequences{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'perl_only_sequences'} = $value;
    }
    return $obj->{'perl_only_sequences'};

}

=head2 perl_only_contigs

 Title   : perl_only_contigs
 Usage   : $obj->perl_only_contigs($newval)
 Function: 
 Returns : value of perl_only_contigs
 Args    : newvalue (optional)


=cut

sub perl_only_contigs{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'perl_only_contigs'} = $value;
    }
    return $obj->{'perl_only_contigs'};

}

=head2 delete_Clone

 Title   : delete_Clone
 Usage   : $obj->delete_Clone($clone_id)
 Function: Deletes clone, including contigs, but not its genes
 Example :
 Returns : 
 Args    :


=cut

sub delete_Clone{
   my ($self,$clone_id) = @_;

   #$self->warn("Obj->delete_Clone is a deprecated method! 
#Calling Clone->delete instead!");
   
   (ref($clone_id)) && $self->throw ("Passing an object reference instead of a variable\n");

   my $clone = new Bio::EnsEMBL::DBOLD::Clone( -id    => $clone_id,
					       -dbobj => $self );
   
   return $clone->delete();
}

=head2 cloneid_to_geneid

 Title   : cloneid_to_geneid
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub cloneid_to_geneid{
   my ($self,$cloneid) = @_;

    $self->warn("Obj->cloneid_to_geneid is a deprecated method! 
Calling Clone->get_all_geneid instead!");

   (ref($cloneid)) && $self->throw ("Passing an object reference instead of a variable!\n");

   my $clone = new Bio::EnsEMBL::DBOLD::Clone( -id    => $cloneid,
					       -dbobj => $self );
   
   return $clone->get_all_my_geneid();
}


=head2 gene_Obj
    
 Title   : gene_Obj
 Usage   : my $geneobj = $db->gene_Obj
 Function: Returns the gene object database handle
 Example : 
 Returns : Bio::EnsEMBL::DB::Gene_ObjI
 Args    : 

=cut

sub gene_Obj {
    my ($self) = @_;

    unless (defined($self->{_gene_obj})) {
	$self->{_gene_obj} = Bio::EnsEMBL::DBOLD::Gene_Obj->new($self);    
    }

    return $self->{_gene_obj};
}

=head2 feature_Obj
    
 Title   : feature_Obj
 Usage   : my $featureobj = $db->feature_Obj
 Function: Returns the feature object database handle
 Example : 
 Returns : Bio::EnsEMBL::DB::Feature_ObjI
 Args    : 

=cut

sub feature_Obj {
    my ($self) = @_;

    unless (defined($self->{_feature_obj})) {
	$self->{_feature_obj} = Bio::EnsEMBL::DBOLD::Feature_Obj->new($self);    
    }

    return $self->{_feature_obj};

}

=head2 get_AnalysisAdaptor

 Title   : get_AnalysisAdaptor
 Usage   : $analysisAdaptor = $dbObj->get_AnalysisAdaptor;
 Function: gives the adaptor to fetch/store Analysis objects.
 Example :
 Returns : the adaptor
 Args    :


=cut

sub get_AnalysisAdaptor {
   my ($self) = @_;
   if( ! defined $self->{_analysisAdaptor} ) {
     $self->{_analysisAdaptor} = 
       Bio::EnsEMBL::DBOLD::AnalysisAdaptor->new($self);
   }
   return $self->{_analysisAdaptor};
}



=head2 get_StaticGoldenPathAdaptor

 Title   : get_StaticGoldenPathAdaptor
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_StaticGoldenPathAdaptor{
   my ($self,@args) = @_;

   return Bio::EnsEMBL::DBOLD::StaticGoldenPathAdaptor->new(-dbobj => $self);
}

=head2 get_KaryotypeBandAdaptor

 Title   : get_KaryotypeBandAdaptor
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_KaryotypeBandAdaptor{
   my ($self,@args) = @_;

   return Bio::EnsEMBL::DBOLD::KaryotypeBandAdaptor->new($self);
}


=head2 find_GenomeHits
    
 Title   : find_GenomeHits
 Usage   : my @features = $self->find_GenomeHits($hid)
 Function: Finds all features in the db that
           are hits to a sequence with id $hid
 Example : 
 Returns : @ Bio::EnsEMBL::FeaturePair
 Args    : string

=cut
 
sub find_GenomeHits {
    my ($self,$arg) = @_;

    return $self->feature_Obj->find_GenomeHits($arg);
}
			     

=head2 deleteObj

    Title   : deleteObj
    Usage   : $dbObj->deleteObj
    Function: Call when you are done with this object. Breaks links between objects. Necessary to clean up memory.
    Example : -
    Returns : -
    Args    : -


=cut

sub deleteObj {

  my  $self=shift;
  my $dummy;

  print STDERR "Destroying DB Obj!\n";       
  $self->DESTROY;
  
  foreach my $name ( keys %{$self} ) {
    eval {
      $dummy = $self->{$name}; 
      $self->{$name}  = undef;
      $dummy->deleteObj;
    };
  }
}



=head2 diff_fh

 Title   : diff_fh
 Usage   : $obj->diff_fh($newval)
 Function: path and name of the file to use for writing the mysql diff dump
 Example : 
 Returns : value of diff_fh
 Args    : newvalue (optional)


=cut

sub diff_fh{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_diff_fh'} = $value;
    }
    return $self->{'_diff_fh'};
    
}


=head2 diffdump

 Title   : diffdump
 Usage   : $obj->diffdump($newval)
 Function: If set to 1 sets $self->_prepare to print the diff sql 
           statementents to the filehandle specified by $self->diff_fh
 Example : 
 Returns : value of diffdump
 Args    : newvalue (optional)


=cut

sub diffdump{
    my ($self,$value) = @_;
    if( defined $value) {
	$self->{'_diffdump'} = $value;
    }
    return $self->{'_diffdump'};
    
}

=head2 contig_overlap_source

 Title   : contig_overlap_source
 Usage   : $obj->contig_overlap_source($newval)
 Function: Gets or sets a subroutine which is used to
           decide which overlap sources are used to
           build virtual contigs.
 Returns : value of contig_overlap_source, or a ref to
           the default subroutine "_default_contig_overlap_source".
 Args    : ref to a subroutine


=cut

sub contig_overlap_source {
    my( $self, $sub ) = @_;
    
    if ($sub) {
        $self->throw("'$sub' is not a CODE reference")
            unless ref($sub) eq 'CODE';
        $self->{'_contig_overlap_source'} = $sub;
    }
    return $self->{'_contig_overlap_source'} || \&_default_contig_overlap_source;
}

=pod

=head2 _default_contig_overlap_source

The default subroutine used by B<contig_overlap_source>.

=cut

sub _default_contig_overlap_source {
    #return $_[0] eq 'phrap';
    return 1;
}


=head2 get_PredictionFeature_as_Transcript

 Title   : get_PredictionFeature_as_Transcript
 Usage   :$obj->get_PredictionFeature_as_Transcript($genscan_id)
 Function:Call get_PredictionFeature_as_Transcript in Feature_obj object,see documentation for this method
 Example :
 Returns : 
 Args    :


=cut

sub get_PredictionFeature_as_Transcript{
   my ($self,$genscan_id) = @_;

   return $self->get_Feature_Obj->get_PredictionFeature_as_Transcript($genscan_id);

}








=head2 overlap_distance_cutoff

 Title   : overlap_distance_cutoff
 Usage   : $obj->overlap_distance_cutoff($int)
 Function: Gets or sets an integer which is used when building
           VirtualContigs.  If the distance in a contig overlap
           is greater than the cutoff, then the overlap will
           not be returned.
 Returns : value of overlap_distance_cutoff, or -1 if it is not set
 Args    : positive integer


=cut


sub overlap_distance_cutoff {
    my( $self, $cutoff ) = @_;
    
    if (defined $cutoff) {
        $self->throw("'$cutoff' is not an positive integer")
            unless $cutoff =~ /^\d+$/;
        $self->{'_overlap_distance_cutoff'} = $cutoff;
    }
    my $ret = $self->{'_overlap_distance_cutoff'};
    if (defined($ret)) {
        return $ret;
    } else {
        return -1;
    }
}


=head2 extension_tables

 Title   : extension_tables
 Usage   : $obj->extension_tables($newval)
 Function: 
 Returns : value of extension_tables
 Args    : newvalue (optional)


=cut

sub extension_tables{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'extension_tables'} = $value;
    }
    return $obj->{'extension_tables'};

}

=head2 static_golden_path_type

 Title   : static_golden_path_type
 Usage   : $obj->static_golden_path_type($newval)
 Function: 
 Example : 
 Returns : value of static_golden_path_type
 Args    : newvalue (optional)


=cut

sub static_golden_path_type{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'static_golden_path_type'} = $value;
    }
    return $obj->{'static_golden_path_type'};

}


=head2 _crossdb

 Title   : _crossdb
 Usage   : $obj->_crossdb($newval)
 Function: 
 Returns : value of _crossdb
 Args    : newvalue (optional)


=cut

sub _crossdb{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'_crossdb'} = $value;
    }
    return $obj->{'_crossdb'};

}


1;
