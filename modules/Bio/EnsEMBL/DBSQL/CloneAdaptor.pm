#
# BioPerl module for DB::Clone
#
# Cared for by EnsEMBL (www.ensembl.org)
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::CloneAdaptor

=head1 SYNOPSIS

    $clone_adaptor = $database_adaptor->get_CloneAdaptor();
    $clone = $clone_adaptor->fetch_by_dbID(1234);
    @contig = $clone->get_all_Contigs();
    @genes    = $clone->get_all_Genes();

=head1 DESCRIPTION

Database adaptor for the creation of clone objects

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal 
methods are usually preceded with a _

=cut

# Let the code begin...


package Bio::EnsEMBL::DBSQL::CloneAdaptor;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Clone;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 _generic_sql_fetch

  Arg [1]    : string $where_clause
               the WHERE clause of the SQL query to be executed
  Example    : $clone = $self->_generic_sql_fetch('embl_acc = AC011082');
  Description: PRIVATE Performs an SQL query and returns an Clone object 
               created from the result.
  Returntype : Bio::EnsEMBL::Clone
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::CloneAdaptor

=cut

sub _generic_sql_fetch {
    my( $self, $where_clause ) = @_;
    
    my $sql = q{
        SELECT clone_id
          , name
          , embl_acc
          , version
          , embl_version
          , htg_phase
          , UNIX_TIMESTAMP(created)
          , UNIX_TIMESTAMP(modified)
	  FROM clone }
        . $where_clause .
        q{ ORDER BY embl_version DESC };

    #print STDERR "issue: $sql\n";

    my $sth = $self->prepare($sql);
    $sth->execute;
    
    if (my @fields = $sth->fetchrow) {
        return Bio::EnsEMBL::Clone->new($self, @fields);
    } else {
        return;
    }
}


=head2 fetch_by_accession

  Arg [1]    : string $acc
               the EMBL accession for the clone
  Example    : $clone = $clone_adaptor->fetch_by_accession('AC011082');
  Description: fetches a clone by its EMBL/GenBank accession number
               It will fetch the highest version in the database 
  Returntype : Bio::EnsEMBL::Clone 
  Exceptions : thrown if $acc not defined or if no clone with accession $acc
               exists in the database
  Caller     : general

=cut

sub fetch_by_accession { 
    my ($self,$acc) = @_;
    
    unless ($acc) {
        $self->throw("Accession not given");
    }

    my $clone = $self->_generic_sql_fetch(
        qq{ WHERE embl_acc = '$acc' }
        );

    if ($clone) {
        return $clone;
    } else {
        $self->throw("No clone with accession '$acc'");
    }
}


=head2 fetch_by_accession_version

  Arg [1]    : string $acc
               the EMBL accession number
  Arg [2]    : int $ver
               the EMBL accession version
  Example    : $clone = $clone_adaptor->fetch_by_accession_version($acc, $ver);
  Description: retrieves a Clone object with accession $acc and version $ver 
               from the database.  
  Returntype : Bio::EnsEMBL::Clone
  Exceptions : thrown if $acc or $ver is not defined or if no clone exists with
               in the database with accession $acc and version $ver
  Caller     : general

=cut

sub fetch_by_accession_version { 
  my ($self, $acc, $ver) = @_;
  
  unless ($acc and $ver) {
    $self->throw("Need both accession (got '$acc') and version (got '$ver')");
  }

  my $clone = 
    $self->_generic_sql_fetch( qq{ WHERE embl_acc = '$acc' 
				   AND embl_version = '$ver' });
  if ($clone) {
    return $clone;
  } else {
    $self->throw("no clone with accession '$acc' and version '$ver'");
  }
}


=head2 fetch_by_name

  Arg [1]    : string $name
               the "name" of the clone.
  Example    : $clone = $clone_adaptor->fetch_by_name($name);
  Description: fetches a single clone by the "name" field. The name
               field is often the EMBL ID of a clone (which itself is often
               the same as the accession number) but is a deliberately open
               slot for another identifier for the clone, eg. the stringified
               tracking identifier in a sequencing centre
  Returntype : Bio::EnsEMBL::Clone
  Exceptions : thrown if name is a non-true value (i.e. undefined or '') or 
               if there is no clone with name $name in the database
  Caller     : general

=cut

sub fetch_by_name {
  my ($self, $name) = @_;
  
  $self->throw("name not given") unless $name;
  
  my $clone = $self->_generic_sql_fetch(qq{ WHERE name = '$name' });
  
  if ($clone) {
    return $clone;
  } else {
    $self->throw("No Clone with name '$name'");
  }
}


=head2 fetch_by_dbID

  Arg [1]    : int $id
               the numeric internal id of the clone table 
  Example    : $clone = $clone_adaptor->fetch_by_dbID(7);
  Description: fetches a clone by the internal id in the database
               of the Clone. Most likely this will be used by
               other adaptors to build objects from references to Clones
  Returntype : Bio::EnsEMBL::Clone
  Exceptions : thrown if $id is not provided
  Caller     : general

=cut

sub fetch_by_dbID {
    my ($self,$id) = @_;
    
    if( !defined $id ) {
	$self->throw("No internal ID provided");
    }

    my $clone = $self->_generic_sql_fetch("WHERE clone_id = $id");

    return $clone;
}


=head2 fetch_all

  Arg [1]    : none
  Example    : @clones = $clone_adaptor->fetch_all();
  Description: Retrieves every clone from the database.  
  Returntype : listref of Bio::EnsEMBL::Clone
  Exceptions : none
  Caller     : none

=cut

sub fetch_all {
  my $self = shift;

  my $sth = $self->prepare("SELECT clone_id, name, embl_acc, version,
                                   embl_version, htg_phase, 
                                   UNIX_TIMESTAMP(created), UNIX_TIMESTAMP(modified)
                            FROM clone");

  $sth->execute();

  my ($clone_id, $name, $embl_acc, $version, $embl_version, $htg_phase,
      $created, $modified);

  $sth->bind_columns(\$clone_id, \$name, \$embl_acc, \$version, \$embl_version,
		     \$htg_phase, \$created, \$modified);

  my @clones;

  while($sth->fetch()) {
    push @clones, new Bio::EnsEMBL::Clone($self, $clone_id, $name, $embl_acc,
					  $version, $embl_version, $htg_phase,
					  $created, $modified);
  }

  return \@clones;
}



=head2 list_embl_version_by_accesssion

  Arg [1]    : string $id
               the EMBL accession of the clone versions to retrieve
  Example    : @vers = $clone->list_embl_version_by_accession($accession) 
  Description: Returns a list of versions for a given EMBL accession
  Returntype : listref of ints
  Exceptions : thrown if $id arg is not defined or if no clone with accession
               $id exists in the database
  Caller     : general

=cut

sub list_embl_version_by_accession {
    my ($self,$id) = @_;
    my (@vers);

    if( !defined $id) {$self->throw("Don't have $id for new adaptor");}

    my $sth = $self->prepare(qq{
	SELECT distinct embl_version
	FROM   clone
	WHERE  embl_acc = '$id'
    });
    my $res = $sth ->execute();

    while( my $rowhash = $sth->fetchrow_hashref) {
	push @vers, $rowhash->{'embl_version'};
    }
    $sth->finish;

    $self->throw("no clone $id") unless scalar @vers > 0;
    
    return \@vers;
}


=head2 remove

  Arg [1]    : Bio::Ensembl::Clone $clone 
  Example    : $clone_adaptor->remove($clone);
  Description: Deletes clone (itself), including contigs and features, 
               but not its genes 
  Returntype : none
  Exceptions : thrown if the clone deletion fails
               throw if attached to the wrong database
  Caller     : ?

=cut

sub remove {
  my ($self, $clone) = @_;
  
  if ($clone->adaptor ne $self) {
    $self->throw("Trying to delete a clone attached to a different database.");
  }

  # Delete all contigs and  related features
  my $contigs = $clone->get_all_Contigs;
  my $rca = $self->db->get_RawContigAdaptor;
  
  foreach my $contig ( @$contigs ) {
    $rca->remove($contig);
  }

  # Delete the row for the clone
  my $sth = $self->prepare("DELETE FROM clone WHERE clone_id = ?");
  $sth->execute($clone->dbID);
  $self->throw("Failed to delete clone for clone_id '$clone->dbID'")
    unless $sth->rows;
}


=head2 get_all_Genes

  Arg [1]    : string $clone_id
               the EMBL accession for the clone from which the genes are to be
               retrieved
  Example    : my @genes = $clone_adaptor->get_all_Genes('AC011082');
  Description: Retrieves a list of Gene objects which are present on a clone.
               It would be better to have this on the gene adaptor but 
               for now it will stay here.
  Returntype : listref of Bio::EnsEMBL::Genes
  Exceptions : thrown if $clone_id is not defined
  Caller     : Clone::get_all_Genes

=cut

sub get_all_Genes {
    my ($self, $clone_id) = @_;

    $self->throw("clone_id not given") unless $clone_id;

    my $sth = $self->prepare(qq{
        SELECT t.gene_id
        FROM transcript t
          , exon_transcript et
          , exon e
          , contig c
        WHERE e.contig_id = c.contig_id
          AND et.exon_id = e.exon_id
          AND t.transcript_id = et.transcript_id
          AND c.clone_id = $clone_id
        });
    $sth->execute();
   
    my $geneAdaptor = $self->db->get_GeneAdaptor();
    my( %got, @genes );
    while (my ($gene_id) = $sth->fetchrow) { 
        unless ($got{$gene_id}) {
            if (my $gene = $geneAdaptor->fetch_by_dbID($gene_id)) {
                push(@genes, $gene);
            }
            $got{$gene_id} = 1;
        }
    }
    return \@genes;
}


=head2 store

  Arg [1]    : Bio::EnsEMBL::Clone $clone 
               the Clone to store in the database 
  Example    : $clone_adaptor->store($clone);
  Description: Stores a clone object and its associated contigs
               in the database and returns the dbID of the new db record
  Returntype : int
  Exceptions : thrown if $clone is not defined or if $clone is not a 
               Bio::EnsEMBL::Clone or if the database insertion fails
  Caller     : general

=cut

sub store{
  my ($self, $clone) = @_;

  unless($clone) {
    $self->throw("trying to write a clone without a clone object : $!\n");
  }

  my @contigs = @{$clone->get_all_Contigs};

  if( !$clone->isa('Bio::EnsEMBL::Clone') ) {
    $self->throw("Clone '$clone' is not a 'Bio::EnsEMBL::Clone'");
  }
  
  my $sql =  "insert into clone(name, 
                                embl_acc, 
                                version, 
                                embl_version, 
                                htg_phase, 
                                created,
                                modified) 
              values( '".$clone->id."' , '".
                         $clone->embl_id."', ".
                         $clone->version.",".
                         $clone->embl_version.", ".
                         $clone->htg_phase.", 
                         FROM_UNIXTIME('".$clone->created."'), 
                         FROM_UNIXTIME('".$clone->modified."'))";

  my $sth = $self->prepare($sql);
  my $rv = $sth->execute();
  $self->throw("Failed to insert clone $clone->id") unless $rv;

  $sth->finish();			   
  $sth = $self->prepare("select last_insert_id()");
  $sth->execute;

  my ($id) = $sth->fetchrow();

  $sth->finish();

  #store the contigs which were on this clone
  my $rca = $self->db->get_RawContigAdaptor();
  foreach my $contig(@contigs){
    $rca->store($contig, $clone);
  }

  #update this clones database identifier
  $clone->dbID($id);
  $clone->adaptor($self);

  return $id;
}


1;
