
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

    # $db is Bio::EnsEMBL::DB::DBAdaptor

    my $da= Bio::EnsEMBL::DBSQL::CloneAdaptor->new($obj);
    my $clone=$da->fetch_by_dbID($id);

    @contig = $clone->get_all_Contigs();
    @genes    = $clone->get_all_Genes();

=head1 DESCRIPTION

Represents information on one Clone

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::CloneAdaptor;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Clone;


@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


sub _generic_sql_fetch {
    my( $self, $where_clause ) = @_;
    
    my $sql = q{
        SELECT clone_id
          , embl_acc
          , name
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

 Function: fetches a clone by its EMBL/GenBank accession number
           It will fetch the highest version in the database

 Args    : accession number, as a string, no version number

 Returns : Bio::EnsEMBL::Clone

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

 Function: fetches a specific version of an accession number
 Args   1: accession number as a string
        2: version as a digit 
 Returns : Bio::EnsEMBL::Clone

=cut

sub fetch_by_accession_version { 
    my ($self, $acc, $ver) = @_;

    unless ($acc and $ver) {
        $self->throw("Need both accession (got '$acc') and version (got '$ver')");
    }

    my $clone = $self->_generic_sql_fetch(
        qq{ WHERE embl_acc = '$acc' AND embl_version = '$ver' }
        );
    if ($clone) {
        return $clone;
    } else {
        $self->throw("no clone with accession '$acc' and version '$ver'");
    }
}

=head2 fetch_by_name

 Function: fetches a single clone by the "name" field. The name
           field is often the EMBL ID of a clone (which itself is often
           the same as the accession number) but is a deliberately open
           slot for another identifier for the clone, eg. the stringified
           tracking identifier in a sequencing centre
 Args    : name as a string
 Returns : Bio::EnsEMBL::Clone

=cut

sub fetch_by_name {
    my ($self, $name) = @_;
    print $name."\n";
    $self->throw("name not given") unless $name;

    my $clone = $self->_generic_sql_fetch(
        qq{ WHERE name = '$name' }
        );

    if ($clone) {
        return $clone;
    } else {
        $self->throw("No Clone with name '$name'");
    }
}

=head2 fetch_by_dbID

 Function: fetches a clone by the internal id in the database
           of the Clone. Most likely this will be used by
           other adaptors to build objects from references to Clones
 Args    : the numeric internal id of the clone table
 Returns : Bio::EnsEMBL::Clone

=cut



sub fetch_by_dbID {
    my ($self,$id) = @_;
    
    if( !defined $id ) {
	$self->throw("No internal ID provided");
    }

    my $clone = $self->_generic_sql_fetch("WHERE clone_id = $id");

    return $clone;
}


sub fetch {
    my ($self,$id) = @_;
    $self->warn("fetch is now deprecated, use fetch_by_accession instead");
    $self->fetch_by_accession($id);
}





=head2 list_embl_version_by_accession

 Title   : list_embl_version_by_accession
 Usage   : @vers = $obj->list_embl_version_by_accession($accession)
 Function:
 Example :
 Returns : @vers
 Args    : $accession


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
    
    return @vers;
}


=head2 delete_by_dbID

 Title   : delete_by_dbID
 Usage   : $clone->delete_by_dbID()
 Function: Deletes clone (itself), including contigs and features, but not its genes
 Example : 
 Returns : nothing
 Args    : none


=cut

sub delete_by_dbID {
    my ($self, $clone_id) = @_;

    my $fadaptor = $self->db->get_FeatureAdaptor;

    # Make a list of all contig and dna entries to delete
    my $sth = $self->prepare(qq{
        SELECT contig_id, dna_id
        FROM contig
        WHERE clone_id = $clone_id
        });
    $sth->execute;

    my( @contigs, @dnas );
    while( my ($c, $d) = $sth->fetchrow_hashref) {
        push(@contigs, $c);
        push(@dnas,    $d);
    }

    # Delete features for each contig, and each contig
    foreach my $contig_id ( @contigs ) {
        $fadaptor->delete_by_RawContig_internal_id($contig_id);
        my $sth = $self->prepare("DELETE FROM contig WHERE contig_id = $contig_id");
        $sth->execute;
        $self->throw("Failed to delete contigs for contig_id '$contig_id'")
            unless $sth->rows;
    }

    # Delete DNA as long as we aren't using a remote DNA database.
    if ($self->db ne $self->db->dnadb) {
        $self->warn("Using a remote dna database - not deleting dna\n");
    } else {
        foreach my $dna_id (@dnas) {
            $sth = $self->prepare("DELETE FROM dna WHERE dna_id = $dna_id");
            $sth->execute;
            $self->throw("Failed to delete dna for dna_id '$dna_id'")
                unless $sth->rows;
        }
    }

    # Delete the row for the clone
    $sth = $self->prepare("DELETE FROM clone WHERE clone_id = $clone_id");
    $sth->execute;
    $self->throw("Failed to delete clone for clone_id '$clone_id'")
        unless $sth->rows;
}


=head2 get_all_Genes

 Title   : get_all_Genes
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

# will contact geneAdaptor when ready
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
    return @genes;
}


=head2 get_Contig

 Title   : get_Contig
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_Contig {
   my ($self,$contigid) = @_;
   
   my $contig = $self->db->get_Contig($contigid);
   
   return $contig->fetch();
}


# creates all tables for this adaptor
# if they exist they are emptied and newly created
sub create_tables {
  my $self = shift;

  my $sth = $self->prepare( "drop table if exists clone" );
  $sth->execute();

  $sth = $self->prepare( qq{
    CREATE TABLE clone (
      internal_id   int(10) unsigned NOT NULL auto_increment,
      id            varchar(40) NOT NULL,
      embl_id       varchar(40) NOT NULL,
      version       int(10) NOT NULL,
      embl_version  int(10) NOT NULL,
      htg_phase     int(10) DEFAULT '-1' NOT NULL,
      created       datetime NOT NULL,
      modified      datetime NOT NULL,
  
      PRIMARY KEY (internal_id),
      UNIQUE embl (embl_id,embl_version),
      UNIQUE id   (id,embl_version)
     )   
  } );
  $sth->execute();
}


sub store{
  my ($self, $clone) = @_;

  $clone || $self->throw("trying to write a clone without a clone object : $!\n");

  if( !$clone->isa('Bio::EnsEMBL::DB::CloneI') ) {
	$self->throw("Clone '$clone' is not a 'Bio::EnsEMBL::DB::CloneI'");
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
                         FROM_UNIXTIME(".$clone->created."), 
                         FROM_UNIXTIME(".$clone->modified."))";
  

  my $sth = $self->prepare($sql);
  my $rv = $sth->execute();
  
  $self->throw("Failed to insert clone $clone->id") unless $rv;			   
  $sth = $self->prepare("select last_insert_id()");
  my $res = $sth->execute;
  my $row = $sth->fetchrow_hashref;
  $sth->finish;
  my $id  = $row->{'last_insert_id()'};

  foreach my $contig($clone->get_all_Contigs()){
    my $rca = $self->db->get_RawContigAdaptor();
    $rca->store($contig, $id);
  }


}


1;
