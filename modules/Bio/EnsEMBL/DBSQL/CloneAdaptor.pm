
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

Bio::EnsEMBL::DBSQL::CloneAdapor

=head1 SYNOPSIS

    # $db is Bio::EnsEMBL::DB::DBAdaptor

    my $da= Bio::EnsEMBL::DBSQL::CloneAdaptor->new($obj);
    my $clone=$da->fetch($id);

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
use Bio::EnsEMBL::DBSQL::Gene_Obj;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



=head2 fetch_by_accession

 Title   : fetch_by_accession
 Usage   :
 Function:
 Example :
 Returns : Bio::EnsEMBL::Clone
 Args    :


=cut


# setup a cache for clone accno's:
my %clone_cache;

sub fetch_by_accession { 
    my ($self,$id) = @_;

    if( !defined $id) {$self->throw("Don't have $id for new adaptor");}

    if ( defined $clone_cache{$id} ) {
        return $clone_cache{$id};
    }

    my $statement="select internal_id,embl_id,version,embl_version,htg_phase,
                          UNIX_TIMESTAMP(created),UNIX_TIMESTAMP(modified),
                          UNIX_TIMESTAMP(stored) 
                   from   clone 
                  where   id = '$id'
                 order by embl_version desc";

   
    my $sth = $self->prepare($statement);    
    my $res = $sth ->execute();

    my ($internal_id,$embl_id,$version,$embl_version,
	$htg_phase,$created,$modified, $stored)= $sth->fetchrow_array;
   
    $self->throw("no clone for $id") unless defined $internal_id;
    
    my @args=($internal_id,$id,$embl_id,$version,$embl_version,$htg_phase,$created,$modified, $stored);

    my $clone = Bio::EnsEMBL::Clone->new($self,@args);

    $clone_cache{$id}= $clone;
    return $clone;
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
	WHERE  id = '$id'
    });
    my $res = $sth ->execute();

    while( my $rowhash = $sth->fetchrow_hashref) {
	push @vers, $rowhash->{'embl_version'};
    }
    $sth->finish;

    $self->throw("no clone $id") unless scalar @vers > 0;
    
    return @vers;
}

=head2 fetch_by_accession_version

 Title   : fetch_by_accession_version
 Usage   :
 Function:
 Example :
 Returns : Bio::EnsEMBL::Clone
 Args    :


=cut

sub fetch_by_accession_version { 
    my ($self,$id,$ver) = @_;

    if( !defined $id) {$self->throw("Don't have $id for new adaptor");}

    my $statement="select internal_id,embl_id,version,embl_version,htg_phase,
                          UNIX_TIMESTAMP(created),UNIX_TIMESTAMP(modified),
                          UNIX_TIMESTAMP(stored) 
                   from   clone 
                   where  id           = '$id'
                   and    embl_version = $ver";

   
    my $sth = $self->prepare($statement);    
    my $res = $sth ->execute();

    my ($internal_id,$embl_id,$version,$embl_version,
	$htg_phase,$created,$modified, $stored)= $sth->fetchrow_array;
   
    $self->throw("no clone $id with version $ver") unless defined $internal_id;
    
    my @args=($internal_id,$id,$embl_id,$version,$embl_version,$htg_phase,$created,$modified, $stored);
    
    return Bio::EnsEMBL::Clone->new($self,@args);
}


sub fetch_by_embl_id {
    my ($self, $embl_id) = @_;

    $self->throw("Can't fetch clone without EMBL id (got '$embl_id')")
        unless $embl_id;

    my $sth = $self->prepare(q{
        SELECT internal_id
          , id
          , embl_id
          , version
          , embl_version
          , htg_phase
          , UNIX_TIMESTAMP(created)
          , UNIX_TIMESTAMP(modified)
          , UNIX_TIMESTAMP(stored)
        FROM clone
        WHERE embl_id = ?
        });
    $sth->execute($embl_id);
    
    my @fields = $sth->fetchrow
        or $self->throw("No Clone with embl_id '$embl_id'");
    
    return Bio::EnsEMBL::Clone->new($self, @fields);
}



sub fetch
{
    my ($self,$id)=@_;
    $self->warn("fetch is now deprecated, use fetch_by_accession instead");
    $self->fetch_by_accession($id);
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
   my ($self,$internal_id) = @_;
   
   my @contigs;
   my @dnas;

   # get a list of contigs to zap

   my $sth = $self->prepare("select internal_id,dna from contig where clone = $internal_id");
   my $res = $sth->execute;

   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@contigs,$rowhash->{'internal_id'});
       push(@dnas,$rowhash->{'dna'});
   }
   
   # Delete from DNA table, Contig table, Clone table
   
   foreach my $contig ( @contigs ) {
       my $sth = $self->prepare("delete from contig where internal_id = $contig");
       my $res = $sth->execute;
   }


   if ($self->db ne $self->db->dnadb) {
     $self->warn("Using a remote dna database - can't delete dna\n");
   } else {
     
     foreach my $dna (@dnas) {
       $sth = $self->prepare("delete from dna where id = $dna");
       $res = $sth->execute;
       
     }

   }

   $sth = $self->prepare("delete from clone where internal_id = $internal_id");
   $res = $sth->execute;
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
   my ($self,$clone_id,$supporting) = @_;
   my @out;

   $self->throw("I need an id") unless $clone_id;
 
   my %got;
   
   my $sth = $self->prepare("
        SELECT t.gene
        FROM transcript t,
             exon_transcript et,
             exon e,
             contig c
        WHERE e.contig = c.internal_id
          AND et.exon = e.id
          AND t.id = et.transcript
          AND c.clone = $clone_id
        ");

    my $res = $sth->execute();
   
    while (my $rowhash = $sth->fetchrow_hashref) { 
            
        if( ! exists $got{$rowhash->{'gene'}}) {  
            
           my $gene_obj = Bio::EnsEMBL::DBSQL::Gene_Obj->new($self->db);             
	   my $gene = $gene_obj->get($rowhash->{'gene'}, $supporting);
           if ($gene) {
	        push(@out, $gene);
           }
	   $got{$rowhash->{'gene'}} = 1;
        }       
    }
   
    if (@out) {
        return @out;
    }
    return;
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

=head2 get_all_geneid

 Title   : get_all_geneid
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_my_geneid {
   my ($self,$cloneid) = @_;

   my $sth = $self->prepare("select count(*),cont.clone ,ex.contig,tran.gene  " .
			    "from   contig          as cont, ".
			    "       transcript      as tran, " .
			    "       exon_transcript as et, " .
			    "       exon            as ex " .
			    "where  ex.id            = et.exon " .
			    "and    tran.id          = et.transcript " .
			    "and    cont.clone       = $cloneid  " .
			    "and    cont.internal_id = ex.contig " .
			    "group by tran.gene");

   my @out;

   $sth->execute;
   while( my $rowhash = $sth->fetchrow_hashref) {
       push(@out,$rowhash->{'gene'});
   }
   return @out;
}

=head2 get_all_Contigs

 Title   : get_Contigs
 Usage   : foreach $contig ( $clone->get_all_Contigs ) 
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_Contigs {
   my ($self,$internal_id,$version) = @_;
   my $sth;
   my @res;
  
   my $sql = "select id,internal_id from contig where clone = $internal_id";

   $sth= $self->prepare($sql);
   my $res  = $sth->execute();
   my $seen = 0;

   my $count   = 0;
   my $total   = 0;
   

   while( my $rowhash = $sth->fetchrow_hashref) {
       my $contig = $self->db->get_Contig( $rowhash->{'id'});
       $contig->internal_id($rowhash->{internal_id});
       $contig->seq_version($version);

       push(@res,$contig);
       $seen = 1;
   }

   if( $seen == 0  ) {
       $self->throw("Clone [$internal_id] has no contigs in the database. Should be impossible, but clearly isn't...");
   }

   return @res;   
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
      stored        datetime NOT NULL,
  
      PRIMARY KEY (internal_id),
      UNIQUE embl (embl_id,embl_version),
      UNIQUE id   (id,embl_version)
     )   
  } );
  $sth->execute();
}



1;
