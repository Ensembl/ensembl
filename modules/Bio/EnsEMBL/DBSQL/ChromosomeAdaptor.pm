#
# Ensembl module for Bio::EnsEMBL::DBSQL::ChromosomeAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::ChromosomeAdaptor - DESCRIPTION of Object

=head1 SYNOPSIS

$chromosome_adaptor = $db_adaptor->get_ChromosomeAdaptor();
$chromosome = $chromosome_adaptor->fetch_by_chr_name('12');

=head1 DESCRIPTION

This is a database adaptor used to retrieve chromosome objects from a database.

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk


=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBSQL::ChromosomeAdaptor;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Chromosome;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

# new inherited from BaseAdaptor



=head2 fetch_by_dbID

  Arg [1]    : int $id 
               unique database identifier for chromosome to retrieve  
  Example    : my $chromosome = $chromosome_adaptor->fetch_by_dbID(1);
  Description: Retrieves a Chromosome object from the database using its
               unique identifier.  Note the the identifier is the dbID and
               does NOT correspond to the chromosome name.  dbID 1 does NOT
               necessarily correspond to chromosome '1'
  Returntype : Bio::EnsEMBL::Chromosome
  Exceptions : thrown if $id not defined
  Caller     : general

=cut

sub fetch_by_dbID {
  my ($self,$id) = @_;

  my $chr = (); 

  unless(defined $id) {
    $self->throw("Chromosome dbID argument required\n");
  }

  unless(defined $self->{'_chr_cache'} ) {
    $self->{'_chr_cache'} = {};
    $self->{'_chr_name_cache'} = {};
  }

  #
  # If there is not already a cached version of this chromosome pull it
  # from the database and add it to the cache.
  #
  unless($chr = $self->{'_chr_cache'}->{$id}) {
    my $sth = $self->prepare( "SELECT name, known_genes, unknown_genes, 
                                      snps, length
                               FROM chromosome
                               WHERE chromosome_id = ?" );
    $sth->execute( $id );
    
    my($name, $known_genes, $unknown_genes, $snps, $length);
    $sth->bind_columns(\$name,\$known_genes,\$unknown_genes,\$snps,\$length); 
    $sth->fetch();
   

    if($@) {
      $self->throw("Could not create chromosome from dbID $id\n" .
		   "Exception: $@\n");
    }

    unless($name) {
      $self->throw("Could determine chromosome name from dbID $id\n");
    }

    $chr = new Bio::EnsEMBL::Chromosome( -adaptor => $self,
                                         -dbID => $id,
					 -chr_name => $name,
					 -known_genes => $known_genes,
					 -unknown_genes => $unknown_genes,
					 -snps => $snps,
					 '-length' => $length );

    $self->{'_chr_cache'}->{$id} = $chr;
    $self->{'_chr_name_cache'}->{$name} = $chr;
    
  }

  return $chr;
}


=head2 fetch_by_chr_name

  Arg [1]    : string $chr_name
               the name of the chromosome to retrieve
  Example    : $chromosome = $chromosome_adaptor->fetch_by_chr_name('X');
  Description: Retrieves a chromosome object from the database using its name.
  Returntype : Bio::EnsEMBL::Chromosome
  Exceptions : none
  Caller     : general

=cut

sub fetch_by_chr_name{
   my ($self,$chr_name) = @_;

  my $chr = undef; 

  unless(defined $chr_name) {
    $self->throw("Chromosome name argument required\n");
  }

  unless(defined $self->{'_chr_cache'} ) {
    $self->{'_chr_cache'} = {};
    $self->{'_chr_name_cache'} = {};
  }

  #
  # If there is not already a cached version of this chromosome pull it
  # from the database and add it to the cache.
  #
  unless($chr = $self->{'_chr_name_cache'}->{$chr_name}) {
    my $sth = $self->prepare( "SELECT chromosome_id, known_genes, unknown_genes, 
                                      snps, length
                               FROM chromosome
                               WHERE name = ?" );
    $sth->execute( $chr_name );
    
    my($dbID, $known_genes, $unknown_genes, $snps, $length);
    $sth->bind_columns(\$dbID,\$known_genes,\$unknown_genes,\$snps,\$length); 

    if ($sth->rows > 0) {
    $sth->fetch();

    $chr = new Bio::EnsEMBL::Chromosome( -adaptor => $self,
					 -dbID => $dbID,
					 -chr_name => $chr_name,
					 -known_genes => $known_genes,
					 -unknown_genes => $unknown_genes,
					 -snps => $snps,
					 '-length' => $length );

    $self->{'_chr_cache'}->{$dbID} = $chr;
    $self->{'_chr_name_cache'}->{$chr_name} = $chr;
  }
  }
  return $chr;
}


=head2 fetch_all

  Args       : none
  Example    : @chromosomes = $chromosome_adaptor->fetch_all(); 
  Description: Retrieves every chromosome object from the database.
  Returntype : listref of Bio::EnsEMBL::Chromosome
  Exceptions : none
  Caller     : general

=cut

sub fetch_all {
  my($self) = @_;
  
  my @chrs = (); 


    my $sth = $self->prepare( "SELECT chromosome_id, name, known_genes, 
                                      unknown_genes, snps, length
                               FROM chromosome" );
    $sth->execute();
    
    my($chromosome_id, $name, $known_genes, $unknown_genes, $snps, $length);
    $sth->bind_columns(\$chromosome_id,\$name,\$known_genes,
                       \$unknown_genes,\$snps,\$length); 

    while($sth->fetch()) {
   
    my $chr = new Bio::EnsEMBL::Chromosome( -adaptor => $self,
					 -chr_name => $name,
					 -dbID => $chromosome_id,
					 -known_genes => $known_genes,
					 -unknown_genes => $unknown_genes,
					 -snps => $snps,
					 '-length' => $length );

    $self->{'_chr_cache'}->{$chromosome_id} = $chr;
    $self->{'_chr_name_cache'}->{$name} = $chr;
    push @chrs, $chr;
  }

  return \@chrs;
}


=head2 get_dbID_by_chr_name

  Arg [1]    : string $chr_name
               the name of the chromosome whose dbID is wanted.
  Example    : $dbID = $chromosome_adaptor->fetch_by_dbID('X') 
  Description: Retrieves a unique database identifier for a chromosome
               using the chromosomes name.  It is not recommended that this
               method be used externally from ChromosomeAdaptor.  It should 
               probably be private and may be made private in the future. A
               better way to obtain a dbID is:
               $dbID = $chromosome_adaptor->fetch_by_chr_name('X')->dbID();
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::ChromosomeAdaptor

=cut

sub get_dbID_by_chr_name {
  my ($self, $chr_name) = @_;

  unless (defined $self->{_chr_name_mapping}) {
    $self->{_chr_name_mapping} = {};

    #get the chromo names and ids from the database
    my $sth = $self->prepare('SELECT name, chromosome_id FROM chromosome');
    $sth->execute();
    
    #Construct the mapping of chromosome name to id
    while(my $a = $sth->fetchrow_arrayref()) {
      $self->{_chr_name_mapping}->{$a->[0]} = $a->[1];
    }
  }

  return $self->{_chr_name_mapping}->{$chr_name};
}    





=head2 get_landmark_MarkerFeatures

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED use Slice::get_landmark_MarkerFeatures instead
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_landmark_MarkerFeatures{
   my ($self,$chr_name,$glob) = @_;

   $self->warn("ChromosomeAdaptor::get_landmark_MarkerFeatures is deprecated " 
	       . "use Slice::get_landmark_MarkerFeatures instead\n");

   if( !defined $glob ) {
       $glob = 500000;
   }

   my $statement= " SELECT  chr_start,
			    chr_end,
			    chr_strand,
			    name 
		    FROM    landmark_marker 
		    WHERE   chr_name = '$chr_name'
		    ORDER BY chr_start
		";
   
   $statement =~ s/\s+/ /g;
   
   my $sth = $self->prepare($statement);
   $sth->execute;
   
   my ($start, $end, $strand, $name);
   
   my $analysis;
   my %analhash;
   
   $sth->bind_columns
       ( undef, \$start, \$end,  \$strand, \$name);
   
   my @out;
   my $prev;
   while( $sth->fetch ) {
       if( defined $prev && $prev->end + $glob > $start  && $prev->id eq $name ) {
           next;
       }

       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->start($start);
       $sf->end($end);
       $sf->strand($strand);
       $sf->id($name);
       push(@out,$sf);
       $prev = $sf;
   } 

   return @out;
}


=head2 get_landmark_MarkerFeatures_old

  Arg [1]    : none
  Example    : none
  Description: DEPRECATED do not use
  Returntype : none
  Exceptions : none
  Caller     : none

=cut

sub get_landmark_MarkerFeatures_old{
   my ($self,$chr_name) = @_;

   my $glob = 1000;
   $self->throw( "Method deprecated. " );

   return ();

#   my $statement= "   SELECT 
#                       IF     (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-1),
#                              (sgp.chr_start+sgp.raw_end-f.seq_end-1)),                                        
#                       IF     (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-1),
#                              (sgp.chr_start+sgp.raw_end-f.seq_start-1)), 
#                              f.score, 
#                       IF     (sgp.raw_ori=1,f.strand,(-f.strand)), 
#                              f.name, f.hstart, f.hend, 
#                              f.hid, f.analysis, c.name 
#                       FROM   contig_landmarkMarker c,
#                              static_golden_path sgp,
#                              feature f
#                       WHERE  f.contig = c.contig
#                       AND    f.hid=c.marker  
#                       AND    sgp.raw_id=f.contig 
#                       AND    sgp.chr_name='$chr_name'";
   
#   $statement =~ s/\s+/ /g;
   
#   my $sth = $self->prepare($statement);
#   $sth->execute;
   
#   my ($start, $end, $score, $strand, $hstart, 
#       $name, $hend, $hid, $analysisid,$synonym);
   
#   my $analysis;
#   my %analhash;
   
#   $sth->bind_columns
#       ( undef, \$start, \$end, \$score, \$strand, \$name, 
#	 \$hstart, \$hend, \$hid, \$analysisid,\$synonym);
   
#   my @out;
#   while( $sth->fetch ) {
#       my $sf = Bio::EnsEMBL::SeqFeature->new();
#       $sf->start($start);
#       $sf->end($end);
#       $sf->strand($strand);
#       $sf->id($synonym);
#       push(@out,$sf);
#   } 

#   return @out;
}


sub store{
  my ($self, $chromosome) = @_;

  my $chr_name = $chromosome->chr_name;
  if(!$chr_name){
    $self->throw("can't store a chromosome without a name");
  }
  my $length = $chromosome->length;
  if(!$length){
    $self->throw("can't store a chromosome without a length");
  }
  my $known_genes = $chromosome->known_genes;
  if(!$known_genes){
    $known_genes = 0;
  } 
  my $unknown_genes = $chromosome->unknown_genes;
  if(!$unknown_genes){
    $unknown_genes = 0;
  } 
  my $snps = $chromosome->snps;
  if(!$snps){
    $snps = 0;
  } 
  my $sql = "insert into chromosome(name,
                                    known_genes,
                                    unknown_genes,
                                    snps,
                                    length)
                             values(?, ?, ?, ?, ?)";

  my $sth = $self->db->prepare($sql);
  $sth->execute($chr_name, $known_genes, $unknown_genes, $snps, $length);
  $chromosome->dbID($sth->{'mysql_insertid'});
  $chromosome->adaptor($self);
}

1;
