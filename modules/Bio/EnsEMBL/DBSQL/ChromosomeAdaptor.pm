

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

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

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

  Arg 1     : txt $chromosome_name
  Function  : return a chromosome object from the dbID of the chromosome.
  Returntype: Bio::EnsEMBL::Chromosome
  Exceptions: if chromosome_name not present in static_golden_path table
  Caller    : ??

=cut

sub fetch_by_dbID {
  my ($self,$id) = @_;

  my $chr = (); 

  unless(defined $id) {
    $self->throw("Chromosome dbID argument required\n");
  }

  unless(defined $self->{'_chr_cache'} ) {
    $self->{'_chr_cache'} = {};
  }

  #
  # If there is not already a cached version of this chromosome pull it
  # from the database and add it to the cache.
  #
  unless($chr = $self->{'_chr_cache'}->{$id}) {
    my $sth = $self->prepare( "SELECT name, known_genes, unknown_genes, 
                                      snps, length
                               FROM chromosome
                               WHERE chromosome_id = $id" );
    $sth->execute();
    
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
					 -chr_name => $name,
					 -chromosome_id => $id,
					 -known_genes => $known_genes,
					 -unknown_genes => $unknown_genes,
					 -snps => $snps,
					 '-length' => $length );

    $self->{'_chr_cache'}->{$id} = $chr;
  }

  return $chr;
}


=head2 fetch_by_chr_name

  Arg  1    : txt $chromosome_name
  Function  : chains to fetch_by_dbID
  Returntype: Bio::EnsEMBL::Chromosome
  Exceptions: see above
  Caller    : ??

=cut


sub fetch_by_chr_name{
   my ($self,$chr_name) = @_;

   #Convert the name to the dbID
   my $dbID = $self->get_dbID_by_chr_name($chr_name);

   unless(defined $dbID) {
     $self->warn("chromosome with name $chr_name not in database");
     return undef;
   } 
   
   return $self->fetch_by_dbID($dbID);
}

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
					 -chromosome_id => $chromosome_id,
					 -known_genes => $known_genes,
					 -unknown_genes => $unknown_genes,
					 -snps => $snps,
					 '-length' => $length );

    $self->{'_chr_cache'}->{$chromosome_id} = $chr;
    push @chrs, $chr;
  }

  return @chrs;
}



=head2 get_landmark_MarkerFeatures_old

  Arg  1    : txt $chromosome_name
  Function  : depracated, defunct
  Returntype: list Bio:EnsEMBL::SeqFeature
  Exceptions: none
  Caller    : ??

=cut


sub get_landmark_MarkerFeatures_old{
   my ($self,$chr_name) = @_;

   my $glob = 1000;
   $self->throw( "Method deprecated. " );

   my $statement= "   SELECT 
                       IF     (sgp.raw_ori=1,(f.seq_start+sgp.chr_start-sgp.raw_start-1),
                              (sgp.chr_start+sgp.raw_end-f.seq_end-1)),                                        
                       IF     (sgp.raw_ori=1,(f.seq_end+sgp.chr_start-sgp.raw_start-1),
                              (sgp.chr_start+sgp.raw_end-f.seq_start-1)), 
                              f.score, 
                       IF     (sgp.raw_ori=1,f.strand,(-f.strand)), 
                              f.name, f.hstart, f.hend, 
                              f.hid, f.analysis, c.name 
                       FROM   contig_landmarkMarker c,
                              static_golden_path sgp,
                              feature f
                       WHERE  f.contig = c.contig
                       AND    f.hid=c.marker  
                       AND    sgp.raw_id=f.contig 
                       AND    sgp.chr_name='$chr_name'";
   
   $statement =~ s/\s+/ /g;
   
   my $sth = $self->prepare($statement);
   $sth->execute;
   
   my ($start, $end, $score, $strand, $hstart, 
       $name, $hend, $hid, $analysisid,$synonym);
   
   my $analysis;
   my %analhash;
   
   $sth->bind_columns
       ( undef, \$start, \$end, \$score, \$strand, \$name, 
	 \$hstart, \$hend, \$hid, \$analysisid,\$synonym);
   
   my @out;
   while( $sth->fetch ) {
       my $sf = Bio::EnsEMBL::SeqFeature->new();
       $sf->start($start);
       $sf->end($end);
       $sf->strand($strand);
       $sf->id($synonym);
       push(@out,$sf);
   } 

   return @out;
}


=head2 get_landmark_MarkerFeatures

  Arg  1    : txt $chromosome_name
  Arg  2    : int $min_feature_distance (glob)
  Function  : landmark marker SeqFeatures from landmark_marker table
              on given chromosome name, ordered ascending. If same marker closer than
              $min_feature_distance on chromosome one is pruned from result.
  Returntype: list Bio::EnsEMBL::SeqFeature
  Exceptions: none
  Caller    : ??

=cut

sub get_landmark_MarkerFeatures{
   my ($self,$chr_name,$glob) = @_;


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


=head2 get_dbID_from_chr_name

  Arg  1    : txt $chr_name
  Function  : Given a chromosome name of the format 'chrN' returns the dbID
              of the given chromosome
  Returntype: int
  Exceptions: none
  Caller    : ??

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


