

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

# Object preamble - inherits from Bio::Root::RootI

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Chromosome;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

# we inheriet new of BaseAdaptor

=head2 fetch_by_dbID

 Title   : fetch_by_dbID
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_dbID{
   my ($self,$id) = @_;

   # should check is correct!
   my $sth= $self->prepare("select chr_name from static_golden_path where chr_name = '$id' limit 1");
   $sth->execute;
   my $a = $sth->fetchrow_arrayref();
   if( !defined $a ) {
       $self->throw("No chromosome of $id");
   }


   my $chr = Bio::EnsEMBL::Chromosome->new( -adaptor => $self,-chrname => $id);

   return $chr;
}

=head2 fetch_by_chrname

 Title   : fetch_by_chrname
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub fetch_by_chrname{
   my ($self,$id) = @_;

   # for the moment, chain to fetch_by_dbID, but sometime this will change
   return $self->fetch_by_dbID($id);
}

=head2 get_landmark_MarkerFeatures

 Title   : get_landmarkMarkers
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_landmark_MarkerFeatures_old{
   my ($self,$chr_name) = @_;

   my $glob = 1000;

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
   #print STDERR "Doing Query ... $statement\n";
   
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

 Title   : get_landmarkMarkers
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_landmark_MarkerFeatures{
   my ($self,$chr_name,$glob) = @_;


   if( !defined $glob ) {
       $glob = 500000;
   }

   my $statement= " SELECT  start,
			    end,
			    strand,
			    name 
		    FROM    contig_landmarkMarker 
		    WHERE   chr_name = '$chr_name'
		    ORDER BY start
		";
   
   $statement =~ s/\s+/ /g;
   #print STDERR "Doing Query ... $statement\n";
   
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




