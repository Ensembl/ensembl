=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::DBSQL::CompressedSequenceAdaptor - Facilitates DB storage and retrieval of compressed sequence

=head1 SYNOPSIS

  $seq_adptr = $database_adaptor->get_SequenceAdaptor();

  $dna =
    ${ $seq_adptr->fetch_by_Slice_start_end_strand( $slice, 1, 1000,
      -1 ) };

=head1 DESCRIPTION

An adaptor for the retrieval of compressed DNA sequence from the EnsEMBL 
database

=head1 METHODS

=cut

package Bio::EnsEMBL::DBSQL::CompressedSequenceAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::SequenceAdaptor;

@ISA = qw(Bio::EnsEMBL::DBSQL::SequenceAdaptor);


sub _fetch_seq {
  my $self          = shift;
  my $seq_region_id = shift;
  my $start         = shift;
  my $len           = shift;

  #calculate the offset and start in the compressed sequence 
  my $comp_start  = ($start-1 >> 2) + 1;
  my $comp_len    = ($len  >> 2) + 2;

  my ($bvector, $nline);

  my $sth = $self->prepare(
               "SELECT SUBSTRING( d.sequence, ?, ?), n_line
                FROM dnac d
                WHERE d.seq_region_id = ?");
  $sth->bind_param(1,$comp_start,SQL_INTEGER);
  $sth->bind_param(2,$comp_len ,SQL_INTEGER);
  $sth->bind_param(3,$seq_region_id,SQL_INTEGER);
  $sth->execute();
  $sth->bind_columns(\$bvector, \$nline);
  $sth->fetch();
  $sth->finish();

  #convert sequence from binary string to 0123 string
  my $bitlen = length($bvector) << 2;
  my $str = '';
  for(my $i=0; $i < $bitlen; $i++) {
    $str .= vec($bvector, $i, 2);
  }

  #convert from 0123 to ACTG
  $str =~ tr/0123/ACTG/;

  $str = substr($str, ($start-1)%4, $len);

  #expand the nlines and place them back in the sequence
  my @nlines = split(/:/, $nline);
  foreach my $nl (@nlines) {
    my ($offset,$char,$nlen) = $nl =~ /(\d+)(\D)(\d+)/;
    
    #skip nlines entirely out of range
    next if(($offset+$nlen-1) < $start || $offset > ($start+$len-1));
    
    #obtain relative offset into requested region
    $offset = $offset -  $start + 1;

    #nlines that partially overlap requested region have to be shrunk
    if($offset < 1) {
      $nlen = $nlen - (1-$offset);
      $offset = 1;
    }
    if($offset + $nlen > $start+$len) {
      $nlen = $len - $offset + 1;
    }

    substr($str,$offset-1,$nlen) = $char x $nlen;    
  }

  return \$str;
}


=head2 store

  Arg [1]    : string $seq_region_id the id of the sequence region this dna
               will be associated with.
  Arg [2]    : string reference  $sequence the dna sequence to be stored in 
               the database
  Example    : $dbID = $seq_adaptor->store(12,\'ACTGGGTACCAAACAAACACAACA');
  Description: stores a dna sequence in the databases dna table and returns the
               database identifier for the new record.
  Returntype : int
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::RawContigAdaptor::store
  Status     : Stable

=cut

sub store {
  my ($self, $seq_region_id, $sequence) = @_;

  if(!$seq_region_id) {
    throw('seq_region_id is required');
  }

  $sequence = uc($sequence);

  my $bvector = '';

  #convert sequence to 0s,1s,2s and 3s
  $sequence =~ tr/ACTG/0123/;

  #nlines cover sequence which is not ACTG such as N
  #nline format is a set of colon delimited int, char, int triplets:
  #<offset><code><length>
  my($nline_char,$nline_len,$nline_off);
  my @nlines;

  my $len    = length($sequence);
  for(my $i=0; $i < $len; $i++) {
    my $char = substr($sequence,$i,1);

    #quickly check if this character was an A,C,T or G (and was converted to
    # a 0,1,2,3)
    if($char =~ /[0-3]/) {
      vec($bvector, $i,2) = $char;
      if($nline_char) {
        #end of an nline
        push @nlines, "$nline_off$nline_char$nline_len";
        $nline_char = undef;
        $nline_len  = 0;
        $nline_off  = 0;
      } 
    } else {
      #this was not an ACTG
      if($nline_char) {
        if($nline_char eq $char) {
          #continuation of an nline
          $nline_len++;
        } else {
          #end of a previous nline and start of a new one
          push @nlines, "$nline_off$nline_char$nline_len";
          $nline_char = $char;
          $nline_len  = 1;
          $nline_off  = $i+1;
        }
      } else {
        #start of a new nline
        $nline_char = $char;
        $nline_len  = 1;
        $nline_off  = $i+1;
      }
      $char = 0; #need to put numeric val into bitvector despite nline
    }

    vec($bvector, $i,2) = $char; 
  }

  my $nline = join(':', @nlines);
  my $statement = $self->prepare(
        "INSERT INTO dnac(seq_region_id, sequence, n_line) VALUES(?,?,?)");

  $statement->bind_param(1,$seq_region_id,SQL_INTEGER);
  $statement->bind_param(2,$bvector,SQL_BLOB);
  $statement->bind_param(3,$nline,SQL_LONGVARCHAR);
  $statement->execute();

  $statement->finish();
  return;
}


1;
