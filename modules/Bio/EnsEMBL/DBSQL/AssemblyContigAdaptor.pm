=head1 NAME

Bio::EnsEMBL::DBSQL::AssemblyContigAdaptor - MySQL Database queries to generate and store assembly information

=head1 SYNOPSIS

=head1 CONTACT

=head1 APPENDIX

=cut

package Bio::EnsEMBL::DBSQL::AssemblyContigAdaptor;

use vars qw( @ISA );
use strict;


use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Virtual::AssemblyContig;

@ISA = qw( Bio::EnsEMBL::DBSQL::BaseAdaptor );

=head2 fetch_by_dbID

 Title   : fetch_by_dbID
 Usage   : $exonAdaptor->fetch_by_dbID($dbid))
 Function: 
 Example : 
 Returns : nothing
 Args    : 

=cut


sub fetch_by_dbID {
  my $self = shift;
  my $dbID = shift;

  my $query = qq {
    SELECT  assembly_contig_id,
            assembly_contig_name,
            assembly_contig_start,
            assembly_contig_end,
            assembly_contig_chr_name,
            assembly_contig_chr_start,
            assembly_contig_chr_end,
            assembly_contig_orientation,
            assembly_contig_type
    FROM assembly_contig
    WHERE assembly_contig_id = $dbID};
      

  my $sth = $self->prepare($query);

  $sth->execute();


  my $hashRef;
  my $contig;

  if( $hashRef = $sth->fetchrow_hashref() ) {
    $contig = $self->_assembly_contig_from_sth( $sth, $hashRef );
  } else {
    $self->throw("Can't fetch AssemblyContig with internal id $dbID\n");
  }

  return $contig;
}

sub fetch_by_chr_start_end {
  my ($self,$chr,$chrstart,$chrend) = @_;

  defined($chr)      || $self->throw("No chr name input");
  defined($chrstart) || $self->throw("No chr start input");
  defined($chrend)   || $self->throw("No chr end input");

  my $query =  qq {
    SELECT  assembly_contig_id,
            assembly_contig_name,
            assembly_contig_start,
            assembly_contig_end,
            assembly_contig_chr_name,
            assembly_contig_chr_start,
            assembly_contig_chr_end,
            assembly_contig_orientation,
            assembly_contig_type
    FROM assembly_contig
    WHERE assembly_contig_chr_name = \'$chr\'
    AND   assembly_contig_chr_end   >= $chrstart
    AND   assembly_contig_chr_start     <= $chrend};
      

  my $sth = $self->prepare($query);

  $sth->execute();
  
  my $hashRef;
  my @contig;

  while ( $hashRef = $sth->fetchrow_hashref) {
    my $contig = $self->_assembly_contig_from_sth( $sth, $hashRef);
    push(@contig,$contig);

  }

  return @contig;

}

sub _assembly_contig_from_sth {
  my ($self,$sth,$hashref) = @_;

  my $id             = $hashref->{assembly_contig_id};
  my $name           = $hashref->{assembly_contig_name};
  my $start          = $hashref->{assembly_contig_start};
  my $end            = $hashref->{assembly_contig_end};
  my $chrname        = $hashref->{assembly_contig_chr_name};
  my $chr_start      = $hashref->{assembly_contig_chr_start};
  my $chr_end        = $hashref->{assembly_contig_chr_end};
  my $orientation    = $hashref->{assembly_contig_orientation};
  my $assembly_type  = $hashref->{assembly_contig_type};


#  print STDERR "$id:$name:$start:$end:$chrname:$chr_start:$chr_end:$orientation:$assembly_type\n";

  my $contig = new Bio::EnsEMBL::Virtual::AssemblyContig(-dbID => $id,
  -display_id => $name,
  -start => $start,
  -end => $end,
  -chr_name => $chrname,
  -chr_start => $chr_start,
  -chr_end  => $chr_end,
  -orientation => $orientation,
  -assembly_type => $assembly_type);

  return $contig;
}

=head2 store

 Title   : store
 Usage   : $assembly_contig_adaptor->store($assembly_contig);
 Function: Stores the assembly_contig
 Example : 
 Returns : nothing
  Args    : Bio::EnsEMBL::Virtual::AssemblyContig

=cut

sub store {
  my ( $self, $contig ) = @_;

  if( ! $contig->isa("Bio::EnsEMBL::Virtual::AssemblyContig")) {
    $self->throw("$contig is not a Bio::EnsEMBL::Virtual::AssemblyContig. Can't store");
  }


  my $sql = q{
    INSERT into assembly_contig (assembly_contig_id,
				 assembly_contig_name,
				 assembly_contig_chr_name,
				 assembly_contig_chr_start,
				 assembly_contig_chr_end,
				 assembly_contig_start,
				 assembly_contig_end,
				 assembly_contig_orientation,
				 assembly_contig_type)
		 VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ? )
		};
  my $sth = $self->prepare($sql);

  $sth->execute(undef,$contig->display_id,
		$contig->chr_name,
		$contig->chr_start,
		$contig->chr_end,
		$contig->start,
		$contig->end,
		$contig->orientation,
		$contig->assembly_type);
		
  my $dbID = $sth->{'mysql_insertid'};
  $contig->dbID($dbID);

}

sub has_AssemblyContigs {
  my ($self,$type) = @_;

  if (!defined($type)) {
    $self->throw("Need to specify the assembly type in has_Assembly_contigs");
  }
  
  my $query = "select * from assembly_contig where assembly_contig_type = '$type' limit 1";
  
  my $sth = $self->prepare($query);
  my $res = $sth->execute;

  if ($sth->rows > 0) {
    return 1;
  } else {
    return 0;
  }
}

1;
