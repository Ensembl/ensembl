#
# BioPerl module for DNA parsing and converting into sql
#
# Cared for by Michele Clamp <michele@sanger.ac.uk>
#
# Copyright Michele Clamp
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

DNA2sql - parses fasta and inserts dna object into ensembl

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Reads a fasta file and inserts a new DNA object into
the ensembl database.

Creation:
   
     my $dna  = new Bio::EnsEMBL::Analysis::DNAparse           # Empty object
     my $dna  = new Bio::EnsEMBL::Analysis::DNAparse($seqio);  # $seqio is Bio::SeqIO
     my $dna  = new Bio::EnsEMBL::Analysis::DNAparse($seq);    # $seqio is Bio::Seq

Manipulation:

     my $sql    = $dna->sql;                # Returns sql needed for entering into the database
     my @sql    = $dna->sql_next            # Returns the sql one sequence at a time
     my $result = $dna->insert($dbh)        # Inserts sql into the database whose handle is $dbh
     my $result = $dna->insert_next($dbh)   # Inserts sql into the database whose handle is $dbh
                                            # a sequence at a time
     
=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Analysis::DNAparse;

use     vars qw(@ISA);
use     strict;

# Object preamble - inherits from Bio::Root::Object

use Bio::EnsEMBL::DB::Obj;
use Bio::Root::Object;
use Bio::SeqIO;

use DBI;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my ($self,@args) = @_;

  my $make = $self->SUPER::_initialize;


  # set stuff in self from @args
  if (defined(@args)) {
    if (ref($args[0]) =~ /Bio::SeqIO/) {
      $self->seqIO($args[0]);

    } elsif ($args[0]->isa("Bio::Seq")) {
      $self->seq($args[0]);

    } else {
      $self->throw("Input argument must be Bio::SeqIO\n");
    }
  } 
  return $self; # success - we hope!
}
      

=pod

=head2 seqIO

  Title   : seqIO
  Usage   : $dna->seqIO($seqio)
  Function: Sets the input stream for the dna parsing
  Returns : nothing
  Args    : Bio::SeqIO

=cut

sub seqIO {
  my ($self,$value) = @_;
  
  if (defined($value)) {
    $self->{'_seqIO'} = $value;
  }
  return $self->{'_seqIO'};
}
=pod

=head2 seq

  Title   : seq
  Usage   : $dna->seq($seq)
  Function: Sets the dna sequence to parse
  Returns : nothing
  Args    : Bio::Seq

=cut

sub seq {
  my ($self,$value) = @_;
  
  if (defined($value)) {
    $self->{'_seq'} = $value;
  }
  return $self->{'_seq'};
}

=pod

=head2 sql

  Title   : sql
  Usage   : $dna->sql
  Function: Returns sql for the ensembl DNA table
  Returns : String
  Args    : none

=cut

sub sql {
  my ($self) = @_;
  
  my @outstr;

  while (my @sql = $self->sql_next) {
    push(@outstr,@sql);
  }
  return @outstr;
}

=pod

=head2 sql_next

  Title   : sql_next
  Usage   : $dna->sql_next
  Function: Returns sql for the ensembl DNA table one sequence at a time
  Returns : String,Array of Strings
  Args    : none

=cut

sub sql_next {
  my ($self,$seq,$contig) = @_;

  if ($seq) {

    my $contig_id =  $seq->id;

    my $clone     = $contig->{clone};

    if (!$clone) {
      $clone = $contig_id;
      $clone =~ s/^(.*)\..*/$1/;
    }

    my $mapbin    = $contig->{mapbin};

    my $seqstr    = $seq->seq;
    my $length    = length($seqstr);

    my $date      = `date '+%Y-%m-%d'`; chomp $date;
    my @str; 

    push(@str, "insert into dna(contig,sequence,created) values('$contig_id','$seqstr','$date')");
    push(@str, "replace into contig(id,dna,length,clone,mapbin) values('$contig_id',LAST_INSERT_ID(),$length,'$clone','$mapbin')");

    return (@str);
  } else {
    return;
  }
}

=pod

=head2 get_seq

  Title   : get_seq
  Usage   : $dna->get_seq
  Function: Returns the next sequence object
  Returns : Bio::Seq
  Args    : none

=cut

sub get_seq {
  my ($self) = @_;
  
  if (defined($self->{_seq})) {
    my $tmp = $self->{_seq};
    delete $self->{_seq};
    return $tmp;

  } 

  if (defined($self->{_seqIO})) {
       return $self->{_seqIO}->next_seq;
  }
}

=pod

=head2 insert

  Title   : insert
  Usage   : $dna->insert($dbh)
  Function: Inserts all dna into database whose handle is $dbh
  Returns : nothing
  Args    : $dbh = DBI->connect('DBI:mysql:pog','root','');

=cut

sub insert {
  my ($self,$dbh) = @_;

  if (!$self->locked) {
    $self->lock($dbh);
  }

  my @sql = $self->sql;
  my $rv;
  foreach my $line (@sql) {
    my $sth = $dbh->prepare($line);
       $rv  = $sth->execute;
  }

  $self->unlock($dbh);

  return $rv;
}

=head2 insert_next

  Title   : insert_next
  Usage   : $dna->insert_next($dbh)
  Function: Inserts one sequence worth of dna into database whose handle is $dbh
  Returns : nothing
  Args    : $dbh = DBI->connect('DBI:mysql:pog','root','');

=cut

sub insert_next {
  my ($self,$dbh) = @_;

  my $seq    = $self->get_seq;

  if ($seq) {
    my $contig = $self->get_contig($dbh,$seq->id);
    my (@sql)  = $self->sql_next  ($seq,$contig);
    
    if (!$self->locked) {
      $self->lock($dbh);
    }
    
    if ($#sql >= 0) {
      my $rv = 0;
      foreach my $line (@sql) {
	my $sth = $dbh->prepare($line);
  	   $rv  = $sth->execute;
      }
      return $seq->id if $rv;
    }
  }
}

=head2 locked 

  Title   : locked
  Usage   : $f = $dna->locked
  Function: Returns whether the contig and dna tables are locked
  Returns : 0,1
  Args    : none

=cut

sub locked {
  my ($self) = @_;
  
  if (!defined($self->{_LOCKED})) {
     $self->{_LOCKED} = 0;
  }
  return $self->{_LOCKED};

}

=head2  lock

  Title   : lock 
  Usage   : $dna->lock($dbh)
  Function: Locks dna and contig tables so we can update them
  Returns : 1,0
  Args    : $dbh = DBI->connect('DBI:mysql:pog','root','');

=cut

sub lock {
  my ($self,$dbh) = @_;

  my $sth = $dbh->prepare('lock tables dna write,contig write');
  my $rv  = $sth->execute;
  
  if ($rv) {
    $self->{_LOCKED} = 1;
  } else {
    $self->{_LOCKED} = 0;
  }
  return $rv;
}

=head2  unlock

  Title   : unlock
  Usage   : $dna->unlock($dbh)
  Function: unlocks tables
  Returns : 1,0
  Args    : $dbh = DBI->connect('DBI:mysql:pog','root','');                     
                                                           
=cut                                                       
    
sub unlock {
  my ($self,$dbh) = @_;
                       
  my $sth = $dbh->prepare('unlock tables');
  my $rv  = $sth->execute;

  return $rv;
}                     


=head2  get_contig

  Title   : get_contig
  Usage   : $dna->get_contig($dbh,$id);
  Function: Gets contig
  Returns : hashref
  Args    : $dbh = DBI->connect('DBI:mysql:pog','root','') , String
                                                           
=cut                                                       
    
sub get_contig {
  my ($self,$dbh,$id) = @_;
                       
  my $sth = $dbh->prepare('select * from contig where id=\'$id\'');
  my $rv  = $sth->execute;

  my $contig = {};

  while(my $rowhash = $sth->fetchrow_hashref) {
    $contig->{id}     = $rowhash->{id};
    $contig->{clone}  = $rowhash->{clone};
    $contig->{mapbin} = $rowhash->{mapbin};
    $contig->{length} = $rowhash->{length};
    $contig->{dna}    = $rowhash->{dna};
  }
  return $contig;
}                         
    
   




1;
