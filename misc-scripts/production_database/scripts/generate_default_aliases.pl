#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw( :config no_ignore_case );
use Pod::Usage;
use POSIX;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Utils::Exception qw(throw);

sub run {
  my ($class) = @_;
  my $self = bless( {}, $class );
  $self->args();
  $self->check_opts();
  my $species = $self->_species();
  foreach my $s ( sort { $a->{production} cmp $b->{production} } values %{$species} ) {
    $self->v( 'Processing %s (species_id: %d)', $s->{production}, $s->{id} );
    my $aliases_to_add = $self->_aliases_to_add($s);
    $self->_write_aliases($aliases_to_add, $s);
    $self->v('Done');
  }
  return;
}

sub args {
  my ($self) = @_;

  my $opts = {

    # Master database location:
    mhost     => 'ens-staging1',
    mport     => 3306,
    mdatabase => 'ensembl_production',
    
    species   => [],
    write     => 0,
    dieonconflict => 1
  };

  my @cmd_opts = qw/
    mhost|mh=s
    mport|mP=i
    muser|mu=s
    mpass|mp=s
    mdatabase|md=s
    species|s=s@
    verbose|v!
    write
    help
    man
    /;
  GetOptions( $opts, @cmd_opts ) or pod2usage( -verbose => 1, -exitval => 1 );
  pod2usage( -verbose => 1, -exitval => 0 ) if $opts->{help};
  pod2usage( -verbose => 2, -exitval => 0 ) if $opts->{man};
  $self->{opts} = $opts;
  return;
}

sub check_opts {
  my ($self) = @_;
  my $o = $self->{opts};

  foreach my $required (qw/mhost muser/) {
    my $msg = "Required parameter --${required} was not given";
    pod2usage( -msg => $msg, -verbose => 1, -exitval => 1 ) if !$o->{$required};
  }
  
  if(! @{$self->{opts}->{species}}){
    my $msg = "Required parameter --species was not given";
    pod2usage( -msg => $msg, -verbose => 1, -exitval => 1 );
  }

  return;
}

sub _write_aliases {
  my ($self, $aliases, $species) = @_;
  my $dbc = $self->_production_dbc();
  $dbc->sql_helper()->transaction(sub {
    my $sql = 'insert ignore into species_alias (species_id, alias, is_current, created_at) values (?,?,1, NOW())';
    my $id = $species->{id};
    $dbc->sql_helper()->batch(-SQL => $sql, -CALLBACK => sub {
      my ($sth) = @_;
      foreach my $a (@{$aliases}) {
        if($self->{opts}->{write}) {
          $sth->execute($id, $a);
        }
        else {
          $self->v('Would have inserted the alias %s for species_id %d', $a, $id);
        }
      }
      return;
    });
  });
  return;
}

sub _species {
  my ($self) = @_;
  my $dbc   = $self->_production_dbc();
  my $h     = $dbc->sql_helper();
  my $sql   = <<'SQL';
select species_id, common_name, web_name, scientific_name, production_name, url_name
from species
where production_name like ?
and is_current = 1
SQL
  my %species;
  foreach my $species (@{$self->{opts}->{species}}) {
    $self->v('Querying production for current species like %s', $species);
    $dbc->sql_helper()->execute_no_return(-SQL => $sql, -PARAMS => [$species], -CALLBACK => sub {
      my ($row) = @_;
      my ($id, $common_name, $web_name, $scientific_name, $production_name, $url_name) = @{$row};
      if (!$common_name) {
        throw("no common name specified, exiting");
      }
      if(!exists $species{$id}) {
        $species{$id} = {
          id => $id,
          production => $production_name, 
          common => $common_name, 
          web => $web_name, 
          scientific => $scientific_name, 
          url => $url_name
        };
      }
      return;
    });
  }
  
  #Enrich after executing
  $self->_enrich(\%species);
  
  return \%species;
}

sub _enrich {
  my ($self, $species) = @_;
  foreach my $id (keys %{$species}) {
    my $s = $species->{$id};
    $s->{automatic_aliases} = $self->_automatic_aliases($s); 
  }
  return;
}

sub _aliases {
  my ($self) = @_;
  if(! exists $self->{_aliases}) {
    my $dbc   = $self->_production_dbc();
    my $aliases = $dbc->sql_helper()->execute_simple(
      -SQL => 'select alias from species_alias where is_current = 1'
    );
    $self->{_aliases} = {map { $_ => 1 } @{$aliases}};
  }
  return $self->{_aliases};
}

sub _automatic_aliases {
  my ($self, $species) = @_;
  my $production_name = $species->{production};
  my $common_name = lc($species->{common});
  my $web_name = lc($species->{web});
    
  my $automatic_aliases = {};
  $automatic_aliases->{$common_name} = 1;
  $automatic_aliases->{$web_name} = 1;

  # *** Assume homo_sapiens ***
  my $alias = $production_name;
  
  #1). homo_sapiens
  $automatic_aliases->{$alias} = 1;

  #2). homo sapiens
  $alias =~ tr [_] [ ];
  $automatic_aliases->{$alias} = 1;

  #3). hsapiens
  $production_name =~ /^(.)[^_]*_(.*)$/;
  $alias = $1 . $2;
  $automatic_aliases->{$alias} = 1;

  #4). hsap
  $production_name =~ /^(.)[^_]*_(...).*$/;
  $alias = $1 . $2;
  $automatic_aliases->{$alias} = 1;

  #5). homosap
  $production_name =~ /^(...)[^_]*_(...).*$/;
  $alias = $1 . $2;
  $automatic_aliases->{$alias} = 1;
  
  return $automatic_aliases;
}

sub _aliases_to_add {
  my ($self, $species) = @_;
  my @aliases_to_add;
  my $aliases = $self->_aliases();
  foreach my $autogenerated (keys %{$species->{automatic_aliases}}) {
    if(exists $species->{aliases}->{$autogenerated}) {
      if($species->{aliases}->{$autogenerated} == 1) {
        $self->v('SKIPPING WARNING: %s as it already registered for another species', $autogenerated);
      }
      elsif($species->{aliases}->{$autogenerated} == 2) {
        $self->v('SKIPPING WARNING: %s as it has been generated for another species in this session', $autogenerated);
      }
    }
    else {
      push(@aliases_to_add, $autogenerated);
      $aliases->{$autogenerated} = 2;
      $self->v('%s is a new alias', $autogenerated);
    }
  }
  return \@aliases_to_add;
}

sub _production_dbc {
  my ($self) = @_;
  my $o      = $self->{opts};
  my %args   = (
    -HOST   => $o->{mhost},
    -PORT   => $o->{mport},
    -DBNAME => $o->{mdatabase},
    -USER   => $o->{muser}
  );
  $args{-PASS} = $o->{mpass} if $o->{mpass};
  return Bio::EnsEMBL::DBSQL::DBConnection->new(%args);
}

sub v {
  my ( $self, $msg, @args ) = @_;
  return unless $self->{opts}->{verbose};
  my $s_msg = sprintf( $msg, @args );
  my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) =
    localtime( time() );
  print sprintf(
    "[%02d-%02d-%04d %02d:%02d:%02d] %s\n",
    $mday, $mon, $year + 1900,
    $hour, $min, $sec, $s_msg
  );
  return;
}

__PACKAGE__->run();

__END__

=pod 

=head1 NAME

generate_default_aliases.pl

=head1 SYNOPSIS

  ./generate_default_aliases.pl 
    -mh host -mp password -mu user [-mP port] \\
    [-md database] \\
    [-th host] [-tP port] \\
    [-tu user] [-tp password] [-td database] \\
    [-species] [-write] \\
    [-v]

=head1 DESCRIPTION

A script used to generate a minimal set of required aliases. Assuming the
production_name I<homo_sapiens> we would generate the following

=over 8

=item B<homo_sapiens>

=item B<homo sapiens>

=item B<hsapiens>

=item B<hsap>

=item B<homsap>

=back

It is up to the user to add more via the admin interface. We do not remove
aliases with this script

=head1 OPTIONS

=over 8

=item B<-mh|--mhost>

Host for the production database

=item B<-mP|--mport>

Port for the production database

=item B<-mu|--muser>

User for the production database

=item B<-mp|--mpass>

Pass for the production database

=item B<-md|--mdatabase>

Name for the production database.

=item B<-s|--species>

Species to generate the names for. Can use a SQL pattern here and multiple
cmd line entries. Please use B<production_names>.

=item B<--write>

Turns writing on. Not on by default 

=item B<--verbose>

Make the script chatty

=item B<--help>

Help message

=item B<--man>

Man page

=back

=cut
