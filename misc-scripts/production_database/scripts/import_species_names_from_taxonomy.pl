#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw( :config no_ignore_case );
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;

sub run {
  my ($class) = @_;
  my $self = bless({}, $class);
  $self->args();
  $self->defaults();
  my $taxons = $self->_get_taxons();
  foreach my $taxon (@{$taxons}) {
    $self->process_taxon($taxon);
  }
  return;
}

sub args {
  my ($self) = @_;
  my $opts = {};
  my @cmd_opts = qw/
    mhost|mh=s
    mport|mP=i
    muser|mu=s
    mpass|mp=s
    mdatabase|md=s
    host|h=s
    port|P=i
    user|u=s
    pass|p=s
    database|d=s
    taxon=i@
    verbose|v!
    help
    man
  /;
  GetOptions($opts, @cmd_opts) or pod2usage(-verbose => 1, -exitval => 1);
  pod2usage(-verbose => 1, -exitval => 0) if $opts->{help};
  pod2usage(-verbose => 2, -exitval => 0) if $opts->{man};
  $self->{opts} = $opts;
  return;
}

sub defaults {
  my ($self) = @_;
  my $o = $self->{opts};

  # Master database location:
  $o->{mhost} = 'ens-staging1' unless $o->{mhost};
  $o->{mport} = 3306 unless $o->{mport};
  
  # User database location (default values):
  $o->{port} = 3306 unless $o->{port};
  
  foreach my $required (qw/host user muser database/) {
    my $msg = "Required parameter --${required} was not given";
    pod2usage(-msg => $msg, -verbose => 1, -exitval => 1) if ! $o->{$required};
  }
  
  return;
}

sub _get_taxons {
  my ($self) = @_;
  my $taxons;
  if(defined $self->{opts}->{taxon} && @{$self->{opts}->{taxon}}) {
    $taxons = $self->{opts}->{taxon};
  }
  else {
    $taxons = $self->_production_dbc()->sql_helper()->execute_simple(-SQL => 'select taxon from species where taxon is not null');
  }
  $self->v('Found %d taxon(s) to process', scalar(@{$taxons}));
  return $taxons;
}

sub process_taxon {
  my ($self, $taxon) = @_;
  $self->v('Processing %s', $taxon);
  $self->v('Getting results from taxonomy DB');
  my $taxonomy = $self->_get_taxonomy($taxon);
  $self->v('Getting results from production DB');
  my $production = $self->_get_production($taxon);
  if(! $production->{species_id}) {
    $self->v("ERROR: %d taxon ID does not exist in the production database", $taxon);
    return;
  }
  $self->v('Updating attributes');
  $self->_update_attributes($taxonomy, $production);
  $self->v('Updating species aliases');
  $self->_update_species_aliases($taxonomy, $production);
  $self->v('Done');
  return;
}

sub _get_taxonomy {
  my ($self, $taxon) = @_;
  my $dbc = $self->_taxonomy_dbc();
  my %aliases;
  my $query_hash = $dbc->sql_helper()->execute_into_hash(
    -SQL => 'select name_class, name from ncbi_taxa_name where taxon_id =?',
    -PARAMS => [$taxon],
    -CALLBACK => sub {
      my ($row, $value) = @_;
      my ($class, $name) = @{$row};
      
      return {} if $class eq 'authority';
      return {} if $class eq 'merged_taxon_id';
      
      $aliases{lc($name)} = 1; 
      if(defined $value) {
        push(@{$value}, $name);
      }
      return [ $name ];
    }
  );
  
  my %return_hash;
  $return_hash{active_aliases} = [sort keys %aliases];
  $return_hash{inactive_aliases} = [];
  
  $return_hash{scientific_name} = $query_hash->{'scientific name'}->[0] if $query_hash->{'scientific name'};
  $return_hash{common_name} = $query_hash->{'genbank common name'}->[0] if $query_hash->{'genbank common name'};
  $return_hash{web_name} = $query_hash->{'common name'}->[0] if $query_hash->{'common name'};
  
  return \%return_hash;
}

sub _get_production {
  my ($self, $taxon_id) = @_;
  my $dbc = $self->_production_dbc();
  my $hash;
  $dbc->sql_helper()->execute_no_return(
    -SQL => 'select species_id, common_name, web_name, scientific_name, taxon from species where taxon =?',
    -PARAMS => [$taxon_id],
    -CALLBACK => sub {
      my ($r) = @_;
      my $i = 0;
      $hash = { 
        species_id => $r->[$i++], 
        common_name => $r->[$i++], 
        web_name => $r->[$i++], 
        scientific_name => $r->[$i++], 
        taxon => $r->[$i++]
      };
    }
  );
  my $sql = 'select alias from species_alias where species_id =? and is_current = ?';
  $hash->{active_aliases} = $dbc->sql_helper()->execute_simple(-SQL => $sql, -PARAMS => [$hash->{species_id}, 1]);
  $hash->{inactive_aliases} = $dbc->sql_helper()->execute_simple(-SQL => $sql, -PARAMS => [$hash->{species_id}, 1]);
  return $hash;
}

sub _decipher_aliases {
  my ($self, $source, $production) = @_;
  
  my %source_active_hash = map { $_ => 1 } @{$source->{active_aliases}};
  my %prod_active_hash = map { $_ => 1 } @{$production->{active_aliases}};
  my %prod_inactive_hash = map { $_ => 1 } @{$production->{inactive_aliases}};
  
  my @new_aliases;
  my @reactivated_aliases;
  my @retired_aliases;
  foreach my $source (keys %source_active_hash) {
    if(! $prod_active_hash{$source}) {
      if(! $prod_inactive_hash{$source}) {
        push(@new_aliases, $source);
      }
      else {
        push(@reactivated_aliases, $source);
      }
    }
  }
  foreach my $target (keys %prod_active_hash) {
    if(! $source_active_hash{$target}) {
      push(@retired_aliases, $target);
    }
  }
  
  return { new => \@new_aliases, retired => \@retired_aliases, reactivated => \@reactivated_aliases };
}

sub _decipher_attributes {
  my ($self, $source, $target) = @_;
  my %final;
  my @keys = qw/common_name web_name scientific_name/;
  foreach my $key (@keys) {
    if($target->{$key}) {
      $final{$key} = $target->{$key};
    }
    else {
      $final{$key} = $source->{$key};
    }
  }
  return \%final;
}

sub _update_attributes {
  my ($self, $source, $production) = @_;
  my $dbc = $self->_production_dbc();
  my $h = $dbc->sql_helper();
  
  my $attributes = $self->_decipher_attributes($source, $production);
    
  my @keys = keys %{$attributes};
  my $set_string = join(q{,}, map { "$_ =?" } @keys );
  my @params = map { $attributes->{$_} } @keys;
  push(@params, $production->{species_id});
  my $sql = sprintf('update species set %s where species_id =?', $set_string);
  $h->execute_update(-SQL => $sql, -PARAMS => \@params);
  
  return;
}

sub _update_species_aliases {
  my ($self, $source, $production) = @_;
  my $dbc = $self->_production_dbc();
  my $h = $dbc->sql_helper();

  my $aliases = $self->_decipher_aliases($source, $production);
  my $species_id = $production->{species_id};
  
  foreach my $alias (@{$aliases->{new}}) {
    my $sql = 'insert into species_alias (species_id, alias, is_current, created_by, created_at) values (?,?,?,?,now())';
    $h->execute_update(-SQL => $sql, -PARAMS => [$species_id, $alias, 1, undef]);
  }
  foreach my $alias (@{$aliases->{reactivated}}) {
    my $sql = 'update species_alias set is_current =?, modified_by =?, modified_at =now() where species_id =? and alias =?';
    $h->execute_update(-SQL => $sql, -PARAMS => [1, undef, $species_id, $alias]);
  }
  foreach my $alias (@{$aliases->{retired}}) {
    my $sql = 'update species_alias set is_current =?, modified_by =?, modified_at =now() where species_id =? and alias =?';
    $h->execute_update(-SQL => $sql, -PARAMS => [0, undef, $species_id, $alias]);
  }
  
  return;  
}

sub _taxonomy_dbc {
  my ($self) = @_;
  my $o = $self->{opts};
  my %args = ( -HOST => $o->{host}, -PORT => $o->{port}, -USER => $o->{user}, -DBNAME => $o->{database} );
  $args{-PASS} = $o->{pass} if $o->{pass};
  return Bio::EnsEMBL::DBSQL::DBConnection->new(%args);
}

sub _production_dbc {
  my ($self) = @_;
  my $o = $self->{opts};
  my %args = ( -HOST => $o->{mhost}, -PORT => $o->{mport}, -DBNAME => $o->{mdatabase}, -USER => $o->{muser} );
  $args{-PASS} = $o->{mpass} if $o->{mpass};
  return Bio::EnsEMBL::DBSQL::DBConnection->new(%args);
}

sub v {
  my ($self, $msg, @args) = @_;
  return unless $self->{opts}->{verbose};
  my $s_msg = sprintf($msg, @args);
  my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) =
    localtime(time());
  print sprintf("[%02d-%02d-%04d %02d:%02d:%02d] %s\n",
                $mday, $mon, $year + 1900,
                $hour, $min, $sec, $s_msg);
  return;
}

__PACKAGE__->run();

__END__

=pod 

=head1 NAME

import_species_names_from_taxonomy.pl

=head1 SYNOPSIS

  ./import_species_names_from_taxonomy.pl -h host [-P port] \\
    -u user [-p password]
    -d database \\
    [-mh host] [-mP port] \\
    [-mu user] [-mp password] [-md database] \\
    [-taxon id]
    [-v]

=head1 DESCRIPTION

TODO

=head1 OPTIONS

=over 8

=item B<-h|--host>

Host of the taxonomy database

=item B<-P|--port>

Port of the taxonomy database

=item B<-u|--user>

Username for the taxonomy database

=item B<-p|--pass>

Password for the taxonomy database

=item B<-mh|--mhost>

Host for the production database

=item B<-mP|--mport>

Port for the production database

=item B<-mu|--muser>

User for the production database

=item B<-mp|--mpass>

Pass for the production database

=item B<--taxon>

Allows for the update of a single species by specifying its taxonomy id. This
is useful for new species updates.

=item B<--verbose>

Make the script chatty

=item B<--help>

Help message

=item B<--man>

Man page

=back

=cut