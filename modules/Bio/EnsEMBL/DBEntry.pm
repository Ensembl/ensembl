#
# EnsEMBL module for DBEntry
#
# Cared for by Arne Stabenau <stabenau@ebi.ac.uk>
#
# Copyright EMBL/EBI 2001
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::DBEntry - Module to collect information about an
 external reference

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 CONTACT

  Arne Stabenau <stabenau@ebi.ac.uk>
  Ewan Birney <birney@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBEntry;
use vars qw(@ISA $AUTOLOAD);
use strict;


@ISA = qw( Bio::Root::RootI Bio::Annotation::DBLink );


sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ( $adaptor, $dbID, $primary_id, $version,
       $dbname, $release  ) = $self->_rearrange
	 ( [ qw { ADAPTOR
		DBID
		PRIMARY_ID
		VERSION
		DBNAME
		RELEASE
	      }], @args );

  if( defined $adaptor ) { $self->adaptor( $adaptor )}
  if( defined $dbID ) { $self->dbID( $dbID ) }
  if( defined $primary_id ) { $self->primary_id( $primary_id ) }
  if( defined $version ) { $self->version( $version ) }
  if( defined $dbname ) { $self->dbname( $dbname ) }
  if( defined $release) { $self->release( $release ) }
  $self->{_synonyms} = [];;

  return $self;
}


sub primary_id {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_primary_id} = $arg;
  } 
  return $self->{_primary_id};
}

sub dbname {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_dbname} = $arg;
  } 
  return $self->{_dbname};
}


sub database {
  my $self = shift;
  return $self->dbname();
}

sub optional_id {
  my $self = shift;
  return "";
}


sub release {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_release} = $arg;
  } 
  return $self->{_release};
}


sub adaptor {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_adaptor} = $arg;
  } 
  return $self->{_adaptor};
}


sub dbID {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_dbID} = $arg;
  } 
  return $self->{_dbID};
}


sub version {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_version} = $arg;
  } 
  return $self->{_version};
}



sub description {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_description} = $arg;
  } 
  return $self->{_description};
}


sub add_synonym {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    push( @{$self->{_synonyms}}, $arg );
  }
}

# get a list of synonym for this db reference
sub get_synonyms {
  my $self = shift;
  return @{$self->{_synonyms}};
}

sub flush_synonyms {
  my $self = shift;
  $self->{_synonyms} = [];
}


# the urlpattern contains a space where the primary id 
# of the refernce can be put. This will lead you to a page
# hopefully giving you more information on the object
# referenced
sub url_pattern {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_urlPattern} = $arg;
  } 
  return $self->{_urlPattern};
}



# returns a url for the object
# you can set it, but stored is only a pattern
sub url {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_url} = $arg;
    return;
  }

  if( defined $self->{_url} ) {
    return $self->{_url};
  }

  my $id = $self->primary_id();
  my $url = $self->urlPattern();
  $url =~ s/ /$id/;
  return $url;
}

sub dump {
  my $self = shift;
  my ($k,$v);

  while( ($k,$v) = each %$self ) {
    print $k," ",$v,"\n";
  }
}

#Cheat to comply with bioperl
sub comment {
    my ($self) = @_;
    if ($self) {
	return $self->description();
    }
}

1;
