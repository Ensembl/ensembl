# EnsEMBL module for Bio::EnsEMBL::DBSQL::Foreign
#
# Cared for Philip lijnzaad@ebi.ac.uk
#
# Copyright EnsEMBL
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBSQL::ObjForeign


=head1 SYNOPSIS

  use Bio::EnsEMBL::DBSQL::Obj;
  use Bio::EnsEMBL::DBSQL::Feature_ObjForeign;

  $db = new Bio::EnsEMBL::DBSQL::Obj( -user => 'root'
                                      , -db => 'pog' 
                                      , -host => 'caldy' 
                                      , -driver => 'mysql' );

  my $translations = { 'static_golden_path' => 'my_golden_path' , 
                       'exon' => 'ens075.exon'};
  my $feature_obj =
   Bio::EnsEMBL::Feature_ObjForeign->new(
                                         -dbobj => $db,
                                         -table_name_translations => $translations);

  # rest as with Feature_Obj. 

=head1 DESCRIPTION

This class is like a utility class for mixing in to other Obj classes
(e.g. Gene_Obj or Feature_Obj) resulting in an adaptor that translate
table names to other names on the fly. The intended use is for an database
that reuses part of another database in order to save on speed (no loading
of all the assembly related data) and disk usage.

This class is typically to be sub-classed.

=head1 CONTACT

Philip lijnzaad@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal
methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::DBSQL::ObjForeign;

use vars qw(@ISA);
use strict;


sub new {
  my ($class, @args) = @_;
  my $self= {};
  bless $self, $class;

  my ($db_obj, $table_name_translations, $read_db_obj) 
    = $self->_rearrange([ qw( DBOBJ TABLE_NAME_TRANSLATIONS READ_DBOBJ) ]
                        ,@args);
  $db_obj || $self->throw("I need a db obj ...");
  $self->_db_obj($db_obj);
  $read_db_obj = $db_obj unless $db_obj;
  $self->_read_db_obj($read_db_obj);
  # table_name_translations are optional  
  $self->_table_name_translations($table_name_translations);
  $self->use_delayed_insert(1);
  return $self; # success - we hope!
}

sub _read_db_obj {
    my ($self, $value) = @_;

    if( defined $value) {
        my $needed = 'Bio::EnsEMBL::DBSQL::Obj';
        if ( ref($value) ne $needed ) {
            $self->throw("expecting a $needed");
        }
      $self->{'_read_db_obj'} = $value;
    }
    return $self->{'_read_db_obj'};
}



=head2 _table_name_tranlations

 Title   : table_name_tranlations
 Usage   :
 Function: sets/gets the hash used for looking up table names
 Example : $self->_table_name_tranlations 
                     { 'static_golden_path' => 'my_golden_path' , # same db
                       'exon' => 'ens075.exon' # different table
                   }
 Returns : 
 Args    :


=cut

sub _table_name_translations {
   my ($self,$value) = @_;
   if( defined $value) {
       if ( ref($value) ne 'HASH') {
           $self->throw('expecting hash ref');
       }
      $self->{'_table_name_translations'} = $value;
    }
    return $self->{'_table_name_translations'};
}

=head2 _lookup_table_name

 Title   :  _lookup_table_name
 Usage   :  
 Function:  translate name to other name
 Example :
 Returns : new name if there is a translation for this table name, 
           the old name otherwise. 
 Args    : the table name.

=cut

sub _lookup_table_name {
   my ($self,$name) = @_;
   
   my $table = $self->_table_name_translations;
   my $newname =  $table->{$name};
   if (defined $table && defined $newname  ) {
       return $newname ;
   }
   return $name;
}

1;
