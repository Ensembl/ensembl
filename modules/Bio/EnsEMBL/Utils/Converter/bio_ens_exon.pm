# Bio::EnsEMBL::Utils::Converter::bio_ens_exon
#
# Created and cared for by Juguang Xiao <juguang@tll.org.sg>
# Created date: 27/3/2003
# 
# Copyright Juguang Xiao
# 
# You may distribute this module under the same terms as perl itself
#
# POD documentation
#

=head1 NAME

Bio::EnsEMBL::Utils::Converter::bio_ens_exon

=head1 SYNOPISIS



=head1 DESCRIPTION


=head1 FEEDBACK

=head2 Mailing Lists

=head2 Reporting Bugs


=head1 AUTHOR Juguang Xiao

Juguang Xiao <juguang@tll.org.sg>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin ...

package Bio::EnsEMBL::Utils::Converter::bio_ens_exon;

use strict;
use vars qw(@ISA);
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Utils::Converter;
@ISA = qw(Bio::EnsEMBL::Utils::Converter);

=head2 convert
  Title   : convert
  Usage   : 
  Function: 
  Return  :
  Args    :
=cut

sub convert {
    my ($self, $arg) = @_;
    $self->throw("needs an array ref") unless(ref($arg) eq 'ARRAY');
    
    my $converter = Bio::EnsEMBL::Utils::Converter(
        -in => 'Bio::SeqFeature::Generic',
        -out => 'Bio::EnsEMBL::Exon',
        -contig => $self->contig
    );

    my $ens_exons = $converter->convert($arg);
    return $ens_exons;
}


sub _convert_single {
    my ($self, @args) = @_;

}


1;
