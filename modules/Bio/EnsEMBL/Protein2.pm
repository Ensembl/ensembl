
#
# BioPerl module for Protein.pm
#
# Cared for by Emmanuel Mongin <mongin@ebi.ac.uk>
#
# Copyright Emmanuel Mongin
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Protein.pm - DESCRIPTION of Object

=head1 SYNOPSIS

my $protein = Bio::EnsEMBL::Protein->new ( -seq =>$sequence,
					      -accession_number  => $id,
					      -display_id => $id,
					      -primary_id => $id,
					      -id => $id,
					      -desc => $desc,
					      -moltype => $moltype,
					      );

=head1 DESCRIPTION

This object inherit from Bio::SeqI for the implementation and from PrimarySeq
A protein object hold the basic PrimarySeq requirements, hold annotation, DBlink, protein_features object.
It will also soon hold family and SNPs objects. The DB adaptor for this object is called Protein_Adaptor 

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Protein2;

use vars qw(@ISA);
use strict;
use Bio::Root::Object;
use Bio::SeqI;
use Bio::Annotation::DBLink;
use Bio::EnsEMBL::Transcript;
use Bio::DBLinkContainerI;
use Bio::SeqIO;
use Bio::Tools::SeqStats;

# Object preamble - inheriets from Bio::Root::Object

@ISA = qw(Bio::PrimarySeq Bio::SeqI);

=head2 SeqI interface features
=cut


=head2 add_date

 Title   : add_date
 Usage   : $self->add_date($ref)
 Function: adds a date
 Example :
 Returns : 
 Args    :

=cut

sub add_date {
   my ($self) = shift;
   foreach my $dt ( @_ ) {
       push(@{$self->{'date'}},$dt);
   }
}

=head2 get_dates

 Title   : get_dates
 Usage   : foreach $dt ( $self->get_dates() )
 Function: gets an array of dates
 Example :
 Returns : 
 Args    :

=cut

sub get_dates {
   my ($self) = @_;
   return @{$self->{'date'}}; 
}


sub each_date { 
   my ($self) = @_;
   return $self->get_dates; 
}

=head2 species

 Title   : species
 Usage   : $obj->species($newval)
 Function: 
 Returns : value of species
 Args    : newvalue (optional)


=cut

sub species{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'species'} = $value;
    }
    return $obj->{'species'};

}

=head2 primary_seq

 Title   : primary_seq
 Usage   : $obj->primary_seq($newval)
 Function: 
 Returns : value of primary_seq
 Args    : newvalue (optional)


=cut

sub primary_seq{
    my ($self) = @_;   
    return $self;
}




=head2 molecule

 Title   : molecule
 Usage   : $obj->molecule($newval)
 Function: 
 Returns : value of molecule
 Args    : newvalue (optional)


=cut

sub molecule{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'molecule'} = $value;
    }
    return $obj->{'molecule'};

}

=head2 moltype

 Title   : moltype
 Usage   : $obj->moltype($newval)
 Function: 
 Returns : value of moltype
 Args    : newvalue (optional)


=cut

sub moltype{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'moltype'} = $value;
    }
    return $obj->{'moltype'};

}


=head2 annotation

 Title   : annotation
 Usage   : $obj->annotation($newval)
 Function: 
 Returns : value of annotation
 Args    : newvalue (optional)


=cut

sub annotation{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'annotation'} = $value;
    }
    return $obj->{'annotation'};

}


=head2 top_SeqFeatures
    
 Title   : top_SeqFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut
    
sub top_SeqFeatures {
    my ($self,@args) = @_;
    my (@f);

    push(@f,$self->get_all_DomainFeatures());
    
    push(@f,$self->get_all_IntronFeatures());
    
    #push(@f,$self->get_all_SnpsFeatures());
    
    push(@f,$self->get_all_SigpFeatures());

    push(@f,$self->get_all_TransmembraneFeatures());
    
    push(@f,$self->get_all_CoilsFeatures());
    
    push(@f,$self->get_all_LowcomplFeatures());
    
    push(@f,$self->get_all_SuperfamilyFeatures());
    
    return @f;
}


=head2 all_SeqFeature

 Title   : all_SeqFeature
 Usage   :
 Function:This method returns the same things than top_SeqFeature, sub SeqFeatures are not currently implemented in the protein object 
 Example :
 Returns : 
 Args    :


=cut

sub all_SeqFeature{
    
    my ($self) = @_;
   
    my @seq_feature = $self->top_SeqFeatures;
    if (@seq_feature) {
	return @seq_feature;
    }

}



=head2 family

 Title   : family
 Usage   : $obj->family($newval)
 Function: This method contains the family object (Not used yet)
 Returns : value of family
 Args    : newvalue (optional)


=cut

sub family{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'family'} = $value;
    }
    return $obj->{'family'};

}

=head2 Protein Specific Features
=cut

=head2 get_all_DomainFeatures

 Title   : get_all_DomainFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_DomainFeatures{
 my ($self) = @_;
 my (@f);
    if (defined ($self->{'_domains'})) {
	return @{$self->{'_domains'}};
    }
   else {
       push(@f,$self->get_all_PrintsFeatures());
       push(@f,$self->get_all_PfamFeatures());
       push(@f,$self->get_all_PrositeFeatures());
       return @f;
    }
}

=head2 add_Domains

 Title   : add_Domains
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Domains{
    my ($self,$value) = @_;
    
    if (!defined $value) {
	$self->throw("Domain is no defined");
    }

   push(@{$self->{'_domains'}},$value); 
}


=head2 get_all_PrintsFeatures

 Title   : get_all_PrintsFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_PrintsFeatures{
    my ($self) = @_;

    if (defined ($self->{'_prints'})) {
	return @{$self->{'_prints'}};
    }
   else {
       my $proteinid = $self->id();
       my @array_features = $self->protfeat_adaptor->fetch_by_feature_and_dbID('PRINTS',$proteinid);
       foreach my $in (@array_features) {
	   $self->add_Prints($in);
       }
       return @{$self->{'_prints'}};
   }

}

=head2 add_Prints

 Title   : add_Prints
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Prints{
    my ($self,$value) = @_;
    
    if (!defined $value) {
	$self->throw("Prints is no defined");
    }

   push(@{$self->{'_prints'}},$value); 

}


=head2 get_all_PfamFeatures

 Title   : get_all_PfamFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_PfamFeatures{
    my ($self) = @_;
    if (defined ($self->{'_pfam'})) {
	return @{$self->{'_pfam'}};
    }
   else {
       my $proteinid = $self->id();
       my @array_features = $self->protfeat_adaptor->fetch_by_feature_and_dbID('Pfam',$proteinid);
       foreach my $in (@array_features) {
	   $self->add_Pfam($in);
       }
       return @{$self->{'_pfam'}};
   }

}

=head2 add_Pfam

 Title   : add_Pfam
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Pfam{
 my ($self,$value) = @_;
    
    if (!defined $value) {
	$self->throw("Pfam is no defined");
    }

   push(@{$self->{'_pfam'}},$value); 
}


=head2 get_all_PrositeFeatures

 Title   : get_all_PrositeFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_PrositeFeatures{
    my ($self) = @_;
   
    if (defined ($self->{'_prosite'})) {
	return @{$self->{'_prosite'}};
    }
   else {
       @{$self->{'_prosite'}} = undef;
       my $proteinid = $self->id();
       my @array_features = $self->protfeat_adaptor->fetch_by_feature_and_dbID('PROSITE',$proteinid);
       foreach my $in (@array_features) {
	   $self->add_Prosite($in);
       }
	   return @{$self->{'_prosite'}};
   }
}

=head2 add_Prosite

 Title   : add_Prosite
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Prosite{
    my ($self,$value) = @_;
    
    if (!defined $value) {
	$self->throw("Prosite is no defined");
    }

   push(@{$self->{'_prosite'}},$value); 

}

=head2 get_all_SigpFeatures

 Title   : get_all_SigpFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SigpFeatures{
    my ($self) = @_;

    if (defined ($self->{'_sigp'})) {
	return @{$self->{'_sigp'}};
    }
    else {
	@{$self->{'_sigp'}} = undef;
       my $proteinid = $self->id();
	my @array_features = $self->protfeat_adaptor->fetch_by_feature_and_dbID('sigp',$proteinid);
	foreach my $in (@array_features) {
	    $self->add_Sigp($in);
	}
	return @{$self->{'_prosite'}};
    }
    
}

=head2 add_Sigp

 Title   : add_Sigp
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Sigp{
    my ($self,$value) = @_;
    if (!defined $value) {
	$self->throw("Sigp is no defined");
    }
    push(@{$self->{'_sigp'}},$value); 
}


=head2 get_all_TransmembraneFeatures

 Title   : get_all_TransmembraneFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_TransmembraneFeatures{
 my ($self) = @_;

    if (defined ($self->{'_transmembrane'})) {
	return @{$self->{'_transmembrane'}};
    }
    else {
	@{$self->{'_transmembrane'}} = undef;
       my $proteinid = $self->id();
	my @array_features = $self->protfeat_adaptor->fetch_by_feature_and_dbID('transmembrane',$proteinid);
	foreach my $in (@array_features) {
	    $self->add_Transmembrane($in);
	}
	return @{$self->{'_transmembrane'}};
    }
    

}

=head2 add_Transmembrane

 Title   : add_Transmembrane
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Transmembrane{
my ($self,$value) = @_;
    if (!defined $value) {
	$self->throw("Transmembrane is no defined");
    }
    push(@{$self->{'_transmembrane'}},$value); 

}

=head2 get_all_CoilsFeatures

 Title   : get_all_CoilsFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_CoilsFeatures{
 my ($self) = @_;

    if (defined ($self->{'_coils'})) {
	return @{$self->{'_coils'}};
    }
    else {
	@{$self->{'_coils'}} = undef;
       my $proteinid = $self->id();
	my @array_features = $self->protfeat_adaptor->fetch_by_feature_and_dbID('coils',$proteinid);
	foreach my $in (@array_features) {
	    $self->add_Coils($in);
	}
	return @{$self->{'_coils'}};
    }
}

=head2 add_Coils

 Title   : add_Coils
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Coils{
my ($self,$value) = @_;
    if (!defined $value) {
	$self->throw("Coils is no defined");
    }
    push(@{$self->{'_coils'}},$value); 

}

=head2 get_all_LowcomplFeatures

 Title   : get_all_LowcomplFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_LowcomplFeatures{
 my ($self) = @_;

    if (defined ($self->{'_lowcompl'})) {
	return @{$self->{'_lowcompl'}};
    }
    else {
	@{$self->{'_lowcompl'}} = undef;
       my $proteinid = $self->id();
	my @array_features = $self->protfeat_adaptor->fetch_by_feature_and_dbID('low_complexity',$proteinid);
	foreach my $in (@array_features) {
	    $self->add_Lowcompl($in);
	}
	return @{$self->{'_lowcompl'}};
    }
}

=head2 add_LowCompl

 Title   : add_LowCompl
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Lowcompl{
my ($self,$value) = @_;
    if (!defined $value) {
	$self->throw("Lowcompl is no defined");
    }
    push(@{$self->{'_lowcompl'}},$value); 

}

=head2 get_all_SuperfamilyFeatures

 Title   : get_all_SuperfamilyFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SuperfamilyFeatures{
 my ($self) = @_;

    if (defined ($self->{'_superfamily'})) {
	return @{$self->{'_superfamily'}};
    }
    else {
	@{$self->{'_superfamily'}} = undef;
       my $proteinid = $self->id();
	my @array_features = $self->protfeat_adaptor->fetch_by_feature_and_dbID('superfamily',$proteinid);
	foreach my $in (@array_features) {
	    $self->add_Lowcompl($in);
	}
	return @{$self->{'_superfamily'}};
    }
}

=head2 add_Superfamily

 Title   : add_Superfamily
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Superfamily{
my ($self,$value) = @_;
    if (!defined $value) {
	$self->throw("Superfamily is no defined");
    }
    push(@{$self->{'_superfamily'}},$value); 

}

=head2 geneac

 Title   : geneac
 Usage   : $obj->geneac($newval)
 Function: 
 Returns : value of geneac
 Args    : newvalue (optional)


=cut

sub geneac{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'geneac'} = $value;
    }
    return $obj->{'geneac'};

}

=head2 transcriptac

 Title   : transcriptac
 Usage   : $obj->transcriptac($newval)
 Function: 
 Returns : value of transcriptac
 Args    : newvalue (optional)


=cut

sub transcriptac{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'transcriptac'} = $value;
    }
    return $obj->{'transcriptac'};

}


=head2 each_Protein_feature

 Title   : each_Protein_feature
 Usage   :my @features = $protein->each_Protein_feature
 Function:Retrieve an array of protein features (FeaturePair objects)
 Example :
 Returns : FeaturePair objects
 Args    :none


=cut

sub each_Protein_feature{
   my ($self,@args) = @_;

   if (defined ($self->{'_prot_feat'})) {
        #print STDERR "ENTERING PROTEIN FEATURE RETURN\n";
  	#foreach my $sf ( @{$self->{'_prot_feat'}} ) {
	#	print STDERR "Protein feature $sf\n";
       ##}
       return @{$self->{'_prot_feat'}};
   }

  return ();
}

=head2 add_intron

 Title   : add_intron
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_intron{
    my ($self,$value) = @_;

    if (!defined $value) {
	$self->throw("The intron is no defined");
    }

   push(@{$self->{'_intron'}},$value);   

}


=head2 get_all_IntronFeatures

 Title   : get_all_IntronFeatures
 Usage   :my @introns_feature = $protein->get_all_IntronFeatures($proteinid)
 Function:Get all of the introns as Protein_FeaturePair for a given peptide and add them to protein features
 Example :
 Returns : Nothing
 Args    :Peptide ID


=cut
 
sub get_all_IntronFeatures{
   my ($self) = @_;
   if (defined ($self->{'_intron'})) {
       return @{$self->{'_intron'}};
   }
   else {
       my $proteinid = $self->id();
       my @array_introns = $self->adaptor->get_Introns($proteinid);
       foreach my $in (@array_introns) {
	   $self->add_intron($in);
       }
       return @{$self->{'_intron'}};
   }
}

=head2 add_snps

 Title   : add_snps
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_snps{
    my ($self,$value) = @_;

    if (!defined $value) {
	$self->throw("The snp is no defined");
    }

   push(@{$self->{'_snp'}},$value);  

}


=head2 get_all_SnpsFeatures

 Title   : get_all_SnpsFeatures
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_all_SnpsFeatures{
   my ($self) = @_;
   if (defined ($self->{'_snp'})) {
       return @{$self->{'_snp'}};
   }
   else {
       my @snps_array = $self->adaptor->get_snps($self);
       foreach my $sn (@snps_array) {
	   $self->add_snps($sn);
       }
       return @{$self->{'_snp'}};
   }

}


=head2 add_Protein_feature

 Title   : add_Protein_feature
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_Protein_feature{
   my ($self,$value) = @_;

   if ((!defined $value) || (!$value->isa('Bio::EnsEMBL::Protein_FeaturePair'))) {
     $self->throw("The Protein Feature added is not defined or is not a protein feature object");
   }

   push(@{$self->{'_prot_feat'}},$value);   

}


=head2 length

 Title   : length
 Usage   :my $lenght = $protein->length
 Function:
 Example :
 Returns :length of the aa sequence 
 Args    :none


=cut

sub length{
   my ($self) = @_;
   my $length = length($self->seq);


   return $length;

}

=head2 start

 Title   : start
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub start{
   my ($self) = @_;
   return 0;

}

=head2 end

 Title   : end
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub end{
   my ($self) = @_;
   my $length = $self->length;
   return $length;

}

=head2 molecular_weight

 Title   : molecular_weight
 Usage   : 
 Function: This method has been placed here for convenience function. This gets and return the molecular weight of the protein
 Example : my $mw = $protein->molecular_weight
 Returns : Moleculat weight of the peptide
 Args    :


=cut

sub molecular_weight{
   my ($self) = @_;
   my $mw = ${Bio::Tools::SeqStats->get_mol_wt($self)}[0];
return $mw
}

=head2 checksum

 Title   : checksum
 Usage   : my $checksum = $protein->checksum()
 Function: Get the crc64 for the sequence of the protein (this method has been copied from _crc64 method which is in swiss.pm in Bioperl and have been placed here for convenience function
 Example : 
 Returns : CRC64
 Args    : Nothing


=cut

sub checksum{
   my ($self) = @_;

   my $sequence = \$self->seq();

   my $POLY64REVh = 0xd8000000;
    my @CRCTableh = 256;
    my @CRCTablel = 256;
    my $initialized;       
    

    my $seq = $$sequence;
      
    my $crcl = 0;
    my $crch = 0;
    if (!$initialized) {
	$initialized = 1;
	for (my $i=0; $i<256; $i++) {
	    my $partl = $i;
	    my $parth = 0;
	    for (my $j=0; $j<8; $j++) {
		my $rflag = $partl & 1;
		$partl >>= 1;
		$partl |= (1 << 31) if $parth & 1;
		$parth >>= 1;
		$parth ^= $POLY64REVh if $rflag;
	    }
	    $CRCTableh[$i] = $parth;
	    $CRCTablel[$i] = $partl;
	}
    }
    
    foreach (split '', $seq) {
	my $shr = ($crch & 0xFF) << 24;
	my $temp1h = $crch >> 8;
	my $temp1l = ($crcl >> 8) | $shr;
	my $tableindex = ($crcl ^ (unpack "C", $_)) & 0xFF;
	$crch = $temp1h ^ $CRCTableh[$tableindex];
	$crcl = $temp1l ^ $CRCTablel[$tableindex];
    }
    my $crc64 = sprintf("%08X%08X", $crch, $crcl);
        
    return $crc64;

}

=head2 get_family

 Title   : get_family
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_family{
   my ($self) = @_;
   my $proteinid = $self->id();
   my $family = $self->adaptor->fetch_Family_by_dbid($proteinid);
   return $family;
}


=head2 adaptor

 Title   : adaptor
 Usage   : $obj->adaptor($newval)
 Function: 
 Returns : value of adaptor
 Args    : newvalue (optional)


=cut

sub adaptor{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;

      $obj->{'adaptor'} = $value;
    }
   return $obj->{'adaptor'};

}

=head2 protfeat_adaptor

 Title   : protfeat_adaptor
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub protfeat_adaptor{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;

      $obj->{'protfeat_adaptor'} = $value;
    }
    return $obj->{'protfeat_adaptor'};
}


1;









