
#
# Config file for SeqContigDraw.pm
#


=head1 NAME

Parameters

=head1 SYNOPSIS

Contains parameters for SeqContigDraw

=head1 DESCRIPTION

Needs more description.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

=cut







package Parameters;

use strict;
use GD;





=head2 contig_image_components

 Title   : contig_image_components
 Usage   : $hashref=&Parameters::contig_image_components
 Function: hashref with image components
 Example :
 Returns : hashref  with image components
 Args    :


=cut






sub contig_image_components
    
{
    my $components={
	sequence    =>{color=>'red',       factor=>0,height=>1,legend=>'contig',},
	exon        =>{color=>'black',     factor=>1,height=>9,legend=>'Ensembl gene',link=>'/perl/geneview?gene=',},
	intron      =>{color=>'black',     factor=>1,height=>1,legend=>'Ensembl gene',link=>'/perl/geneview?gene=',},
	swir        =>{color=>'darkorange',factor=>3,height=>6,legend=>'swissprot',   link=>'http://www.ebi.ac.uk/cgi-bin/swissfetch?',},
	dbest       =>{color=>'red',       factor=>5,height=>6,legend=>'EST',       link=>"http://www.ebi.ac.uk/cgi-bin/emblfetch?",},
	vert        =>{color=>'pine',      factor=>6,height=>6,legend=>'EMBL',        link=>'http://www.ebi.ac.uk/cgi-bin/emblfetch?',},
	PfamFrag    =>{color=>'flora',     factor=>4,height=>6,legend=>'pfam',        link=>'/perl/getpfamacc?pfamentry=',},
	RepeatMasker=>{color=>'darkgreen', factor=>7,height=>6,legend=>'repeat',},
	genscan     =>{color=>'blue2',     factor=>2,height=>6,legend=>'genscan',},
	scale       =>{color=>'black',     factor=>8,height=>0,legend=>'scale (bp)',length=>1000,step=>10,times=>5,},		       
    };   
    
    return $components;       
}


=head2 gene_image_components

 Title   : gene_image_components
 Usage   : $hashref=&Parameters::gene_image_components
 Function: hashref with image components
 Example :
 Returns : hashref  with image components
 Args    :


=cut




sub gene_image_components
    
{
    my $components={	
	exon     =>{color=>'black',factor=>1,height=>30,legend=>'Ensembl gene',link=>'/perl/geneview?gene=',legend=>'Ensembl gene',},
	intron   =>{color=>'black',factor=>1,height=>2, legend=>'Ensembl gene',link=>'/perl/geneview?gene=',legend=>'Ensembl gene',},	

	scale_bar=>{color=>'blue' ,legend_color=>'black',height=>1,descr=>'bp',length=>1000,},
    };   
    
    return $components;    
}




=head2 contig_image_par

 Title   : contig_image_par
 Usage   : $hashref=&Parameters::contig_image_par
 Function: hashref with image parameters
 Example :
 Returns : hashref  with image parameters
 Args    :


=cut




sub contig_image_par
{
    
    my $image_param_ref={};
    
    $image_param_ref->{x_img_len}     = 600;
    $image_param_ref->{y_img_len}     = 220;
    $image_param_ref->{left_margin}   = 120;
    $image_param_ref->{right_margin}  = 30;
    $image_param_ref->{top_margin}    = 0;
    $image_param_ref->{bottom_margin} = 0;
    $image_param_ref->{legend_margin} = 10;
    $image_param_ref->{y_cent}        = 0.5 * $image_param_ref->{y_img_len};
    $image_param_ref->{gap}           = 11;
       
    return $image_param_ref;
}



=head2 gene_image_par

 Title   : contig_image_par
 Usage   : $hashref=&Parameters::contig_image_par
 Function: hashref with image parameters
 Example :
 Returns : hashref  with image parameters
 Args    :


=cut





sub gene_image_par
{
        
    my $image_param_ref={};
    
    $image_param_ref->{x_img_len}     = 550;
    $image_param_ref->{y_img_len}     = 100;
    $image_param_ref->{left_margin}   = 30;
    $image_param_ref->{right_margin}  = 30;
    $image_param_ref->{top_margin}    = 0;
    $image_param_ref->{bottom_margin} = 0;
    $image_param_ref->{y_cent}        = 0.5 * $image_param_ref->{y_img_len};
    $image_param_ref->{y_scale}       = $image_param_ref->{y_cent}+40;
    

    return $image_param_ref;
}



=head2 colors

 Title   : colors
 Usage   : $hashref=&Parameters::colors($im)
 Function: hashref with GD colors
 Example :
 Returns : hashref  with GD colors
 Args    :


=cut



sub colors
{

    my ($im)=@_;
    
    my $colors_ref={};
    
    $colors_ref->{white}      = $im->colorAllocate(255,255,255);
    $colors_ref->{grey}       = $im->colorAllocate(192,192,192);
    $colors_ref->{lightgrey}  = $im->colorAllocate(140,140,140);
    $colors_ref->{darkgrey}   = $im->colorAllocate(80,80,80);  
    $colors_ref->{black}      = $im->colorAllocate(0,0,0);         
    $colors_ref->{red}        = $im->colorAllocate(238,0,0);
    $colors_ref->{green}      = $im->colorAllocate(0,255,0);
    $colors_ref->{blue}       = $im->colorAllocate(0,0,255);
    $colors_ref->{blue2}      = $im->colorAllocate(0,4,182);
    $colors_ref->{darkgreen}  = $im->colorAllocate(50,80,50);
    $colors_ref->{chocolate}  = $im->colorAllocate(210,105,30);
    $colors_ref->{brown}      = $im->colorAllocate(165,42,42);
    $colors_ref->{blueviolet} = $im->colorAllocate(138,43,226);
    $colors_ref->{darkorange} = $im->colorAllocate(255,127,0);
    $colors_ref->{mustard}    = $im->colorAllocate(205,207,33);
    $colors_ref->{violet3}    = $im->colorAllocate(226,226,255);
    $colors_ref->{pine}       = $im->colorAllocate(0,150,28);
    $colors_ref->{flora}      = $im->colorAllocate(202,209,13);
    
    return $colors_ref;
}


1;



