package Parameters;

use strict;
use GD;


sub contig_image_components
    
{
    my $components={
	sequence=>{color=>'red',factor=>0,height=>1,legend=>'contig'},
	exon=>{color=>'black',factor=>1,height=>6,link=>'http://www.ensembl.org/cgi-bin/geneview.pl?gene=',legend=>'Ensembl gene',},
	intron=>{color=>'black',factor=>1,height=>1,link=>'http://www.ensembl.org/cgi-bin/geneview.pl?gene=',legend=>'Ensembl gene',},
	RepeatMasker=>{color=>'darkgreen',factor=>7,height=>6,legend=>'repeat'},
	swir=>{color=>'darkorange',factor=>3,height=>6,link=>'http://www.ebi.ac.uk/cgi-bin/swissfetch?',legend=>'swissprot',},
	PfamFrag=>{color=>'black',factor=>4,height=>6,legend=>'pfam',
		   link=>'http://www.sanger.ac.uk/cgi-bin/Pfam/querypfam.pl?loc=http%3A%2F%2Fwww.sanger.ac.uk%2Fcgi-bin%2FPfam&db=pfam&db=prosite&terms=',},
	genscan=>{color=>'blue',factor=>2,height=>6,legend=>'genscan',},
	dbest=>{color=>'red',factor=>5,height=>6,legend=>'dbEST',
		link=>"http://www.ncbi.nlm.nih.gov/htbin-post/Entrez/query?form=6&dopt=g&db=n&uid=",},
	vert=>{color=>'black',factor=>6,height=>6,link=>'http://www.ebi.ac.uk/cgi-bin/emblfetch?',legend=>'vert',},
	scale=>{color=>'black',factor=>8, height=>0,legend=>'scale (bp)',length=>1000,step=>10,times=>5,},		       
		   };   
    
    return $components;       
}


sub gene_image_components
    
{
    my $components={	
	exon=>{color=>'black',factor=>1,height=>30,link=>'http://www.ensembl.org/cgi-bin/geneview.pl?gene=',legend=>'Ensembl gene',},
	intron=>{color=>'black',factor=>1,height=>2,link=>'http://www.ensembl.org/cgi-bin/geneview.pl?gene=',legend=>'Ensembl gene',},	
	scale_bar=>{color=>'blue',legend_color=>'black',height=>1,descr=>'bp',length=>1000,},
    };   
    
    return $components;    
}



sub contig_image_par
{
    
    my $image_param_ref={};
    
    $image_param_ref->{x_img_len} = 600;
    $image_param_ref->{y_img_len} = 220;
    $image_param_ref->{left_margin} = 120;
    $image_param_ref->{right_margin} = 30;
    $image_param_ref->{top_margin} = 0;
    $image_param_ref->{bottom_margin} =0;
    $image_param_ref->{legend_margin} =10;
    $image_param_ref->{y_cent} =0.5*$image_param_ref->{y_img_len};
    $image_param_ref->{gap} =11;
       
    return $image_param_ref;
}



sub gene_image_par
{
        
    my $image_param_ref={};
    
    $image_param_ref->{x_img_len} = 550;
    $image_param_ref->{y_img_len} = 100;
    $image_param_ref->{left_margin} = 30;
    $image_param_ref->{right_margin} = 30;
    $image_param_ref->{top_margin} = 0;
    $image_param_ref->{bottom_margin} =0;
    $image_param_ref->{y_cent} =0.5*$image_param_ref->{y_img_len};
    $image_param_ref->{y_scale} =$image_param_ref->{y_cent}+40;
    

    return $image_param_ref;
}



sub colors
{

    my ($im)=@_;
    
    my $colors_ref={};
    
    $colors_ref->{white} = $im->colorAllocate(255,255,255);
    $colors_ref->{grey} = $im->colorAllocate(192,192,192);
    $colors_ref->{lightgrey}=$im->colorAllocate(140,140,140);
    $colors_ref->{darkgrey}= $im->colorAllocate(80,80,80);  
    $colors_ref->{black} = $im->colorAllocate(0,0,0);         
    $colors_ref->{red} = $im->colorAllocate(238,0,0);
    $colors_ref->{green}= $im->colorAllocate(0,255,0);
    $colors_ref->{blue} = $im->colorAllocate(0,0,255);
    $colors_ref->{blue2} = $im->colorAllocate(0,4,182);
    $colors_ref->{darkgreen} = $im->colorAllocate(50,80,50);
    $colors_ref->{chocolate} = $im->colorAllocate(210,105,30);
    $colors_ref->{brown} = $im->colorAllocate(165,42,42);
    $colors_ref->{blueviolet} = $im->colorAllocate(138,43,226);
    $colors_ref->{darkorange} = $im->colorAllocate(255,127,0);
    $colors_ref->{mustard} = $im->colorAllocate(205,207,33);
    $colors_ref->{violet3} = $im->colorAllocate(226,226,255);
    $colors_ref->{pine} = $im->colorAllocate(0,150,28);
    $colors_ref->{flora} = $im->colorAllocate(202,209,13);
    
    return $colors_ref;
}


1;



