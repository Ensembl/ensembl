

#
# Module for drawing contig gif
#

=head1 NAME

SeqContigDraw

=head1 SYNOPSIS

Contains subroutines for drawing contigs

=head1 DESCRIPTION

Needs more description.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

=cut



package SeqContigDraw;
use strict;
use GD;
use Parameters;




=head2 draw_contig_image

 Title   : draw_contig_image
 Usage   : &SeqContigDraw::draw_contig_image($im,$contig)
 Function: draws an image of a contig with all sequence features
 Example :
 Returns : 
 Args    :


=cut


sub draw_contig_image
{    
    my ($im,$contig)=@_;
    my $image_comp_ref=&Parameters::contig_image_components;

    my $gd_col_ref=&Parameters::colors;
    my $bcg_color=$gd_col_ref->{white};

    # print legend   
    foreach my $type(keys %{$image_comp_ref})
    {	
	my @legend_args=($im,$type);
	&print_legend (@legend_args);	    	
    }

    # draw sequence
    &draw_seq($im);

    # draw genes
    my $seq_len=$contig->seq->seq_len;
    my @genes = $contig->get_all_Genes();    
    foreach my $gene(@genes){&draw_gene($im,$gene,$seq_len,$contig->id);}
   
    # draw sequence features
    my @features= $contig->get_all_SeqFeatures();
    foreach my $ft (@features)
    {
	my $type;
        
	if ($ft->analysis && $ft->analysis->db ne ""){$type=$ft->analysis->db;}
	else{$type=$ft->analysis->gff_source;}

	my $color=$gd_col_ref->{$image_comp_ref->{$type}{color}};	

        if ($ft->sub_SeqFeature) {
           foreach my $f ($ft->sub_SeqFeature) {
	      &draw_feature($im,$f,$seq_len,$type,$color);	   	
           }
        } else {
	   &draw_feature($im,$ft,$seq_len,$type,$color);	   	
        }
    }
     
    # draw scale
    my @args=($im,$contig);
    &draw_scale (@args);
}




=head2 print_legend

 Title   : print_legend
 Usage   : &SeqContigDraw::print_legend($im,$type)
 Function: prints a legend for a contig image
 Example :
 Returns : 
 Args    :


=cut




sub print_legend
{

    my ($im,$type)=@_;
    my $gd_col_ref=&Parameters::colors;    
    my $image_param_ref=&Parameters::contig_image_par;
    my $image_comp_ref=&Parameters::contig_image_components;
    my $x_start=$image_param_ref->{legend_margin}; 
    my $color=$gd_col_ref->{$image_comp_ref->{scale}{color}};
    my $name=$image_comp_ref->{$type}{legend};
    my $width=0;

    my @strands=(1,-1);
    foreach my $strand (@strands){
	my @arg=($type,$width,$strand);
	my @y_coord=&calc_y_coord(@arg);
	my $y_start=(shift @y_coord)-8;
	$im->string(gdSmallFont,$x_start,$y_start,$name,$color);
    }
    
}



=head2 draw_seq

 Title   : draw_seq
 Usage   : &SeqContigDraw::draw_seq($im)
 Function: draws a sequence
 Example :
 Returns : 
 Args    :


=cut




sub draw_seq
{
    my ($im)=@_;
    
    my $type="sequence";
    my $image_param_ref=&Parameters::contig_image_par;    
    my $image_comp_ref=&Parameters::contig_image_components;
    my $gd_col_ref=&Parameters::colors;
    my $width=1;
    my $strand=1;
    my $color=$gd_col_ref->{$image_comp_ref->{sequence}{color}};
    my @arg=($type,$width,$strand);
    my ($y_start,$y_end)=&calc_y_coord(@arg);
    my $x_start=$image_param_ref->{left_margin};
    my $x_end=$image_param_ref->{x_img_len}-$image_param_ref->{right_margin};
    
    $im->filledRectangle($x_start,$y_start,$x_end,$y_end,$color);        
}




=head2 draw_gene

 Title   : draws_gene
 Usage   : &SeqContigDraw::draw_gene($im,$gene,$seq_len,$contig_id)
 Function: draws a gene 
 Example :
 Returns : 
 Args    :


=cut



sub draw_gene
{
    my ($im,$gene,$seq_len,$contig_id)=@_;
    
    my @types=('exon','intron');
    
    my $image_param_ref;
    my $image_comp_ref;

    # if it is a gene image
    if (!defined $seq_len){ $image_param_ref=&Parameters::gene_image_par;$image_comp_ref=&Parameters::gene_image_components;}

    # if it is a contig image
    else{ $image_param_ref=&Parameters::contig_image_par;$image_comp_ref=&Parameters::contig_image_components;}   
      
    my $gd_col_ref=&Parameters::colors;  	    
    my $new_gene_status=1;        
    my $keep_coord=0;
    my $len;
    my $fixed;  
    my $substr;

  
    my @exons = sort { $a->start <=> $b->start } $gene->each_unique_Exon();
    foreach my $exon ( @exons )
    {	
	
	if( defined $seq_len ){$len=$seq_len;}
	    else {$len=&gene_length($gene);$fixed=1;}
	if (!defined $seq_len && $new_gene_status==1){$substr=$exon->start;}
		
	my $width=1;
	my @x_args;
	
	foreach my $type(@types){

	    # calculate x coord 
	    if ($type eq 'exon'){@x_args=($len,$exon->start-$substr,$exon->end-$substr,$exon->end-$exon->start,$fixed)};
	    if ($type eq 'intron'){@x_args=($len,$keep_coord-$substr,$exon->start-$substr,$exon->start-$keep_coord,$fixed);}
	    my  ($x_start,$x_end)=&calc_x_coord(@x_args);
	 
	    # calculate y coord
	    if (defined $seq_len){$image_param_ref=undef;}
	    my @arg=($type,$width,$exon->strand,$image_param_ref);
	    my ($y_start,$y_end)=&calc_y_coord(@arg);

	    # start drawing from first exon 
	    if 	($type eq 'exon' || $new_gene_status !=1){  	
	   
		# do not draw unless exons and introns are on the contig (contig image) or it is a gene image
		if (!defined $seq_len || $exon->contig_id eq $contig_id){
		$im->filledRectangle($x_start,$y_start,$x_end,$y_end,$gd_col_ref->{$image_comp_ref->{$type}{color}});
	    }	    
		# do not print map for a gene image
	    if (defined $seq_len){&print_map($x_start,$y_start,$x_end,$y_end,$gene->id,$image_comp_ref->{$type}{link},$gene->id);}
	    }
	}
	# do not reset variables unless exons and introns are on the contig (contig image) or it is a gene image
	if (!defined $seq_len || $exon->contig_id eq $contig_id){	
	    $new_gene_status=0;
	    $keep_coord=$exon->end;  
	}	    
    }
    # draw a scale bar for a gene image
    if (! defined $seq_len){ &draw_scale_bar($im,$len);}    

}




=head2 draw_feature

 Title   : draws_feature
 Usage   : &SeqContigDraw::draw_feature($im,$ft,$seq_len,$type,$color)
 Function: draws a sequence feature
 Example :
 Returns : 
 Args    :


=cut


sub draw_feature
{
 
    my ($im,$ft,$seq_len,$type,$color)=@_;

    my $image_comp_ref=&Parameters::contig_image_components;  
    my $url;
    if ($ft->analysis->db ne ""){ $url=$image_comp_ref->{$ft->analysis->db}{link};}
				  else {$url=$image_comp_ref->{$ft->analysis->gff_source}{link};}

    my $y_start;
    my $y_end;
    my $strand;
    my ($x_start,$x_end)=&calc_x_coord($seq_len,$ft->start,$ft->end,$ft->length);  
    my $width=1;
    my @arg=($type,$width,$ft->strand);
    ($y_start,$y_end)=&calc_y_coord(@arg);

    my $name;
    my $db;
    if ($ft->isa("Bio::EnsEMBL::FeaturePair"))
    {
	if ($ft->analysis->db eq "swir"){($db,$name)=split /:/,$ft->hseqname;}
	if ($ft->analysis->db eq "dbest"){($db,$name)=split /\|/,$ft->hseqname;}
	else {$name=$ft->hseqname;}
    }
    if ($ft->isa("Bio::EnsEMBL::SeqFeature")){$name=$ft->analysis->gff_source}

    unless ($url eq ""){&print_map($x_start,$y_start,$x_end,$y_end,$name,$url,$name);}
    $im->filledRectangle($x_start,$y_start,$x_end,$y_end,$color);    

}





=head2 draw_scale

 Title   : draw_scale
 Usage   : &SeqContigDraw::draw_scale($im,$contig)
 Function: draws an image scale
 Example :
 Returns : 
 Args    :


=cut



sub draw_scale
{

    my ($im,$contig)=@_;
    
    my $type="scale";
    my $image_param_ref=&Parameters::contig_image_par;
    my $image_comp_ref=&Parameters::contig_image_components;     
    my $gd_col_ref=&Parameters::colors;
    my $color=$gd_col_ref->{$image_comp_ref->{scale}{color}};
    my $x_start=$image_param_ref->{left_margin};
    my $x_end=$image_param_ref->{x_img_len}-$image_param_ref->{right_margin};       
    my $width=0;
     
    my @strands=(1,-1);
    foreach my $strand(@strands){
	my $x_start=$image_param_ref->{left_margin};
	my $x_end=$image_param_ref->{x_img_len}-$image_param_ref->{right_margin};
		
	my @arg=($type,$width,$strand);
	my ($y_start,$y_end)=&calc_y_coord(@arg);
	$im->line($x_start,$y_start,$x_end,$y_end,$color);    
			
	my $seq_len=$contig->seq->seq_len;
	my $x_img=$image_param_ref->{x_img_len};
	my $times;
	my $scale_len=$image_comp_ref->{$type}{length};

	
	$times= $seq_len/$scale_len;
	
	# calculate the no of bp coordinates that can fit on a picture
	until ($times<$image_comp_ref->{$type}{step})
	{$scale_len=$image_comp_ref->{$type}{times}*$scale_len;  $times= $seq_len/$scale_len;}
	
	my $pixel_bases=($x_img-$image_param_ref->{left_margin}-$image_param_ref->{right_margin})/$seq_len;
	
	# draw scale
	my $y_start_new;
	if ($strand==1){$y_start_new=$y_start-15;}
	if ($strand==-1){$y_start_new=$y_start+5;}
	my $i;   
	for ($i=0;$i<$times-1;$i++)
	{
	    $x_start=$image_param_ref->{left_margin}+($i*$pixel_bases*$scale_len);
	    $x_end=$x_start;
	    my $bases=$i*$scale_len;
	    my $print_start=$x_start-15;
	    if ($i==0){$print_start=$x_start-2;$bases=1;}
	    
	    $im->string(gdSmallFont,$print_start,$y_start_new,$bases,$color);
	    $im->line($x_start,$y_start-2,$x_end,$y_end+2,$color);	    
	}
	
	$x_start=$x_end=$image_param_ref->{x_img_len}-$image_param_ref->{right_margin};
	
	my $print_start=$x_start-15;
	$im->string(gdSmallFont,$print_start,$y_start_new,$seq_len,$color);
	$im->line($x_start,$y_start-2,$x_end,$y_end+2,$color);
		
    }
    
}



=head2 draw_scale-bar

 Title   : draws_scale_bar
 Usage   : &SeqContigDraw::draw_scale_bar($im,$len)
 Function: draws a scale bar on a gene image
 Example :
 Returns : 
 Args    :


=cut




sub draw_scale_bar
{

    my ($im,$len)=@_;

    my $type='scale_bar';
    my $image_param_ref=&Parameters::gene_image_par;
    my $image_comp_ref=&Parameters::gene_image_components;    
    my $gd_col_ref=&Parameters::colors;
    my $scale_bar=$image_comp_ref->{$type}{length};
    if ($len <=1000){$scale_bar=100;}
    my $fixed=1;
    my ($scale_start,$scale_end)=&calc_x_coord($len,0,$scale_bar,$scale_bar,$fixed);
    
    $im->filledRectangle($scale_start,$image_param_ref->{y_scale},$scale_end,$image_param_ref->{y_scale}+$image_param_ref->{scale_width},$gd_col_ref->{$image_comp_ref->{$type}{color}});
    $im->string(gdSmallFont,$scale_end+8,$image_param_ref->{y_scale}-5,"$scale_bar bp",$gd_col_ref->{$image_comp_ref->{$type}{legend_color}});

}


=head2 gene_length

 Title   : gene_length
 Usage   : &SeqContigDraw::gene_length($gene)
 Function: returns gene length
 Example :
 Returns : gene length
 Args    :


=cut



sub gene_length
{
    
    my ($gene) = @_;
    
    my @exons = sort { $a->start <=> $b->start } $gene->each_unique_Exon();    
    my $f= shift @exons;       
    my $first=$f->start; 
    my $l=pop @exons;
    my $last=$l->end;
    my $gene_len=$last-$first;
    
    return $gene_len;    
}


=head2 calc_x_coord

 Title   : calc_x_coord
 Usage   : &SeqContigDraw::calc_x_coord($seq_len,$start,$end,$length,$fixed)
 Function: calculates x coordinates
 Example :
 Returns : start and end coordinates
 Args    :


=cut




sub calc_x_coord

{
    my ($seq_len,$start,$end,$length,$fixed)=@_;

    my $image_param_ref=&Parameters::contig_image_par;
    if ($fixed==1){$image_param_ref=&Parameters::gene_image_par;}
    my $x_start;
    my $x_end;

    my	$img_len=$image_param_ref->{x_img_len}-$image_param_ref->{left_margin}-$image_param_ref->{right_margin};
    my	$pixel_bases=$img_len/$seq_len;	    
    $x_start=$start*$pixel_bases+$image_param_ref->{left_margin};   
    my $feature_len=$length*$pixel_bases;
    $x_end=$x_start+$feature_len;

    return my @x_coord=($x_start,$x_end);
}



=head2 calc_y_coord

 Title   : calc_y_coord
 Usage   : &SeqContigDraw::calc_y_coord($type,$width,$strand,$image_param_ref)
 Function: calculates y coordinates
 Example :
 Returns : start and end coordinates
 Args    :


=cut



sub calc_y_coord
    
{    
    my ($type,$width,$strand,$image_param_ref)= @_;
    
    my $image_comp_ref;
    if (! defined $image_param_ref){$image_param_ref=&Parameters::contig_image_par;$image_comp_ref=&Parameters::contig_image_components;}
    else{$image_comp_ref=&Parameters::gene_image_components;}    
    my $gap=$image_param_ref->{gap};
    my $factor=$image_comp_ref->{$type}{factor}; 
    my $y_start;
    my $y_end;
      
    $y_start=$image_param_ref->{y_cent} -$strand*$gap*$factor;
    $y_end=$y_start;
   
    my $height;
    unless ($width==0){
	$y_start=$y_start-0.5*$image_comp_ref->{$type}{height};
	$y_end=$y_start+$image_comp_ref->{$type}{height};
    }
 
    return my @y_coord=($y_start,$y_end);    
}


 
=head2 print_map

 Title   : print map
 Usage   : &SeqContigDraw::print_map($type,$width,$strand,$image_param_ref)
 Function: prints image map
 Example :
 Returns : 
 Args    :


=cut



sub print_map
{
    my ($x_start,$y_start,$x_end,$y_end,$name,$url,$alt_text)=@_;
    
    print  "<AREA Shape=\"Rect\" coords = \"",int $x_start," ",int $y_start,
    " ",int $x_end," ",int $y_end," \"  HREF=\"",$url,$name,"\"","alt=\"",$alt_text,"\">","\n";
    
}





1;





