package GenDraw;

use strict;
use GD;
#use Bio::EnsEMBL::DBSQL::Obj;
use ImageMap;
use DrawingRoutines;

 

sub get_contig_image_parameters
{  

    my($contig_obj)=@_;
 
 my $y_img_len=65;
 my $x_img_len=500;
 my $margin=30;
 my $im = new GD::Image($x_img_len+2*$margin,$y_img_len);
 my $GAP=25;
 
 my @pm = &GenDraw::get_seq_len($x_img_len,$contig_obj);
 my $s_len=@pm[0];
 my $fixed_bases_per_pixel=$s_len/($x_img_len-1);

 my @param=($im,$GAP,$x_img_len,$margin,$fixed_bases_per_pixel);

    return @param;

# return my $arr_point=\@param;

}



sub get_gene_image_parameters
{
    my($bases_per_pixel)=@_;
    my $y_img_len=80;
    my $x_img_len=500;
    my $margin=30;
    my $im = new GD::Image($x_img_len+2*$margin,$y_img_len);
    my @param=($im,$y_img_len,$x_img_len,$margin);
    return my $arr_point=\@param; 
}




sub build_image
{
    my($bases_per_pixel,$contig_obj)=@_;
    
    my $GAP=25;
    my $x_img_len=500;
    my @pm = &GenDraw::scale_image($bases_per_pixel,$x_img_len,$contig_obj);
    my $repeats=@pm[0];
    my $y_img_len=$repeats*$GAP;

#    my $fixed_y_img_len=65;
 #   $y_img_len=$fixed_y_img_len;

    my $margin=30;
    my $im = new GD::Image($x_img_len+2*$margin,$y_img_len);
   # my $s_len=@pm[1];
  #  my $fixed_bases_per_pixel=$s_len/($x_img_len-1);

    return my @param= ($im,$GAP,$x_img_len,$margin,$bases_per_pixel);
}



sub scale_image 
{
    my ($bases_per_pixel,$x_img_len,$contig_obj)=@_;
    
    my $seq = $contig_obj->seq();
    my $repeats =$seq->seq_len/($x_img_len*$bases_per_pixel)+2;    

    # fixed length adjustment
    my $s_len=$seq->seq_len;
    return my @param=($repeats,$s_len);
}



sub get_seq_len 
{
    my ($x_img_len,$contig_obj)=@_;
    
    my $seq = $contig_obj->seq();
    #my $repeats =$seq->seq_len/($x_img_len*$bases_per_pixel)+2;    

    # fixed length adjustment
    my $s_len=$seq->seq_len;
    return my @param=($s_len);
}





sub draw_genes 
{
    my ($contig_obj,$arr_point,$gene_array,$img_map_fh,$gif_url) = @_;

# get image parameters
 
    my $im_fh=${$arr_point}[0];
my $GAP= ${$arr_point}[1];
my $x_img_len= ${$arr_point}[2];
my $margin= ${$arr_point}[3];
my $bases_per_pixel=${$arr_point}[4];    

# set colours
   my $white = $im_fh->colorAllocate(255,255,255);
    my $grey = $im_fh->colorAllocate(192,192,192);# background
	my $lightgrey=$im_fh->colorAllocate(140,140,140);
    my $darkgrey= $im_fh->colorAllocate(80,80,80);
    my $white = $im_fh->colorAllocate(255,255,255);
    my $black = $im_fh->colorAllocate(0,0,0);
    my $red = $im_fh->colorAllocate(255,0,0);
    my $blue = $im_fh->colorAllocate(0,0,255);
 my $green= $im_fh->colorAllocate(0,255,0);

# set the scale bar (bp)
    my $scale_bar =1000;


# set the the vertical distance between each drawn object and the sequence
    my $center= ($GAP)/4;
    my $gap=0.1*$GAP;
    my $adjusted_gap=0;
    my $y_cent_seq=$center;
    my $n=4;
    my $y_cent_exon=$center+$n*$gap;
    my $y_cent_baseno=$center-12*$gap;
    my $y_cent_contig_name=-6*$gap;
    my $y_scale_bar=12*$gap;

# set width of each object
    my $seq_width=0.04*$GAP;
    my $exon_width=0.3*$GAP;
    my $gene_width=0.03*$GAP;
    my $scale_bar_width=0.04*$GAP;
    my $y_start_seq=$y_cent_seq-0.5*$seq_width;
    my $y_start_exon=$y_cent_exon-0.5*$exon_width;
    my $y_start_baseno=$y_cent_baseno;
    my $y_end_seq=$y_start_seq+$seq_width;
    my $y_end_exon=$y_start_exon+$exon_width;
    my $y_end_scale_bar=$y_scale_bar+$scale_bar_width;    
              
# bases count parameters
    my $no_of_bases=0;
    my $last_no_of_bases=0;
    my $first_no_of_bases=0;
    my $remaining_bases=0;

# x coordinate parameters
    my $x_start=$margin;
    my $x_end=$margin;
    my $x_exon_start=0;
    my $x_exon_end=0;
    my $keep_gene_start=0;
    my $keep_gene_end=0;    
    
# drawing parameters
    my $gene_status=0;
    my $drawing_status=1;
    my $drawing_counter=0;  
    my $factor=1/$bases_per_pixel;
  





#
# contig object drawer	
#

    # start of an image map    
    &ImageMap::print_map_start($img_map_fh,$gif_url);

    # create contig object and get seq
    my $contig = $contig_obj;
    my $seq = $contig->seq();    
    $remaining_bases=$seq->seq_len;

    # print contig id
#    $im_fh->string(gdMediumBoldFont,$x_img_len/2-50 ,20,$contig->id,$black);
   

    # adjust drawing parameters according to scale   
    while ($drawing_status==1)
    {
	if ($remaining_bases>=$x_img_len*$bases_per_pixel)
	{    
	    $no_of_bases = $x_img_len*$bases_per_pixel; 
	    $x_end=$x_img_len+$margin;
	    $drawing_counter++;
	    $first_no_of_bases=$last_no_of_bases+1;
	    $last_no_of_bases=$no_of_bases +$last_no_of_bases;
	    $remaining_bases=$seq->seq_len-$last_no_of_bases;		 	    
	}
	else
	{
	    $drawing_status=0;	 
	    $x_end=$margin+$remaining_bases*$factor;
	    $no_of_bases = $remaining_bases;
	    $drawing_counter++;
	    $first_no_of_bases=$last_no_of_bases+1;
	    $last_no_of_bases=$no_of_bases +$last_no_of_bases;     
	}
	
	
	#
	# loop through all genes/transcripts/exons
	#

	
	
	foreach my $gene ( @{$gene_array} ) {
	    my $seen = 0;
	    foreach my $c ( $gene->unique_contig_ids() ) {
		
		if( $c eq $contig->id() ) {
                  $seen = 1;
                  last;
                 }
            }
            if( $seen == 0 ) { next; }
	    $gene_status=1;
	    foreach my $trans ( $gene->each_Transcript() ) {
		foreach my $exon ( $trans->each_Exon() ) {
		    
		    #
		    # draw exons
		    #

		    # adjust the gap according to strand
		    if($exon->strand == 1 && $adjusted_gap == 0){$adjusted_gap=-2*$n*$gap;}
		    if($exon->strand == -1 && $adjusted_gap<0){$adjusted_gap=0;}

		    # calculate x exon coordinates
		    $x_exon_start=$margin+($exon->start)*$factor-($drawing_counter-1)*$x_img_len;
		    $x_exon_end=$margin+($exon->end)*$factor-($drawing_counter-1)*$x_img_len;
		    
		    # don't go beyond the margin
		    if ($x_exon_start <$margin){$x_exon_start=$margin;}
		    if ($x_exon_end >$x_img_len+$margin){$x_exon_end=$x_img_len+$margin;}
		    
		    # calculate y exon coordinates
		    my    $y_ex_start= ($y_start_exon+$adjusted_gap)+$drawing_counter*$GAP;
		    my    $y_ex_end= ($y_end_exon+$adjusted_gap)+$drawing_counter*$GAP;
		    
		    # and finally draw exons
		    &DrawingRoutines::draw_a_rectangle($im_fh,$x_exon_start,$y_ex_start,$x_exon_end,$y_ex_end,$black);


		    # prepare for drawing a gene ...		    
		    # find the coordinates of the first and the last gene exon
		    
		    # if it's a first exon of a gene 
		    if ($gene_status==1){
			$keep_gene_start=$x_exon_start;
			$keep_gene_end=$x_exon_end;
			$gene_status=2;
		    }
		    
		    # if it's a subsequent exon
		    if ($gene_status==2){
			if($keep_gene_end<$x_exon_end){$keep_gene_end=$x_exon_end;}
			if($keep_gene_start>$x_exon_end){$keep_gene_start=$x_exon_start;}
		    }
		}
		
		# calculate gene image map coordinates
		my	$y_ex_img_start= ($y_start_exon+$adjusted_gap)+$drawing_counter*$GAP;
		my	$y_ex_img_end= ($y_end_exon+$adjusted_gap)+$drawing_counter*$GAP;
		
		# calculate gene coordinates 
		my $x_start=$keep_gene_start;
		my $x_end=$keep_gene_end;
		my $y_gene_line_start=($y_cent_exon-0.5*$gene_width+$adjusted_gap)+$drawing_counter*$GAP;
		my $y_gene_line_end = ($y_cent_exon+0.5*$gene_width+$adjusted_gap)+$drawing_counter*$GAP;
		
		# draw genes
		&DrawingRoutines::draw_a_rectangle($im_fh,$x_start,$y_gene_line_start,$x_end,$y_gene_line_end,$black);
		
		# draw gene image map
		&ImageMap::print_map($im_fh,$x_start,$y_ex_img_start,$x_end,$y_ex_img_end,$gene->id,"geneview.pl?gene_id=",$img_map_fh);	
	    }    
	}
	
	# draw seqeunce
	&DrawingRoutines::draw_a_rectangle($im_fh,$x_start,$y_start_seq +$drawing_counter*$GAP,$x_end,$y_end_seq+$drawing_counter*$GAP,$red);	
	


	# draw a scale bar
	&DrawingRoutines::draw_a_rectangle($im_fh,$x_start,$y_scale_bar +$drawing_counter*$GAP,$x_start+($scale_bar/$bases_per_pixel),$y_end_scale_bar+$drawing_counter*$GAP,$green);	

	# print scale bar bp
	$im_fh->string(gdSmallFont,$x_start+($scale_bar/$bases_per_pixel)+10,$y_scale_bar-5+$drawing_counter*$GAP,$scale_bar,$black);
	$im_fh->string(gdSmallFont,$x_start+($scale_bar/$bases_per_pixel)+40,$y_scale_bar-5+$drawing_counter*$GAP,"bp",$black);


	# print no of bases
	$im_fh->string(gdSmallFont,$x_start,$y_start_baseno+$drawing_counter*$GAP,$first_no_of_bases,$black);
	$im_fh->string(gdSmallFont,$x_img_len+$margin-30,$y_start_baseno+$drawing_counter*$GAP,$last_no_of_bases,$black);
	
	#print $last_no_of_bases," bp ...  done\n";
    }
    
    # end of an image map
    &ImageMap::print_map_end($img_map_fh);   
}



sub draw_one_gene 
{
    my ($gene,$arr_point,$img_map_fh,$gif_url) = @_;


   
# get image parameters
    my $im_fh=${$arr_point}[0];
my $y_img_len=${$arr_point}[1];
my $x_img_len=${$arr_point}[2];
my $margin=${$arr_point}[3];



# set colors
my $white = $im_fh->colorAllocate(255,255,255);
my $grey = $im_fh->colorAllocate(192,192,192);
my $lightgrey=$im_fh->colorAllocate(140,140,140);
my $darkgrey= $im_fh->colorAllocate(80,80,80);
my $white = $im_fh->colorAllocate(255,255,255);
my $black = $im_fh->colorAllocate(0,0,0);
my $red = $im_fh->colorAllocate(255,0,0);
my $blue = $im_fh->colorAllocate(0,0,255);
my $green= $im_fh->colorAllocate(0,255,0);

# get gene parameters
my @param=&get_gene_length($gene);
my $gene_len=@param[0];
my $gene_start=@param[1];
my $gene_end=@param[2];


# set drawing parameters
my $keep_gene_end;
my $keep_gene_start;
my $keep_end=0;
my $first_exon=1;
my $dist=0;
my $x_start=0;
my $exon_len=0;
my $x_end=0;
my $keep_end_coord=0;
my $half_dist;
my $line_break;

my $y_start=$y_img_len*0.5-15;
my $y_end=$y_img_len*0.5+15;

my $scale_bar=1000;
my $y_scale_bar=$y_img_len*0.5+30;

my $bases_pixel=$gene_len/$x_img_len;





    my @exons = sort { $a->start <=> $b->start } $gene->each_unique_Exon();

    foreach my $exon ( @exons ){

	my $strand= $exon->strand;


	# if it is a first exon
	if ($first_exon==1){
	    $dist=0;$x_start=$margin; $first_exon=0;
	}


	
	# if it IS NOT a first exon
	elsif ($first_exon==0){
	    $dist=$exon->start-$keep_end_coord;
	    $x_start=$keep_end+$dist/$bases_pixel;
	    $line_break=$keep_end+$half_dist/$bases_pixel;
	}
	

	$exon_len=($exon->end-$exon->start);
	
	$x_end=$x_start+$exon_len/$bases_pixel;


		

	$keep_end_coord=$exon->end;
	$keep_end=$x_end;


	# draw exons
	&DrawingRoutines::draw_a_rectangle($im_fh,$x_start,$y_start,$x_end,$y_end,$black);
	&ImageMap::print_map($im_fh,$x_start,$y_start,$x_end,$y_end,$exon->id,$gif_url,$img_map_fh);
	
	
    }

    my $line_start=$y_img_len*0.5-1;
    my $line_end=$y_img_len*0.5+1;

    &DrawingRoutines::draw_a_rectangle($im_fh,$margin,$line_start,$margin+($gene_len/$bases_pixel),$line_end,$black);


    # draw a scale bar
    &DrawingRoutines::draw_a_rectangle($im_fh,$margin,$y_scale_bar,$margin+($scale_bar/$bases_pixel),$y_scale_bar+1,$green);

    # print scale bar bp
    $im_fh->string(gdSmallFont,$margin+($scale_bar/$bases_pixel)+10,$y_scale_bar-5,$scale_bar,$black);
    $im_fh->string(gdSmallFont,$margin+($scale_bar/$bases_pixel)+40,$y_scale_bar-5,"bp",$black);


}






sub get_gene_length 
{
    my ($gene) = @_;

    my $gene_status=1;
    my $keep_gene_end;
    my $keep_gene_start;
    my $gene_len;

    foreach my $exon ( $gene->each_unique_Exon ){
	

	# if it's a first exon of a gene 
	if ($gene_status==1){
	    $keep_gene_start= $exon->start;
	    $keep_gene_end= $exon->end;
	    $gene_status=2;
	}
		    
	# if it's a subsequent exon
	if ($gene_status==2){
	    if($keep_gene_end< $exon->end){$keep_gene_end= $exon->end;}
	    if($keep_gene_start> $exon->end){$keep_gene_start= $exon->start;}
	}
	

	 $gene_len= $keep_gene_end-$keep_gene_start;
	

    }


    return my @param=($gene_len,$keep_gene_start,$keep_gene_end);

}








1;















