package GenDraw;

use strict;
use GD;
use Bio::EnsEMBL::DBSQL::Obj;
use ImageMap;
use DrawingRoutines;

 



sub build_image
{
    my($bases_per_pixel,$contig_obj)=@_;
    
    my $GAP=50;
    my $x_img_len=600;
    my $repeats = &GenDraw::adjust_image_size($bases_per_pixel,$x_img_len,$contig_obj);
    my $y_img_len=$repeats*$GAP;
    my $margin=30;
    my $im = new GD::Image($x_img_len+2*$margin,$y_img_len);
    
    return my @param= ($im,$GAP,$x_img_len,$margin);
}



sub adjust_image_size 
{
    my ($bases_per_pixel,$x_img_len,$contig_obj)=@_;
    
    my $dbuser = 'ensro';
    my $dbname = 'ensdev';
    my $host = 'obi-wan';
    my $dbpass = undef;

    my $seq = $contig_obj->seq();
    my $repeats =$seq->seq_len/($x_img_len*$bases_per_pixel)+2;    

    return $repeats;
}





sub draw_genes 
{
    my ($im_fh,$contig_obj,$bases_per_pixel,$img_map_fh,$gif_url,$genearray) = @_;

# get image parameters
    my @parms=&build_image($bases_per_pixel,$contig_obj);
    my $GAP=@parms[1];
    my $x_img_len=@parms[2];
    my $margin=@parms[3];
    
# set colours
    my $white = $im_fh->colorAllocate(255,255,255);
    my $grey = $im_fh->colorAllocate(192,192,192);# background
	my $lightgrey=$im_fh->colorAllocate(140,140,140);
    my $darkgrey= $im_fh->colorAllocate(80,80,80);
    my $white = $im_fh->colorAllocate(255,255,255);
    my $black = $im_fh->colorAllocate(0,0,0);
    my $red = $im_fh->colorAllocate(255,0,0);
    my $blue = $im_fh->colorAllocate(0,0,255);


# set the the vertical distance between each drawn object and the sequence
    my $center= ($GAP)/2;
    my $gap=0.1*$GAP;
    my $adjusted_gap=0;
    my $y_cent_seq=$center;
    my $n=2;
    my $y_cent_exon=$center+$n*$gap;
    my $y_cent_baseno=$center-6*$gap;
    my $y_cent_contig_name=-6*$gap;

# set width of each object
    my $seq_width=0.01*$GAP;
    my $exon_width=0.2*$GAP;
    my $gene_width=0.01*$GAP;
    my $y_start_seq=$y_cent_seq-0.5*$seq_width;
    my $y_start_exon=$y_cent_exon-0.5*$exon_width;
    my $y_start_baseno=$y_cent_baseno;
    my $y_end_seq=$y_start_seq+$seq_width;
    my $y_end_exon=$y_start_exon+$exon_width;
              
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
    $im_fh->string(gdMediumBoldFont,$x_img_len/2-50 ,20,$contig->id,$black);
   

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
	
	foreach my $gene ( @{$genearray} ) { 
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
	
	# print no of bases
	$im_fh->string(gdSmallFont,$x_start,$y_start_baseno+$drawing_counter*$GAP,$first_no_of_bases,$black);
	$im_fh->string(gdSmallFont,$x_img_len+$margin-30,$y_start_baseno+$drawing_counter*$GAP,$last_no_of_bases,$black);
	
	#print $last_no_of_bases," bp ...  done\n";
    }
    
    # end of an image map
    &ImageMap::print_map_end($img_map_fh);   
}

1;















