use strict;
use GD;
use ImageMap;
use DrawingRoutines;
use GenDraw;


# read in arguments
my $contig_name = shift @ARGV;
my $bases_per_pixel = shift @ARGV;
if( !$contig_name || !$bases_per_pixel ){die"\n","*****  USAGE: contig_drawer.pl contig_name (e.g.AC000001.00010) bases_per_pixel (e.g. 100)*****","\n","\n";}

# open output HTML file
open(HTML,">output.html");

# build image 
my @parms=&GenDraw::build_image($bases_per_pixel,$contig_name);
my $im=shift @parms;

# draw genes
&GenDraw::draw_genes($im,$contig_name,$bases_per_pixel,\*HTML);

# draw image
open(GIF,">output.gif");
print GIF $im->gif;
close(GIF);






