use strict;
use GD;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use SeqContigDraw;
use Parameters;


my $image_param_ref=&Parameters::contig_image_par;
my $im=new GD::Image($image_param_ref->{x_img_len},$image_param_ref->{y_img_len});

my $dbuser = 'ensembl';
my $dbname = 'test500';
my $host = 'sol28';
my $dbpass = undef; 
my $locator = "Bio::EnsEMBL::DBSQL::Obj/host=$host;dbname=$dbname;user=$dbuser;pass=$dbpass";
my $db =  Bio::EnsEMBL::DBLoader->new($locator);

#if ($db){print "connect to the db\n";}

my $contig_id=shift @ARGV; $contig_id || die 'I need contig id e.g. AC004471.00001  ';
my $contig=$db->get_Contig($contig_id);



my $vc = Bio::EnsEMBL::DB::VirtualContig->new( -focuscontig => $contig,
                                              -focusposition => 25000,
                                              -ori => 25000,
                                              -left => 16000,
                                              -right => 1000
                                              );




#print "vcontig ",$contig->rawcontig_ids,"\n";

print "<body bgcolor=\"#FFFFFF\">\n";
print  "<IMG SRC=\"contig.gif\" BORDER=\"0\" USEMAP=\"#Ngs22\"><BR>\n";
print  "<MAP Name=\"Ngs22\">\n";


&SeqContigDraw::draw_contig_image($im,$vc);

open(GIF,">contig.gif");
print GIF $im->gif;
close(GIF);

print  "</map>\n";



