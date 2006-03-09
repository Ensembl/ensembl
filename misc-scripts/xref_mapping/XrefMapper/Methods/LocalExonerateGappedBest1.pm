package XrefMapper::Methods::LocalExonerateGappedBest1;

use strict;
use File::Basename;

use XrefMapper::Methods::ExonerateBasic;

use vars '@ISA';

@ISA = qw{XrefMapper::Methods::ExonerateBasic};

our $exonerate_path = 'exonerate'; # Make sure exonerate binary is in
                                   # $PATH

sub options {

  return ('--model', 'affine:local', '--subopt', 'no', '--bestn', '1');

}

sub query_identity_threshold {

  #return 90;
  return 50;
}

sub target_identity_threshold {

  #return 90;
  return 50;
}

sub submit_exonerate{
  my ($self, $query, $target, $root_dir, @options) = @_;

  # Exonerate can run individual chunks of a large job.
  # Determine the number of chunks that will be used, and add to total
  my $num_jobs 
      = XrefMapper::Methods::ExonerateBasic::calculate_num_jobs($query);
  $self->jobcount($self->jobcount()+$num_jobs);

  # Determine the output filename
  my $queryfile = basename($query);
  my $targetfile = basename($target);
  my $query_alphabet = $queryfile =~ /peptide/ ? 'peptide' : 'dna';
  my $outfile = join( '_',
                      $self->get_class_name() , 
                      $query_alphabet, 
                      q[%s.map] # Substitute the chunk ID
                      );

  # Build a template representing the exonerate command
  my $options_str = join(" ", @options);
  my $exonerate_command = join
      ( " \\\n ", 
        $exonerate_path, 
        qq[--target $target],
        qq[--query $query], 
        qq[--querychunktotal $num_jobs],
        q [--querychunkid %s], # Substitute the chunk ID
        q [--showvulgar false],
        q [--showalignment FALSE],
        q [--ryo "%s"],
        $options_str );
  my $ryo_format = 'xref:%qi:%ti:%ei:%ql:%tl:%qab:%qae:%tab:%tae:%C:%s\n';
  my $shell_command = qq[%s \\\n| grep '^xref' \\\n> %s/%s];

  # Run exonerate for each job chunk
  foreach( my $i=1; $i<=$num_jobs; $i++ ){
    my $command = sprintf( $shell_command, 
                           sprintf($exonerate_command, $i, $ryo_format),
                           $root_dir,
                           sprintf($outfile, $i) );
    my $status = system( $command );
    unless( $status == 0 ){ warn( "Bad exit status ($?) from:\n $command" )}
  }
}

sub submit_depend_job{
  die;
}

1;
