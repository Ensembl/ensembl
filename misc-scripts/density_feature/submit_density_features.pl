#!/usr/bin/env perl

# The script lists databases/species which should have density features updated at the specified stage in the release cycle. There's an option for submitting a selected script. 

use strict;
use Getopt::Long;
use DBI qw( :sql_types );
use Switch;

use FindBin qw($Bin);
use vars qw($SERVERROOT);

BEGIN {
    $SERVERROOT = "$Bin";
}

my $gene_gc_path = "$SERVERROOT/..";

my $getdbs;
my $submit_script;
my $outdir;
my @host;
my @user;
my @pass;
my @port;

GetOptions( "getdbs|g",  \$getdbs,
	    "submit|s=s", \$submit_script,
            "outdir|o=s", \$outdir,
	    "host|h=s",\@host,
	    "user|u=s",\@user,
	    "pass|p=s",\@pass,
	    "port=s",\@port,
	    "help" ,     \&usage
	  );

my $host_count = @host;
my $user_count = @user;
my $pass_count = @pass;

usage() if (!defined @host || !defined @user || !defined @pass );  


my $port_count = @port;
my $host_string = join("\|",@host);

# if we have fewer user names specified than hosts copy user name from the first -u parameter
if ($user_count < $host_count) {
    for (my $i = $user_count; $i < $host_count; $i++) {
	push(@user,$user[0]);
    }
}

# if we have fewer passwords specified than hosts copy password from the first -p parameter
if ($pass_count < $host_count) {
    for (my $i = $pass_count; $i < $host_count; $i++) {
	push(@pass,$pass[0]);
    }
}

if ( (!defined @port) || ($port_count < $host_count) ) { 

    if (!defined $port[0]) {
	$port[0] = 3306;
    }

    for (my $i=1; $i<$host_count;$i++) {

	if (!defined $port[$i]) {
	    push(@port,$port[0]);
	}
    } 
}

my @hosts = qw(ens-staging1, ens-staging2);

if (!$outdir) { 
    $outdir = $ENV{'PWD'}; 
} else {
  #strip the final /
  $outdir =~ s/\/$//;  
  #test if the directory exists
  if (!-d $outdir) {
      die("Directory $outdir does not exist."); 
  }

}

# production/master database location:
my ( $mhost, $mport ) = ( 'ens-staging1', '3306' );
my ( $muser, $mpass ) = ( 'ensro',        undef );
my $mdbname = 'ensembl_production';

#connect to the production database
my $prod_dsn = sprintf( 'DBI:mysql:host=%s;port=%d;database=%s',
                     $mhost, $mport, $mdbname );
my $prod_dbh = DBI->connect( $prod_dsn, $muser, $mpass,
                          { 'PrintError' => 1, 'RaiseError' => 1 } );

#get the current release number

my ($current_release ) = $prod_dbh->selectrow_array('select max(db_release) from db where is_current = 1');
 


if ( !defined $submit_script ) {
    $getdbs = 1;
}

if ( defined $getdbs ) { 
# ask the user which stage of the release cycle we're in
  print <<CYCLE; 

The output will list density features scripts which can be run at a specified point in the release.
Select option 3 if variation dbs are handed over. The script will assume that no density features scripts
have been run yet and list all of them.
In case you would like to run some scripts earlier in the release, for instance when genebuild and genebuild
xrefs are complete, select option 1 and the output will only contain density features scripts which can be
run at this early point in the release. Later, when variation dbs are handed over, run this script again with
option 3 and bare in mind that it will list all density feature scripts which can be run, inluding
those that you may have already run before. Make sure that you keep the output from the density feature 
script runs so you don\'t submit those jobs again.  


Where in the release cycle are we (0,1,2,3)?

0 - None of the below - exit program
1 - Genebuild and genebuild xrefs complete (excluding projected xrefs), all gene healthchecks cleared
2 - Compara homologies handed over and core xref projections complete
3 - Variation dbs handed over

CYCLE

my $response = <>;

if ( !defined $response or $response !~ /[0-4]/ ) {
 
  print <<CYCLE2;

Please specify a valid option: 0, 1, 2 or 3:

Where in the release cycle are we (0,1,2,3)?

0 - None of the below - exit program
1 - Genebuild and genebuild xrefs complete (excluding projected xrefs), all gene healthchecks cleared
2 - Compara homologies handed over and core xref projections complete
3 - Variation dbs hand

CYCLE2

   $response = <>;
}


if ($response == 0) {
  $prod_dbh->disconnect;
  exit(0);
}


#dbs for new species or changed assembly
my @new_sp_assem;

#dbs with changed gene sequence

my @chg_seq;

#dbs with changed repeats
my @chg_repeats;

#dbs with new variation positions
my @chg_variation;

#core dbs which have a corresponding variation db
my @core_with_variation;

if ($response >= 1) {
    #get new dbs, or changed assembly  
    @new_sp_assem =  map { $_->[0] }  @{ $prod_dbh->selectall_arrayref("select distinct concat(full_db_name,'|',db_host) from db_list dl join db d using (db_id) where db_release = $current_release and db_type = 'core' and species_id not in (select distinct species_id from db where db_release <> $current_release and db_type = 'core') union
select distinct concat(full_db_name,'|',db_host) from db_list dl join db d using (db_id) where db_type = 'core' and is_current = 1 and species_id in (select distinct species_id from changelog_species cs join changelog c using (changelog_id) where release_id = $current_release and status not in ('cancelled','postponed') and assembly = 'Y')") };
    #get dbs with changed sequence
    @chg_seq = map { $_->[0] }  @{ $prod_dbh->selectall_arrayref("select distinct concat(full_db_name,'|',db_host) from db_list dl join db d using (db_id) where db_type = 'core' and is_current = 1 and species_id in (select distinct species_id from changelog_species cs join changelog c using (changelog_id) where release_id = $current_release and status not in ('cancelled','postponed') and gene_set = 'Y')") };

    #get dbs with changed repeats
    @chg_repeats =  map { $_->[0] }  @{ $prod_dbh->selectall_arrayref("select distinct concat(full_db_name,'|',db_host) from db_list dl join db d using (db_id) where db_type = 'core' and is_current = 1 and species_id in (select distinct species_id from changelog_species cs join changelog c using (changelog_id) where release_id = $current_release and status not in ('cancelled','postponed') and repeat_masking = 'Y')") };


    print "1. Density features scripts which can be run when qenebuild and genebuild xrefs (excluding projected xrefs) are complete and all gene healthchecks are cleared: \n\n";

    print "gene_gc.pl - run on all core databases (use the commands below or script submit_density_features.pl -submit gene_gc):\n";

    for (my $i=0; $i<$host_count;$i++) {
	print "\nbsub -q normal -J genegc_stats -M2000000 -R'select[mem>2000] rusage[mem=2000]' -oo $outdir/core_dbs_$current_release"."_".$host[$i]."_genegc.out -eo $outdir/core_dbs_$current_release" ."_".$host[$i]."_genegc.err perl $gene_gc_path/gene_gc.pl -h ".$host[$i]." -port ".$port[$i]." -u ".$user[$i]." -p ".$pass[$i]." -pattern 'core_$current_release'\n";
    }

    print "\npercent_gc_calc.pl – run on core databases for new species, or where sequence or assembly have changed (db names will be stored in file $outdir/percent_gc_data.txt, to submit run submit_density_features.pl -submit percent_gc -h ens-staging1 -h ens-staging2 -u ensadmin -p xxxx): \n";
 
    my %array_union = ();
    foreach my $element (@new_sp_assem, @chg_seq) { $array_union{$element}++ }
    my @dbnames_hosts = sort(keys %array_union); 

    my $file_path = "$outdir/percent_gc_data.txt";
 
    open(DATAFILE, ">$file_path") or die("Failed to open file $file_path for writing\n"); 
    foreach my $dbname_host (@dbnames_hosts) {
	my ($db_name, $host) = split(/\|/,$dbname_host);
	if ( ($db_name =~ /core_/) && ( $host_string =~ /$host/) ) {	    
	    print DATAFILE $db_name."\t".$host."\n";
	    print $db_name."\t".$host."\n";
	}	
    }
    close DATAFILE;

    print "\n\nrepeat_coverage_calc.pl – run on core databases for new species, or where sequence, assembly or repeats have changed (db names will be stored in file $outdir/repeat_coverage_data.txt, to submit run submit_density_features.pl -submit repeat_coverage  -h ens-staging1 -h ens-staging2 -u ensadmin -p xxxx): \n";   

    foreach my $element (@chg_repeats) { $array_union{$element}++ }
    @dbnames_hosts = sort(keys %array_union); 

    my $file_path = "$outdir/repeat_coverage_data.txt";
    open(DATAFILE, ">$file_path") or die("Failed to open file $file_path for writing\n"); 
    foreach my $dbname_host (@dbnames_hosts) {
	my ($db_name, $host) = split(/\|/,$dbname_host);
	if ( ($db_name =~ /core_/) && ( $host_string =~ /$host/) ) {	    
	    print DATAFILE $db_name."\t".$host."\n";
	    print $db_name."\t".$host."\n";
	}	
    }
    close DATAFILE;

} 


if ($response >= 2) {
    print "\n\n2. Density features scripts which can be run when Compara homologies are handed over and core xref projections are complete:\n\n";
    print "gene_density_calc.pl - run on all core dbs (use the commands below or script submit_density_features.pl -submit gene_density  -h ens-staging1 -h ens-staging2 -u ensadmin -p xxxx)\n";
    
    for (my $i=0; $i<$host_count;$i++) {
	print "\nbsub -q normal -J gene_density -M2000000 -R'select[mem>2000] rusage[mem=2000]' -oo $outdir/core_dbs_$current_release"."_".$host[$i]."_gene.out -eo $outdir/core_dbs_$current_release"."_".$host[$i]."_gene.err perl $SERVERROOT/gene_density_calc.pl -h ".$host[$i]." -port ".$port[$i]." -u ".$user[$i]." -p ".$pass[$i]." -pattern 'core_$current_release'\n";
    }

    print "\n\nseq_region_stats.pl (gene stats option only) - run on all core databases (use the commands below or script submit_density_features.pl -submit seq_region_stats_gene  -h ens-staging1 -h ens-staging2 -u ensadmin -p xxxx)\n";
    for (my $i=0; $i<$host_count;$i++) {
	print "\nbsub -q normal -J seqreg_stats_gene -M2000000 -R'select[mem>2000] rusage[mem=2000]' -oo $outdir/core_dbs_$current_release"."_".$host[$i]."_seqreg_gene.out -eo $outdir/core_dbs_$current_release"."_".$host[$i]."_seqreg_gene.err perl $SERVERROOT/seq_region_stats.pl -h ".$host[$i]." -port ".$port[$i]." -u ".$user[$i]." -p ".$pass[$i]." -pattern 'core_$current_release' -s gene\n";
    }
}


if ($response == 3) {

    print "\n\n3. Density features scripts which can be run when Variation dbs are handed over:\n";

    print "\nvariation_density.pl - run for new species or where the core assembly has changed, or if there are any changes to variation positions in the variation database (species will be stored in file $outdir/variation_density_data.txt, to submit run submit_density_features.pl -submit variation_density  -h ens-staging1 -h ens-staging2 -u ensadmin -p xxxx):\n";

    #get species for new dbs or changed assembly or where variation positions have changed
    @core_with_variation =  map { $_->[0] }  @{ $prod_dbh->selectall_arrayref("select distinct concat(full_db_name,'|',db_host) from db_list dl join db d using (db_id) where db_release = $current_release and db_type = 'core' and species_id in (select distinct species_id from db where db_release = $current_release and db_type = 'variation');") };

    @chg_variation =  map { $_->[0] }  @{ $prod_dbh->selectall_arrayref("select distinct concat(full_db_name,'|',db_host) from db_list dl join db d using (db_id) where db_type = 'core' and is_current = 1 and species_id in (select distinct species_id from changelog_species cs join changelog c using (changelog_id) where release_id = $current_release and status not in ('cancelled','postponed') and variation_pos_changed = 'Y')") };

    my %core_with_variation = map { $_ => 1 } @core_with_variation;


    my @new_sp_assem_var;

    foreach my $element (@new_sp_assem) {

	if ( exists($core_with_variation{$element}) ) {
	    push (@new_sp_assem_var, $element);
	}
    }

    my %array_union = ();
    foreach my $element (@new_sp_assem_var, @chg_variation) { $array_union{$element}++ }
    my @dbnames_hosts = sort(keys %array_union); 
 
    my $file_path = "$outdir/variation_density_data.txt";
 
    open(DATAFILE, ">$file_path") or die("Failed to open file $file_path for writing\n"); 
    foreach my $dbname_host (@dbnames_hosts) {
	my ($db_name, $host) = split(/\|/,$dbname_host);
	if ( ($db_name =~ /([^\s]+)_core/) && ( $host_string =~ /$host/) ) {
	    $db_name =~ /([^\s]+)_core/;	    
	    print DATAFILE $1 ."\t".$host."\n";
	    print $1 ."\t".$host."\n";
	}	
    }
    close DATAFILE;
 
    print "\n\nseq_region_stats.pl (snp stats option only) - run on core databases for new species or if the assembly changed, or if the variation positions have changed in the corresponding variation db (db names will be stored in file $outdir/seq_region_stats_snp_data.txt, to submit run submit_density_features.pl -submit seq_region_stats_snp  -h ens-staging1 -h ens-staging2 -u ensadmin -p xxxx):\n";

    my $file_path = "$outdir/seq_region_stats_snp_data.txt";

    open(DATAFILE, ">$file_path") or die("Failed to open file $file_path for writing\n"); 
    foreach my $dbname_host (@dbnames_hosts) {
	my ($db_name, $host) = split(/\|/,$dbname_host);
	if ( ($db_name =~ /core_/) && ( $host_string =~ /$host/) ) {	    
	    print DATAFILE $db_name."\t".$host."\n";
	    print $db_name."\t".$host."\n";
	}	
    }
    close DATAFILE;

}

$prod_dbh->disconnect;

} else {

    $prod_dbh->disconnect;

#submit selected script
    my @cmd;
    my $data_file;
    my @print_message;
    my $error;
    my %error_message;
    my $queue;
    my $file_name_end;
    my $script;
    my $script_title;
    my $job_name;
    my $option;
    switch ($submit_script) {
	case 'gene_gc' {
	    for (my $i=0; $i<$host_count;$i++) {
		push(@cmd, "bsub -q normal -J genegc_stats -M2000000 -R'select[mem>2000] rusage[mem=2000]' -oo $outdir/core_dbs_$current_release"."_". $host[$i] ."_genegc.out -eo $outdir/core_dbs_$current_release"."_" .$host[$i]. "_genegc.err perl $gene_gc_path/gene_gc.pl -h ".$host[$i]." -port ".$port[$i]." -u ".$user[$i]." -p ".$pass[$i]." -pattern 'core_$current_release'");
		push(@print_message,"Submitting gene GC calculation for host ".$host[$i]." to queue 'normal'. The output from this job goes to the file $outdir/core_dbs_$current_release"."_".$host[$i]."_genegc.out\n");
	    }

	}
	case 'percent_gc' {
	    $data_file = "$outdir/percent_gc_data.txt";
	    $queue = "normal";
	    $job_name = "gc_calc";
	    $file_name_end = "_gc";
	    $script = "$SERVERROOT/percent_gc_calc.pl";
	    $script_title = "percent GC calculation";
	    $option = " -d ";
	}
	case 'repeat_coverage' {
	    $data_file = "$outdir/repeat_coverage_data.txt";
	    $queue = "long";
	    $job_name = "repeat_cov";
	    $file_name_end = "_repeat";
	    $script = "$SERVERROOT/repeat_coverage_calc.pl";
	    $script_title = "repeat coverage calculation";
	    $option = " -d ";
	}
	case 'gene_density' {
	    for (my $i=0; $i<$host_count;$i++) {
		push(@cmd, "bsub -q normal -J gene_density -M2000000 -R'select[mem>2000] rusage[mem=2000]' -oo $outdir/core_dbs_$current_release"."_".$host[$i]."_gene.out -eo $outdir/core_dbs_$current_release"."_".$host[$i]."_gene.err perl $SERVERROOT/gene_density_calc.pl -h ".$host[$i]." -port ".$port[$i]." -u ".$user[$i]." -p ".$pass[$i]." -pattern 'core_$current_release'");
		push(@print_message,"Submitting gene density calculation for host ".$host[$i]." to queue 'normal'. The output from this job goes to the file ".$outdir."/core_dbs_".$current_release."_".$host[$i]."_gene.out\n");
	    }
	}
	case 'seq_region_stats_gene' {
	    for (my $i=0; $i<$host_count;$i++) {
		push(@cmd, "bsub -q normal -J seqreg_stats -M2000000 -R'select[mem>2000] rusage[mem=2000]' -oo $outdir/core_dbs_$current_release"."_".$host[$i]."_seqreg_gene.out -eo $outdir/core_dbs_$current_release"."_".$host[$i]."_seqreg_gene.err perl $SERVERROOT/seq_region_stats.pl -h ".$host[$i]." -port ".$port[$i]." -u ".$user[$i]." -p ".$pass[$i]." -pattern 'core_$current_release' -s gene");
		push(@print_message,"Submitting seq region gene stats for host ".$host[$i]." to queue 'normal'. The output from this job goes to the file ".$outdir."/core_dbs_".$current_release."_".$host[$i]."_seqreg_gene.out\n");
	    }
	}
	case 'variation_density' {
	    $data_file = "$outdir/variation_density_data.txt";
	    $queue = "normal";
	    $job_name = "var_density";
	    $file_name_end = "_var";
	    $script = "$SERVERROOT/variation_density.pl";
	    $script_title = "variation density calculation";
	    $option = " -s ";
	}
	case 'seq_region_stats_snp' {
	    $data_file = "$outdir/seq_region_stats_snp_data.txt";
	    $queue = "normal";
	    $job_name = "seqreg_stats_snp";
	    $file_name_end = "_seqreg_snp";
	    $script = "$SERVERROOT/seq_region_stats.pl";
	    $script_title = "seq region snp stats";
	    $option = " -s snp -d ";
	}
	else { usage(); }
    }

    if (defined $data_file) {
	    open(DATAFILE, "<$data_file") or die("Failed to open file $data_file for reading\n"); 
	    while( my $line = <DATAFILE> ) {
		chomp $line;
		my ($db_name, $host_name) = split(/\t/,$line);
		if ( $host_string =~ /$host_name/) {
		    #get user and password for host
		    my ( $index )= grep { $host[$_] =~ /$host_name/ } 0..$#host;
		    push(@cmd,  "bsub -q $queue -J $job_name -M2000000 -R'select[mem>2000] rusage[mem=2000]' -oo $outdir/$db_name$file_name_end.out -eo $outdir/$db_name$file_name_end.err perl $script -h $host_name -port ".$port[$index]." -u ".$user[$index]." -p ".$pass[$index].$option. $db_name);
		    push(@print_message,"Submitting ".$script_title." for ".$db_name ." on host ".$host_name." to queue '".$queue."'. The output from this job goes to the file $outdir/$db_name$file_name_end.out\n");
		}
		else {
		    $error = 1;
		    $error_message{$host_name} = "Host info for $host_name not found in the script's arguments. Please run the script again providing the user and password for host $host_name.\n";
		}
	
	    }
	    close DATAFILE;

    } 
  
    #submit jobs to the farm
    if ($error) {
	foreach my $host_in_error (sort(keys %error_message)) {
	    print $error_message{$host_in_error};
	} 

    } else {
        my $cmd_count = 0;
        foreach my $cmd (@cmd) {	
	   print $print_message[$cmd_count];
	   #for testing
	   #print "\n\n". $cmd . "\n\n";
	   system($cmd);
	   $cmd_count++;
        }
    }
}

sub usage {
  my $indent = ' ' x length($0);
  print <<EOF; exit(0);


Options -h -u -p are mandatory and need to be specified for at least one host. When using more than 
one host it\'s possible to leave out the user name and password for the second host and they will
be copied from the first host: e.g. -h ens-staging1 -h ens-staging2 -u ensadmin -p xxxx (user ensadmin
and password xxxx will be used for both ens-staging1 and ens-staging2).
The script lists databases/species which should have density features updated at the specified stage in the release cycle when using option -g, e.g. 

submit_density_features.pl -g -h ens-staging1 -h ens-staging2 -u ensadmin -p xxxx

There\'s an option for submitting a selected script, e.g.

submit_density_features.pl -s gene_gc -h ens-staging1 -h ens-staging2 -u ensadmin -p xxxx


Usage:

  $0 -h host [-h host] -u user [-u user] -p password [-p password] 
  $indent -port port_number [-port port_number]
  $indent [-g] [-s script_name]
  $indent [-o output_directory_path]
  $indent [-help]  


  -h|host              Database host (multiple hosts can be specified) 

  -u|user              Database user (each host needs a user specified, if multiple -h|host options
                       are given and fewer -u|user options are specified, the first user name will be used
                       for the hosts where no user name was given)

  -p|pass              User password (each host needs a password specified, if multiple -h|host options
                       are given and fewer -p|pass options are specified, the first password will be used
                       for the hosts where no password was given))

  -port                Database port (default 3306)  

  -g|getdbs            Use this option to generate input files for the -submit option

  -s|submit            Use this option to submit a density feature script to the farm:

                       gene_gc - the script will run on all core databases
		       
		       percent_gc - the script will run on dbs listed in [outdir]/percent_gc_data.txt
		       
                       repeat_coverage - the script will run on dbs listed in [outdir]/repeat_coverage_data.txt

		       gene_density - the script will run on all core databases

                       seq_region_stats_gene - the script will run on all core databases 

		       variation_density - the script will run for species listed in [outdir]/variation_density_data.txt

                       seq_region_stats_snp - the script will run on dbs listed in [outdir]/seq_region_stats_snp_data.txt

  -o|outdir            Path for farm job output and error files as well as data files generated by the script using                       
                       option -g(-getdbs) (current path if not specified)

  -help                This message


EOF

}

  
