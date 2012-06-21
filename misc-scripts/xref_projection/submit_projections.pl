use strict;

use Data::Dumper;
use Bio::EnsEMBL::ApiVersion qw/software_version/;

$Data::Dumper::Useqq=1;
$Data::Dumper::Terse = 1;
$Data::Dumper::Indent = 0;

# Submits the display name and GO term projections as farm jobs
# Remember to check/set the various config optons

# ------------------------------ config -------------------------------
my $release = software_version();


my $base_dir = "mydir";

my $conf = "release_${release}.ini"; # registry config file, specifies Compara location

# location of other databases

my @config = ( {
    '-host'       => 'HOST',
    '-port'       => 'PORT',
    '-user'       => 'USER',
    '-pass'       => 'PASS',
    '-db_version' => $release
  },
  {
    '-host'       => 'HOST',
    '-port'       => 'PORT',
    '-user'       => 'USER',
    '-pass'       => 'PASS',
    '-db_version' => $release
  } );

my $registryconf = Dumper(\@config);

# load limit for ens-staging MySQL instance above which jobs won't be started
my $limit = 200;

# -------------------------- end of config ----------------------------

# check that base directory exists
die ("Cannot find base directory $base_dir") if (! -e $base_dir);

# create release subdir if necessary
my $dir = $base_dir. $release;
if (! -e $dir) {
  mkdir $dir;
  print "Created $dir\n";
} else {
  print "Cleaning and re-using $dir\n";
  unlink <$dir/*.out>, <$dir/*.err>;
}

# common options
my $script_opts = "-conf '$conf' -registryconf '$registryconf' -version '$release' -release '$release' -quiet";

my $bsub_opts = "";
$bsub_opts .= "-M2000000 -R'select[mem>2000] rusage[mem=2000]'";

my %names_1_1;

######
# When editing xref projection lists below, remember to check the species is in 
# the execution order array that follows.
######
$names_1_1{'human'} =  [qw(
    alpaca
    anolis
    armadillo
    bushbaby
    cat
    chicken
    chimp
    coelacanth
    cow
    dolphin
    elephant
    gibbon
    gorilla
    ground_shrew
    guinea_pig
    horse
    hyrax
    macaque
    marmoset
    megabat
    microbat
    mouse_lemur
    opossum
    orang_utan
    panda
    pig
    pika
    platypus
    sloth
    squirrel
    rabbit
    tarsier
    tasmanian_devil
    tenrec
    tree_shrew
    turkey
    wallaby
    western_european_hedgehog
    xenopus
    zebrafinch
    )];

$names_1_1{'mouse'} = [qw(
    kangaroo_rat
    rat    
)];

my %names_1_many;
$names_1_many{'human'} = [qw(
    cod
    fugu
    lamprey
    medaka
    stickleback
    tetraodon
    tilapia
    zebrafish
)];

my %go_terms;
$go_terms{'human'} = [qw(
    alpaca
    anolis
    armadillo
    bushbaby
    cat
    chicken
    chimp
    cow
    dog
    dolphin
    elephant
    gibbon
    gorilla
    ground_shrew
    guinea_pig
    horse
    hyrax
    kangaroo_rat
    macaque
    marmoset
    megabat
    microbat
    mouse
    mouse_lemur
    opossum
    orang_utan
    panda
    pig
    platypus
    sloth
    squirrel
    pika
    rabbit
    rat
    tarsier
    tasmanian_devil
    tenrec
    tree_shrew
    turkey
    wallaby
    western_european_hedgehog
    zebrafinch
)];
$go_terms{'mouse'} = [qw(
    alpaca
    anolis
    armadillo
    bushbaby
    cat
    chicken
    chimp
    cow
    dog
    dolphin
    elephant
    gorilla
    ground_shrew
    guinea_pig
    horse
    hyrax
    human
    kangaroo_rat
    macaque
    marmoset
    megabat
    microbat
    mouse_lemur
    opossum
    orang_utan
    panda
    pig
    pika
    platypus
    rabbit
    rat
    sloth
    squirrel
    tarsier
    tasmanian_devil
    tenrec
    tree_shrew
    turkey
    wallaby
    western_european_hedgehog
    zebrafinch
)];
$go_terms{'rat'} = [qw(
    human
    mouse
)];
$go_terms{'zebrafish'} = [qw(
    cod
    coelacanth
    fugu
    lamprey
    stickleback
    tetraodon
    tilapia
    xenopus
)];
$go_terms{'xenopus'} = [qw(zebrafish)];

# order to run projections in, just in case they are order-sensitive.
my @execution_order = ( 'human', 'mouse', 'rat', 'zebrafish');
# except of course order is irrelevant to the job queue. Consider provisional for when
# someone desires jobs that wait for others to finish.


# ----------------------------------------
# Display names

print "Deleting projected names (one to one)\n";
foreach my $species (keys %names_1_1) {
    foreach my $to (@{$names_1_1{$species}}) {
        system "perl project_display_xrefs.pl $script_opts -to $to -delete_names -delete_only\n";
    };
}

# 1:1
foreach my $from (@execution_order) {
    if (not exists($names_1_1{$from})) {next;}
    foreach my $to (@{$names_1_1{$from}}) {
        my $o = "$dir/names_${from}_$to.out";
        my $e = "$dir/names_${from}_$to.err";
        my $n = substr("n_${from}_$to", 0, 10); # job name display limited to 10 chars
        my $all = ($from eq "human") ? "" : "--all_sources"; # non-human from species -> use all sources
        print "Submitting name projection from $from to $to\n";
        system "bsub $bsub_opts -o $o -e $e -J $n perl project_display_xrefs.pl $script_opts -from $from -to $to -names -no_database $all\n";
    }
}

print "Deleting projected names (one to many)\n";
foreach my $from (keys %names_1_many) {
    foreach my $to (@{$names_1_many{$from}}) {
        system "perl project_display_xrefs.pl $script_opts -to $to -delete_names -delete_only\n";
    }
}

# 1:many
foreach my $from (@execution_order) {
    if (not exists($names_1_many{$from})) {next;}
    foreach my $to (@{$names_1_many{$from}}) {
        my $o = "$dir/names_${from}_$to.out";        
        my $e = "$dir/names_${from}_$to.err";
        my $n = substr("n_${from}_$to", 0, 10);
        print "Submitting name projection from $from to $to (1:many)\n";
        system "bsub $bsub_opts -o $o -e $e -J $n perl project_display_xrefs.pl $script_opts -from $from -to $to -names -no_database -one_to_many\n";
    }
}

# ----------------------------------------
# GO terms

$script_opts .= " -nobackup";

print "Deleting projected GO terms\n";
foreach my $from (keys %go_terms) {
    foreach my $to (@{$go_terms{$from}}) {
        system "perl project_display_xrefs.pl $script_opts -to $to -delete_go_terms -delete_only\n";
    }
}



foreach my $from (@execution_order) {
    if (not exists($go_terms{$from})) {next;}
    foreach my $to (@{$go_terms{$from}}) {
        my $o = "$dir/go_${from}_$to.out";
        my $e = "$dir/go_${from}_$to.err";
        my $n = substr("g_${from}_$to", 0, 10);
        print "Submitting GO term projection from $from to $to\n";
        system "bsub $bsub_opts -q long -o $o -e $e -J $n perl project_display_xrefs.pl $script_opts -from $from -to $to -go_terms\n";
    }
}


# ----------------------------------------



