use strict;

# Submits the display name and GO term projections as farm jobs
# Remember to check/set the various config optons

# ------------------------------ config -------------------------------
my $release = 47;

my $base_dir = "/lustre/work1/ensembl/gp1/projections/";

my $conf = "release_47.ini"; # registry config file, specifies Compara location

# location of other databases
my $host = "ens-staging";
my $port = 3306;
my $user = "ensadmin";
my $pass = "ensembl";

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
my $opts = "-conf $conf -host $host -user $user -port $port -pass $pass -version $release -release $release -quiet -nobackup";

my @names_1_1 = (["human", "chimp"            ],
		 ["human", "opossum"          ],
		 ["human", "dog"              ],
		 ["human", "cow"              ],
		 ["human", "macaque"          ],
		 ["human", "chicken"          ],
		 ["human", "xenopus"          ],
		 ["human", "guinea_pig"       ],
		 ["human", "armadillo"        ],
		 ["human", "small_hedgehog"   ],
		 ["human", "european_hedgehog"],
		 ["human", "cat"              ],
		 ["human", "elephant"         ],
		 ["human", "bat"              ],
		 ["human", "platypus"         ],
		 ["human", "rabbit"           ],
		 ["human", "galago"           ],
		 ["human", "european_shrew"   ],
		 ["human", "squirrel"         ],
		 ["human", "ground_shrew"     ],
		 ["mouse", "rat"              ]);

my @names_1_many = (["human", "zebrafish"  ],
		    ["human", "medaka"     ],
		    ["human", "tetraodon"  ],
		    ["human", "fugu"       ],
		    ["human", "stickleback"]);


my @go_terms = (["human",      "mouse"     ],
		["human",      "rat"       ],
		["human",      "dog"       ],
		["human",      "chicken"   ],
		["human",      "cow"       ],
		["human",      "chimp"     ],
		["human",      "macaque"   ],
		["human",      "guinea_pig"],
		["drosophila", "anopheles" ],
		["drosophila", "aedes"     ],
		["mouse",      "human"     ],
		["mouse",      "rat"       ],
		["mouse",      "dog"       ],
		["mouse",      "chicken"   ],
		["mouse",      "cow"       ],
		["rat",        "human"     ],
		["rat",        "mouse"     ],
		["danio",      "xenopus"   ],
		["danio",      "fugu"      ],
		["danio",      "tetraodon" ],
		["xenopus",    "danio"     ]);

my ($from, $to, $o, $e, $n);

# ----------------------------------------
# Display names

# 1:1
foreach my $pair (@names_1_1) {
  ($from, $to) = @$pair;
  $o = "$dir/names_${from}_$to.out";
  $e = "$dir/names_${from}_$to.err";
  $n = substr("n_${from}_$to", 0, 10); # job name display limited to 10 chars
  print "Submitting name projection from $from to $to\n";
  system "bsub -o $o -e $e -J $n perl project_display_xrefs.pl $opts -from $from -to $to -names -delete_names -no_database";
}

# 1:many
foreach my $pair (@names_1_many) {
  ($from, $to) = @$pair;
  $o = "$dir/names_${from}_$to.out";
  $e = "$dir/names_${from}_$to.err";
  $n = substr("n_${from}_$to", 0, 10);
  print "Submitting name projection from $from to $to (1:many)\n";
  system "bsub -o $o -e $e -J $n perl project_display_xrefs.pl $opts -from from -to $to -names -delete_names -no_database -one_to_many";
}

# ----------------------------------------
# GO terms

$opts .= " -nobackup";

foreach my $pair (@go_terms) {
  ($from, $to) = @$pair;
  $o = "$dir/go_${from}_$to.out";
  $e = "$dir/go_${from}_$to.err";
  $n = substr("g_${from}_$to", 0, 10);
  print "Submitting GO term projection from $from to $to\n";
  system "bsub -o $o -e $e -J $n perl project_display_xrefs.pl $opts -from $from -to $to -go_terms -delete_go_terms";
}

# ----------------------------------------



