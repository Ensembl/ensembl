use strict;

# Submits the display name and GO term projections as farm jobs
# Remember to check/set the various config optons

# ------------------------------ config -------------------------------
my $release = 51;

my $base_dir = "/lustre/scratch1/ensembl/gp1/projections/";

my $conf = "release_51.ini"; # registry config file, specifies Compara location

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
my $script_opts = "-conf $conf -host $host -user $user -port $port -pass $pass -version $release -release $release -quiet";

#my $bsub_opts = "-R'select[myens-staging<200]'";
my $bsub_opts = "";

my @names_and_go = (["human", "alpaca"           ],
		    ["human", "armadillo"        ],
		    ["human", "bushbaby"         ],
		    ["human", "cat"              ],
		    ["human", "chicken"          ],
		    ["human", "chimp"            ],
		    ["human", "cow"              ],
		    ["human", "dog"              ],
		    ["human", "dolphin"          ],
		    ["human", "elephant"         ],
		    ["human", "ground_shrew"     ],
		    ["human", "guinea_pig"       ],
		    ["human", "horse"            ],
		    ["human", "hyrax"            ],
		    ["human", "kangaroo_rat"     ],
		    ["human", "squirrel"         ],
		    ["human", "tarsier"          ],
		    ["human", "tenrec"           ],
		    ["human", "tree_shrew"       ],
		    ["human", "western_european_hedgehog"],
		    ["human", "xenopus"          ],
		    ["mouse", "kangaroo_rat"     ],
		    ["mouse", "rat"              ]);

my @names_1_many = (["human", "zebrafish"  ],
		    ["human", "medaka"     ],
		    ["human", "tetraodon"  ],
		    ["human", "fugu"       ],
		    ["human", "stickleback"]);

my @go_only = (["human",      "mouse"          ],
	       ["human",      "rat"            ],
	       ["mouse",      "alpaca"         ],
	       ["mouse",      "armadillo"      ],
	       ["mouse",      "bushbaby"       ],
	       ["mouse",      "cat"            ],
	       ["mouse",      "chicken"        ],
	       ["mouse",      "chimp"          ],
	       ["mouse",      "cow"            ],
	       ["mouse",      "dog"            ],
	       ["mouse",      "dolphin"        ],
	       ["mouse",      "elephant"       ],
	       ["mouse",      "ground_shrew"   ],
	       ["mouse",      "guinea_pig"     ],
	       ["mouse",      "human"          ],
	       ["mouse",      "hyrax"          ],
	       ["mouse",      "kangaroo_rat"   ],
	       ["mouse",      "macaque"        ],
	       ["mouse",      "megabat"        ],
	       ["mouse",      "microbat"       ],
	       ["mouse",      "mouse_lemur"    ],
	       ["mouse",      "opossum"        ],
	       ["mouse",      "orang_utan"     ],
	       ["mouse",      "pika"           ],
	       ["mouse",      "platypus"       ],
	       ["mouse",      "rabbit"         ],
	       ["mouse",      "rat"            ],
	       ["mouse",      "squirrel"       ],
	       ["mouse",      "tarsier"        ],
	       ["mouse",      "tenrec"         ],
	       ["mouse",      "tree_shrew"     ],
	       ["mouse",      "western_european_hedgehog"],
	       ["mouse",      "horse"          ],
	       ["rat",        "human"          ],
	       ["rat",        "mouse"          ],
	       ["drosophila", "anopheles"      ],
	       ["drosophila", "aedes"          ],
	       ["danio",      "xenopus"        ],
	       ["danio",      "fugu"           ],
	       ["danio",      "tetraodon"      ],
	       ["xenopus",    "danio"          ]);

my ($from, $to, $o, $e, $n);

# ----------------------------------------
# Display names and GO terms

foreach my $pair (@names_and_go) {
  ($from, $to) = @$pair;
  $o = "$dir/names_${from}_$to.out";
  $e = "$dir/names_${from}_$to.err";
  $n = substr("n_${from}_$to", 0, 10); # job name display limited to 10 chars
  my $all = ($from eq "human") ? "" : "--all_sources"; # non-human from species -> use all sources
  print "Submitting name and GO term projection from $from to $to\n";
  system "bsub $bsub_opts -o $o -e $e -J $n perl project_display_xrefs.pl $script_opts -from $from -to $to -names -delete_names -go_terms -delete_go_terms -no_database $all";
}

# 1:many name projections
foreach my $pair (@names_1_many) {
  ($from, $to) = @$pair;
  $o = "$dir/names_${from}_$to.out";
  $e = "$dir/names_${from}_$to.err";
  $n = substr("n_${from}_$to", 0, 10);
  print "Submitting name projection from $from to $to (1:many)\n";
  system "bsub $bsub_opts -o $o -e $e -J $n perl project_display_xrefs.pl $script_opts -from $from -to $to -names -delete_names -no_database -one_to_many";
}

# ----------------------------------------
# GO terms only

$script_opts .= " -nobackup";

foreach my $pair (@go_only) {
  ($from, $to) = @$pair;
  $o = "$dir/go_${from}_$to.out";
  $e = "$dir/go_${from}_$to.err";
  $n = substr("g_${from}_$to", 0, 10);
  print "Submitting GO term projection from $from to $to\n";
  system "bsub $bsub_opts -q long -o $o -e $e -J $n perl project_display_xrefs.pl $script_opts -from $from -to $to -go_terms -delete_go_terms";
}

# ----------------------------------------



