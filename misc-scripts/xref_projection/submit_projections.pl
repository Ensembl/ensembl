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

my @names_1_1 = (["human", "chimp"            ],
		 ["human", "opossum"          ],
		 ["human", "dog"              ],
		 ["human", "cow"              ],
		 ["human", "macaque"          ],
		 ["human", "chicken"          ],
		 ["human", "xenopus"          ],
		 ["human", "guinea_pig"       ],
		 ["human", "armadillo"        ],
		 ["human", "tenrec"           ],
		 ["human", "western_european_hedgehog"],
		 ["human", "cat"              ],
		 ["human", "elephant"         ],
		 ["human", "microbat"         ],
		 ["human", "platypus"         ],
		 ["human", "rabbit"           ],
		 ["human", "bushbaby"         ],
		 ["human", "ground_shrew"     ],
		 ["human", "squirrel"         ],
		 ["human", "tree_shrew"       ],
		 ["human", "pika"             ],
		 ["human", "mouse_lemur"      ],
		 ["human", "horse"            ],
		 ["human", "orang_utan"       ],
		 ["human", "dolphin"          ],
		 ["human", "hyrax"            ],
		 ["human", "megabat"          ],
		 ["human", "tarsier"          ],
		 ["human", "alpaca"           ],
		 ["mouse", "kangaroo_rat"     ],
		 ["mouse", "rat"              ]);

my @names_1_many = (["human", "zebrafish"  ],
		    ["human", "medaka"     ],
		    ["human", "tetraodon"  ],
		    ["human", "fugu"       ],
		    ["human", "stickleback"]);

my @go_terms = (["human",      "mouse"          ],
		["human",      "rat"            ],
		["human",      "dog"            ],
		["human",      "chicken"        ],
		["human",      "cow"            ],
		["human",      "chimp"          ],
		["human",      "macaque"        ],
		["human",      "guinea_pig"     ],
		["human",      "bushbaby"       ],
		["human",      "rabbit"         ],
		["human",      "cat"            ],
		["human",      "ground_shrew"   ],
		["human",      "western_european_hedgehog"],
		["human",      "microbat"       ],
		["human",      "armadillo"      ],
		["human",      "elephant"       ],
		["human",      "tenrec"         ],
		["human",      "opossum"        ],
		["human",      "platypus"       ],
		["human",      "squirrel"       ],
		["human",      "tree_shrew"     ],
		["human",      "pika"           ],
		["human",      "mouse_lemur"    ],
		["human",      "horse"          ],
		["human",      "orang_utan"     ],
		["human",      "dolphin"        ],
		["human",      "hyrax"          ],
		["human",      "megabat"        ],
		["human",      "tarsier"        ],
		["human",      "alpaca"         ],
		["human",      "kangaroo_rat"   ],
		["mouse",      "human"          ],
		["mouse",      "rat"            ],
		["mouse",      "dog"            ],
		["mouse",      "chicken"        ],
		["mouse",      "cow"            ],
		["mouse",      "chimp"          ],
		["mouse",      "macaque"        ],
		["mouse",      "guinea_pig"     ],
		["mouse",      "bushbaby"       ],
		["mouse",      "rabbit"         ],
		["mouse",      "cat"            ],
		["mouse",      "ground_shrew"   ],
		["mouse",      "western_european_hedgehog"],
		["mouse",      "microbat"       ],
		["mouse",      "armadillo"      ],
		["mouse",      "elephant"       ],
		["mouse",      "tenrec"         ],
		["mouse",      "opossum"        ],
		["mouse",      "platypus"       ],
		["mouse",      "squirrel"       ],
		["mouse",      "tree_shrew"     ],
		["mouse",      "pika"           ],
        	["mouse",      "horse"          ],
		["mouse",      "orang_utan"     ],
		["mouse",      "mouse_lemur"    ],
		["mouse",      "dolphin"        ],
		["mouse",      "hyrax"          ],
		["mouse",      "megabat"        ],
		["mouse",      "tarsier"        ],
		["mouse",      "alpaca"         ],
		["mouse",      "kangaroo_rat"   ],
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
# Display names

# 1:1
foreach my $pair (@names_1_1) {
  ($from, $to) = @$pair;
  $o = "$dir/names_${from}_$to.out";
  $e = "$dir/names_${from}_$to.err";
  $n = substr("n_${from}_$to", 0, 10); # job name display limited to 10 chars
  my $all = ($from eq "human") ? "" : "--all_sources"; # non-human from species -> use all sources
  print "Submitting name projection from $from to $to\n";
  system "bsub $bsub_opts -o $o -e $e -J $n perl project_display_xrefs.pl $script_opts -from $from -to $to -names -delete_names -no_database $all";
}

# 1:many
foreach my $pair (@names_1_many) {
  ($from, $to) = @$pair;
  $o = "$dir/names_${from}_$to.out";
  $e = "$dir/names_${from}_$to.err";
  $n = substr("n_${from}_$to", 0, 10);
  print "Submitting name projection from $from to $to (1:many)\n";
  system "bsub $bsub_opts -o $o -e $e -J $n perl project_display_xrefs.pl $script_opts -from $from -to $to -names -delete_names -no_database -one_to_many";
}

# ----------------------------------------
# GO terms

$script_opts .= " -nobackup";

foreach my $pair (@go_terms) {
  ($from, $to) = @$pair;
  $o = "$dir/go_${from}_$to.out";
  $e = "$dir/go_${from}_$to.err";
  $n = substr("g_${from}_$to", 0, 10);
  print "Submitting GO term projection from $from to $to\n";
  system "bsub $bsub_opts -q long -o $o -e $e -J $n perl project_display_xrefs.pl $script_opts -from $from -to $to -go_terms -delete_go_terms";
}

# ----------------------------------------



