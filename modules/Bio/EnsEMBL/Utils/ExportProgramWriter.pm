package ExportProgramWriter;

use CreateTable;

sub create_export_program {
  (my $package_label,
   my $new_create_table,
   my $old_create_table) = @_;

  my $new_table_name = $new_create_table->get_table_name();
  my $old_table_name = $old_create_table->get_table_name();

  my $template;
  open(TEMPLATE, "ExportProgramTemplate.txt") || die "can not open file: $!";;
  while(<TEMPLATE>) {
    $template .= $_;
  }

  my @old_fields = @{$old_create_table->get_fields()};
  my $old_fields = join(", ", @old_fields);

  my @new_fields = @{$new_create_table->get_fields()};
  my $new_fields = join(", ", @new_fields);
  my $new_table_create_statement = $new_create_table->get_create_statement();


  my $old_max = scalar(@old_fields);
  my $new_max = scalar(@new_fields);
  my $max = ($old_max>$new_max) ? $old_max : $new_max;
  
  my $mapping = ($old_max!=$new_max) ? "# WARNING INCOMPATIBLE NUMBER OF FIELDS\n": "\n";
  $mapping .= "my \@selection_list;\n";
  for (my $i=0; $i<$max; ++$i) {
    my $old_field = @old_fields[$i];
    my $new_field = @new_fields[$i];
    my $too_many = ($i>=$new_max) ? "**WARNNING TOO MANY OLD FIELDS**" : "";
    my $hide = ($i>=$new_max) ? "#" : "";
    $mapping .= "$hide\@selection_list[$i] \t= \"" . $old_field  . "\"; # $new_field $too_many\n";
  }

  # search and replace KEYWORDS in template
  my $filename = $new_table_name . ".txt";
  $template=~ s/OUT_FILE_PATH/$filename/s;
  $template=~ s/OLD_TABLE_NAME/$old_table_name/s;
  $template=~ s/NEW_TABLE_CREATE_STATEMENT/$new_table_create_statement/s;
  $template=~ s/INPUT_ARRAY/$mapping/s;
  $template=~ s/CONNECTION_STRING/"DBI:mysql:homo_sapiens_core_110:ensrv3:3306",/s;
  $template=~ s/USER/"ensro"/s;
  $template=~ s/PASSWORD/""/s;
  my $table_statment =  $new_table_create_statement;
  $table_statment =~ s/^|\n/\n#/gs;
  $template=~ s/OUTPUT_TABLE_CREATE_STATEMENT/$table_statment/s;

  return $template;
}

1;
