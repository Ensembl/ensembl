use strict;
use warnings;

use StatMsg;

my $tr_split_count = 0;
my $tr_translates_count = 0;
my $tr_partial_translates_count = 0;
my $tr_entire_translates_count = 0;
my $tr_no_sequence_left_count = 0;
my $tr_no_cds_left_count = 0;
my $tr_doesnt_translate_count = 0;
my $tr_partial_doesnt_translate_count = 0;
my $tr_entire_doesnt_translate_count = 0;

my $ex_delete_count = 0;
my $ex_short_delete_count  = 0;
my $ex_medium_delete_count = 0;
my $ex_long_delete_count   = 0;
my $ex_entire_delete_count = 0;
my $ex_cds_delete_count    = 0;
my $ex_utr_delete_count    = 0;
my $ex_frameshift_delete_count = 0;
my $ex_long_frameshift_delete_count = 0;
my $ex_medium_frameshift_delete_count = 0;
my $ex_short_frameshift_delete_count = 0;

my $ex_insert_count        = 0;
my $ex_short_insert_count  = 0;
my $ex_medium_insert_count = 0;
my $ex_long_insert_count   = 0;
my $ex_cds_insert_count    = 0;
my $ex_utr_insert_count    = 0;
my $ex_frameshift_insert_count = 0;
my $ex_long_frameshift_insert_count = 0;
my $ex_medium_frameshift_insert_count = 0;
my $ex_short_frameshift_insert_count = 0;

my $ex_split_count    = 0;
my $ex_confused_count = 0;

while(my $code = <>) {
  chomp($code);

  if($code & StatMsg::TRANSCRIPT) {
    if($code & StatMsg::TRANSLATES) {
      $tr_translates_count++;

      if($code & StatMsg::PARTIAL) {
        $tr_partial_translates_count++;
      }
      if($code & StatMsg::ENTIRE) {
        $tr_entire_translates_count++;
      }
    }

    if($code & StatMsg::DOESNT_TRANSLATE) {
      $tr_doesnt_translate_count++;

      if($code & StatMsg::PARTIAL) {
        $tr_partial_doesnt_translate_count++;
      }
      if($code & StatMsg::ENTIRE) {
        $tr_entire_doesnt_translate_count++;
      }
    }


    if($code & StatMsg::NO_SEQUENCE_LEFT) {
      $tr_no_sequence_left_count++;
    }

    if($code & StatMsg::NO_CDS_LEFT) {
      $tr_no_cds_left_count++;
    }

    if($code & StatMsg::SPLIT) {
      $tr_split_count++;
    }
  }

  elsif($code & StatMsg::EXON) {
    if($code & StatMsg::DELETE) {
      $ex_delete_count++;

      if($code & StatMsg::LONG) {
        $ex_long_delete_count++;
      }
      elsif($code & StatMsg::MEDIUM) {
        $ex_medium_delete_count++;
      }
      elsif($code & StatMsg::SHORT) {
        $ex_short_delete_count++;
      }

      if($code & StatMsg::ENTIRE) {
        $ex_entire_delete_count++;
      }

      if($code & StatMsg::CDS) {
        $ex_cds_delete_count++;
      }

      if($code & StatMsg::UTR) {
        $ex_utr_delete_count++;
      }

      if($code & StatMsg::FRAMESHIFT) {
        $ex_frameshift_delete_count++;

        if($code & StatMsg::SHORT) {
          $ex_short_frameshift_delete_count++;
        }
        elsif($code & StatMsg::MEDIUM) {
          $ex_medium_frameshift_delete_count++;
        }
        elsif($code & StatMsg::LONG) {
          $ex_long_frameshift_delete_count++;
        }
      }
    }
    elsif($code & StatMsg::INSERT) {
      $ex_insert_count++;

      if($code & StatMsg::LONG) {
        $ex_long_insert_count++;
      }
      elsif($code & StatMsg::MEDIUM) {
        $ex_medium_insert_count++;
      }
      elsif($code & StatMsg::SHORT) {
        $ex_short_insert_count++;
      }

      if($code & StatMsg::CDS) {
        $ex_cds_insert_count++;
      }

      if($code & StatMsg::UTR) {
        $ex_utr_insert_count++;
      }

      if($code & StatMsg::FRAMESHIFT) {
        $ex_frameshift_insert_count++;

        if($code & StatMsg::SHORT) {
          $ex_short_frameshift_insert_count++;
        }
        elsif($code & StatMsg::MEDIUM) {
          $ex_medium_frameshift_insert_count++;
        }
        elsif($code & StatMsg::LONG) {
          $ex_long_frameshift_insert_count++;
        }
      }
    }

    if($code & StatMsg::CONFUSED) {
      $ex_confused_count++;
    }

    if($code & StatMsg::SPLIT) {
      $ex_split_count++;
    }

  }

}

print "Transcript Summary:\n";
print "  Split Events = $tr_split_count\n";
print "  Translates (Partial/Entire) = $tr_translates_count\n";
print "    Partial = $tr_partial_translates_count\n";
print "    Entire = $tr_entire_translates_count\n";
print "  Doesn't Translate = $tr_doesnt_translate_count\n";
print "    Partial = $tr_partial_doesnt_translate_count\n";
print "    Entire = $tr_entire_doesnt_translate_count\n";
print "  No Sequence Left = $tr_no_sequence_left_count\n";
print "  No CDS Left = $tr_no_cds_left_count\n";

print "\nExon Summary:\n";
print "  Split Events = $ex_split_count\n";
print "  Inserts = $ex_insert_count\n";
print "    Long = $ex_long_insert_count\n";
print "    Medium = $ex_medium_insert_count\n";
print "    Short = $ex_short_insert_count\n";
print "    CDS = $ex_cds_insert_count\n";
print "      Frameshifting = $ex_frameshift_insert_count\n";
print "        Long = $ex_long_frameshift_insert_count\n";
print "        Medium = $ex_medium_frameshift_insert_count\n";
print "        Short = $ex_short_frameshift_insert_count\n";
print "    UTR = $ex_utr_insert_count\n";
print "  Deletes = $ex_delete_count\n";
print "    Entire = $ex_entire_delete_count\n";
print "    Long = $ex_long_delete_count\n";
print "    Medium = $ex_medium_delete_count\n";
print "    Short = $ex_short_delete_count\n";
print "    CDS = $ex_cds_delete_count\n";
print "      Frameshifting = $ex_frameshift_delete_count\n";
print "        Long = $ex_long_frameshift_delete_count\n";
print "        Medium = $ex_medium_frameshift_delete_count\n";
print "        Short = $ex_short_frameshift_delete_count\n";
print "    UTR = $ex_utr_delete_count\n";
print "  Confused Events = $ex_confused_count\n";
