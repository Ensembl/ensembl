use Bio::SeqIO;
my $ntfile = shift(@ARGV);
my $seqin  = Bio::SeqIO->new(-file => "$ntfile" , '-format' => 'Fasta');
while ( my $seq = $seqin->next_seq() ) {
    $seq->id =~ /ref\|(\S+)\|(\S+)/;
    my $id=$1;
    my $name = "$id.fa";
    my $seqout = Bio::SeqIO->new(-file => ">$name" , '-format' => 'Fasta');
    $seqout->write_seq($seq);
}
