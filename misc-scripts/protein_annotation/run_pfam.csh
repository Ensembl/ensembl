#!/bin/tcsh

#test_tmp1 (you may change it) is the name of the file containing your peptides on wolverine (the computer connected to teh paracel box) I have harcoded it on purpose...see bellow
ssh -l paracel wolverine " cd mongin; cat > test_tmp1.in " < $1 || exit 1

#This command will load the protein file to the paracel box itself
ssh -l paracel wolverine " cd mongin; btkload src=test_tmp1.in dst=/fdf/build109/gm0/0/mongin seqtype=protein" || exit 1

#This command run btk hmm, this produces a blast like output. NB: test_tmp1.in is called test_tmp1.prot on the paracel box
ssh -l paracel wolverine "cd mongin; btk hmm database=/fdf/build109/gm0/0/mongin/test_tmp1.prot mhsp=data query=pfam_query_1.tbl -invert format=blast " > $1.out || exit 1

