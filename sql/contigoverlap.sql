

create table contigoverlap (
       contig_a varchar(40) DEFAULT '' NOT NULL,
       contig_b varchar(40) DEFAULT '' NOT NULL,
       contig_a_position int(10) unsigned,
       contig_b_position int(10) unsigned,
       overlap_type set('right2left','left2right','left2left','right2right'),
       PRIMARY_KEY(contig_a,contig_b)
  );
       
