

create table contigoverlap (
       contig_a varchar(40) DEFAULT '' NOT NULL,
       contig_b varchar(40) DEFAULT '' NOT NULL,
       contig_a_position int(10) unsigned,
       contig_b_position int(10) unsigned,
       contig_a_version int(10) unsigned,
       contig_b_version int(10) unsigned,
       type varchar(40) DEFAULT '' NOT NULL,
       overlap_type set('right2left','left2right','left2left','right2right'),
       PRIMARY KEY(contig_a,contig_a_version,contig_b,contig_b_version,type)
  );
       
