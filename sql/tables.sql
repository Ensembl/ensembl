use pog;

create table exon(id VARCHAR(40) not null,primary key(id), contig VARCHAR(40) not null, version VARCHAR(10) not null, created DATE not null, modified DATE not null, start INT(10) not null, end INT(10) not null, strand VARCHAR(1) not null, key id_contig (id,contig));

create table transcript(id VARCHAR(40) not null,primary key(id), geneid VARCHAR(40) not null, key id_geneid (id,geneid));

create table homolhit(id VARCHAR(40) not null,primary key(id), name VARCHAR(40) not null, start INT(10) not null, end INT(10) not null,strand VARCHAR(1) not null, type VARCHAR(40) not null);

create table homolset (homolid VARCHAR(40) not null, setid VARCHAR(40) not null, key homol_set(homolid,setid));

create table homol(id VARCHAR(40) not null, contig VARCHAR(40) not null, start INT(10) not null, end INT(10) not null, hstart INT(10) not null, hend INT(10) not null,score INT(10) not null, strand VARCHAR(1),analysis VARCHAR(40) not null, homolset VARCHAR(40), key id_contig (id,contig),key(start));

create table gene(id VARCHAR(40) not null, primary key(id), version VARCHAR(40) not null, created DATE not null, modified DATE not null);

create table contig(id VARCHAR(40) not null, clone VARCHAR(40) not null,primary key(id), map VARCHAR(40) not null, start INT(10), end INT(10));

create table mapbin(id VARCHAR(40) not null,primary key(id), chromosome VARCHAR(2) not null);

create table analysis(id VARCHAR(40) not null, primary key(id),db VARCHAR(40),db_version VARCHAR(5), program VARCHAR(40) not null, program_version VARCHAR(5), method VARCHAR(40) not null);

create table exon_bridge(exon VARCHAR(40) not null, transcript VARCHAR(40) not null, rank INT(10) not null, key exon_transcript (exon,transcript));

create table homol_bridge(homol VARCHAR(40) not null, homolhit VARCHAR(40) not null, rank INT(10) not null);

create table test(id int  not null auto_increment, primary key(id), text VARCHAR(20));
