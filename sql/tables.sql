use pog;

create table exon(
    id         VARCHAR(40) not null,  primary key(id), 
    contig     VARCHAR(40) not null, 
    version    VARCHAR(10) not null, 
    created    DATE        not null, 
    modified   DATE        not null, 
    start      INT(10)     not null, 
    end        INT(10)     not null, 
    strand     VARCHAR(1)  not null, 
    key id_contig (id,contig));

create table transcript(
    id         VARCHAR(40) not null, primary key(id), 
    geneid     VARCHAR(40) not null, 
    key id_geneid (id,geneid));

create table featurehit(
    id         VARCHAR(40) not null, primary key(id), 
    name       VARCHAR(40) not null, 
    start      INT(10)     not null, 
    end        INT(10)     not null,
    strand     VARCHAR(1)  not null, 
    type       VARCHAR(40) not null);

create table feature_bridge(
    feature    VARCHAR(40) not null,
    transcript VARCHAR(40) not null,
    rank       INT,
    type       VARCHAR(40));

create table featureset(
    id         VARCHAR(40) not null,
    feature    VARCHAR(40) not null, 
    rank       INT,
    key feature_id (feature,id));

create table feature(
    id         VARCHAR(40) not null, 
    contig     VARCHAR(40) not null, 
    start      INT(10)     not null, 
    end        INT(10)     not null, 
    hstart     INT(10), 
    hend       INT(10),
    score      INT(10)     not null, 
    strand     VARCHAR(1), 
    analysis   VARCHAR(40) not null, 
    featureset   VARCHAR(40), 
    key id_contig (id,contig),key(start));

create table gene(
    id         VARCHAR(40) not null, primary key(id), 
    version    VARCHAR(40) not null, 
    created    DATE        not null, 
    modified   DATE        not null);

create table contig(
    id         VARCHAR(40) not null, primary key(id),
    clone      VARCHAR(40) not null,
    mapbin     VARCHAR(40) not null, 
    start      INT(10), 
    end INT(10));

create table mapbin(
    id         VARCHAR(40) not null, primary key(id), 
    chromosome VARCHAR(2)  not null);

create table analysis(
    id         VARCHAR(40) not null, primary key(id),
    db         VARCHAR(40),
    db_version VARCHAR(5), 
    program    VARCHAR(40) not null, 
    program_version VARCHAR(5), 
    method     VARCHAR(40) not null);

create table exon_bridge(
    exon       VARCHAR(40) not null, 
    transcript VARCHAR(40) not null, 
    rank       INT(10)     not null, 
    key exon_transcript (exon,transcript));
